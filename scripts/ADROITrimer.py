#!/usr/bin/env python3
"""
ADROITrimer - Python 3 conversion of ADROITrimer1_0-3.R

Usage:
    python3 ADROITrimer.py <repair_sequence.fasta> <subtype> <muscle_exe>

Dependencies:
    pip install biopython pandas numpy
    MUSCLE aligner binary (v3/v5) must be available.
"""

import sys
import os
import re
import time
import string
import subprocess
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter


# ==================================================
# Helper functions
# ==================================================

def is_letter(ch: str) -> bool:
    """Return True if ch is an alphabetic character."""
    return ch.isalpha()


def read_fasta(path: str) -> tuple[list[str], list[list[str]]]:
    """
    Read a FASTA file and return (ids, alignment_matrix).
    alignment_matrix is a list of lists of single characters (uppercase).
    """
    records = list(SeqIO.parse(path, "fasta"))
    ids = [r.id for r in records]
    seqs = [list(str(r.seq).upper()) for r in records]
    return ids, seqs


def write_fasta(ids: list[str], seqs: list[list[str]], path: str) -> None:
    """Write a FASTA file from ids and per-sequence character lists."""
    records = []
    for sid, seq in zip(ids, seqs):
        records.append(SeqRecord(Seq("".join(seq)), id=sid, description=""))
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fasta")


def write_fasta_single(sid: str, seq: list[str], path: str) -> None:
    """Write a single-sequence FASTA file."""
    write_fasta([sid], [seq], path)


def consensus(ali: list[list[str]]) -> tuple[list[str], list[dict]]:
    """
    Compute a consensus sequence (freq >= 0.6) and per-column AA frequencies.
    Returns (consensus_seq_as_char_list, list_of_freq_dicts).
    """
    if not ali:
        return [], []
    n_cols = len(ali[0])
    n_rows = len(ali)
    cons_seq = []
    freq_list = []
    for col in range(n_cols):
        col_chars = [ali[row][col] for row in range(n_rows)]
        counter = Counter(col_chars)
        freq = {aa: cnt / n_rows for aa, cnt in counter.items()}
        freq_list.append(freq)
        best_aa, best_freq = max(freq.items(), key=lambda x: x[1])
        cons_seq.append(best_aa if best_freq >= 0.6 else "-")
    return cons_seq, freq_list


def parse_rsa(path: str) -> pd.DataFrame:
    """
    Parse an RSA file produced by NACCESS or similar tool.
    Returns a DataFrame with columns: atom, aa, chain, resno, sasa.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("RES ") or line.startswith("HEM "):
                chars = list(line)
                try:
                    atom = "".join(chars[0:3]).strip()
                    aa = "".join(chars[4:7]).strip()
                    chain = chars[8].strip()
                    resno_str = "".join(chars[9:13]).strip()
                    sasa_str = "".join(chars[36:41]).strip()
                    resno = int(resno_str) if resno_str else None
                    sasa = float(sasa_str) if sasa_str else None
                    rows.append({"atom": atom, "aa": aa, "chain": chain,
                                 "resno": resno, "sasa": sasa})
                except (IndexError, ValueError):
                    pass
    return pd.DataFrame(rows)


def parse_pdb_ca(path: str) -> pd.DataFrame:
    """
    Read a PDB file and extract CA atoms with chain, resno, and insert code.
    Returns DataFrame with columns: chain, resno, insert.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21]
                resno = int(line[22:26].strip())
                insert = line[26].strip() or None
                rows.append({"chain": chain, "resno": resno, "insert": insert})
    return pd.DataFrame(rows)


# ==================================================
# Main ADROITrimer function
# ==================================================

def adroi_trimer(
    file_reference_alignment: str,
    file_repair_sequence: str,
    subtype: str,
    exe_muscle: str,
    input_prefix: str = None,
    output_prefix: str = None,
) -> None:

    if input_prefix is None:
        input_prefix = f"input_{round(time.time())}"
    if output_prefix is None:
        output_prefix = f"output_{round(time.time())}"

    # --------------------------------------------------
    # Set variables
    print("Set variables")

    base = file_reference_alignment.replace(".fasta", "")
    base_subtype = f"{base}_subtype{subtype}.fasta"
    base_subtype_repair = base_subtype.replace(".fasta", f"_{input_prefix}.fasta")
    base_subtype_repair_muscle = base_subtype.replace(".fasta", f"_{input_prefix}_muscle.fasta")
    base_subtype_repair_muscle_no_hxb2 = base_subtype.replace(".fasta", f"_{input_prefix}_muscle_noHXB2.fasta")
    file_new_aln = base_subtype.replace(".fasta", f"_{output_prefix}.fasta")

    # --------------------------------------------------
    # Load raw HIV alignment and select for subtype of interest
    print("Load raw HIV alignment and select for subtype of interest")

    aln_ref_id, aln_ref_ali = read_fasta(file_reference_alignment)

    # Extract first component of the dot-separated ID for subtype label
    raw_subtypes = [sid.split(".")[0] for sid in aln_ref_id]

    # Keep only uppercase alphabetic characters and deduplicate
    vec_subtypes_clean = []
    for t in raw_subtypes:
        upper_chars = list(dict.fromkeys(ch for ch in t if is_letter(ch) and ch == ch.upper()))
        vec_subtypes_clean.append("".join(upper_chars))

    # --------------------------------------------------
    # Identify if sufficient number (>=50) of subtypes available
    print("Identify if sufficient number (>=50) of subtypes available")

    subtype_indices = [i for i, s in enumerate(vec_subtypes_clean) if s == subtype]

    if len(subtype_indices) >= 50:
        # Always regenerate subtype file (avoid stale cache)
        keep = [0] + subtype_indices
        write_fasta(
            [aln_ref_id[i] for i in keep],
            [aln_ref_ali[i] for i in keep],
            base_subtype,
        )
        file_ref_aln_subtype = base_subtype
    else:
        file_ref_aln_subtype = file_reference_alignment

    # Load subtype-only sequences
    aln_subtype_id, aln_subtype_ali = read_fasta(file_ref_aln_subtype)

    # --------------------------------------------------
    # Load test sequence
    print("Load test sequence")

    aln_test_id, aln_test_ali = read_fasta(file_repair_sequence)
    aln_test_ali_flat = aln_test_ali[0]  # single sequence

    # --------------------------------------------------
    # Align test sequence with subtype alignment
    print("Align test sequence with subtype alignment")

    n_cols = len(aln_subtype_ali[0])
    padded_test = aln_test_ali_flat + ["-"] * (n_cols - len(aln_test_ali_flat))
    aln_subtype_test_ali = aln_subtype_ali + [padded_test]
    aln_subtype_test_id = aln_subtype_id + aln_test_id

    write_fasta(aln_subtype_test_id, aln_subtype_test_ali, base_subtype_repair)

    # Align using MUSCLE with 2 iterations
    subprocess.run(
        [exe_muscle, "-in", base_subtype_repair,
         "-out", base_subtype_repair_muscle, "-maxiters", "2"],
        check=True,
    )

    # --------------------------------------------------
    # Load new subtype alignment including HXB2 and test sequence
    print("Load new subtype alignment including HXB2 and test sequence")

    aln_muscle_id, aln_muscle_ali = read_fasta(base_subtype_repair_muscle)

    # Remove the first HXB2 sequence; save it separately
    hxb2_indices = [i for i, sid in enumerate(aln_muscle_id) if "HXB2" in sid]
    hxb2_idx = hxb2_indices[0]
    vec_seq_hxb2 = aln_muscle_ali[hxb2_idx]

    ali_no_hxb2 = [s for i, s in enumerate(aln_muscle_ali) if i != hxb2_idx]
    id_no_hxb2 = [sid for i, sid in enumerate(aln_muscle_id) if i != hxb2_idx]

    write_fasta(id_no_hxb2, ali_no_hxb2, base_subtype_repair_muscle_no_hxb2)

    # Reload to stay consistent with R
    aln_no_hxb2_id, aln_no_hxb2_ali = read_fasta(base_subtype_repair_muscle_no_hxb2)

    # Extract aligned test sequence
    test_id_match = aln_test_id[0]
    test_idx_in_muscle = next(i for i, sid in enumerate(aln_muscle_id) if sid == test_id_match)
    vec_seq_test_aligned = aln_muscle_ali[test_idx_in_muscle]

    # --------------------------------------------------
    # Generate initial data frame
    print("Generate initial data frame")

    inserts = list(string.ascii_uppercase) + list(string.ascii_lowercase)
    rows = []
    ind_hxb2 = 0
    ind_letters = 0

    for i, hxb2_char in enumerate(vec_seq_hxb2):
        if hxb2_char != "-":
            ind_hxb2 += 1
            insert = " "
            ind_letters = 0
        else:
            ind_letters += 1
            insert = inserts[ind_letters - 1]

        rows.append({
            "no_aln": i + 1,          # 1-based column index
            "no_hxb2": ind_hxb2,
            "insert_hxb2": insert,
            "aa_hxb2": hxb2_char,
            "aa_newSequence": vec_seq_test_aligned[i],
        })

    df = pd.DataFrame(rows)

    # --------------------------------------------------
    # Calculate consensus sequence with default (0.6) cutoff
    print("Calculate consensus sequence with default (0.6) cutoff")

    cons_seq, freq_list = consensus(aln_no_hxb2_ali)

    # Build AA frequency dataframe (one row per alignment column)
    all_aas = sorted({aa for fd in freq_list for aa in fd})
    freq_rows = [{aa: fd.get(aa, 0.0) for aa in all_aas} for fd in freq_list]
    df_aa_freq = pd.DataFrame(freq_rows)

    # --------------------------------------------------
    # Data frame containing consensus frequencies
    print("Data frame containing consensus frequencies")

    df["aa_con"] = cons_seq
    for col in df_aa_freq.columns:
        df[col] = df_aa_freq[col].values

    # --------------------------------------------------
    # Upload PDB
    pdb_df = parse_pdb_ca("5fyl_trimer.pdb")
    pdb_protomer = pdb_df[pdb_df["chain"].isin(["A", "B"])].copy().reset_index(drop=True)

    # Upload RSAs
    rsa_df = parse_rsa("5fyl_trimer.rsa")
    rsa_protomer = rsa_df[rsa_df["chain"].isin(["A", "B"])].copy().reset_index(drop=True)
    rsa_protomer["sasa"] = pd.to_numeric(rsa_protomer["sasa"], errors="coerce")

    # --------------------------------------------------
    # Merge RSA with data frame
    df["RSA"] = np.nan
    df["no_hxb2"] = pd.to_numeric(df["no_hxb2"], errors="coerce")

    for _, atom_row in pdb_protomer.iterrows():
        resno = atom_row["resno"]
        insert = atom_row["insert"]

        if resno in (321, 322):
            continue

        # Find matching RSA row (same position in rsa_protomer)
        rsa_row = rsa_protomer[rsa_protomer["resno"] == resno]
        if rsa_row.empty:
            continue
        sasa_val = rsa_row.iloc[0]["sasa"]

        if insert is None or (isinstance(insert, float) and np.isnan(insert)):
            mask = (df["no_hxb2"] == resno) & (df["insert_hxb2"] == " ")
        else:
            mask = (df["no_hxb2"] == resno) & (df["insert_hxb2"] == insert)

        df.loc[mask, "RSA"] = sasa_val

    # --------------------------------------------------
    # Repair sequences
    print("Repair sequence")

    var_rsa_30pct = 0.3 * df["RSA"].max(skipna=True)
    vec_seq_repair = []
    vec_aa_new_seq_freq = [np.nan] * len(df)

    for i, row in df.iterrows():
        aa_new = row["aa_newSequence"]

        # Manually ignore positions 321 and 322
        if row["no_hxb2"] in range(321, 323):
            vec_seq_repair.append(aa_new)
            continue

        if row["no_hxb2"] in range(508, 512):
            vec_seq_repair.append(aa_new)
        elif aa_new == "-":
            vec_seq_repair.append(aa_new)
        else:
            aa_freq = row.get(aa_new, 0.0)
            vec_aa_new_seq_freq[i] = aa_freq

            if aa_freq < 0.02:
                aa_con = row["aa_con"]
                vec_seq_repair.append(aa_new if aa_con == "-" else aa_con)
            elif 0.02 <= aa_freq < 0.075:
                rsa = row["RSA"]
                if pd.isna(rsa):
                    vec_seq_repair.append(aa_new)
                else:
                    if rsa < var_rsa_30pct:
                        aa_con = row["aa_con"]
                        vec_seq_repair.append(aa_new if aa_con == "-" else aa_con)
                    else:
                        vec_seq_repair.append(aa_new)
            else:
                vec_seq_repair.append(aa_new)

    df["aa_new_seq_freq"] = vec_aa_new_seq_freq
    df["repaired_aa"] = vec_seq_repair
    df["mutation"] = np.nan

    # --------------------------------------------------
    # Record mutations
    for i in range(len(df)):
        orig = df.at[i, "aa_newSequence"]
        rep = df.at[i, "repaired_aa"]
        if orig != rep:
            no_hxb2 = df.at[i, "no_hxb2"]
            ins = df.at[i, "insert_hxb2"]
            if ins == " ":
                df.at[i, "mutation"] = f"{orig}{int(no_hxb2)}{rep}"
            else:
                df.at[i, "mutation"] = f"{orig}{int(no_hxb2)}{ins}{rep}"

    mutation_list = [v for v in df["mutation"].dropna().tolist()]

    # --------------------------------------------------
    # Save final CSV file
    df.to_csv(f"{output_prefix}_FullAnalysis.csv", index=False, quoting=csv.QUOTE_NONE)

    # Save mutation list
    with open(f"{output_prefix}_repaired_mutList.txt", "w") as fh:
        fh.write("\n".join(mutation_list) + "\n")

    # --------------------------------------------------
    # Save repaired FASTA (no gaps)
    seq_repair_no_gaps = [aa for aa in vec_seq_repair if aa != "-"]
    write_fasta_single(output_prefix, seq_repair_no_gaps, f"{output_prefix}_repaired.fasta")

    # ==================================================
    # Additional mutations  (DS-SOSIP-RnS with Nt/Ct tags)
    print("Additional mutations")

    sub_cols = list(df.columns[:3]) + ["repaired_aa", "mutation"]
    df_repaired = df[sub_cols].copy().reset_index(drop=True)

    vec_final = []
    df_tmp_rows = []
    index = 4  # starts at 4 (after Nt tag "AENL" would be 1-4)

    SIX_R_POSITIONS = {508: "R", 509: "R", 510: "R"}
    DS_POSITIONS    = {201: "C", 433: "C"}
    SOSIP_POSITIONS = {501: "C", 605: "C"}
    STABILIZATION   = {
        559: "P",
        535: "N",
        556: "P",
        588: "E",
        589: "V",
        651: "F",
        655: "I",
        658: "V",
    }

    for i, row in df_repaired.iterrows():
        no_hxb2 = row["no_hxb2"]
        ins = row["insert_hxb2"]

        if no_hxb2 in range(1, 35):
            continue
        if no_hxb2 == 665:
            break

        # 6R at 508–511
        if no_hxb2 in SIX_R_POSITIONS and ins == " ":
            index += 1
            vec_final.append(SIX_R_POSITIONS[no_hxb2])
            df_tmp_rows.append([index, int(no_hxb2), "6R"])
        elif no_hxb2 == 511 and ins == " ":
            for k, pos_label in enumerate(["511", "511a", "511b"]):
                index += 1
                vec_final.append("R")
                df_tmp_rows.append([index, pos_label, "6R"])
        elif no_hxb2 in DS_POSITIONS and ins == " ":
            index += 1
            vec_final.append("C")
            df_tmp_rows.append([index, int(no_hxb2), "DS"])
        elif no_hxb2 in SOSIP_POSITIONS and ins == " ":
            index += 1
            vec_final.append("C")
            df_tmp_rows.append([index, int(no_hxb2), "SOSIP"])
        elif no_hxb2 in STABILIZATION and ins == " ":
            index += 1
            vec_final.append(STABILIZATION[no_hxb2])
            df_tmp_rows.append([index, int(no_hxb2), "Stabilization"])
        elif row["repaired_aa"] != "-":
            index += 1
            vec_final.append(row["repaired_aa"])
            mut = row["mutation"]
            if not pd.isna(mut):
                if ins == " ":
                    df_tmp_rows.append([index, int(no_hxb2), f"Repair_{mut}"])
                else:
                    df_tmp_rows.append([index, f"{int(no_hxb2)}_{ins}", f"Repair_{mut}"])

    # Add N- and C-terminal tags
    nt_tag = list("MPMGSLQPLATLYLLGMLVASVLAQQAENL")
    ct_tag = list("GSAPTKAKRRVVQREKR")
    final_with_tags = nt_tag + vec_final + ct_tag
    final_with_tags_no_gaps = [aa for aa in final_with_tags if aa != "-"]

    write_fasta_single(output_prefix, final_with_tags_no_gaps,
                       f"{output_prefix}_DS-SOSIP-RnS_NtCt-tags.fasta")

    # Save mutations CSV
    df_mutations = pd.DataFrame(df_tmp_rows, columns=["Position", "HXB2", "Mutation"])
    df_mutations.to_csv("Mutations.csv", index=False, quoting=csv.QUOTE_NONE)

    # ==================================================
    # DS-SOSIP-RnS (no Janssen / with AENL prefix only)
    print("DS-SOSIP-RnS (no Janssen)")

    sub_cols2 = list(df.columns[:3]) + ["repaired_aa"]
    df_repaired2 = df[sub_cols2].copy().reset_index(drop=True)

    vec_final2 = []

    for i, row in df_repaired2.iterrows():
        no_hxb2 = row["no_hxb2"]
        ins = row["insert_hxb2"]

        if no_hxb2 in range(1, 35):
            continue
        if no_hxb2 == 665:
            break

        if no_hxb2 in SIX_R_POSITIONS and ins == " ":
            vec_final2.append(SIX_R_POSITIONS[no_hxb2])
        elif no_hxb2 == 511 and ins == " ":
            vec_final2.append("RRR")
        elif no_hxb2 in DS_POSITIONS and ins == " ":
            vec_final2.append("C")
        elif no_hxb2 in SOSIP_POSITIONS and ins == " ":
            vec_final2.append("C")
        elif no_hxb2 == 559 and ins == " ":
            vec_final2.append("P")
        elif no_hxb2 in STABILIZATION and ins == " ":
            vec_final2.append(STABILIZATION[no_hxb2])
        else:
            vec_final2.append(row["repaired_aa"])

    nt_tag2 = list("AENL")
    final2 = nt_tag2 + vec_final2
    final2_no_gaps = [aa for aa in final2 if aa != "-"]
    write_fasta_single(output_prefix, final2_no_gaps, f"{output_prefix}_DS-SOSIP-RnS.fasta")

    # ==================================================
    # DS-SOSIP only (no 6R, no Janssen, with Nt/Ct tags)
    print("DS-SOSIP only")

    sub_cols3 = list(df.columns[:3]) + ["aa_newSequence"]
    df_repaired3 = df[sub_cols3].copy().reset_index(drop=True)

    vec_final3 = []

    for i, row in df_repaired3.iterrows():
        no_hxb2 = row["no_hxb2"]
        ins = row["insert_hxb2"]

        if no_hxb2 in range(1, 35):
            continue
        if no_hxb2 == 665:
            break

        if no_hxb2 in SIX_R_POSITIONS and ins == " ":
            vec_final3.append(SIX_R_POSITIONS[no_hxb2])
        elif no_hxb2 == 511 and ins == " ":
            vec_final3.append("RRR")
        elif no_hxb2 in {201, 433} and ins == " ":
            vec_final3.append("C")
        elif no_hxb2 in {501, 605} and ins == " ":
            vec_final3.append("C")
        elif no_hxb2 == 559 and ins == " ":
            vec_final3.append("P")
        else:
            vec_final3.append(row["aa_newSequence"])

    nt_tag3 = list("MPMGSLQPLATLYLLGMLVASVLAQQAENL")
    ct_tag3 = list("GSAPTKAKRRVVQREKR")
    final3 = nt_tag3 + vec_final3 + ct_tag3
    final3_no_gaps = [aa for aa in final3 if aa != "-"]
    write_fasta_single(output_prefix, final3_no_gaps, f"{output_prefix}_DS-SOSIP.fasta")


# ==================================================
# Main
# ==================================================

if __name__ == "__main__":
    start = time.time()

    if len(sys.argv) < 4:
        print("Usage: python3 ADROITrimer.py <repair_sequence.fasta> <subtype> <muscle_exe>")
        sys.exit(1)

    file_repair_sequence = sys.argv[1]
    subtype = sys.argv[2]
    exe_muscle = sys.argv[3]

    # Resolve reference alignment relative to the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ref_name = "HIV1_FLT_2016_env_PRO_cleanLastCol_final.fasta"
    candidates = [
        os.path.join(script_dir, "required_files", ref_name),
        os.path.join(script_dir, ref_name),
        ref_name,
    ]
    file_reference_alignment = next(
        (p for p in candidates if os.path.exists(p)), ref_name
    )

    adroi_trimer(file_reference_alignment, file_repair_sequence, subtype, exe_muscle)

    duration = (time.time() - start) / 60
    print(f"Script duration: {round(duration, 2)} min")
