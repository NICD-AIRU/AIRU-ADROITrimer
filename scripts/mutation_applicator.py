#!/usr/bin/env python3
"""
mutation_applicator.py
======================
Batch-processes one or more ADROITrimer *_DS-SOSIP-RnS_NtCt-tags.fasta files,
applies user-defined named mutation sets (e.g. 3mut, 2G), and for every
output writes:
  • A new FASTA with the mutations applied
  • A Mutations.csv that exactly matches the format produced by the R script
    (Position, HXB2, Mutation columns — Position values mirror R's `index`
     variable starting at 4, so they are identical to the original R output)

For each FASTA the script auto-locates the matching *_FullAnalysis.csv in
the same directory.  Both files are required; missing pairs are skipped with
a warning.

Usage
-----
  # All FASTA files in a directory, all mutation-set combinations:
  python mutation_applicator.py /path/to/fastas/ --outdir results/

  # Explicit file list (mix files and directories freely):
  python mutation_applicator.py sample1_DS-SOSIP-RnS_NtCt-tags.fasta \\
                                sample2_DS-SOSIP-RnS_NtCt-tags.fasta \\
                                --outdir results/

  # Shell glob (let the shell expand it):
  python mutation_applicator.py *_DS-SOSIP-RnS_NtCt-tags.fasta --outdir results/

  # Explicitly choose the sets to combine (default = all sets in MUTATION_SETS):
  python mutation_applicator.py /path/to/fastas/ --sets 3mut 2G --outdir results/

Only the FULLY COMBINED variant is written (all requested sets applied
together).  No base, 3mut-only or 2G-only files are produced.

Output files per sample (example prefix "H703"):
  H703_3mut_2G_DS-SOSIP-RnS_NtCt.fasta + H703_3mut_2G_Mutation.csv

Dependencies: standard library only — no pip installs needed.
"""

import argparse
import csv
import glob as _glob
import sys
from pathlib import Path


# =============================================================================
# USER-CONFIGURABLE MUTATION SETS
# Keys   = set name (used in --sets and in output filenames)
# Values = {HXB2_position (int): target_amino_acid (str)}
#
# To add a new set, simply insert another entry here — no other changes needed.
# =============================================================================
MUTATION_SETS = {
    "3mut": {302: "M", 320: "L", 329: "P"},
    "2G":   {569: "G", 636: "G"},
}


# =============================================================================
# TERMINAL TAGS  (mirror the R script for _DS-SOSIP-RnS_NtCt-tags.fasta)
# =============================================================================
NT_TAG = "MPMGSLQPLATLYLLGMLVASVLAQQAENL"   # 30 aa N-terminal leader
CT_TAG = "GSAPTKAKRRVVQREKR"                 # 17 aa C-terminal tag


# =============================================================================
# HARD-WIRED STRUCTURAL MUTATIONS
# Mirrors the R "Additional mutations" block exactly.
# Format: {hxb2_pos: (forced_aa_string, Mutation_csv_label)}
# The "RRR" at position 511 is handled specially: it adds THREE residues
# and produces three rows in the CSV (511, 511a, 511b).
# =============================================================================
STRUCTURAL = {
    201: ("C",   "DS"),
    433: ("C",   "DS"),
    501: ("C",   "SOSIP"),
    559: ("P",   "SOSIP"),
    605: ("C",   "SOSIP"),
    508: ("R",   "6R"),
    509: ("R",   "6R"),
    510: ("R",   "6R"),
    511: ("RRR", "6R"),
    535: ("N",   "Stabilization"),
    556: ("P",   "Stabilization"),
    588: ("E",   "Stabilization"),
    589: ("V",   "Stabilization"),
    651: ("F",   "Stabilization"),
    655: ("I",   "Stabilization"),
    658: ("V",   "Stabilization"),
}

SKIP_POSITIONS = set(range(1, 35))   # HXB2 1–34: N-terminal signal peptide
TRUNCATE_AT    = 665                  # stop before this HXB2 position
FASTA_SUFFIX   = "_DS-SOSIP-RnS_NtCt-tags.fasta"
CSV_SUFFIX     = "_FullAnalysis.csv"


# =============================================================================
# FASTA I/O
# =============================================================================

def write_fasta(path: Path, seq_id: str, sequence: str, wrap: int = 60):
    """Write a single-record FASTA with lines wrapped at `wrap` characters."""
    with open(path, "w") as fh:
        fh.write(f">{seq_id}\n")
        for i in range(0, len(sequence), wrap):
            fh.write(sequence[i : i + wrap] + "\n")
    print(f"      FASTA  → {path.name}  ({len(sequence)} aa)")


# =============================================================================
# FullAnalysis CSV reader
# =============================================================================

def read_full_analysis(csv_path: Path) -> list:
    """Parse *_FullAnalysis.csv; return a list of row dicts."""
    with open(csv_path, newline="") as fh:
        return list(csv.DictReader(fh))


# =============================================================================
# SEQUENCE BUILDER
# Mirrors R's "Additional mutations" loop for _DS-SOSIP-RnS_NtCt-tags.fasta.
#
# R initialises  `index <- 4`  because the intermediate "_DS-SOSIP-RnS.fasta"
# prepends "AENL" (4 aa), so the first HIV residue (HXB2 35) lands at
# position 5.  We replicate this so the Position column in the CSV is
# byte-for-byte identical to what the R script produces.
# =============================================================================

def build_base(rows: list) -> tuple:
    """
    Walk the FullAnalysis rows and produce:

      hiv_seq       : str   — HIV-only sequence (no Nt/Ct tags, gaps removed)
      hxb2_to_seqidx: dict  — {hxb2_pos_int → 0-based index in hiv_seq}
                              Only non-insert, non-structural positions are
                              stored here; used to apply user mutations later.
      hxb2_to_pos   : dict  — {hxb2_pos_int → R index value (= Position in CSV)}
                              Same keys as hxb2_to_seqidx.
      base_log      : list  — [(Position, HXB2_label, Mutation_label)] ready
                              for writing to Mutations.csv.
    """
    hiv_chars      = []     # grows one character at a time
    hxb2_to_seqidx = {}     # {int_pos → 0-based index in hiv_chars}
    hxb2_to_pos    = {}     # {int_pos → R index value = Position in CSV}
    base_log       = []     # rows for Mutations.csv

    index = 4               # ← mirrors  `index <- 4`  in the R script

    for row in rows:
        # no_hxb2 is always an integer in FullAnalysis.csv;
        # insert_hxb2 carries the insert letter (" " = regular position)
        try:
            pos = int(str(row["no_hxb2"]).strip())
        except ValueError:
            continue   # safety net only

        insert       = row["insert_hxb2"]                   # " " or letter
        repaired     = str(row["repaired_aa"]).strip()
        mutation_col = str(row.get("mutation", "NA")).strip()  # "W85V" or "NA"

        # ── 1. Skip N-terminal signal peptide ────────────────────────────────
        if pos in SKIP_POSITIONS:
            continue

        # ── 2. Truncate after gp41 ───────────────────────────────────────────
        if pos >= TRUNCATE_AT:
            break

        # ── 3. Structural / hard-wired mutations (non-insert positions only) ─
        if pos in STRUCTURAL and insert == " ":
            aa_str, label = STRUCTURAL[pos]
            if aa_str == "RRR":
                # Special case: 6R at position 511 adds THREE R residues and
                # produces three rows in the CSV: 511, 511a, 511b
                for k, ch in enumerate(aa_str):
                    index += 1
                    hiv_chars.append(ch)
                    hxb2_label = str(pos) if k == 0 else f"{pos}{'ab'[k - 1]}"
                    base_log.append((index, hxb2_label, label))
            else:
                index += 1
                hxb2_to_seqidx[pos] = len(hiv_chars)
                hxb2_to_pos[pos]    = index
                hiv_chars.append(aa_str)
                base_log.append((index, str(pos), label))
            continue

        # ── 4. Skip alignment gaps ────────────────────────────────────────────
        if repaired == "-":
            continue

        # ── 5. Normal residue (includes insert positions relative to HXB2) ───
        index += 1
        if insert == " ":
            # Only non-insert positions are tracked for user point-mutations
            hxb2_to_seqidx[pos] = len(hiv_chars)
            hxb2_to_pos[pos]    = index
        hiv_chars.append(repaired)

        # Record repair mutations (non-NA mutation column = residue was changed)
        if mutation_col not in ("NA", "", "na"):
            if insert == " ":
                hxb2_label = str(pos)
            else:
                hxb2_label = f"{pos}_{insert}"
            base_log.append((index, hxb2_label, f"Repair_{mutation_col}"))

    return "".join(hiv_chars), hxb2_to_seqidx, hxb2_to_pos, base_log


# =============================================================================
# POINT-MUTATION APPLICATOR
# =============================================================================

def apply_mutations(
    hiv_seq:       str,
    hxb2_to_seqidx: dict,
    hxb2_to_pos:   dict,
    mutations:     dict,
    base_log:      list,
    set_label:     str,
) -> tuple:
    """
    Overlay user-defined point mutations on hiv_seq.

    Returns (mutated_hiv_seq, updated_log) where updated_log = base_log plus
    new entries for each applied mutation, sorted by Position.

    The Mutation label for user mutations follows the same notation as repair
    mutations:  Mut_{original_aa}{hxb2_pos}{new_aa}_{set_label}
    e.g.  Mut_T302M_3mut
    """
    seq_list  = list(hiv_seq)
    extra_log = []

    for pos, new_aa in sorted(mutations.items()):
        if pos not in hxb2_to_seqidx:
            print(f"        [WARN] HXB2 {pos} not in sequence "
                  f"(gap, insert, or out of range) — skipped.", file=sys.stderr)
            continue

        idx    = hxb2_to_seqidx[pos]
        old_aa = seq_list[idx]

        if old_aa == new_aa:
            print(f"        [INFO] HXB2 {pos}: already '{new_aa}', no change.")
            continue

        seq_list[idx] = new_aa
        position      = hxb2_to_pos[pos]   # exact R index value
        extra_log.append((position, str(pos), f"Mut_{old_aa}{pos}{new_aa}_{set_label}"))
        print(f"        HXB2 {pos}: {old_aa} → {new_aa}  "
              f"(seq idx {idx},  Position {position})")

    combined = sorted(base_log + extra_log, key=lambda x: x[0])
    return "".join(seq_list), combined


# =============================================================================
# Mutations.csv writer  (replicates R's write.csv with quote=FALSE)
# =============================================================================

def write_mutations_csv(path: Path, log: list):
    """
    Write Mutations.csv with columns  Position,HXB2,Mutation
    in the exact same plain-text format as R's write.csv(quote=FALSE).
    """
    with open(path, "w", newline="\n") as fh:
        fh.write("Position,HXB2,Mutation\n")
        for position, hxb2, mutation in log:
            fh.write(f"{position},{hxb2},{mutation}\n")
    print(f"      CSV    → {path.name}  ({len(log)} entries)")


# =============================================================================
# FILE COLLECTION
# =============================================================================

def collect_fastas(inputs: list) -> list:
    """
    Resolve a mixed list of file paths, directories and shell globs into a
    deduplicated, sorted list of *_DS-SOSIP-RnS_NtCt-tags.fasta Path objects.
    """
    collected = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            collected.extend(sorted(p.glob(f"*{FASTA_SUFFIX}")))
        else:
            matched = sorted(_glob.glob(inp))
            if matched:
                collected.extend(Path(m) for m in matched)
            elif p.exists():
                collected.append(p)
            else:
                print(f"  [WARN] No files matched: {inp}", file=sys.stderr)

    # Deduplicate while preserving order
    seen, unique = set(), []
    for p in collected:
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            unique.append(p)
    return unique


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch-process ADROITrimer *_DS-SOSIP-RnS_NtCt-tags.fasta files "
            "with named mutation sets. Writes one FASTA + one Mutations.csv "
            "per combination, all in a single output directory."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "inputs", nargs="+", metavar="FASTA_OR_DIR",
        help=(
            f"One or more *{FASTA_SUFFIX} files, directories containing them, "
            "or shell globs. Each FASTA must have a co-located *_FullAnalysis.csv."
        ),
    )
    parser.add_argument(
        "--sets", nargs="*", default=list(MUTATION_SETS.keys()), metavar="SET",
        help=(
            f"Named mutation sets to apply. Available: {list(MUTATION_SETS.keys())}. "
            "All requested sets are applied together in ONE combined output "
            "per sample (at least 2 sets required). Default: all."
        ),
    )
    parser.add_argument(
        "--outdir", default="mutation_outputs",
        help=(
            "Single directory for ALL output FASTAs and CSVs "
            "(created if absent). Default: ./mutation_outputs/"
        ),
    )
    args = parser.parse_args()

    # ── Validate requested sets ───────────────────────────────────────────────
    for name in (args.sets or []):
        if name not in MUTATION_SETS:
            sys.exit(
                f"ERROR: Unknown mutation set '{name}'. "
                f"Available: {list(MUTATION_SETS.keys())}"
            )

    if len(set(args.sets or [])) < 2:
        sys.exit("ERROR: Provide at least two mutation sets via --sets "
                 "(only combined-mutation outputs are produced).")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_files = collect_fastas(args.inputs)
    if not fasta_files:
        sys.exit("ERROR: No matching FASTA files found.")

    # Pre-build the combination list (same for every sample)
    requested = args.sets or []
    combos = [tuple(requested)]   # single combination: every set together

    print(f"\n{'='*62}")
    print(f"  Samples found : {len(fasta_files)}")
    print(f"  Mutation sets : {requested}")
    print(f"  Combinations  : {len(combos)}  (combined only)")
    print(f"  Output dir    : {outdir.resolve()}")
    print(f"{'='*62}")

    ok = skipped = 0

    for fasta_path in fasta_files:
        print(f"\n{'─'*62}")
        print(f"  {fasta_path.name}")

        # Derive prefix and locate the matching FullAnalysis CSV
        if not fasta_path.name.endswith(FASTA_SUFFIX):
            print(f"    [SKIP] Does not end with '{FASTA_SUFFIX}'.", file=sys.stderr)
            skipped += 1
            continue

        prefix   = fasta_path.name[: -len(FASTA_SUFFIX)]
        csv_path = fasta_path.parent / f"{prefix}{CSV_SUFFIX}"

        if not csv_path.exists():
            print(
                f"    [SKIP] Matching FullAnalysis CSV not found:\n"
                f"           {csv_path}",
                file=sys.stderr,
            )
            skipped += 1
            continue

        # ── Read FullAnalysis CSV and build base sequence ─────────────────────
        print(f"    Reading  : {csv_path.name}")
        rows = read_full_analysis(csv_path)

        print(f"    Building base sequence ...")
        hiv_seq, hxb2_to_seqidx, hxb2_to_pos, base_log = build_base(rows)
        base_full = NT_TAG + hiv_seq + CT_TAG
        print(
            f"    Sequence : {len(NT_TAG)} (Nt) + {len(hiv_seq)} (HIV) "
            f"+ {len(CT_TAG)} (Ct) = {len(base_full)} aa  |  "
            f"{len(base_log)} mutations logged"
        )

        # ── Per-combination outputs ───────────────────────────────────────────
        for combo in combos:
            label = "_".join(combo)
            print(f"    ── {label} ──")

            # Merge mutations from all sets in this combination
            merged: dict = {}
            for name in combo:
                for pos, aa in MUTATION_SETS[name].items():
                    if pos in merged and merged[pos] != aa:
                        print(
                            f"      [WARN] HXB2 {pos} conflict — "
                            f"'{merged[pos]}' overridden by '{aa}' (set '{name}').",
                            file=sys.stderr,
                        )
                    merged[pos] = aa

            mutated_hiv, ext_log = apply_mutations(
                hiv_seq, hxb2_to_seqidx, hxb2_to_pos, merged, base_log, label
            )
            mutated_full = NT_TAG + mutated_hiv + CT_TAG

            write_fasta(
                outdir / f"{prefix}_{label}_DS-SOSIP-RnS_NtCt.fasta",
                f"{prefix}_{label}",
                mutated_full,
            )
            write_mutations_csv(
                outdir / f"{prefix}_{label}_Mutation.csv",
                ext_log,
            )

        ok += 1

    print(f"\n{'='*62}")
    print(f"  ✓ Done.  Processed: {ok}   Skipped: {skipped}")
    print(f"  All outputs in: {outdir.resolve()}")
    print(f"{'='*62}\n")


if __name__ == "__main__":
    main()
