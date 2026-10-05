#!/usr/bin/env Rscript
## =============================================================================
## ADROITrimer RnS-DS-SOSIP Mutation QC Dashboard
##
## Reads the output files produced by ADROITrimer1.0-3.R
## (https://github.com/RedaRawi/ADROITrimer; Rawi et al. 2020, Cell Reports,
## https://doi.org/10.1016/j.celrep.2020.108432) for one or more HIV-1 Env
## strains, and confirms which of the fixed RnS-DS-SOSIP structural
## mutations were successfully added by the pipeline, which are missing,
## and lists the strain-specific consensus "repair" substitutions alongside
## them.
##
## Primary input per strain: Mutations.csv
##   (columns: Position, HXB2, Mutation -- written by ADROITrimer's
##   "Additional mutations" block; Mutation is one of "SOSIP", "DS", "6R",
##   "Stabilization", "3mut", "2G", or "Repair_<aa>")
## Multiple Mutations.csv files (one per strain) can be uploaded in a single
## batch on the Upload tab.
##
## The dashboard also ships an offline "AI Assistant" tab that talks to a local
## Ollama server (default model: dolphin3:8b). See the OLLAMA_* settings below.
##
## NOTE: ADROITrimer1.0-3.R writes Mutations.csv under a fixed filename
## (not prefixed with output.prefix), so it gets overwritten on every run.
## Rename/copy it per strain (e.g. CAP256_Mutations.csv) before uploading
## here for multi-strain / batch comparison -- the dashboard uses the
## uploaded file name as the strain label.
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

## -----------------------------------------------------------------------------
## Offline AI assistant settings (Ollama running on this machine)
##   Override without editing this file via environment variables, e.g.
##   ADROIT_OLLAMA_URL=http://127.0.0.1:11434  ADROIT_OLLAMA_MODEL=dolphin3:8b
## The assistant needs: curl, jsonlite, promises, later (all installed with
## shiny in a normal setup). If any is missing the rest of the app still
## works and the assistant tab explains what to install.
## -----------------------------------------------------------------------------
OLLAMA_URL   <- sub("/+$", "", Sys.getenv("ADROIT_OLLAMA_URL", "http://127.0.0.1:11434"))
OLLAMA_MODEL <- Sys.getenv("ADROIT_OLLAMA_MODEL", "dolphin3:8b")
llm_pkgs_ok  <- all(vapply(c("curl", "jsonlite", "promises", "later"),
                            requireNamespace, logical(1), quietly = TRUE))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## -----------------------------------------------------------------------------
## Canonical RnS-DS-SOSIP mutation checklist (HXB2 numbering)
##
## Source: Rawi et al. 2020 Cell Reports, Figure 1 (structure-based
## stabilization set: 535N, 556P, 588E, 589V, 651F, 655I, 658V) plus the
## fixed DS / SOSIP / 6R substitutions hardcoded in the "Additional
## mutations" block of ADROITrimer1.0-3.R (positions 201, 433, 501, 605,
## 559, and 508-511/511a/511b), plus the two RnS-3mut-2G add-ons:
##   3mut = N302M / T320L / A329P  (apex stabilization; Ou et al. 2020 JVI)
##   2G   = 569G / 636G            (glycine helix-breaking substitutions)
## These 23 positions are constant across every ADROITrimer run regardless
## of input strain; everything else the pipeline changes (the "Repair_*"
## rows) is consensus-driven and strain-specific, so it is tracked
## separately rather than checked off here.
## NOTE: wild-type residues for 3mut are the BG505 residues from the
## literature; the checklist matches on HXB2 position + category only, so a
## different wild-type residue in your strain does not affect Present/Missing.
## -----------------------------------------------------------------------------
canonical_mutations <- tibble::tribble(
  ~category,       ~hxb2_position, ~wt_aa, ~mut_aa, ~label,                    ~note,
  "DS",             "201",          "I",    "C",    "I201C",                   "Kwong-lab DS disulfide (gp120-gp41)",
  "DS",             "433",          "A",    "C",    "A433C",                   "Kwong-lab DS disulfide (gp120-gp41)",
  "SOSIP",          "501",          "A",    "C",    "A501C",                   "SOS disulfide (gp120-gp41)",
  "SOSIP",          "605",          "T",    "C",    "T605C",                   "SOS disulfide (gp120-gp41)",
  "SOSIP",          "559",          "I",    "P",    "I559P",                   "SOSIP proline (gp41 helix)",
  "6R",             "508",          NA,     "R",    "508R",                    "Hexa-arginine furin cleavage site",
  "6R",             "509",          NA,     "R",    "509R",                    "Hexa-arginine furin cleavage site",
  "6R",             "510",          NA,     "R",    "510R",                    "Hexa-arginine furin cleavage site",
  "6R",             "511",          NA,     "R",    "511R",                    "Hexa-arginine furin cleavage site",
  "6R",             "511a",         NA,     "R",    "511aR",                   "Hexa-arginine furin cleavage site (insertion)",
  "6R",             "511b",         NA,     "R",    "511bR",                   "Hexa-arginine furin cleavage site (insertion)",
  "Stabilization",  "535",          "M",    "N",    "M535N",                   "Structure-based stabilization",
  "Stabilization",  "556",          "L",    "P",    "L556P",                   "Structure-based stabilization",
  "Stabilization",  "588",          "K",    "E",    "K588E",                   "Structure-based stabilization",
  "Stabilization",  "589",          "D",    "V",    "D589V",                   "Structure-based stabilization",
  "Stabilization",  "651",          "N",    "F",    "N651F",                   "Structure-based stabilization",
  "Stabilization",  "655",          "K",    "I",    "K655I",                   "Structure-based stabilization",
  "Stabilization",  "658",          "K",    "V",    "K658V",                   "Structure-based stabilization",
  "3mut",           "302",          "N",    "M",    "N302M",                   "3mut apex stabilization (wt shown = BG505)",
  "3mut",           "320",          "T",    "L",    "T320L",                   "3mut apex stabilization (wt shown = BG505)",
  "3mut",           "329",          "A",    "P",    "A329P",                   "3mut apex stabilization (wt shown = BG505)",
  "2G",             "569",          NA,     "G",    "569G",                    "2G glycine helix-breaking substitution (gp41)",
  "2G",             "636",          NA,     "G",    "636G",                    "2G glycine helix-breaking substitution (gp41)"
) |>
  mutate(
    hxb2_sort = suppressWarnings(as.numeric(str_extract(hxb2_position, "^[0-9]+"))),
    label = if_else(!is.na(wt_aa), paste0(wt_aa, hxb2_position, mut_aa), label)
  )

category_order <- c("DS", "SOSIP", "6R", "Stabilization", "3mut", "2G")
category_colors <- c(
  "DS"            = "#2E86AB",
  "SOSIP"         = "#5B2C6F",
  "6R"            = "#B9770E",
  "Stabilization" = "#1E8449",
  "3mut"          = "#E056A0",
  "2G"            = "#48C9B0",
  "Repair"        = "#7B8894"
)

## e.g. "23 fixed structural positions tracked: 2 DS + 3 SOSIP + 6 6R + ..."
checklist_summary <- local({
  n <- table(factor(canonical_mutations$category, levels = category_order))
  paste0(nrow(canonical_mutations), " fixed structural positions tracked: ",
         paste(n, names(n), collapse = " + "))
})
status_colors <- c("Present" = "#1E8449", "Missing" = "#C0392B")

## -----------------------------------------------------------------------------
## Parsing helpers
## -----------------------------------------------------------------------------

## Parses one Mutations.csv (ADROITrimer1.0-3.R output). Returns a tibble
## with an added `category` (DS/SOSIP/6R/Stabilization/Repair) and, for
## Repair rows, the amino acid the position was repaired to.
parse_mutations_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE,
                         col_types = readr::cols(.default = "c"))
  req <- c("Position", "HXB2", "Mutation")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols)) {
    stop("This file doesn't look like an ADROITrimer Mutations.csv -- missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  df |>
    mutate(
      HXB2 = str_trim(HXB2),
      Mutation = str_trim(Mutation),
      category = if_else(str_starts(Mutation, "Repair"), "Repair",
                          str_split_i(Mutation, "_", 1)),
      ## tolerate case differences for the known labels (e.g. "3Mut", "2g")
      category = if_else(str_to_lower(category) %in% str_to_lower(category_order),
                          category_order[match(str_to_lower(category), str_to_lower(category_order))],
                          category),
      repaired_to_aa = if_else(category == "Repair",
                                str_split_i(Mutation, "_", 2), NA_character_)
    )
}

## Compares the canonical checklist against one strain's detected fixed
## (non-Repair) mutations and returns Present/Missing per canonical row.
build_checklist_status <- function(canonical, detected) {
  detected_fixed <- detected |> filter(category != "Repair")
  canonical |>
    left_join(
      detected_fixed |> select(HXB2, category, Position),
      by = c("hxb2_position" = "HXB2", "category" = "category")
    ) |>
    mutate(status = if_else(is.na(Position), "Missing", "Present")) |>
    select(-Position) |>
    arrange(hxb2_sort)
}


## -----------------------------------------------------------------------------
## Offline AI assistant helpers (Ollama, local only)
## -----------------------------------------------------------------------------

## Names of models installed in the local Ollama server, or NULL if the server
## is unreachable. Short timeout so the UI never hangs on a dead server.
ollama_list_models <- function(host = OLLAMA_URL, timeout = 3) {
  if (!llm_pkgs_ok) return(NULL)
  h <- curl::new_handle(timeout = timeout, connecttimeout = timeout)
  res <- tryCatch(curl::curl_fetch_memory(paste0(host, "/api/tags"), handle = h),
                  error = function(e) NULL)
  if (is.null(res) || res$status_code != 200) return(NULL)
  parsed <- tryCatch(jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed)) return(NULL)
  vapply(parsed$models %||% list(), function(m) m$name %||% "", character(1))
}

model_installed <- function(model, installed) {
  model %in% installed || (!grepl(":", model, fixed = TRUE) && paste0(model, ":latest") %in% installed)
}

## Non-blocking call to Ollama's /api/chat. Returns a promise that resolves to
## the assistant's reply text (character) or rejects with an error. Built on
## curl's multi interface + later, so the Shiny session stays responsive while
## the model generates.
ollama_chat_async <- function(messages, host = OLLAMA_URL, model = OLLAMA_MODEL,
                              temperature = 0.2, num_ctx = 8192, timeout = 900) {
  body <- jsonlite::toJSON(
    list(model = model, messages = messages, stream = FALSE,
         options = list(temperature = temperature, num_ctx = num_ctx)),
    auto_unbox = TRUE
  )
  promises::promise(function(resolve, reject) {
    h <- curl::new_handle()
    curl::handle_setopt(h, post = TRUE, postfields = charToRaw(enc2utf8(as.character(body))),
                        timeout = timeout, connecttimeout = 5)
    curl::handle_setheaders(h, `Content-Type` = "application/json")
    pool <- curl::new_pool()
    curl::curl_fetch_multi(
      paste0(host, "/api/chat"), handle = h, pool = pool,
      done = function(res) {
        tryCatch({
          txt <- rawToChar(res$content); Encoding(txt) <- "UTF-8"
          parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
          if (res$status_code != 200 || !is.null(parsed$error)) {
            stop(parsed$error %||% paste("HTTP", res$status_code))
          }
          resolve(trimws(parsed$message$content %||% ""))
        }, error = function(e) reject(e))
      },
      fail = function(msg) {
        reject(simpleError(paste0("Could not reach Ollama at ", host, " (", msg, ")")))
      }
    )
    poll <- function() {
      st <- curl::multi_run(timeout = 0, pool = pool)
      if (st$pending > 0) later::later(poll, 0.15)
    }
    poll()
  })
}

## Plain-text summary of the loaded data that is handed to the model as
## context. Active strain gets detail; other strains get one line each.
build_llm_context <- function(strains, active, max_other = 12) {
  if (length(strains) == 0) return("No strains are loaded in the dashboard yet.")
  describe <- function(lbl, detailed) {
    d <- strains[[lbl]]$detected
    st <- build_checklist_status(canonical_mutations, d)
    miss <- st[st$status == "Missing", ]
    line <- paste0(lbl, ": ", sum(st$status == "Present"), "/", nrow(st), " checklist positions present; ",
                   if (nrow(miss)) paste0("missing ", paste0(miss$label, " [", miss$category, "]", collapse = ", "))
                   else "none missing",
                   "; ", sum(d$category == "Repair"), " consensus repair substitutions.")
    if (!detailed) return(line)
    by_cat <- st |> group_by(category) |>
      summarise(txt = paste0(category[1], " ", sum(status == "Present"), "/", n()), .groups = "drop")
    by_cat <- by_cat[match(category_order, by_cat$category), ]
    rep_rows <- d[d$category == "Repair", ]
    rep_txt <- if (nrow(rep_rows)) {
      shown <- utils::head(rep_rows, 60)
      paste0(paste0(shown$HXB2, "->", shown$repaired_to_aa, collapse = ", "),
             if (nrow(rep_rows) > 60) paste0(" (+", nrow(rep_rows) - 60, " more)") else "")
    } else "none"
    paste0(line, "\n  By category: ", paste(stats::na.omit(by_cat$txt), collapse = "; "),
           "\n  Repair substitutions (HXB2 position -> residue): ", rep_txt)
  }
  act <- if (!is.null(active) && active %in% names(strains)) active else names(strains)[1]
  others <- setdiff(names(strains), act)
  paste0(
    "ACTIVE STRAIN\n", describe(act, TRUE),
    if (length(others)) paste0(
      "\n\nOTHER LOADED STRAINS\n",
      paste(vapply(utils::head(others, max_other), describe, character(1), detailed = FALSE), collapse = "\n"),
      if (length(others) > max_other) paste0("\n(+", length(others) - max_other, " more strains not shown)") else ""
    ) else ""
  )
}

build_llm_system_prompt <- function(context = NULL) {
  ck <- canonical_mutations |>
    group_by(category) |>
    summarise(txt = paste0(category[1], ": ", paste(label, collapse = ", ")), .groups = "drop")
  ck <- ck[match(category_order, ck$category), ]
  paste0(
    "You are the QC assistant built into the ADROITrimer RnS-DS-SOSIP dashboard, used by a ",
    "bioinformatics/immunology lab. You run fully offline on the user's own machine.\n\n",
    "BACKGROUND\n",
    "ADROITrimer (Rawi et al. 2020, Cell Reports) designs stabilized prefusion-closed HIV-1 Env trimers. ",
    "It writes Mutations.csv listing (a) a fixed set of structural mutations added to every strain and ",
    "(b) strain-specific consensus 'Repair' substitutions. Positions use HXB2 numbering.\n",
    "The dashboard checks these ", nrow(canonical_mutations), " fixed positions:\n",
    paste0("- ", ck$txt, collapse = "\n"), "\n",
    "'6R' is the hexa-arginine furin cleavage site (508-511 plus insertions 511a/511b). ",
    "'3mut' (N302M/T320L/A329P) stabilizes the trimer apex. '2G' (569G/636G) are glycine helix-breaking ",
    "substitutions in gp41. 'DS' is the I201C-A433C disulfide; 'SOSIP' is A501C-T605C plus I559P.\n",
    "'Present' means the position appears in that strain's Mutations.csv; 'Missing' means it does not.\n\n",
    "RULES\n",
    "- Statements about specific strains must come only from the DATA section. If something is not in the DATA, ",
    "say you do not have it. Never invent positions, residues, counts or literature results.\n",
    "- The dashboard only reads Mutations.csv; it does not re-read the FASTA sequence. A Missing position usually means ",
    "the pipeline could not place the mutation (for example an indel or alignment problem near that HXB2 position), ",
    "so suggest checking the alignment/input sequence and re-running ADROITrimer, and say this is a likely explanation, not a certainty.\n",
    "- You are a helper, not a substitute for lab verification. Keep answers concise (short paragraphs or bullets) and use HXB2 numbers.\n\n",
    "DATA\n", context %||% "(The user chose not to share the loaded strain data with you.)"
  )
}

## -----------------------------------------------------------------------------
## UI
## -----------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "ADROITrimer RnS-DS-SOSIP QC", titleWidth = 320),

  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "sidebar",
      menuItem("Upload", tabName = "upload", icon = icon("upload")),
      menuItem("RnS Checklist", tabName = "checklist", icon = icon("list-check")),
      menuItem("Repair Mutations", tabName = "repairs", icon = icon("wrench")),
      menuItem("Sequence Map", tabName = "seqmap", icon = icon("dna")),
      menuItem("Batch Compare", tabName = "batch", icon = icon("layer-group")),
      menuItem("AI Assistant", tabName = "assistant", icon = icon("robot")),
      menuItem("About", tabName = "about", icon = icon("circle-info"))
    ),
    br(),
    div(
      style = "padding: 0 15px; color: #b8c7ce; font-size: 12px;",
      "Reference: Rawi et al. 2020, Cell Reports (RnS-DS-SOSIP / ADROITrimer). ",
      paste0(checklist_summary, ".")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #0f1b24; }
      .box { border-top-color: #1F6FB2; }
      .small-box h3, .small-box p { color: #ffffff; }
      .status-present { color: #1E8449; font-weight: 600; }
      .status-missing { color: #C0392B; font-weight: 700; }
      .strain-pill {
        display: inline-block; padding: 2px 10px; border-radius: 12px;
        background: #1F6FB2; color: #fff; font-size: 12px; margin-right: 6px;
      }
      /* AI assistant */
      #chat_log_wrap { height: 470px; overflow-y: auto; padding: 10px; background: #0f1b24;
                       border: 1px solid #1c2b36; border-radius: 4px; }
      .chat-msg { margin-bottom: 12px; max-width: 92%; }
      .chat-msg.user { margin-left: auto; }
      .chat-role { font-size: 11px; color: #9fb0bc; margin-bottom: 2px; }
      .chat-msg.user .chat-role { text-align: right; }
      .chat-text { padding: 8px 12px; border-radius: 10px; white-space: pre-wrap;
                   line-height: 1.45; color: #e6edf3; }
      .chat-msg.user .chat-text { background: #1F6FB2; }
      .chat-msg.assistant .chat-text { background: #1c2b36; }
      .chat-msg.error .chat-text { background: #3a1414; color: #F08080; }
      .chat-thinking { font-style: italic; color: #9fb0bc !important; }
      .chat-empty { color: #9fb0bc; padding: 20px; }
      .chat-chip { margin: 0 6px 6px 0; }
      #chat_input { background: #0f1b24; color: #e6edf3; border-color: #1c2b36; resize: none; }
      .llm-dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; }
      .llm-ok { background: #1E8449; } .llm-warn { background: #D68910; } .llm-bad { background: #C0392B; }
    ")),
    tags$script(HTML("
      $(function() {
        function submitChat(text) {
          text = (text || '').trim();
          if (!text) return;
          Shiny.setInputValue('chat_submit', {text: text, nonce: Math.random()}, {priority: 'event'});
        }
        $(document).on('click', '#chat_send', function() {
          var el = $('#chat_input'); submitChat(el.val()); el.val('');
        });
        $(document).on('keydown', '#chat_input', function(e) {
          if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); $('#chat_send').click(); }
        });
        $(document).on('click', '.chat-chip', function() { submitChat($(this).attr('data-prompt')); });
      });
      $(document).on('shiny:value', function(e) {
        if (e.name === 'chat_log') {
          setTimeout(function() {
            var w = document.getElementById('chat_log_wrap');
            if (w) w.scrollTop = w.scrollHeight;
          }, 60);
        }
      });
    "))),

    tabItems(

      ## ---------------- Upload ----------------
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "Load Mutations.csv", status = "primary", solidHeader = TRUE,
            width = 6,
            p("Upload one or more ", code("Mutations.csv"), " files written by ",
              code("ADROITrimer1.0-3.R"), " (one per strain/run). Select multiple ",
              "files at once to batch-load them. Rename each file per strain before ",
              "uploading if you're comparing multiple runs (the pipeline writes the ",
              "same filename every time) -- the dashboard uses each file's name as ",
              "its strain label."),
            fileInput("mutations_csv", "Mutations.csv (one or more files)",
                      accept = ".csv", multiple = TRUE),
            textInput("strain_label", "Strain label (optional, single-file uploads only)",
                      placeholder = "e.g. CAP256.SU"),
            actionButton("load_active", "Load strain(s)", icon = icon("check"), class = "btn-primary")
          ),
          box(
            title = "Active strain", status = "info", solidHeader = TRUE, width = 6,
            uiOutput("active_strain_summary")
          )
        ),
        fluidRow(
          box(
            title = "Loaded strains (for Batch Compare)", status = "primary", width = 12,
            DTOutput("loaded_strains_table"),
            br(),
            actionButton("clear_all", "Clear all loaded strains", icon = icon("trash"), class = "btn-danger")
          )
        )
      ),

      ## ---------------- RnS Checklist ----------------
      tabItem(
        tabName = "checklist",
        fluidRow(
          valueBoxOutput("vb_total"),
          valueBoxOutput("vb_present"),
          valueBoxOutput("vb_missing")
        ),
        fluidRow(
          box(
            title = "Fixed RnS-DS-SOSIP mutation checklist", status = "primary",
            solidHeader = TRUE, width = 12,
            p("Every canonical position ADROITrimer is designed to always introduce ",
              "(DS / SOSIP / 6R / Stabilization / 3mut / 2G), and whether it was found in this strain's ",
              "Mutations.csv."),
            DTOutput("checklist_table")
          )
        )
      ),

      ## ---------------- Repair Mutations ----------------
      tabItem(
        tabName = "repairs",
        fluidRow(
          valueBoxOutput("vb_nrepairs", width = 4)
        ),
        fluidRow(
          box(
            title = "Consensus repair substitutions (strain-specific)", status = "warning",
            solidHeader = TRUE, width = 12,
            p("These are not part of the fixed checklist -- they're rare-residue repairs ",
              "ADROITrimer made against the subtype consensus for this particular input ",
              "sequence, so the count and positions vary strain to strain."),
            DTOutput("repairs_table")
          )
        )
      ),

      ## ---------------- Sequence Map ----------------
      tabItem(
        tabName = "seqmap",
        fluidRow(
          box(
            title = "Mutation map along HXB2 Env numbering", status = "primary",
            solidHeader = TRUE, width = 12,
            p("Env drawn as a gp120 / gp41 gene track. Detected mutations are shown as ",
              "lollipops above the backbone, positioned by HXB2 coordinate and colored by ",
              "category. Missing checklist positions are shown as open red circles below ",
              "the backbone at their expected location."),
            plotOutput("seqmap_plot", height = 440)
          )
        )
      ),

      ## ---------------- Batch Compare ----------------
      tabItem(
        tabName = "batch",
        fluidRow(
          box(
            title = "Checklist status across all loaded strains", status = "primary",
            solidHeader = TRUE, width = 12,
            p("Upload additional strains from the Upload tab (each with a distinct filename) ",
              "to build this comparison. Green = present, red = missing."),
            plotOutput("batch_heatmap", height = 420)
          )
        ),
        fluidRow(
          box(
            title = "Batch summary table", status = "primary", width = 12,
            DTOutput("batch_table")
          )
        )
      ),

      ## ---------------- AI Assistant (offline, Ollama) ----------------
      tabItem(
        tabName = "assistant",
        fluidRow(
          column(
            width = 8,
            box(
              title = "Ask about your strains", status = "primary", solidHeader = TRUE, width = 12,
              div(id = "chat_log_wrap", uiOutput("chat_log")),
              br(),
              tags$textarea(id = "chat_input", class = "form-control", rows = 2,
                            placeholder = "Ask a question... (Enter to send, Shift+Enter for a new line)"),
              br(),
              tags$button(id = "chat_send", type = "button", class = "btn btn-primary",
                          icon("paper-plane"), " Send"),
              actionButton("chat_clear", "Clear conversation", icon = icon("eraser"))
            )
          ),
          column(
            width = 4,
            box(
              title = "Local model", status = "info", solidHeader = TRUE, width = 12,
              uiOutput("llm_status"),
              actionButton("chat_check", "Re-check connection", icon = icon("rotate")),
              tags$hr(),
              checkboxInput("chat_use_context", "Share loaded-strain data with the assistant", value = TRUE),
              sliderInput("chat_temp", "Temperature (lower = more literal)", min = 0, max = 1,
                          value = 0.2, step = 0.05),
              tags$small("Runs entirely on this machine through Ollama; nothing is sent to the internet. ",
                         "An 8B model can still make mistakes -- treat answers as a QC aid, not a verdict.")
            ),
            box(
              title = "Quick questions", status = "primary", width = 12,
              tags$button(type = "button", class = "btn btn-default btn-sm chat-chip",
                          `data-prompt` = "Summarise the QC result for the active strain.",
                          "Summarise active strain"),
              tags$button(type = "button", class = "btn btn-default btn-sm chat-chip",
                          `data-prompt` = "Which checklist mutations are missing in the active strain, and what could explain that?",
                          "Why is something missing?"),
              tags$button(type = "button", class = "btn btn-default btn-sm chat-chip",
                          `data-prompt` = "Compare all loaded strains and point out any that need attention.",
                          "Compare loaded strains"),
              tags$button(type = "button", class = "btn btn-default btn-sm chat-chip",
                          `data-prompt` = "Explain what the 3mut and 2G mutations are and why they are in the checklist.",
                          "What are 3mut and 2G?")
            )
          )
        )
      ),

      ## ---------------- About ----------------
      tabItem(
        tabName = "about",
        fluidRow(
          box(
            title = "About this dashboard", status = "primary", solidHeader = TRUE, width = 12,
            HTML("
              <p>This dashboard QCs the output of <b>ADROITrimer</b>, an automated
              structure-based stabilization + consensus-repair pipeline for designing
              prefusion-closed HIV-1 Env trimers.</p>
              <ul>
                <li>Paper: Rawi <i>et al.</i> 2020, <i>Cell Reports</i> --
                  <a href='https://doi.org/10.1016/j.celrep.2020.108432' target='_blank'>
                  doi:10.1016/j.celrep.2020.108432</a></li>
                <li>Code: <a href='https://github.com/RedaRawi/ADROITrimer' target='_blank'>
                  github.com/RedaRawi/ADROITrimer</a></li>
              </ul>
              <p>The 23-position checklist tracked here is 2 DS, 3 SOSIP, 6 6R and 7 structure-based
              Stabilization mutations (hardcoded identically in every run of
              <code>ADROITrimer1.0-3.R</code>'s 'Additional mutations' block), plus the RnS-3mut-2G
              add-ons: <b>3mut</b> (N302M, T320L, A329P; apex stabilization) and <b>2G</b>
              (569G, 636G; glycine helix-breaking substitutions). These are independent of
              input strain -- so unlike the consensus 'Repair' substitutions, a missing entry
              here indicates the pipeline could not introduce that mutation for this
              particular sequence (e.g. an indel disrupted the expected position) rather than
              a strain difference. 3mut and 2G rows are matched on the <code>Mutation</code>
              labels <code>3mut</code> and <code>2G</code> in Mutations.csv (case-insensitive).</p>
              <p>The <b>AI Assistant</b> tab uses a local Ollama server (model
              <code>dolphin3:8b</code> by default). Nothing is sent to the internet.</p>
            ")
          )
        )
      )
    )
  )
)

## -----------------------------------------------------------------------------
## Server
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    strains = list(),       # named list: label -> list(detected=tibble)
    active = NULL           # currently viewed strain label
  )

  ## ---- Load one or more strains from the Upload tab ----
  observeEvent(input$load_active, {
    req(input$mutations_csv)
    files <- input$mutations_csv   ## data.frame: one row per uploaded file
    label_input <- str_trim(input$strain_label)
    multi <- nrow(files) > 1

    loaded_labels <- character(0)
    for (i in seq_len(nrow(files))) {
      detected <- tryCatch(
        parse_mutations_csv(files$datapath[i]),
        error = function(e) {
          showNotification(paste0(files$name[i], ": ", conditionMessage(e)),
                            type = "error", duration = 8)
          NULL
        }
      )
      if (is.null(detected)) next

      label <- if (!multi && label_input != "") label_input else tools::file_path_sans_ext(files$name[i])
      ## avoid silent overwrite of a same-named strain already loaded
      orig_label <- label
      n <- 1
      while (label %in% names(rv$strains)) {
        n <- n + 1
        label <- paste0(orig_label, " (", n, ")")
      }

      rv$strains[[label]] <- list(detected = detected)
      loaded_labels <- c(loaded_labels, label)
    }

    req(length(loaded_labels) > 0)
    rv$active <- loaded_labels[length(loaded_labels)]
    updateTextInput(session, "strain_label", value = "")
    if (length(loaded_labels) == 1) {
      showNotification(paste0("Loaded strain '", loaded_labels, "' as active."), type = "message")
    } else {
      showNotification(paste0("Loaded ", length(loaded_labels), " strains; '",
                                rv$active, "' set as active."), type = "message")
    }
  })

  observeEvent(input$clear_all, {
    rv$strains <- list()
    rv$active <- NULL
  })

  ## ---- Loaded strains table (Upload tab) ----
  output$loaded_strains_table <- renderDT({
    if (length(rv$strains) == 0) {
      return(datatable(tibble(Message = "No strains loaded yet."), rownames = FALSE))
    }
    summary_df <- purrr_map_df_strains(rv$strains)
    datatable(
      summary_df, rownames = FALSE, selection = "single",
      options = list(dom = "t", paging = FALSE)
    )
  })

  ## click a row in loaded_strains_table to switch active strain
  observeEvent(input$loaded_strains_table_rows_selected, {
    sel <- input$loaded_strains_table_rows_selected
    req(sel)
    labels <- names(rv$strains)
    if (length(labels) >= sel) rv$active <- labels[sel]
  })

  output$active_strain_summary <- renderUI({
    if (is.null(rv$active)) {
      return(tags$p("No active strain yet -- upload one or more Mutations.csv files and click ",
                     tags$b("Load strain(s)"), "."))
    }
    d <- rv$strains[[rv$active]]$detected
    status <- build_checklist_status(canonical_mutations, d)
    n_present <- sum(status$status == "Present")
    n_total <- nrow(status)
    n_repair <- sum(d$category == "Repair")
    tagList(
      tags$p(tags$span(class = "strain-pill", rv$active)),
      tags$p(tags$b(n_present), " / ", n_total, " fixed RnS-DS-SOSIP mutations present"),
      tags$p(tags$b(n_repair), " consensus repair substitutions recorded"),
      if (n_present < n_total)
        tags$p(class = "status-missing",
               paste0(n_total - n_present, " checklist position(s) missing -- see RnS Checklist tab"))
      else
        tags$p(class = "status-present", "All fixed checklist positions present.")
    )
  })

  ## ---- helper: build a one-row-per-strain summary table ----
  purrr_map_df_strains <- function(strains) {
    rows <- lapply(names(strains), function(lbl) {
      d <- strains[[lbl]]$detected
      status <- build_checklist_status(canonical_mutations, d)
      tibble(
        Strain = lbl,
        `Checklist present` = sum(status$status == "Present"),
        `Checklist missing` = sum(status$status == "Missing"),
        `Repair mutations` = sum(d$category == "Repair"),
        Active = if (identical(lbl, rv$active)) "\u2713" else ""
      )
    })
    bind_rows(rows)
  }

  ## ---- Active strain's checklist status (reactive) ----
  active_status <- reactive({
    req(rv$active)
    build_checklist_status(canonical_mutations, rv$strains[[rv$active]]$detected)
  })

  active_detected <- reactive({
    req(rv$active)
    rv$strains[[rv$active]]$detected
  })

  ## ---- RnS Checklist tab ----
  output$vb_total <- renderValueBox({
    valueBox(nrow(canonical_mutations), "Fixed checklist positions", icon = icon("list"), color = "blue")
  })
  output$vb_present <- renderValueBox({
    st <- tryCatch(active_status(), error = function(e) NULL)
    n <- if (is.null(st)) 0 else sum(st$status == "Present")
    valueBox(n, "Present in active strain", icon = icon("check"), color = "green")
  })
  output$vb_missing <- renderValueBox({
    st <- tryCatch(active_status(), error = function(e) NULL)
    n <- if (is.null(st)) nrow(canonical_mutations) else sum(st$status == "Missing")
    valueBox(n, "Missing in active strain", icon = icon("triangle-exclamation"),
              color = if (n > 0) "red" else "green")
  })

  output$checklist_table <- renderDT({
    st <- tryCatch(active_status(), error = function(e) NULL)
    validate(need(!is.null(st), "Load and select an active strain on the Upload tab first."))
    disp <- st |>
      select(Category = category, `HXB2 position` = hxb2_position, Mutation = label,
             Note = note, Status = status)
    datatable(disp, rownames = FALSE,
              options = list(pageLength = nrow(canonical_mutations), dom = "ftip")) |>
      formatStyle("Status", target = "row",
                  backgroundColor = styleEqual(c("Present", "Missing"),
                                                c("#0f2e1a", "#3a1414")),
                  color = styleEqual(c("Present", "Missing"),
                                      c("#5CD68A", "#F08080")))
  })

  ## ---- Repair Mutations tab ----
  output$vb_nrepairs <- renderValueBox({
    d <- tryCatch(active_detected(), error = function(e) NULL)
    n <- if (is.null(d)) 0 else sum(d$category == "Repair")
    valueBox(n, "Consensus repair substitutions", icon = icon("wrench"), color = "yellow")
  })
  output$repairs_table <- renderDT({
    d <- tryCatch(active_detected(), error = function(e) NULL)
    validate(need(!is.null(d), "Load and select an active strain on the Upload tab first."))
    rep_tbl <- d |>
      filter(category == "Repair") |>
      transmute(Position = Position, `HXB2 position` = HXB2, `Repaired to` = repaired_to_aa) |>
      arrange(suppressWarnings(as.numeric(str_extract(`HXB2 position`, "^[0-9]+"))))
    validate(need(nrow(rep_tbl) > 0, "No consensus repair substitutions recorded for this strain."))
    datatable(rep_tbl, rownames = FALSE, options = list(pageLength = 15, dom = "ftip"))
  })

  ## ---- Sequence Map tab ----
  ## Genome/gene-track style: Env drawn as a gp120/gp41 backbone, detected
  ## mutations shown as lollipops (stem + point) above it colored by
  ## category, missing checklist positions shown as open circles below it.
  output$seqmap_plot <- renderPlot({
    d <- tryCatch(active_detected(), error = function(e) NULL)
    st <- tryCatch(active_status(), error = function(e) NULL)
    validate(need(!is.null(d) && !is.null(st), "Load and select an active strain on the Upload tab first."))

    env_len <- 856   ## HXB2 Env (gp160) length used for the gene-track backbone

    d_plot <- d |>
      mutate(pos_num = suppressWarnings(as.numeric(str_extract(HXB2, "^[0-9]+"))),
             cat_plot = category) |>
      filter(!is.na(pos_num)) |>
      arrange(pos_num) |>
      mutate(
        lane = (row_number() - 1) %% 4,
        lollipop_h = 0.55 + lane * 0.42
      )

    missing_plot <- st |>
      filter(status == "Missing") |>
      arrange(hxb2_sort) |>
      mutate(
        lane = (row_number() - 1) %% 3,
        stem_bottom = -0.45 - lane * 0.35
      )

    domains <- tibble::tibble(
      domain = c("gp120", "gp41"),
      xmin   = c(1, 512),
      xmax   = c(511, env_len),
      xmid   = c(1 + (511 - 1) / 2, 512 + (env_len - 512) / 2)
    )

    ggplot() +
      ## gene-track backbone (gp120 / gp41)
      geom_rect(data = domains, aes(xmin = xmin, xmax = xmax, ymin = -0.16, ymax = 0.16, fill = domain),
                color = "#0a141b", linewidth = 0.6) +
      geom_text(data = domains, aes(x = xmid, y = 0, label = domain),
                color = "#e6edf3", size = 4.2, fontface = "bold") +
      scale_fill_manual(values = c(gp120 = "#1F6FB2", gp41 = "#5B2C6F"), guide = "none") +
      ## lollipop stems + heads for detected mutations
      geom_segment(data = d_plot, aes(x = pos_num, xend = pos_num, y = 0.16, yend = lollipop_h, color = cat_plot),
                   linewidth = 0.8, alpha = 0.9) +
      geom_point(data = d_plot, aes(x = pos_num, y = lollipop_h, color = cat_plot),
                 size = 3.6, alpha = 0.95) +
      ## missing checklist positions, dropped below the backbone
      geom_segment(data = missing_plot, aes(x = hxb2_sort, xend = hxb2_sort, y = -0.16, yend = stem_bottom),
                   color = "#C0392B", linetype = "22", linewidth = 0.9) +
      geom_point(data = missing_plot, aes(x = hxb2_sort, y = stem_bottom),
                 shape = 21, size = 3.8, color = "#C0392B", fill = "#0f1b24", stroke = 1.4) +
      scale_color_manual(values = category_colors, name = "Category (detected)") +
      scale_x_continuous(limits = c(1, env_len), breaks = seq(0, 800, 100),
                          expand = expansion(mult = 0.015), name = "HXB2 Env position") +
      scale_y_continuous(limits = c(-1.6, 2.4), expand = c(0, 0)) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "#1c2b36", linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "#0f1b24", color = NA),
        panel.background = element_rect(fill = "#0f1b24", color = NA),
        text = element_text(color = "#e6edf3"),
        axis.text.x = element_text(color = "#e6edf3"),
        legend.background = element_rect(fill = "#0f1b24"),
        legend.key = element_rect(fill = "#0f1b24"),
        legend.text = element_text(color = "#e6edf3"),
        legend.title = element_text(color = "#e6edf3"),
        plot.title = element_text(color = "#e6edf3", face = "bold"),
        plot.subtitle = element_text(color = "#9fb0bc", size = 10),
        plot.caption = element_text(color = "#9fb0bc", size = 9, hjust = 0)
      ) +
      labs(
        y = NULL,
        title = paste0("Mutation map -- ", rv$active),
        subtitle = "Lollipops above the backbone = detected mutations, colored by category",
        caption = "Open red circles below the backbone (dashed stem) = missing checklist positions"
      )
  })

  ## ---- Batch Compare tab ----
  batch_matrix <- reactive({
    validate(need(length(rv$strains) > 0, "No strains loaded yet."))
    rows <- lapply(names(rv$strains), function(lbl) {
      st <- build_checklist_status(canonical_mutations, rv$strains[[lbl]]$detected)
      st |> mutate(Strain = lbl)
    })
    bind_rows(rows)
  })

  output$batch_heatmap <- renderPlot({
    m <- tryCatch(batch_matrix(), error = function(e) NULL)
    validate(need(!is.null(m) && nrow(m) > 0, "Load at least one strain on the Upload tab first."))

    m <- m |> mutate(label = factor(label, levels = canonical_mutations$label[order(canonical_mutations$hxb2_sort)]))

    ggplot(m, aes(x = Strain, y = label, fill = status)) +
      geom_tile(color = "#0f1b24", linewidth = 1) +
      scale_fill_manual(values = status_colors, name = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "#e6edf3"),
        axis.text.y = element_text(color = "#e6edf3"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "#0f1b24", color = NA),
        panel.background = element_rect(fill = "#0f1b24", color = NA),
        legend.text = element_text(color = "#e6edf3"),
        text = element_text(color = "#e6edf3")
      ) +
      labs(x = NULL, y = NULL, title = "RnS-DS-SOSIP checklist across loaded strains")
  })

  output$batch_table <- renderDT({
    m <- tryCatch(batch_matrix(), error = function(e) NULL)
    validate(need(!is.null(m) && nrow(m) > 0, "Load at least one strain on the Upload tab first."))
    wide <- m |>
      select(Strain, label, status) |>
      pivot_wider(names_from = Strain, values_from = status)
    datatable(wide, rownames = FALSE, options = list(pageLength = nrow(canonical_mutations), dom = "ftip")) |>
      formatStyle(names(wide)[-1],
                  backgroundColor = styleEqual(c("Present", "Missing"), c("#0f2e1a", "#3a1414")),
                  color = styleEqual(c("Present", "Missing"), c("#5CD68A", "#F08080")))
  })

  ## ---------------------------------------------------------------------------
  ## AI Assistant (offline; Ollama on this machine)
  ## ---------------------------------------------------------------------------
  rv$chat <- list()                                  # list of list(role, content)
  rv$busy <- FALSE
  rv$llm  <- list(state = "unknown", msg = "Not checked yet.")

  check_llm <- function() {
    if (!llm_pkgs_ok) {
      rv$llm <- list(state = "nopkg",
                     msg = "Missing R package(s). Install: install.packages(c('curl','jsonlite','promises','later'))")
      return(invisible())
    }
    models <- ollama_list_models()
    rv$llm <- if (is.null(models)) {
      list(state = "down",
           msg = paste0("Ollama is not reachable at ", OLLAMA_URL, ". Start it with: ollama serve"))
    } else if (!model_installed(OLLAMA_MODEL, models)) {
      list(state = "nomodel",
           msg = paste0("Ollama is running, but '", OLLAMA_MODEL, "' is not installed. Run once (needs internet): ollama pull ", OLLAMA_MODEL))
    } else {
      list(state = "ok", msg = paste0("Ready: ", OLLAMA_MODEL, " at ", OLLAMA_URL))
    }
    invisible()
  }
  observeEvent(input$chat_check, check_llm())        # "Re-check connection" button
  observeEvent(TRUE, check_llm(), once = TRUE)       # and once when the session starts

  output$llm_status <- renderUI({
    st <- rv$llm
    cls <- switch(st$state, ok = "llm-ok", unknown = "llm-warn", nomodel = "llm-warn", "llm-bad")
    tags$p(tags$span(class = paste("llm-dot", cls)), st$msg)
  })

  output$chat_log <- renderUI({
    msgs <- rv$chat
    if (length(msgs) == 0 && !rv$busy) {
      return(div(class = "chat-empty",
                 "Ask about the active strain's checklist, missing positions, repair mutations, ",
                 "or what 3mut / 2G mean. Try a quick question on the right."))
    }
    bubbles <- lapply(msgs, function(m) {
      who <- switch(m$role, user = "You", assistant = "Assistant", "Error")
      div(class = paste("chat-msg", m$role),
          div(class = "chat-role", who),
          div(class = "chat-text", m$content))       # text child => HTML-escaped
    })
    if (isTRUE(rv$busy)) {
      bubbles <- c(bubbles, list(
        div(class = "chat-msg assistant",
            div(class = "chat-role", "Assistant"),
            div(class = "chat-text chat-thinking", "Thinking..."))))
    }
    tagList(bubbles)
  })

  observeEvent(input$chat_clear, {
    rv$chat <- list()
  })

  observeEvent(input$chat_submit, {
    txt <- str_trim(as.character(input$chat_submit$text %||% ""))
    req(nzchar(txt))
    if (isTRUE(rv$busy)) {
      showNotification("Still answering the previous question...", type = "warning")
      return(invisible())
    }
    if (!llm_pkgs_ok) {
      rv$chat <- c(rv$chat, list(list(role = "error", content = rv$llm$msg)))
      return(invisible())
    }

    rv$chat <- c(rv$chat, list(list(role = "user", content = txt)))
    rv$busy <- TRUE

    history <- Filter(function(m) m$role %in% c("user", "assistant"), rv$chat)
    history <- utils::tail(history, 12)              # keep the prompt small for an 8B model
    ctx <- if (isTRUE(input$chat_use_context)) build_llm_context(rv$strains, rv$active) else NULL
    msgs <- c(list(list(role = "system", content = build_llm_system_prompt(ctx))), history)

    p <- ollama_chat_async(msgs, temperature = input$chat_temp %||% 0.2)
    promises::then(
      p,
      onFulfilled = function(reply) {
        if (!nzchar(reply)) reply <- "(The model returned an empty reply.)"
        rv$chat <- c(isolate(rv$chat), list(list(role = "assistant", content = reply)))
        rv$busy <- FALSE
      },
      onRejected = function(e) {
        rv$chat <- c(isolate(rv$chat), list(list(role = "error", content = conditionMessage(e))))
        rv$busy <- FALSE
        check_llm()
      }
    )
    invisible(NULL)                                  # don't hold the session flush on the promise
  })
}

shinyApp(ui, server)
