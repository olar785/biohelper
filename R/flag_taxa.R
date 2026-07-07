#' Flag taxa against expected ecological and geographic context
#'
#' @description
#' Builds an evidence-first ecological interpretation table for checking
#' whether taxa make sense in an expected environment, habitat, and region.
#' By default, `prompt_only = TRUE` returns the deterministic evidence-first
#' table, writes a manual-review prompt to a `.txt` file, and does not call an
#' LLM. The normal LLM-backed workflow is `prompt_only = FALSE` with an ellmer
#' chat object supplied via `chat`. Externally generated structured results can
#' be integrated with `llm_result` or `apply_flag_taxa_llm_result()`;
#' `mock_llm_result` remains available for tests and backward compatibility.
#'
#' @param x A taxonomy table as a data frame, matrix, phyloseq taxonomy table,
#'   or a phyloseq/speedyseq object containing a taxonomy table.
#' @param expected_environment Character scalar describing the expected broad
#'   environment, for example `"marine"`, `"freshwater"`, or `"terrestrial"`.
#' @param expected_habitat Optional character scalar describing the expected
#'   habitat, for example `"estuary"` or `"kelp forest"`.
#' @param expected_region Optional character scalar describing the expected
#'   geographic region.
#' @param tax_ranks Character vector of taxonomic ranks to inspect, ordered from
#'   most specific to broadest.
#' @param evidence_sources Optional character vector of on-demand evidence
#'   sources to fetch before prompt generation. Currently only `"worms"` is
#'   supported. Use `NULL` to avoid automatic evidence retrieval.
#' @param taxon_evidence Optional data frame containing user-supplied evidence
#'   about taxon environment, habitat, and/or region, including output from
#'   `fetch_worms_evidence()`. Required columns are `taxon_name`, `taxon_rank`,
#'   `source`, `evidence_type`, `evidence_summary`, and `reference`. Optional
#'   columns are `accepted_name`, `accepted_rank`, `source_taxon_id`,
#'   `source_record_id`, `environment`, `habitat`, `region`, `locality`,
#'   `decimal_latitude`, `decimal_longitude`, `basis_of_record`,
#'   `occurrence_count`, `reference_url`, `doi`, `checked_at`, and lineage
#'   columns such as `kingdom`, `phylum`, `class`, `order`, `family`, and
#'   `genus`.
#' @param ccz_evidence_path Optional path to a locally downloaded WoRMS
#'   Clarion-Clipperton Zone checklist file. When supplied, the file is read
#'   with `read_worms_ccz_evidence()` and combined with `taxon_evidence` and
#'   any on-demand WoRMS evidence. This never downloads data.
#' @param chat Optional ellmer chat object. When `prompt_only = FALSE`,
#'   `llm_result` and `mock_llm_result` are `NULL`, `flag_taxa()` sends the
#'   generated prompt to `chat$chat_structured()` using an ellmer
#'   structured-output schema. Users are responsible for configuring their own
#'   ellmer provider and API keys.
#' @param allow_llm_tools Logical scalar. If `TRUE`, prompt instructions allow
#'   an ellmer chat object to use tools that the user has already registered on
#'   that chat object, for example Scite MCP literature tools. Tool use is
#'   governed by `tool_requirement`: by default, every taxon selected for
#'   LLM/manual judgement must report a tool attempt; with
#'   `tool_requirement = "optional"`, prompts keep looser selective tool
#'   guidance. `biohelper` does not configure or require these tools.
#' @param prompt_tools Optional character vector naming tools that should be
#'   mentioned explicitly in manual/LLM prompts, for example `"Scite"`. This
#'   only changes prompt instructions; `biohelper` does not configure or call
#'   these tools unless they are already available on the supplied chat object.
#' @param tool_requirement Character scalar controlling how tool traceability is
#'   handled when `allow_llm_tools = TRUE`. `"required_for_llm"` is the default:
#'   every taxon selected for LLM/manual judgement is instructed to use the
#'   registered tool workflow, and missing tool attempts are conservatively
#'   repaired as failed/no-evidence tool attempts before validation.
#'   `"optional"` keeps the looser tool-use guidance for exploratory prompts.
#' @param judgement_mode Character scalar controlling how far `flag_taxa()`
#'   escalates beyond deterministic local evidence. `"evidence_only"` performs
#'   no LLM or tool calls and returns the conservative preliminary judgement.
#'   `"llm_missing_evidence"` sends only taxa with no informative matched local
#'   evidence to the LLM/tools. `"llm_flagged"` sends taxa whose deterministic
#'   action is `flag_for_review` or `flag_possible_exclusion`.
#'   `"llm_possible_exclusion"` sends only taxa whose deterministic action is
#'   `flag_possible_exclusion`. `"llm_all"` sends all unique query taxa. The
#'   older `"llm_review"` value is accepted as a deprecated alias for
#'   `"llm_flagged"`. If omitted, no-chat calls default to `"evidence_only"`
#'   and LLM-backed calls keep the previous `"llm_all"` behaviour for
#'   compatibility.
#' @param prompt_only Logical scalar. When `TRUE`, run deterministic
#'   evidence-first classification, write one prompt `.txt` file for manual LLM
#'   review, and return the evidence-first result table. No ellmer object is
#'   required. When `FALSE`, provide `llm_result`, `mock_llm_result`, or an
#'   ellmer chat object via `chat`.
#' @param prompt_path Optional prompt destination used only when
#'   `prompt_only = TRUE`. If `NULL`, `flag_taxa_prompt.txt` is written in the
#'   working directory. If a directory is supplied, that file is written inside
#'   it. If a file stem/path is supplied, a `.txt` extension is added unless the
#'   path already ends in `.txt`. Multiple prompt chunks are written into this
#'   one file with clear chunk separators.
#' @param verbose Logical scalar. If `TRUE`, report the number of unique taxa
#'   assessed and the number of feature rows used or dropped during taxonomy
#'   preparation.
#' @param worms_batch_size Positive integer scalar. Maximum number of unique
#'   taxa sent to WoRMS per live query batch when on-demand WoRMS evidence is
#'   requested. Defaults to `25`.
#' @param worms_max_tries Positive integer scalar. Maximum number of attempts
#'   for transient WoRMS/API failures per batch. Defaults to `3`.
#' @param worms_retry_sleep Non-negative numeric scalar. Seconds to wait between
#'   WoRMS retry attempts. Defaults to `5`.
#' @param max_evidence_rows_per_query Positive integer scalar. Maximum number
#'   of evidence rows shown in the prompt for each unique query taxon before
#'   matching evidence is summarised. Defaults to `5`.
#' @param max_prompt_chars Optional positive integer scalar. If supplied and the
#'   compact prompt is longer than this value, `flag_taxa()` errors before
#'   calling an LLM. Oversized internal chunks are split automatically where
#'   possible; a one-taxon prompt that is still too large errors clearly.
#' @param max_tool_taxa Non-negative integer scalar or `Inf`. Maximum number of
#'   unique query taxa to send to the optional tool-evidence pass. Defaults to
#'   `25`.
#' @param tool_batch_size Positive integer scalar. Number of candidate taxa per
#'   optional tool-evidence batch. Defaults to `10`.
#' @param llm_result Optional externally generated structured LLM result. This
#'   supports the manual workflow where prompts are generated with
#'   `prompt_only = TRUE`, answered outside `biohelper`, read back into R, and
#'   then validated and joined to feature rows by `flag_taxa()`. The preferred
#'   form is one row per unique query taxon keyed by `query_id`.
#' @param mock_llm_result Optional data frame used for tests and backward
#'   compatibility. It follows the same validation path as `llm_result`. Supply
#'   only one of `llm_result` and `mock_llm_result`.
#'
#' @return A data frame with one row per original input feature. When
#'   `prompt_only = TRUE`, the table is the deterministic evidence-first result
#'   and has a `prompt_path` attribute pointing to the written prompt file.
#'   When `prompt_only = FALSE`, the table includes any supplied or LLM-backed
#'   judgements merged back to all original feature rows.
#'
#' @section Taxon evidence table:
#' `taxon_evidence` can be used to provide curated or previously fetched
#' evidence that should be considered before any future online evidence lookup.
#' The manual workflow is to call `extract_taxa_for_evidence()`, then
#' `fetch_worms_evidence()`, then pass the result as `taxon_evidence`.
#'
#' The convenience workflow is to set `evidence_sources = "worms"`. This performs
#' the same on-demand WoRMS evidence retrieval internally. Live WoRMS queries
#' require the suggested `worrms` package and depend on external service
#' availability. For cached workflows, call `fetch_worms_evidence()` yourself and
#' pass the resulting table through `taxon_evidence`.
#'
#' Regional checklist evidence, such as a local WoRMS Clarion-Clipperton Zone
#' taxlist read with `read_worms_ccz_evidence()`, can be supplied through
#' `taxon_evidence` or `ccz_evidence_path`. Checklist evidence is kept as a
#' separate row from generic WoRMS evidence so context-specific regional and
#' deep-sea information remains visible in the prompt.
#'
#' Before evidence is added to the prompt, `flag_taxa()` filters
#' `taxon_evidence` to rows relevant to the unique input query taxa. Evidence
#' can match exact `taxon_name` or `accepted_name` values, matching taxon/source
#' IDs, or lineage columns such as `family`, `class`, and `phylum` for
#' higher-rank query taxa. Large sets of lineage matches are summarised to avoid
#' prompt overload. Full references, URLs, source IDs, and coordinates remain in
#' the evidence object, but the prompt receives compact evidence fields by
#' default.
#'
#' LLM assessment is done once per unique highest-resolution taxon with a useful
#' assignment, not once per duplicated feature/ASV row. `flag_taxa()` keeps an
#' internal feature-to-taxon map and joins taxon-level assessments back to
#' `feature_id` in the final result. The final table includes all original input
#' feature rows. Features without useful taxonomy are not sent to the LLM; they
#' are returned as deterministic `flag_for_review` rows with
#' `insufficient_taxonomic_resolution` statuses. The compact lineage sent to the
#' LLM contains available ranks from broadest to most specific, for example
#' `Phylum=Porifera; Class=Demospongiae; Order=NA; Family=NA; Genus=Plenaster;
#' Species=Plenaster craigi`. Missing useful ranks are represented as `NA`.
#'
#' The prompt stays compact, but the returned table includes richer deterministic
#' context and provenance columns: `expected_environment`, `expected_habitat`,
#' `expected_region`, `lineage`, `worms_environment`, `evidence_sources`, and
#' `evidence_summary`. Full reference URLs and source identifiers remain in
#' evidence objects rather than being pasted into the LLM prompt by default.
#'
#' @section Ecological interpretation columns:
#' The final output separates component evidence statuses from user action.
#' `ecological_status` describes ecological fit as `compatible`, `plausible`,
#' `unlikely_resident`, `incompatible_resident`, `unknown`, or
#' `insufficient_taxonomic_resolution`. `occurrence_interpretation` describes
#' how to interpret the sequence in the sample: `expected_resident`,
#' `plausible_resident`, `possible_transient_or_allochthonous`,
#' `possible_contaminant_or_misassignment`,
#' `likely_contaminant_or_misassignment`, `uncertain`, or `not_assessed`.
#' `recommended_action` is one of `retain`, `flag_for_review`,
#' `flag_possible_exclusion`, or `exclude`. `flag_possible_exclusion` means the
#' taxon is likely inappropriate for resident-community interpretation but
#' should not be removed automatically without user review.
#'
#' Absence of records is not ecological incompatibility. Known ecology or
#' defining biology can, however, support `unlikely_resident` or
#' `incompatible_resident`. For example, photosynthetic taxa in aphotic
#' deep-sea benthic samples may be resident-incompatible while still plausible
#' as transported or allochthonous DNA. Potential contamination or taxonomic
#' misassignment is represented separately from allochthonous DNA.
#'
#' @section Evidence-first judgement:
#' `flag_taxa()` first builds a deterministic preliminary judgement from
#' matched local evidence for each unique query taxon. This is deliberately
#' conservative: missing evidence, WoRMS no-matches, lookup failures, absent
#' CCZ rows, absent Scite/tool records, and absent deep-sea records are treated
#' as uncertainty rather than incompatibility. A taxon with local marine
#' evidence and an expected marine habitat subtype such as `"deep sea"` is not
#' excluded only because no explicit deep-sea record was found. For broad ranks
#' above species, lack of habitat or region records normally remains
#' `flag_for_review`, often with `mixed_within_rank`,
#' `insufficient_taxonomic_resolution`,
#' `no_known_habitat_evidence`, `possible_likely`, `possible_unlikely`, or
#' `unknown`.
#'
#' `judgement_mode = "evidence_only"` returns this preliminary classification
#' without requiring `chat`. Use it for fast conservative screening.
#' `judgement_mode = "llm_flagged"` sends taxa with deterministic
#' `flag_for_review` or `flag_possible_exclusion` actions for LLM/manual
#' escalation. Use `"llm_possible_exclusion"` when you only want LLM/tools to
#' check taxa with a positive ecological concern that may affect exclusion or
#' resident-community interpretation. Use `"llm_missing_evidence"` when local
#' evidence should be trusted and the LLM should only fill gaps. Use
#' `"llm_all"` for exploratory or audit runs; it is slower and requires careful
#' interpretation because models can over-read absence of evidence despite
#' conservative prompts. `judgement_mode` selects taxa for LLM/manual
#' judgement; `allow_llm_tools` only controls whether those selected taxa may
#' use registered tools such as Scite.
#'
#' `expected_habitat_status` includes `no_known_habitat_evidence`,
#' `possible_likely`, and `possible_unlikely`. These values represent missing
#' or plausibility-level habitat evidence, not positive habitat
#' incompatibility. `exclude` should be used only when explicit positive
#' evidence demonstrates strong incompatibility, such as species-level local
#' evidence showing restriction to freshwater/terrestrial environments when the
#' expected environment is marine. For broad-rank or allochthonous-DNA cases,
#' prefer `flag_possible_exclusion` over `exclude`.
#'
#' `prompt_only = TRUE` supports prompt inspection, prompt-size debugging, and
#' manual/external LLM use when API access is not available. It writes one
#' `.txt` prompt file and returns the evidence-first result table. The prompt
#' asks the external LLM to create a downloadable JSON file named
#' `flag_taxa_llm_result.json`, with one object for each listed `query_id`. If
#' the interface cannot create a file, it should output only the raw JSON array
#' so the user can save it. Use `apply_flag_taxa_llm_result()` to read and
#' validate that JSON and merge it into the evidence-first table.
#' Deterministic rows not selected for LLM/manual judgement are preserved
#' unchanged and are marked with `llm_prompt_selected = FALSE`.
#'
#' If `allow_llm_tools = TRUE`, the prompt tells the LLM that registered tools
#' may be available through the supplied chat object and may be used only to
#' gather explicit supporting evidence. With the default
#' `tool_requirement = "required_for_llm"`, every taxon selected by
#' `judgement_mode` for LLM/manual judgement is instructed to use the registered
#' tool workflow before returning final JSON and to report `llm_tool_name`,
#' `llm_tool_query`, `llm_tool_evidence_summary`, and
#' `llm_tool_references`. If no evidence is found, the tool attempt is recorded
#' explicitly as missing evidence, not incompatibility. With
#' `tool_requirement = "optional"`, the prompt keeps looser selective tool
#' guidance for weak, ambiguous, contradictory, flagged, or possible-exclusion
#' cases.
#'
#' Tools should be used to find positive evidence, not to prove absence. If
#' tool evidence does not resolve the case, the result should remain
#' `flag_for_review`. `biohelper` does not call
#' `mcptools::mcp_tools()`, read MCP configuration, or depend on Scite; Scite MCP
#' is only an example of a tool the user may register on an ellmer chat object
#' before calling `flag_taxa()`. Local `taxon_evidence` remains the first source
#' of evidence. If tools are unavailable or no explicit evidence is found, the
#' result should remain conservative.
#'
#' When `prompt_only = FALSE`, `allow_llm_tools = TRUE`, and `chat` is supplied,
#' `flag_taxa()` can run an explicit, targeted tool-evidence pass before the
#' final structured-output judgement. `biohelper` still does not configure or
#' call Scite directly; it only asks the user-supplied ellmer chat object to use
#' tools that the user already registered. Candidate taxa follow
#' `judgement_mode`: no tools for `"evidence_only"`, missing-evidence taxa for
#' `"llm_missing_evidence"`, `flag_for_review` and
#' `flag_possible_exclusion` taxa for `"llm_flagged"`, only
#' `flag_possible_exclusion` taxa for `"llm_possible_exclusion"`, and all
#' unique taxa for `"llm_all"`. Tools are not used for deterministic retain rows
#' except in `"llm_all"`. Tool evidence gathering is batched separately with
#' `tool_batch_size` and may substantially increase runtime. If a required tool
#' attempt fails or is missing from the LLM result, the default repair path marks
#' the row with conservative failed/no-evidence tool metadata rather than
#' allowing fabricated references or stopping the whole run.
#'
#' When `verbose = TRUE`, progress messages report elapsed time and distinguish
#' original feature rows, feature rows with useful taxonomy, feature rows
#' returned as not assessed, unique taxa assessed by the LLM, local evidence
#' rows and coverage, prompt sizes, optional tool-evidence batches, final LLM
#' calls, and final returned rows.
#'
#' HTTP 413 errors from an LLM backend indicate that the final judgement prompt
#' or request payload is too large for the selected model route, not provider
#' quota. HTTP 429 errors indicate LLM provider request quota/rate limiting; they
#' do not necessarily mean Scite/tool quota and do not necessarily indicate a
#' prompt-size problem. Premature EOF, parse errors, or truncated JSON from a
#' structured-output call usually mean the prompt/output was too large for the
#' model or output limit. Small model routes, including some GitHub model
#' routes, can be too constrained for tool-enabled analyses; use a
#' larger-context model where possible.
#'
#' Absence of evidence is not evidence of incompatibility. Tool searches that
#' find no explicit evidence should normally lead to `unknown`,
#' `no_known_habitat_evidence`, `no_distribution_evidence`, and
#' `flag_for_review`, not exclusion. `exclude` requires explicit positive
#' evidence of incompatibility.
#' Scite or other registered tools are useful for finding explicit supporting
#' evidence; they should not be used as absence checkers that prove a taxon
#' cannot occur in a habitat or region. Broad-rank no-record cases should
#' usually remain `flag_for_review`, `mixed_within_rank`, or
#' `insufficient_taxonomic_resolution`.
#'
#' LLM/tool traceability is reported through self-reported structured-output
#' columns: `evidence_basis`, `llm_tool_used`, `llm_tool_name`,
#' `llm_tool_query`, `llm_tool_evidence_summary`, and
#' `llm_tool_references`. These columns describe whether the LLM relied on
#' local evidence, registered tools, both, or conservative reasoning only. When
#' a tool is used, `references` and `llm_tool_references` should contain
#' verifiable titles, short citations, DOIs, or URLs when available; if the tool
#' returned no usable references, `llm_tool_references` should say so
#' explicitly. AphiaIDs, DOIs, URLs, article titles, author names, and formal
#' citations must come from matched local evidence or returned tool evidence;
#' general biological judgement belongs in `rationale` and must not be turned
#' into fabricated references. These columns are distinct from the deterministic local evidence
#' provenance columns `evidence_sources` and `evidence_summary`, which summarise
#' matched `taxon_evidence` rows passed into the prompt.
#'
#' The required columns are: `taxon_name`, `taxon_rank`, `source`,
#' `evidence_type`, `evidence_summary`, and `reference`.
#'
#' Optional columns are: `accepted_name`, `accepted_rank`, `source_taxon_id`,
#' `source_record_id`, `environment`, `habitat`, `region`, `locality`,
#' `decimal_latitude`, `decimal_longitude`, `basis_of_record`,
#' `occurrence_count`, `reference_url`, `doi`, `checked_at`, `kingdom`,
#' `phylum`, `class`, `order`, `family`, and `genus`.
#' @export
#'
#' @examples
#' tax <- data.frame(
#'   kingdom = "Animalia",
#'   phylum = "Chordata",
#'   family = "Salmonidae",
#'   genus = "Salmo",
#'   species = "Salmo salar"
#' )
#'
#' evidence_first <- flag_taxa(
#'   tax,
#'   expected_environment = "marine",
#'   expected_habitat = "coastal water",
#'   expected_region = "North Atlantic",
#'   prompt_path = tempdir()
#' )
#'
#' \dontrun{
#' data("ps_test_data_euk")
#'
#' # Manual evidence workflow
#' taxa_to_query <- extract_taxa_for_evidence(ps_test_data_euk)
#' worms_evidence <- fetch_worms_evidence(
#'   taxa_to_query,
#'   by = "name",
#'   cache_path = "worms_cache.rds"
#' )
#' evidence_first <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   taxon_evidence = worms_evidence,
#'   prompt_only = TRUE
#' )
#' attr(evidence_first, "prompt_path")
#' # final <- apply_flag_taxa_llm_result(
#' #   evidence_first,
#' #   llm_json_path = "flag_taxa_llm_result.json"
#' # )
#'
#' # Convenience workflow
#' evidence_first <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   evidence_sources = "worms",
#'   prompt_only = TRUE
#' )
#'
#' # Local CCZ checklist evidence can be combined with generic WoRMS evidence.
#' ccz_evidence <- read_worms_ccz_evidence("CCZ_taxlist.csv")
#' taxon_evidence <- dplyr::bind_rows(worms_evidence, ccz_evidence)
#' evidence_first <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   expected_habitat = "deep sea",
#'   expected_region = "Clarion-Clipperton Zone",
#'   taxon_evidence = taxon_evidence,
#'   prompt_only = TRUE
#' )
#'
#' # Optional ellmer-backed workflow. Configure the provider/API key yourself.
#' chat <- ellmer::chat_openai(model = "gpt-4.1-mini")
#' result <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   expected_habitat = "deep sea",
#'   taxon_evidence = worms_evidence,
#'   chat = chat,
#'   prompt_only = FALSE
#' )
#'
#' # Optional tool-enabled ellmer workflow. Tools are configured outside
#' # biohelper; Scite MCP is an example, not a required dependency.
#' chat <- ellmer::chat_google_gemini(api_key = apikey_gemini)
#' scite_tools <- mcptools::mcp_tools(
#'   config = "~/.config/Code/mcpServers.json"
#' )
#' chat$set_tools(scite_tools)
#' result <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   expected_habitat = "deep sea",
#'   expected_region = "Clarion-Clipperton Zone",
#'   taxon_evidence = taxon_evidence,
#'   chat = chat,
#'   allow_llm_tools = TRUE,
#'   prompt_only = FALSE
#' )
#'
#' # Manual prompt/debug workflow
#' prompt <- flag_taxa(
#'   ps_test_data_euk,
#'   expected_environment = "marine",
#'   taxon_evidence = taxon_evidence,
#'   prompt_only = TRUE,
#'   prompt_path = "flag_taxa_prompt.txt"
#' )
#' # external_result <- jsonlite::fromJSON("external_flag_taxa_result.json")
#' # result <- flag_taxa(
#' #   ps_test_data_euk,
#' #   expected_environment = "marine",
#' #   taxon_evidence = taxon_evidence,
#' #   prompt_only = FALSE,
#' #   llm_result = external_result
#' # )
#' }
flag_taxa <- function(
  x,
  expected_environment,
  expected_habitat = NULL,
  expected_region = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  evidence_sources = NULL,
  taxon_evidence = NULL,
  ccz_evidence_path = NULL,
  chat = NULL,
  allow_llm_tools = FALSE,
  prompt_tools = NULL,
  tool_requirement = c("required_for_llm", "optional"),
  judgement_mode = c(
    "evidence_only",
    "llm_missing_evidence",
    "llm_flagged",
    "llm_possible_exclusion",
    "llm_all"
  ),
  prompt_only = TRUE,
  prompt_path = NULL,
  verbose = FALSE,
  worms_batch_size = 25,
  worms_max_tries = 3,
  worms_retry_sleep = 5,
  max_evidence_rows_per_query = 5,
  max_prompt_chars = NULL,
  max_tool_taxa = 25,
  tool_batch_size = 10,
  llm_result = NULL,
  mock_llm_result = NULL
) {
  start_time <- Sys.time()
  judgement_mode_missing <- missing(judgement_mode)
  llm_max_tries <- 3L
  llm_retry_sleep <- 5
  repair_invalid_llm <- TRUE
  worms_cache_path <- NULL
  worms_sleep <- 0.2
  continue_on_worms_error <- TRUE
  worms_verbose <- verbose
  warn_prompt_chars <- 30000L
  max_evidence_summary_chars <- 200L
  max_tool_evidence_summary_chars <- 300L
  max_lineage_chars <- 300L
  auto_reduce_llm_chunk_size <- TRUE
  tool_use_policy <- if (isTRUE(allow_llm_tools)) "review_or_exclude" else "never"
  run_tools_in_prompt_only <- FALSE
  llm_chunk_size <- Inf
  expected_environment <- .validate_required_scalar_character(
    expected_environment,
    "expected_environment"
  )
  expected_habitat <- .validate_optional_scalar_character(
    expected_habitat,
    "expected_habitat"
  )
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_flag_taxa_evidence_sources(evidence_sources)
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  ccz_evidence_path <- .validate_optional_path(ccz_evidence_path, "ccz_evidence_path")
  allow_llm_tools <- .validate_logical_scalar(allow_llm_tools, "allow_llm_tools")
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  prompt_only <- .validate_logical_scalar(prompt_only, "prompt_only")
  prompt_path <- .validate_optional_path(prompt_path, "prompt_path")
  if (!is.null(llm_result) && !is.null(mock_llm_result)) {
    stop(
      "Provide only one of `llm_result` and `mock_llm_result`.",
      call. = FALSE
    )
  }
  supplied_llm_result <- if (!is.null(llm_result)) llm_result else mock_llm_result
  judgement_mode <- .validate_flag_taxa_judgement_mode(
    judgement_mode = judgement_mode,
    use_default = judgement_mode_missing,
    prompt_only = prompt_only,
    chat = chat,
    supplied_llm_result = supplied_llm_result
  )
  if (identical(judgement_mode, "evidence_only")) {
    allow_llm_tools <- FALSE
  }
  verbose <- .validate_logical_scalar(verbose, "verbose")
  worms_cache_path <- .validate_optional_cache_path(worms_cache_path)
  worms_sleep <- .validate_non_negative_numeric_scalar(worms_sleep, "worms_sleep")
  worms_batch_size <- .validate_positive_integer_scalar(worms_batch_size, "worms_batch_size")
  worms_max_tries <- .validate_positive_integer_scalar(worms_max_tries, "worms_max_tries")
  worms_retry_sleep <- .validate_non_negative_numeric_scalar(worms_retry_sleep, "worms_retry_sleep")
  worms_verbose <- .validate_logical_scalar(worms_verbose, "worms_verbose")
  max_evidence_rows_per_query <- .validate_positive_integer_scalar(
    max_evidence_rows_per_query,
    "max_evidence_rows_per_query"
  )
  max_prompt_chars <- .validate_optional_positive_integer_scalar(
    max_prompt_chars,
    "max_prompt_chars"
  )
  max_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_evidence_summary_chars,
    "max_evidence_summary_chars"
  )
  max_tool_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_tool_evidence_summary_chars,
    "max_tool_evidence_summary_chars"
  )
  max_lineage_chars <- .validate_positive_integer_scalar(
    max_lineage_chars,
    "max_lineage_chars"
  )
  max_tool_taxa <- .validate_non_negative_count_or_inf(max_tool_taxa, "max_tool_taxa")
  tool_batch_size <- .validate_positive_integer_scalar(tool_batch_size, "tool_batch_size")
  tax_table <- extract_tax_table(x)
  .flag_taxa_progress(
    verbose,
    start_time,
    "input contains ",
    nrow(tax_table),
    " feature rows."
  )
  all_feature_map <- .prepare_flag_taxa_all_feature_map(
    tax_table = tax_table,
    tax_ranks = tax_ranks
  )
  feature_query_map <- .retained_flag_taxa_feature_query_map(all_feature_map)
  not_assessed_rows <- nrow(all_feature_map) - nrow(feature_query_map)
  .flag_taxa_progress(
    verbose,
    start_time,
    "retained ",
    nrow(feature_query_map),
    " feature rows with useful taxonomy; ",
    not_assessed_rows,
    " will be returned as not assessed."
  )
  unique_query_table <- .unique_flag_taxa_query_table(feature_query_map)
  feature_query_map <- .assign_flag_taxa_query_ids_to_features(
    feature_query_map,
    unique_query_table
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "assessing ",
    nrow(unique_query_table),
    " unique taxa from ",
    nrow(feature_query_map),
    " retained feature rows."
  )

  if (is.null(supplied_llm_result)) {
    taxon_evidence <- .collect_flag_taxa_requested_evidence(
      x = x,
      evidence_sources = evidence_sources,
      taxon_evidence = taxon_evidence,
      ccz_evidence_path = ccz_evidence_path,
      worms_cache_path = worms_cache_path,
      worms_sleep = worms_sleep,
      worms_batch_size = worms_batch_size,
      worms_max_tries = worms_max_tries,
      worms_retry_sleep = worms_retry_sleep,
      continue_on_worms_error = continue_on_worms_error,
      worms_verbose = worms_verbose
    )
  }

  evidence_stats <- .flag_taxa_evidence_match_stats(
    unique_query_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "matched ",
    evidence_stats$evidence_rows,
    " local evidence rows covering ",
    evidence_stats$query_taxa_with_evidence,
    " / ",
    nrow(unique_query_table),
    " unique taxa."
  )

  preliminary_judgement <- .build_flag_taxa_preliminary_judgement(
    unique_query_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    tax_ranks = tax_ranks
  )
  llm_query_table <- .select_flag_taxa_llm_query_table(
    unique_query_table = unique_query_table,
    preliminary_judgement = preliminary_judgement,
    judgement_mode = judgement_mode
  )
  prompt_query_table <- llm_query_table
  .flag_taxa_progress(
    verbose,
    start_time,
    "judgement_mode = '",
    judgement_mode,
    "' selected ",
    nrow(llm_query_table),
    " / ",
    nrow(unique_query_table),
    " unique taxa for LLM escalation."
  )
  if (isTRUE(prompt_only)) {
    resolved_taxa <- max(nrow(unique_query_table) - nrow(prompt_query_table), 0L)
    .flag_taxa_progress(
      verbose,
      start_time,
      "prompt includes ",
      nrow(prompt_query_table),
      " / ",
      nrow(unique_query_table),
      " unique taxa selected for LLM/manual judgement."
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      resolved_taxa,
      " / ",
      nrow(unique_query_table),
      " unique taxa were resolved deterministically and are not included in the prompt."
    )
  }

  needs_llm <- !isTRUE(prompt_only) &&
    is.null(supplied_llm_result) &&
    !identical(judgement_mode, "evidence_only") &&
    nrow(llm_query_table) > 0
  needs_tool_pass <- .should_run_flag_taxa_tool_pass(
    allow_llm_tools = allow_llm_tools,
    tool_use_policy = tool_use_policy,
    prompt_only = prompt_only,
    run_tools_in_prompt_only = run_tools_in_prompt_only
  ) && is.null(supplied_llm_result) && nrow(llm_query_table) > 0
  if ((needs_llm || needs_tool_pass) && is.null(chat)) {
    stop(
      "No LLM backend supplied for judgement_mode = '",
      judgement_mode,
      "'. Use judgement_mode = 'evidence_only', use prompt_only = TRUE, ",
      "provide `llm_result` for an external structured result, provide ",
      "`mock_llm_result` for testing, or pass an ellmer chat object via `chat`.",
      call. = FALSE
    )
  }

  tool_evidence <- NULL
  if (needs_tool_pass && nrow(llm_query_table) > 0) {
    .flag_taxa_progress(
      verbose,
      start_time,
      "starting optional tool evidence pass using policy '",
      tool_use_policy,
      "'."
    )
    tool_evidence <- tryCatch(
      .run_flag_taxa_tool_evidence_pass(
        unique_query_table = llm_query_table,
        preliminary_judgement = preliminary_judgement,
        judgement_mode = judgement_mode,
        taxon_evidence = taxon_evidence,
        expected_environment = expected_environment,
        expected_habitat = expected_habitat,
        expected_region = expected_region,
        tax_ranks = tax_ranks,
        max_evidence_rows_per_query = max_evidence_rows_per_query,
        allow_llm_tools = allow_llm_tools,
        tool_use_policy = tool_use_policy,
        max_tool_taxa = max_tool_taxa,
        tool_batch_size = tool_batch_size,
        chat = chat,
        llm_max_tries = llm_max_tries,
        llm_retry_sleep = llm_retry_sleep,
        verbose = verbose,
        start_time = start_time
      ),
      error = function(error) {
        if (!identical(tool_requirement, "required_for_llm")) {
          stop(error)
        }
        .flag_taxa_progress(
          verbose,
          start_time,
          "required tool evidence pass failed; recording conservative failed-tool rows for ",
          nrow(llm_query_table),
          " selected taxa."
        )
        .flag_taxa_required_tool_failure_evidence(
          query_table = llm_query_table,
          error_message = conditionMessage(error),
          prompt_tools = prompt_tools,
          expected_environment = expected_environment,
          expected_habitat = expected_habitat,
          expected_region = expected_region
        )
      }
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed optional tool evidence pass; returned ",
      nrow(tool_evidence),
      " tool evidence rows."
    )
  } else {
    .flag_taxa_progress(
      verbose,
      start_time,
      "skipping optional tool evidence pass."
    )
  }

  prompt_chunks <- .build_flag_taxa_prompt_chunks(
    tax_table = tax_table,
    unique_query_table = prompt_query_table,
    preliminary_judgement = preliminary_judgement,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    tax_ranks = tax_ranks,
    evidence_sources = evidence_sources,
    taxon_evidence = taxon_evidence,
    tool_evidence = tool_evidence,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars,
    max_tool_evidence_summary_chars = max_tool_evidence_summary_chars,
    max_lineage_chars = max_lineage_chars,
    allow_llm_tools = allow_llm_tools,
    prompt_tools = prompt_tools,
    tool_requirement = tool_requirement,
    llm_chunk_size = llm_chunk_size,
    max_prompt_chars = max_prompt_chars,
    auto_reduce_llm_chunk_size = auto_reduce_llm_chunk_size
  )
  query_chunks <- attr(prompt_chunks, "query_chunks", exact = TRUE)
  .report_flag_taxa_prompt_chunks(
    prompt_chunks = prompt_chunks,
    query_chunks = query_chunks,
    verbose = verbose,
    start_time = start_time
  )
  if (!isTRUE(prompt_only)) {
    .warn_flag_taxa_large_prompt_chunks(
      prompt_chunks = prompt_chunks,
      query_chunks = query_chunks,
      warn_prompt_chars = warn_prompt_chars,
      llm_chunk_size = llm_chunk_size
    )
  }

  if (isTRUE(prompt_only)) {
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed prompt generation for ",
      nrow(unique_query_table),
      " unique taxa."
    )
    resolved_prompt_path <- .write_flag_taxa_prompt_file(
      prompt_chunks = prompt_chunks,
      prompt_path = prompt_path
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "wrote prompt to ",
      normalizePath(resolved_prompt_path, mustWork = FALSE),
      "."
    )
    out <- .build_flag_taxa_evidence_first_output(
      preliminary_judgement = preliminary_judgement,
      feature_query_map = feature_query_map,
      unique_query_table = unique_query_table,
      selected_query_table = prompt_query_table,
      all_feature_map = all_feature_map,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      taxon_evidence = taxon_evidence,
      tool_evidence = tool_evidence,
      tax_ranks = tax_ranks,
      max_evidence_rows_per_query = max_evidence_rows_per_query,
      max_evidence_summary_chars = max_evidence_summary_chars
    )
    attr(out, "prompt_path") <- resolved_prompt_path
    attr(out, "prompt_chunks") <- names(prompt_chunks)
    attr(out, "allow_llm_tools") <- allow_llm_tools
    attr(out, "tool_requirement") <- tool_requirement
    attr(out, "prompt_tools") <- prompt_tools
    return(out)
  }

  if (!is.null(supplied_llm_result)) {
    assessed <- .process_flag_taxa_llm_result(
      result = supplied_llm_result,
      feature_query_map = feature_query_map,
      unique_query_table = unique_query_table,
      llm_query_table = .resolve_flag_taxa_supplied_result_query_table(
        result = supplied_llm_result,
        unique_query_table = unique_query_table,
        llm_query_table = llm_query_table
      ),
      preliminary_judgement = preliminary_judgement,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      taxon_evidence = taxon_evidence,
      tool_evidence = tool_evidence,
      tax_ranks = tax_ranks,
      max_evidence_rows_per_query = max_evidence_rows_per_query,
      max_evidence_summary_chars = max_evidence_summary_chars,
      allow_llm_tools = allow_llm_tools,
      tool_requirement = tool_requirement,
      prompt_tools = prompt_tools,
      repair_invalid_llm = repair_invalid_llm,
      verbose = verbose,
      start_time = start_time
    )
    out <- .restore_unassessed_flag_taxa_rows(
      assessed_output = assessed,
      all_feature_map = all_feature_map,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed. Returned ",
      nrow(out),
      " feature-level rows, including ",
      not_assessed_rows,
      " not assessed rows.",
      total = TRUE
    )
    return(out)
  }

  if (nrow(unique_query_table) == 0) {
    out <- .restore_unassessed_flag_taxa_rows(
      assessed_output = NULL,
      all_feature_map = all_feature_map,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed. Returned ",
      nrow(out),
      " feature-level rows, including ",
      not_assessed_rows,
      " not assessed rows.",
      total = TRUE
    )
    return(out)
  }

  if (!needs_llm) {
    out <- .build_flag_taxa_evidence_first_output(
      preliminary_judgement = preliminary_judgement,
      feature_query_map = feature_query_map,
      unique_query_table = unique_query_table,
      selected_query_table = llm_query_table,
      all_feature_map = all_feature_map,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      taxon_evidence = taxon_evidence,
      tool_evidence = tool_evidence,
      tax_ranks = tax_ranks,
      max_evidence_rows_per_query = max_evidence_rows_per_query,
      max_evidence_summary_chars = max_evidence_summary_chars
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed evidence-first deterministic assessment. Returned ",
      nrow(out),
      " feature-level rows, including ",
      not_assessed_rows,
      " not assessed rows.",
      total = TRUE
    )
    return(out)
  }

  .flag_taxa_progress(
    verbose,
    start_time,
    .format_flag_taxa_llm_start_message(
      unique_taxa = nrow(llm_query_table),
      chunks = length(prompt_chunks),
      llm_chunk_size = llm_chunk_size
    )
  )
  llm_result <- .call_flag_taxa_ellmer_chunks(
    prompt_chunks = prompt_chunks,
    query_chunks = query_chunks,
    chat = chat,
    expected_region = expected_region,
    max_tries = llm_max_tries,
    retry_sleep = llm_retry_sleep,
    allow_llm_tools = allow_llm_tools,
    tool_requirement = tool_requirement,
    prompt_tools = prompt_tools,
    repair_invalid_llm = repair_invalid_llm,
    verbose = verbose,
    start_time = start_time
  )
  .flag_taxa_progress(verbose, start_time, "LLM call completed.")

  assessed <- .process_flag_taxa_llm_result(
    result = llm_result,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table,
    llm_query_table = llm_query_table,
    preliminary_judgement = preliminary_judgement,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    taxon_evidence = taxon_evidence,
    tool_evidence = tool_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars,
    allow_llm_tools = allow_llm_tools,
    repair_invalid_llm = repair_invalid_llm,
    verbose = verbose,
    start_time = start_time
  )
  out <- .restore_unassessed_flag_taxa_rows(
    assessed_output = assessed,
    all_feature_map = all_feature_map,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "validated and joined LLM output to feature rows."
  )
  .flag_taxa_progress(
    verbose,
    start_time,
    "completed. Returned ",
    nrow(out),
    " feature-level rows, including ",
    not_assessed_rows,
    " not assessed rows.",
    total = TRUE
  )
  out
}

#' Apply a manual flag_taxa LLM result
#'
#' @description
#' Integrates a structured JSON/data-frame result from a manually answered
#' `flag_taxa(prompt_only = TRUE)` prompt into the evidence-first table returned
#' by that call.
#'
#' @param result Evidence-first result table returned by
#'   `flag_taxa(prompt_only = TRUE)`.
#' @param llm_result Optional data frame or list parsed from JSON. Supply either
#'   `llm_result` or `llm_json_path`.
#' @param llm_json_path Optional path to a JSON file containing an array of LLM
#'   result objects. Read with `jsonlite::fromJSON()`.
#'
#' @return A `data.frame` in the same final output format as
#'   `flag_taxa(prompt_only = FALSE)`, with LLM/manual judgements merged into
#'   rows where `llm_selected = TRUE` and deterministic rows preserved.
#' @export
#'
#' @examples
#' tax <- data.frame(
#'   phylum = "Chordata",
#'   genus = "Salmo",
#'   species = "Salmo salar"
#' )
#'
#' evidence_first <- flag_taxa(
#'   tax,
#'   expected_environment = "marine",
#'   prompt_only = TRUE,
#'   prompt_path = tempdir()
#' )
#'
#' # After manually saving a valid JSON response:
#' # final <- apply_flag_taxa_llm_result(
#' #   evidence_first,
#' #   llm_json_path = "flag_taxa_llm_result.json"
#' # )
apply_flag_taxa_llm_result <- function(
  result,
  llm_result = NULL,
  llm_json_path = NULL
) {
  if (!inherits(result, "data.frame")) {
    stop("`result` must be the data.frame returned by `flag_taxa(prompt_only = TRUE)`.", call. = FALSE)
  }
  if (!is.null(llm_result) && !is.null(llm_json_path)) {
    stop("Supply only one of `llm_result` and `llm_json_path`.", call. = FALSE)
  }
  if (is.null(llm_result) && is.null(llm_json_path)) {
    stop("Supply `llm_result` or `llm_json_path`.", call. = FALSE)
  }
  if (!("llm_selected" %in% colnames(result)) && "llm_prompt_selected" %in% colnames(result)) {
    result$llm_selected <- result$llm_prompt_selected
  }
  required_result_columns <- c("query_id", "taxon_name", "taxon_rank", "lineage", "llm_selected")
  missing_result_columns <- setdiff(required_result_columns, colnames(result))
  if (length(missing_result_columns) > 0) {
    stop(
      "`result` is missing required manual integration columns: ",
      paste(missing_result_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(llm_json_path)) {
    llm_result <- .read_flag_taxa_llm_json(llm_json_path)
  }
  llm_result <- .as_flag_taxa_result_table(llm_result)

  selected <- as.logical(result$llm_selected)
  selected[is.na(selected)] <- FALSE
  selected <- selected & !is.na(result$query_id)
  query_rows <- result[selected, c("query_id", "taxon_name", "taxon_rank", "lineage"), drop = FALSE]
  query_rows$query_name <- query_rows$taxon_name
  query_rows$query_rank <- query_rows$taxon_rank
  query_rows <- query_rows[!duplicated(as.character(query_rows$query_id)), , drop = FALSE]
  query_taxonomy <- query_rows[, c("query_id", "query_name", "query_rank", "lineage"), drop = FALSE]
  rownames(query_taxonomy) <- NULL

  expected_region <- .flag_taxa_result_expected_region(result)
  allow_llm_tools <- isTRUE(attr(result, "allow_llm_tools", exact = TRUE))
  tool_requirement <- attr(result, "tool_requirement", exact = TRUE)
  if (is.null(tool_requirement)) {
    tool_requirement <- "optional"
  }
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  prompt_tools <- attr(result, "prompt_tools", exact = TRUE)
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  llm_result <- .repair_flag_taxa_required_tool_attempts(
    result = llm_result,
    tool_requirement = tool_requirement,
    allow_llm_tools = allow_llm_tools,
    repair_invalid_llm = TRUE,
    query_taxonomy = query_taxonomy,
    prompt_tools = prompt_tools,
    verbose = FALSE,
    start_time = NULL
  )
  llm_result <- .repair_flag_taxa_invalid_llm_output(
    result = llm_result,
    repair_invalid_llm = TRUE,
    expected_region = expected_region
  )
  validated <- validate_flag_taxa_taxon_output(
    result = llm_result,
    query_taxonomy = query_taxonomy,
    expected_region = expected_region,
    allow_llm_tools = allow_llm_tools ||
      any(as.logical(.optional_result_column(result, "llm_tool_used", nrow(result))), na.rm = TRUE)
  )

  out <- result
  update_columns <- intersect(
    setdiff(.flag_taxa_structured_output_columns(include_feature_id = FALSE, include_query_id = TRUE), "query_id"),
    colnames(out)
  )
  for (row_index in seq_len(nrow(validated))) {
    query_id <- as.character(validated$query_id[[row_index]])
    target <- which(as.character(out$query_id) == query_id)
    if (length(target) == 0) {
      next
    }
    for (column_name in update_columns) {
      if (column_name %in% colnames(validated)) {
        out[[column_name]][target] <- validated[[column_name]][[row_index]]
      }
    }
  }

  .finalise_flag_taxa_output(out, expected_region = expected_region)
}

.read_flag_taxa_llm_json <- function(llm_json_path) {
  llm_json_path <- .validate_optional_path(llm_json_path, "llm_json_path")
  if (is.null(llm_json_path) || !file.exists(llm_json_path)) {
    stop("`llm_json_path` must point to an existing JSON file.", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "The jsonlite package is required to read `llm_json_path`. Install it with install.packages('jsonlite') or pass `llm_result` directly.",
      call. = FALSE
    )
  }
  jsonlite::fromJSON(llm_json_path, simplifyDataFrame = TRUE)
}

.flag_taxa_result_expected_region <- function(result) {
  if (!("expected_region" %in% colnames(result))) {
    return(NULL)
  }
  values <- as.character(result$expected_region)
  values <- values[.is_non_empty_value(values)]
  values <- values[!tolower(values) %in% c("na", "null")]
  if (length(values) == 0) {
    return(NULL)
  }
  unique(values)[[1]]
}

#' Validate a flag_taxa result table
#'
#' @param result Data frame returned by an LLM-like review.
#' @param original_taxonomy Original taxonomy table used to build the prompt.
#' @param expected_region Optional expected region supplied to `flag_taxa()`.
#'
#' @return `result`, unchanged.
#' @noRd
validate_flag_taxa_output <- function(
  result,
  original_taxonomy,
  expected_region = NULL,
  allow_llm_tools = FALSE
) {
  if (!inherits(result, "data.frame")) {
    stop("`result` must be a data.frame.", call. = FALSE)
  }
  if (!inherits(original_taxonomy, "data.frame")) {
    stop("`original_taxonomy` must be a data.frame.", call. = FALSE)
  }
  allow_llm_tools <- .validate_logical_scalar(allow_llm_tools, "allow_llm_tools")
  result <- .normalise_flag_taxa_tool_columns(result)
  result <- .normalise_flag_taxa_interpretation_columns(result)
  result <- .normalise_flag_taxa_categorical_columns(
    result,
    expected_region = expected_region
  )

  include_feature_id <- "feature_id" %in% colnames(original_taxonomy)
  required_columns <- .flag_taxa_result_columns(
    result = result,
    include_feature_id = include_feature_id
  )
  missing_columns <- setdiff(required_columns, colnames(result))
  if (length(missing_columns) > 0) {
    stop(
      "`result` is missing required flag_taxa columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(result) != nrow(original_taxonomy)) {
    stop(
      "`result` must have the same number of rows as `original_taxonomy`. ",
      "`result` has ",
      nrow(result),
      " row(s); `original_taxonomy` has ",
      nrow(original_taxonomy),
      " row(s).",
      call. = FALSE
    )
  }

  if (
    "feature_id" %in% colnames(result) &&
      "feature_id" %in% colnames(original_taxonomy)
  ) {
    result_feature_id <- as.character(result$feature_id)
    original_feature_id <- as.character(original_taxonomy$feature_id)
    if (!identical(result_feature_id, original_feature_id)) {
      stop(
        "`feature_id` values in `result` must match `original_taxonomy` ",
        "in the same row order.",
        call. = FALSE
      )
    }
  }

  allowed_values <- .flag_taxa_allowed_values(for_validation = TRUE)
  for (column_name in names(allowed_values)) {
    .validate_flag_taxa_enum_column(
      result = result,
      column_name = column_name,
      allowed_values = allowed_values[[column_name]]
    )
  }
  .validate_flag_taxa_tool_traceability(result)
  .validate_flag_taxa_exclusion_evidence(result, allow_llm_tools)
  .validate_flag_taxa_region_assessment(result, expected_region)

  result
}

.flag_taxa_result_columns <- function(result, include_feature_id) {
  structured_columns <- .flag_taxa_structured_output_columns(
    include_feature_id,
    include_query_id = TRUE
  )
  structured_columns_without_query_id <- .flag_taxa_structured_output_columns(
    include_feature_id,
    include_query_id = FALSE
  )
  legacy_columns <- .flag_taxa_legacy_output_columns(include_feature_id)

  if (all(structured_columns %in% colnames(result))) {
    return(structured_columns)
  }
  if (all(structured_columns_without_query_id %in% colnames(result))) {
    return(structured_columns_without_query_id)
  }
  if (all(legacy_columns %in% colnames(result))) {
    return(legacy_columns)
  }

  legacy_markers <- c(
    "query_rank",
    "query_name",
    "lineage",
    "expected_environment",
    "expected_habitat",
    "expected_region"
  )
  if (any(legacy_markers %in% colnames(result))) {
    return(legacy_columns)
  }

  structured_columns
}

.prepare_flag_taxa_all_feature_map <- function(tax_table, tax_ranks) {
  if (!inherits(tax_table, "data.frame")) {
    stop("`tax_table` must be a data.frame.", call. = FALSE)
  }
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")

  reserved_columns <- .flag_taxa_reserved_input_columns()
  collisions <- colnames(tax_table)[tolower(colnames(tax_table)) %in% reserved_columns]
  if (length(collisions) > 0) {
    stop(
      "`x` already contains columns reserved by `flag_taxa()`: ",
      paste(collisions, collapse = ", "),
      call. = FALSE
    )
  }

  matched_ranks <- tolower(tax_ranks) %in% tolower(colnames(tax_table))
  if (!any(matched_ranks)) {
    stop(
      "`x` must contain at least one requested useful taxonomic rank column. ",
      "Expected one of: ",
      paste(tax_ranks, collapse = ", "),
      call. = FALSE
    )
  }

  query_taxa <- lapply(seq_len(nrow(tax_table)), function(row_index) {
    choose_query_taxon(tax_table[row_index, , drop = FALSE], tax_ranks = tax_ranks)
  })
  feature_id <- .flag_taxa_feature_id(tax_table)
  out <- data.frame(
    original_row = seq_len(nrow(tax_table)),
    feature_id = feature_id,
    query_name = vapply(
      query_taxa,
      function(query_taxon) query_taxon$query_name,
      character(1)
    ),
    query_rank = vapply(
      query_taxa,
      function(query_taxon) query_taxon$query_rank,
      character(1)
    ),
    lineage = vapply(
      seq_len(nrow(tax_table)),
      function(row_index) .build_compact_lineage(tax_table[row_index, , drop = FALSE], tax_ranks),
      character(1)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out$useful_taxonomy <- vapply(out$query_name, .is_informative_taxon, logical(1)) &
    vapply(out$query_rank, .is_informative_taxon, logical(1))
  rownames(out) <- NULL
  out
}

.retained_flag_taxa_feature_query_map <- function(all_feature_map) {
  if (!inherits(all_feature_map, "data.frame")) {
    stop("`all_feature_map` must be a data.frame.", call. = FALSE)
  }
  if (!("useful_taxonomy" %in% colnames(all_feature_map))) {
    stop("`all_feature_map` must contain useful_taxonomy.", call. = FALSE)
  }

  out <- all_feature_map[all_feature_map$useful_taxonomy, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.prepare_flag_taxa_feature_query_map <- function(tax_table, tax_ranks) {
  out <- .retained_flag_taxa_feature_query_map(
    .prepare_flag_taxa_all_feature_map(tax_table, tax_ranks)
  )
  if (nrow(out) == 0) {
    stop(
      "`x` does not contain any useful taxonomic assignments at the requested ranks.",
      call. = FALSE
    )
  }

  out
}

.unique_flag_taxa_query_table <- function(feature_query_map) {
  required_columns <- c("query_name", "query_rank", "lineage")
  if (!all(required_columns %in% colnames(feature_query_map))) {
    stop(
      "`feature_query_map` must contain query_name, query_rank, and lineage.",
      call. = FALSE
    )
  }

  unique_key <- .flag_taxa_query_key(feature_query_map)
  out <- feature_query_map[!duplicated(unique_key), required_columns, drop = FALSE]
  out <- cbind(
    data.frame(
      query_id = .flag_taxa_query_ids(nrow(out)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    out,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

.flag_taxa_query_ids <- function(n) {
  if (n < 1) {
    return(character())
  }
  sprintf("q%04d", seq_len(n))
}

.assign_flag_taxa_query_ids_to_features <- function(feature_query_map, unique_query_table) {
  if (nrow(feature_query_map) == 0) {
    feature_query_map$query_id <- character()
    return(feature_query_map)
  }
  match_index <- match(
    .flag_taxa_query_key(feature_query_map),
    .flag_taxa_query_key(unique_query_table)
  )
  if (any(is.na(match_index))) {
    stop(
      "Internal error: could not assign query IDs to feature rows.",
      call. = FALSE
    )
  }
  feature_query_map$query_id <- unique_query_table$query_id[match_index]
  feature_query_map
}

.build_flag_taxa_preliminary_judgement <- function(
  unique_query_table,
  taxon_evidence,
  expected_environment,
  expected_habitat,
  expected_region,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum")
) {
  if (nrow(unique_query_table) == 0) {
    return(.empty_flag_taxa_preliminary_judgement())
  }

  rows <- lapply(seq_len(nrow(unique_query_table)), function(row_index) {
    query_row <- unique_query_table[row_index, , drop = FALSE]
    evidence_rows <- .raw_flag_taxa_relevant_evidence(
      query_table = unique_query_table,
      taxon_evidence = taxon_evidence,
      row_index = row_index
    )
    .classify_flag_taxa_preliminary_row(
      query_row = query_row,
      evidence_rows = evidence_rows,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region
    )
  })

  out <- .bind_rows_fill(rows)
  rownames(out) <- NULL
  out
}

.empty_flag_taxa_preliminary_judgement <- function() {
  data.frame(
    query_id = character(),
    query_name = character(),
    query_rank = character(),
    lineage = character(),
    prelim_environment_status = character(),
    prelim_habitat_status = character(),
    prelim_region_status = character(),
    prelim_ecological_status = character(),
    prelim_occurrence_interpretation = character(),
    prelim_recommended_action = character(),
    prelim_rationale = character(),
    prelim_evidence_basis = character(),
    prelim_references = character(),
    has_matched_local_evidence = logical(),
    has_informative_local_evidence = logical(),
    needs_llm = logical(),
    needs_tools = logical(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.raw_flag_taxa_relevant_evidence <- function(query_table, taxon_evidence, row_index) {
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  if (is.null(taxon_evidence) || nrow(taxon_evidence) == 0) {
    return(NULL)
  }

  matched <- .taxon_evidence_match_indices(
    tax_table = query_table,
    taxon_evidence = taxon_evidence,
    row_index = row_index
  )
  if (length(matched) == 0) {
    return(NULL)
  }

  taxon_evidence[matched, , drop = FALSE]
}

.classify_flag_taxa_preliminary_row <- function(
  query_row,
  evidence_rows,
  expected_environment,
  expected_habitat,
  expected_region
) {
  query_rank <- as.character(query_row$query_rank[[1]])
  query_name <- as.character(query_row$query_name[[1]])
  is_species <- identical(tolower(query_rank), "species")
  informative_evidence <- .flag_taxa_informative_evidence_rows(evidence_rows)
  has_matched <- !is.null(evidence_rows) && nrow(evidence_rows) > 0
  has_informative <- !is.null(informative_evidence) && nrow(informative_evidence) > 0

  if (!has_informative) {
    env_status <- "unknown"
    habitat_status <- if (is.null(expected_habitat)) "unknown" else "no_known_habitat_evidence"
    region_status <- if (is.null(expected_region)) "not_assessed" else "no_distribution_evidence"
    interpretation <- .flag_taxa_prelim_interpretation(
      env_status = env_status,
      habitat_status = habitat_status,
      region_status = region_status,
      query_rank = query_rank,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      evidence_rows = evidence_rows,
      is_species = is_species,
      has_positive_environment_incompatibility = FALSE
    )
    action <- interpretation$recommended_action
    basis <- "conservative_reasoning_only"
    rationale <- paste0(
      "No informative matched local evidence was available for ",
      query_name,
      ". Missing evidence, failed lookups, no matches, or absent records are treated as uncertainty rather than incompatibility."
    )
    references <- NA_character_
  } else {
    env_values <- .flag_taxa_evidence_values(informative_evidence, "environment")
    habitat_values <- .flag_taxa_evidence_values(informative_evidence, "habitat")
    region_values <- .flag_taxa_evidence_values(informative_evidence, "region")
    supports_environment <- .flag_taxa_values_contain_context(env_values, expected_environment)
    contradicts_environment <- .flag_taxa_environment_contradicts_expected(
      env_values = env_values,
      expected_environment = expected_environment
    )

    env_status <- .flag_taxa_prelim_environment_status(
      env_values = env_values,
      expected_environment = expected_environment,
      supports_environment = supports_environment,
      contradicts_environment = contradicts_environment,
      is_species = is_species
    )
    habitat_status <- .flag_taxa_prelim_habitat_status(
      habitat_values = habitat_values,
      expected_habitat = expected_habitat,
      env_status = env_status,
      supports_environment = supports_environment
    )
    if (.flag_taxa_photosynthetic_aphotic_case(
      ecology_text = .flag_taxa_evidence_ecology_text(informative_evidence),
      expected_habitat = expected_habitat
    )) {
      habitat_status <- "incompatible"
    }
    region_status <- .flag_taxa_prelim_region_status(
      region_values = region_values,
      expected_region = expected_region
    )
    action <- .flag_taxa_prelim_recommended_action(
      env_status = env_status,
      habitat_status = habitat_status,
      region_status = region_status,
      is_species = is_species,
      has_positive_environment_incompatibility = isTRUE(contradicts_environment),
      evidence_rows = informative_evidence,
      query_rank = query_rank,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat
    )
    interpretation <- .flag_taxa_prelim_interpretation(
      env_status = env_status,
      habitat_status = habitat_status,
      region_status = region_status,
      query_rank = query_rank,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      evidence_rows = informative_evidence,
      is_species = is_species,
      has_positive_environment_incompatibility = isTRUE(contradicts_environment)
    )
    action <- interpretation$recommended_action
    basis <- "local_evidence"
    rationale <- .flag_taxa_prelim_rationale(
      query_name = query_name,
      query_rank = query_rank,
      env_status = env_status,
      habitat_status = habitat_status,
      region_status = region_status,
      action = action,
      env_values = env_values,
      habitat_values = habitat_values,
      region_values = region_values,
      is_species = is_species,
      evidence_rows = informative_evidence,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      ecological_status = interpretation$ecological_status,
      occurrence_interpretation = interpretation$occurrence_interpretation
    )
    references <- .collapse_distinct_values(.optional_evidence_column(informative_evidence, "reference"), max_values = 6)
  }

  data.frame(
    query_id = as.character(query_row$query_id[[1]]),
    query_name = query_name,
    query_rank = query_rank,
    lineage = as.character(query_row$lineage[[1]]),
    prelim_environment_status = env_status,
    prelim_habitat_status = habitat_status,
    prelim_region_status = region_status,
    prelim_ecological_status = interpretation$ecological_status,
    prelim_occurrence_interpretation = interpretation$occurrence_interpretation,
    prelim_recommended_action = action,
    prelim_rationale = rationale,
    prelim_evidence_basis = basis,
    prelim_references = references,
    has_matched_local_evidence = has_matched,
    has_informative_local_evidence = has_informative,
    needs_llm = action %in% c("flag_for_review", "flag_possible_exclusion"),
    needs_tools = action %in% c("flag_for_review", "flag_possible_exclusion"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.flag_taxa_informative_evidence_rows <- function(evidence_rows) {
  if (is.null(evidence_rows) || nrow(evidence_rows) == 0) {
    return(NULL)
  }

  evidence_type <- tolower(as.character(.optional_evidence_column(evidence_rows, "evidence_type")))
  summary <- tolower(as.character(.optional_evidence_column(evidence_rows, "evidence_summary")))
  non_informative_type <- evidence_type %in% c("not_queried", "taxonomy_lookup_failed")
  no_match_summary <- grepl(
    "no worms record found|no record found|no match|lookup failed|not queried|no_evidence_found",
    summary
  )
  has_context <- .row_has_any_non_empty_evidence_value(
    evidence_rows,
    c("environment", "habitat", "region")
  )
  keep <- has_context & !(non_informative_type | no_match_summary)
  out <- evidence_rows[keep, , drop = FALSE]
  if (nrow(out) == 0) {
    return(NULL)
  }
  out
}

.row_has_any_non_empty_evidence_value <- function(evidence_rows, columns) {
  available <- intersect(columns, colnames(evidence_rows))
  if (length(available) == 0) {
    return(rep(FALSE, nrow(evidence_rows)))
  }
  values <- evidence_rows[available]
  Reduce(`|`, lapply(values, .is_non_empty_value))
}

.flag_taxa_evidence_values <- function(evidence_rows, column_name) {
  if (is.null(evidence_rows) || !(column_name %in% colnames(evidence_rows))) {
    return(character())
  }
  values <- as.character(evidence_rows[[column_name]])
  values <- values[.is_non_empty_value(values)]
  if (length(values) == 0) {
    return(character())
  }
  split_values <- unlist(strsplit(values, "\\s*;\\s*|\\s*,\\s*|\\s*/\\s*"), use.names = FALSE)
  split_values <- trimws(split_values)
  unique(tolower(split_values[.is_non_empty_value(split_values)]))
}

.flag_taxa_values_contain_context <- function(values, context) {
  if (is.null(context) || !.is_non_empty_value(context) || length(values) == 0) {
    return(FALSE)
  }
  context <- tolower(trimws(as.character(context)))
  any(
    values == context |
      grepl(context, values, fixed = TRUE) |
      vapply(values, function(value) grepl(value, context, fixed = TRUE), logical(1))
  )
}

.flag_taxa_environment_contradicts_expected <- function(env_values, expected_environment) {
  if (length(env_values) == 0 || is.null(expected_environment)) {
    return(FALSE)
  }
  expected <- tolower(trimws(expected_environment))
  supports <- .flag_taxa_values_contain_context(env_values, expected)
  if (supports) {
    return(FALSE)
  }

  incompatible_sets <- list(
    marine = c("freshwater", "terrestrial"),
    freshwater = c("marine", "terrestrial"),
    terrestrial = c("marine", "freshwater", "brackish")
  )
  if (!(expected %in% names(incompatible_sets))) {
    return(FALSE)
  }
  any(env_values %in% incompatible_sets[[expected]])
}

.flag_taxa_prelim_environment_status <- function(
  env_values,
  expected_environment,
  supports_environment,
  contradicts_environment,
  is_species
) {
  if (length(env_values) == 0) {
    return("unknown")
  }
  if (isTRUE(supports_environment)) {
    if (length(env_values) > 1 && !isTRUE(is_species)) {
      return("mixed_within_rank")
    }
    return("compatible")
  }
  if (isTRUE(contradicts_environment)) {
    if (isTRUE(is_species)) {
      return("incompatible")
    }
    return("mixed_within_rank")
  }
  "unknown"
}

.flag_taxa_prelim_habitat_status <- function(
  habitat_values,
  expected_habitat,
  env_status,
  supports_environment
) {
  if (is.null(expected_habitat)) {
    return("unknown")
  }
  if (.flag_taxa_values_contain_context(habitat_values, expected_habitat)) {
    return("compatible")
  }
  if (identical(env_status, "incompatible")) {
    return("unknown")
  }
  if (isTRUE(supports_environment) || env_status %in% c("compatible", "mixed_within_rank")) {
    return("no_known_habitat_evidence")
  }
  "unknown"
}

.flag_taxa_prelim_region_status <- function(region_values, expected_region) {
  if (is.null(expected_region)) {
    return("not_assessed")
  }
  if (.flag_taxa_values_contain_context(region_values, expected_region)) {
    return("known_in_region")
  }
  if (length(region_values) > 0) {
    return("known_elsewhere_only")
  }
  "no_distribution_evidence"
}

.flag_taxa_prelim_recommended_action <- function(
  env_status,
  habitat_status,
  region_status,
  is_species,
  has_positive_environment_incompatibility,
  evidence_rows = NULL,
  query_rank = NA_character_,
  expected_environment = NULL,
  expected_habitat = NULL
) {
  interpretation <- .flag_taxa_prelim_interpretation(
    env_status = env_status,
    habitat_status = habitat_status,
    region_status = region_status,
    query_rank = query_rank,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    evidence_rows = evidence_rows,
    is_species = is_species,
    has_positive_environment_incompatibility = has_positive_environment_incompatibility
  )
  interpretation$recommended_action
}

.flag_taxa_prelim_interpretation <- function(
  env_status,
  habitat_status,
  region_status,
  query_rank,
  expected_environment,
  expected_habitat,
  evidence_rows = NULL,
  is_species = FALSE,
  has_positive_environment_incompatibility = FALSE
) {
  ecology_text <- .flag_taxa_evidence_ecology_text(evidence_rows)
  photosynthetic_deep_sea <- .flag_taxa_photosynthetic_aphotic_case(
    ecology_text = ecology_text,
    expected_habitat = expected_habitat
  )
  host_associated_mismatch <- .flag_taxa_host_associated_habitat_mismatch_case(
    ecology_text = ecology_text,
    expected_habitat = expected_habitat
  )
  pelagic_benthic_mismatch <- .flag_taxa_pelagic_benthic_mismatch_case(
    ecology_text = ecology_text,
    expected_habitat = expected_habitat
  )
  questionable_marine_broad <- .flag_taxa_questionable_marine_broad_case(
    ecology_text = ecology_text,
    expected_environment = expected_environment,
    is_species = is_species
  )

  if (identical(env_status, "insufficient_taxonomic_resolution") ||
      identical(habitat_status, "insufficient_taxonomic_resolution")) {
    return(list(
      ecological_status = "insufficient_taxonomic_resolution",
      occurrence_interpretation = "not_assessed",
      recommended_action = "flag_for_review"
    ))
  }

  if (
    isTRUE(is_species) &&
      identical(env_status, "incompatible") &&
      isTRUE(has_positive_environment_incompatibility)
  ) {
    return(list(
      ecological_status = "incompatible_resident",
      occurrence_interpretation = "likely_contaminant_or_misassignment",
      recommended_action = "exclude"
    ))
  }

  if (isTRUE(photosynthetic_deep_sea)) {
    return(list(
      ecological_status = "incompatible_resident",
      occurrence_interpretation = "possible_transient_or_allochthonous",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (isTRUE(host_associated_mismatch)) {
    return(list(
      ecological_status = "unlikely_resident",
      occurrence_interpretation = "possible_transient_or_allochthonous",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (isTRUE(pelagic_benthic_mismatch)) {
    return(list(
      ecological_status = "unlikely_resident",
      occurrence_interpretation = "possible_transient_or_allochthonous",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (isTRUE(questionable_marine_broad)) {
    return(list(
      ecological_status = "unlikely_resident",
      occurrence_interpretation = "possible_contaminant_or_misassignment",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (
    !isTRUE(has_positive_environment_incompatibility) &&
      habitat_status == "compatible" &&
      region_status %in% c("known_in_region", "known_near_region") &&
      env_status %in% c("compatible", "mixed_within_rank", "possible_likely")
  ) {
    return(list(
      ecological_status = "compatible",
      occurrence_interpretation = "expected_resident",
      recommended_action = "retain"
    ))
  }

  if (identical(habitat_status, "incompatible")) {
    return(list(
      ecological_status = "incompatible_resident",
      occurrence_interpretation = "possible_transient_or_allochthonous",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  if (
    identical(env_status, "compatible") &&
      habitat_status %in% c("compatible", "possible_likely", "unknown") &&
      region_status %in% c("known_in_region", "known_near_region", "not_assessed", "unknown")
  ) {
    return(list(
      ecological_status = "compatible",
      occurrence_interpretation = if (identical(habitat_status, "compatible")) "expected_resident" else "plausible_resident",
      recommended_action = "retain"
    ))
  }

  if (
    env_status %in% c("compatible", "mixed_within_rank", "possible_likely") &&
      habitat_status %in% c("no_known_habitat_evidence", "possible_likely", "unknown")
  ) {
    return(list(
      ecological_status = "plausible",
      occurrence_interpretation = "plausible_resident",
      recommended_action = "flag_for_review"
    ))
  }

  if (identical(env_status, "incompatible")) {
    return(list(
      ecological_status = "incompatible_resident",
      occurrence_interpretation = "possible_contaminant_or_misassignment",
      recommended_action = "flag_possible_exclusion"
    ))
  }

  list(
    ecological_status = "unknown",
    occurrence_interpretation = "uncertain",
    recommended_action = "flag_for_review"
  )
}

.flag_taxa_evidence_ecology_text <- function(evidence_rows) {
  if (is.null(evidence_rows) || nrow(evidence_rows) == 0) {
    return("")
  }
  columns <- intersect(
    c("taxon_name", "taxon_rank", "environment", "habitat", "region", "evidence_summary", "reference"),
    colnames(evidence_rows)
  )
  text <- unlist(evidence_rows[columns], use.names = FALSE)
  tolower(paste(as.character(text), collapse = " "))
}

.flag_taxa_photosynthetic_aphotic_case <- function(ecology_text, expected_habitat) {
  if (is.null(expected_habitat) || !.is_non_empty_value(expected_habitat)) {
    return(FALSE)
  }
  habitat <- tolower(as.character(expected_habitat))
  aphotic <- grepl("deep sea|deep-sea|aphotic|hadal|abyssal|seafloor|benthic", habitat)
  photosynthetic <- grepl(
    "photosynthetic|phototroph|alga|algal|red algal|florideophyceae|chlorophy|diatom|cyanobacter",
    ecology_text
  )
  isTRUE(aphotic && photosynthetic)
}

.flag_taxa_host_associated_habitat_mismatch_case <- function(ecology_text, expected_habitat) {
  if (is.null(expected_habitat) || !.is_non_empty_value(expected_habitat)) {
    return(FALSE)
  }
  habitat <- tolower(as.character(expected_habitat))
  resident_habitat <- grepl("deep sea|deep-sea|benthic|sediment|seafloor|abyssal|hadal", habitat)
  host_associated <- grepl(
    "host[ -]?associated|parasite|parasitic|pathogen|pathogenic|obligate host|host dependent|host-dependent",
    ecology_text
  )
  isTRUE(resident_habitat && host_associated)
}

.flag_taxa_pelagic_benthic_mismatch_case <- function(ecology_text, expected_habitat) {
  if (is.null(expected_habitat) || !.is_non_empty_value(expected_habitat)) {
    return(FALSE)
  }
  habitat <- tolower(as.character(expected_habitat))
  benthic <- grepl("benthic|sediment|seafloor|bottom|abyssal|hadal", habitat)
  pelagic <- grepl("pelagic|planktonic|plankton|water column", ecology_text)
  isTRUE(benthic && pelagic)
}

.flag_taxa_questionable_marine_broad_case <- function(ecology_text, expected_environment, is_species) {
  if (isTRUE(is_species) || is.null(expected_environment)) {
    return(FALSE)
  }
  expected <- tolower(as.character(expected_environment))
  if (!identical(expected, "marine")) {
    return(FALSE)
  }
  has_nonmarine <- grepl("freshwater|terrestrial|soil", ecology_text)
  has_questionable <- grepl("questionable|mostly|primarily|true .*freshwater|marine records? may be|marine records? questionable", ecology_text)
  has_nonmarine && has_questionable
}

.flag_taxa_prelim_rationale <- function(
  query_name,
  query_rank,
  env_status,
  habitat_status,
  region_status,
  action,
  env_values,
  habitat_values,
  region_values,
  is_species,
  evidence_rows = NULL,
  expected_environment = NULL,
  expected_habitat = NULL,
  ecological_status = NA_character_,
  occurrence_interpretation = NA_character_
) {
  ecology_text <- .flag_taxa_evidence_ecology_text(evidence_rows)
  pieces <- c(
    paste0("Matched local evidence for ", query_name, " at rank ", query_rank, "."),
    paste0("Environment status: ", env_status, "."),
    paste0("Habitat status: ", habitat_status, "."),
    paste0("Region status: ", region_status, ".")
  )
  if (length(env_values) > 0) {
    pieces <- c(pieces, paste0("Environment evidence: ", paste(env_values, collapse = "; "), "."))
  }
  if (length(habitat_values) > 0) {
    pieces <- c(pieces, paste0("Habitat evidence: ", paste(habitat_values, collapse = "; "), "."))
  }
  if (length(region_values) > 0) {
    pieces <- c(pieces, paste0("Region evidence: ", paste(region_values, collapse = "; "), "."))
  }
  if (!isTRUE(is_species) && identical(action, "flag_for_review")) {
    pieces <- c(pieces, "Because the query rank is broad, missing habitat or region evidence is treated cautiously and does not justify exclusion.")
  }
  if (identical(habitat_status, "no_known_habitat_evidence")) {
    pieces <- c(pieces, "No direct expected-habitat evidence was matched; this is uncertainty, not habitat incompatibility.")
  }
  if (.flag_taxa_photosynthetic_aphotic_case(ecology_text, expected_habitat)) {
    pieces <- c(
      pieces,
      "Biological reason: local ecology evidence describes photosynthetic or algal biology, while the expected habitat is aphotic/deep-sea/benthic; resident photosynthesis is not plausible there, although transported or allochthonous DNA remains possible."
    )
  }
  if (.flag_taxa_host_associated_habitat_mismatch_case(ecology_text, expected_habitat)) {
    pieces <- c(
      pieces,
      "Biological reason: local ecology evidence indicates host-associated, parasitic, or pathogenic biology, which makes interpretation as a free-living resident of the expected habitat uncertain or unlikely without explicit host/resident evidence."
    )
  }
  if (.flag_taxa_pelagic_benthic_mismatch_case(ecology_text, expected_habitat)) {
    pieces <- c(
      pieces,
      "Biological reason: local ecology evidence indicates pelagic or planktonic ecology, so detection in a benthic or sediment habitat may reflect transient or allochthonous material rather than a resident benthic organism."
    )
  }
  if (.flag_taxa_questionable_marine_broad_case(ecology_text, expected_environment, is_species)) {
    pieces <- c(
      pieces,
      "Biological reason: local evidence indicates predominantly non-marine ecology or questionable marine records for this broad taxon; marine flags prevent a hard exclusion, but resident marine/deep-sea interpretation remains a positive ecological concern."
    )
  }
  if (
    as.character(ecological_status) %in% c("unlikely_resident", "incompatible_resident") ||
      as.character(occurrence_interpretation) %in% c(
        "possible_transient_or_allochthonous",
        "possible_contaminant_or_misassignment",
        "likely_contaminant_or_misassignment"
      ) ||
      action %in% c("flag_possible_exclusion", "exclude")
  ) {
    pieces <- c(
      pieces,
      paste0(
        "Interpretation reason: ecological_status = ",
        ecological_status,
        " and occurrence_interpretation = ",
        occurrence_interpretation,
        " are based on positive ecological context above, not on missing records alone."
      )
    )
  }
  if (identical(action, "exclude")) {
    pieces <- c(pieces, "Exclusion is used only because species-level local evidence positively contradicts the expected environment.")
  }

  paste(pieces, collapse = " ")
}

.flag_taxa_taxon_result_from_preliminary <- function(preliminary_judgement) {
  data.frame(
    query_id = as.character(preliminary_judgement$query_id),
    taxon_name = as.character(preliminary_judgement$query_name),
    taxon_rank = as.character(preliminary_judgement$query_rank),
    expected_environment_status = as.character(preliminary_judgement$prelim_environment_status),
    expected_habitat_status = as.character(preliminary_judgement$prelim_habitat_status),
    expected_region_status = as.character(preliminary_judgement$prelim_region_status),
    ecological_status = as.character(preliminary_judgement$prelim_ecological_status),
    occurrence_interpretation = as.character(preliminary_judgement$prelim_occurrence_interpretation),
    recommended_action = as.character(preliminary_judgement$prelim_recommended_action),
    rationale = as.character(preliminary_judgement$prelim_rationale),
    references = as.character(preliminary_judgement$prelim_references),
    evidence_basis = as.character(preliminary_judgement$prelim_evidence_basis),
    llm_tool_used = rep(FALSE, nrow(preliminary_judgement)),
    llm_tool_name = rep(NA_character_, nrow(preliminary_judgement)),
    llm_tool_query = rep(NA_character_, nrow(preliminary_judgement)),
    llm_tool_evidence_summary = rep(NA_character_, nrow(preliminary_judgement)),
    llm_tool_references = rep(NA_character_, nrow(preliminary_judgement)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.select_flag_taxa_llm_query_table <- function(
  unique_query_table,
  preliminary_judgement,
  judgement_mode
) {
  if (nrow(unique_query_table) == 0 || identical(judgement_mode, "evidence_only")) {
    return(unique_query_table[0, , drop = FALSE])
  }

  selected <- .flag_taxa_judgement_mode_selection(
    preliminary_judgement = preliminary_judgement,
    judgement_mode = judgement_mode
  )
  out <- unique_query_table[selected, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.flag_taxa_judgement_mode_selection <- function(preliminary_judgement, judgement_mode) {
  if (is.null(preliminary_judgement) || nrow(preliminary_judgement) == 0) {
    return(logical())
  }

  action <- as.character(preliminary_judgement$prelim_recommended_action)
  has_informative <- as.logical(preliminary_judgement$has_informative_local_evidence)
  has_informative[is.na(has_informative)] <- FALSE

  switch(
    judgement_mode,
    evidence_only = rep(FALSE, nrow(preliminary_judgement)),
    llm_missing_evidence = !has_informative,
    llm_flagged = action %in% c("flag_for_review", "flag_possible_exclusion"),
    llm_possible_exclusion = action == "flag_possible_exclusion",
    llm_all = rep(TRUE, nrow(preliminary_judgement)),
    stop("Unsupported judgement_mode.", call. = FALSE)
  )
}

.merge_flag_taxa_llm_result_with_preliminary <- function(
  llm_taxon_result,
  preliminary_judgement
) {
  out <- .flag_taxa_taxon_result_from_preliminary(preliminary_judgement)
  if (nrow(llm_taxon_result) == 0) {
    return(out)
  }
  match_index <- match(as.character(llm_taxon_result$query_id), as.character(out$query_id))
  if (any(is.na(match_index))) {
    stop("Internal error: LLM result query IDs could not be merged with preliminary judgements.", call. = FALSE)
  }
  out[match_index, colnames(llm_taxon_result)] <- llm_taxon_result
  out
}

.resolve_flag_taxa_supplied_result_query_table <- function(
  result,
  unique_query_table,
  llm_query_table
) {
  if (!inherits(result, "data.frame") || "feature_id" %in% colnames(result)) {
    return(unique_query_table)
  }
  if ("query_id" %in% colnames(result)) {
    result_ids <- as.character(result$query_id)
    if (setequal(result_ids, as.character(unique_query_table$query_id))) {
      return(unique_query_table)
    }
    if (setequal(result_ids, as.character(llm_query_table$query_id))) {
      return(llm_query_table)
    }
  } else if (nrow(result) == nrow(unique_query_table)) {
    return(unique_query_table)
  } else if (nrow(result) == nrow(llm_query_table)) {
    return(llm_query_table)
  }
  llm_query_table
}

.process_flag_taxa_llm_result <- function(
  result,
  feature_query_map,
  unique_query_table,
  llm_query_table = unique_query_table,
  preliminary_judgement = NULL,
  expected_environment,
  expected_habitat,
  expected_region,
  taxon_evidence = NULL,
  tool_evidence = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  max_evidence_rows_per_query = 5,
  max_evidence_summary_chars = 160,
  allow_llm_tools = FALSE,
  tool_requirement = "optional",
  prompt_tools = NULL,
  repair_invalid_llm = TRUE,
  verbose = FALSE,
  start_time = NULL
) {
  repair_invalid_llm <- .validate_logical_scalar(repair_invalid_llm, "repair_invalid_llm")
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  if (
    inherits(result, "data.frame") &&
      "feature_id" %in% colnames(result)
  ) {
    result <- .repair_flag_taxa_required_tool_attempts(
      result = result,
      tool_requirement = tool_requirement,
      allow_llm_tools = allow_llm_tools,
      repair_invalid_llm = repair_invalid_llm,
      query_taxonomy = feature_query_map,
      prompt_tools = prompt_tools,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      verbose = verbose,
      start_time = start_time
    )
    result <- .repair_flag_taxa_invalid_llm_output(
      result = result,
      repair_invalid_llm = repair_invalid_llm,
      expected_region = expected_region,
      verbose = verbose,
      start_time = start_time
    )
    feature_result <- validate_flag_taxa_output(
      result = result,
      original_taxonomy = feature_query_map,
      expected_region = expected_region,
      allow_llm_tools = allow_llm_tools
    )
    .validate_feature_level_taxon_identity(feature_result, feature_query_map)
    return(.enrich_flag_taxa_final_output(
      result = feature_result,
      feature_query_map = feature_query_map,
      unique_query_table = unique_query_table,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      taxon_evidence = taxon_evidence,
      tool_evidence = tool_evidence,
      tax_ranks = tax_ranks,
      max_evidence_rows_per_query = max_evidence_rows_per_query,
      max_evidence_summary_chars = max_evidence_summary_chars
    ))
  }

  result <- .augment_flag_taxa_result_with_tool_evidence(
    result = result,
    query_taxonomy = llm_query_table,
    taxon_evidence = taxon_evidence,
    tool_evidence = tool_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  result <- .repair_flag_taxa_required_tool_attempts(
    result = result,
    tool_requirement = tool_requirement,
    allow_llm_tools = allow_llm_tools,
    repair_invalid_llm = repair_invalid_llm,
    query_taxonomy = llm_query_table,
    prompt_tools = prompt_tools,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    verbose = verbose,
    start_time = start_time
  )
  result <- .repair_flag_taxa_invalid_llm_output(
    result = result,
    repair_invalid_llm = repair_invalid_llm,
    expected_region = expected_region,
    verbose = verbose,
    start_time = start_time
  )
  taxon_result <- validate_flag_taxa_taxon_output(
    result = result,
    query_taxonomy = llm_query_table,
    expected_region = expected_region,
    allow_llm_tools = allow_llm_tools
  )
  if (!is.null(preliminary_judgement)) {
    taxon_result <- .merge_flag_taxa_llm_result_with_preliminary(
      llm_taxon_result = taxon_result,
      preliminary_judgement = preliminary_judgement
    )
  }
  feature_result <- .join_flag_taxa_taxon_result_to_features(
    taxon_result = taxon_result,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table
  )
  .enrich_flag_taxa_final_output(
    result = feature_result,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    taxon_evidence = taxon_evidence,
    tool_evidence = tool_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
}

validate_flag_taxa_taxon_output <- function(
  result,
  query_taxonomy,
  expected_region = NULL,
  allow_llm_tools = FALSE
) {
  if (!inherits(result, "data.frame")) {
    stop("`result` must be a data.frame.", call. = FALSE)
  }
  if (!inherits(query_taxonomy, "data.frame")) {
    stop("`query_taxonomy` must be a data.frame.", call. = FALSE)
  }
  allow_llm_tools <- .validate_logical_scalar(allow_llm_tools, "allow_llm_tools")
  result <- .normalise_flag_taxa_tool_columns(result)
  result <- .normalise_flag_taxa_interpretation_columns(result)
  result <- .normalise_flag_taxa_categorical_columns(
    result,
    expected_region = expected_region
  )

  required_columns <- .flag_taxa_taxon_result_columns(result)
  missing_columns <- setdiff(required_columns, colnames(result))
  if (length(missing_columns) > 0) {
    stop(
      "`result` is missing required flag_taxa columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(result) != nrow(query_taxonomy)) {
    stop(
      "`result` must have one row per unique query taxon. `result` has ",
      nrow(result),
      " row(s); the unique query table has ",
      nrow(query_taxonomy),
      " row(s).",
      call. = FALSE
    )
  }

  if ("query_id" %in% colnames(result)) {
    .validate_flag_taxa_result_query_ids(result$query_id, query_taxonomy$query_id)
    order_index <- match(query_taxonomy$query_id, as.character(result$query_id))
    result <- result[order_index, , drop = FALSE]
  } else {
    expected_key <- .flag_taxa_result_query_key(
      query_taxonomy$query_name,
      query_taxonomy$query_rank
    )
    result_key <- .flag_taxa_result_query_key(
      result$taxon_name,
      result$taxon_rank
    )
    if (!identical(result_key, expected_key)) {
      stop(
        "`result` taxon_name/taxon_rank values must match the unique query taxa ",
        "in the same row order.",
        call. = FALSE
      )
    }
    result$query_id <- query_taxonomy$query_id
  }

  allowed_values <- .flag_taxa_allowed_values(for_validation = TRUE)
  for (column_name in names(allowed_values)) {
    .validate_flag_taxa_enum_column(
      result = result,
      column_name = column_name,
      allowed_values = allowed_values[[column_name]]
    )
  }
  .validate_flag_taxa_tool_traceability(result)
  .validate_flag_taxa_exclusion_evidence(result, allow_llm_tools)
  .validate_flag_taxa_region_assessment(result, expected_region)

  result[, .flag_taxa_structured_output_columns(
    include_feature_id = FALSE,
    include_query_id = TRUE
  ), drop = FALSE]
}

.join_flag_taxa_taxon_result_to_features <- function(
  taxon_result,
  feature_query_map,
  unique_query_table
) {
  if ("query_id" %in% colnames(taxon_result) && "query_id" %in% colnames(feature_query_map)) {
    match_index <- match(
      as.character(feature_query_map$query_id),
      as.character(taxon_result$query_id)
    )
  } else {
    query_key <- .flag_taxa_query_key(unique_query_table)
    feature_key <- .flag_taxa_query_key(feature_query_map)
    match_index <- match(feature_key, query_key)
  }
  if (any(is.na(match_index))) {
    stop(
      "Could not join taxon-level flag_taxa results back to all feature rows.",
      call. = FALSE
    )
  }

  out <- cbind(
    data.frame(
      feature_id = as.character(feature_query_map$feature_id),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    taxon_result[match_index, , drop = FALSE]
  )
  rownames(out) <- NULL
  out
}

.build_flag_taxa_evidence_first_output <- function(
  preliminary_judgement,
  feature_query_map,
  unique_query_table,
  selected_query_table,
  all_feature_map,
  expected_environment,
  expected_habitat,
  expected_region,
  taxon_evidence,
  tool_evidence,
  tax_ranks,
  max_evidence_rows_per_query,
  max_evidence_summary_chars
) {
  preliminary_result <- .flag_taxa_taxon_result_from_preliminary(preliminary_judgement)
  feature_result <- .join_flag_taxa_taxon_result_to_features(
    taxon_result = preliminary_result,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table
  )
  assessed <- .enrich_flag_taxa_final_output(
    result = feature_result,
    feature_query_map = feature_query_map,
    unique_query_table = unique_query_table,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    taxon_evidence = taxon_evidence,
    tool_evidence = tool_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  selected_ids <- if (is.null(selected_query_table) || nrow(selected_query_table) == 0) {
    character()
  } else {
    as.character(selected_query_table$query_id)
  }
  if ("query_id" %in% colnames(assessed)) {
    assessed$llm_selected <- as.character(assessed$query_id) %in% selected_ids
    assessed$llm_prompt_selected <- assessed$llm_selected
  }
  .restore_unassessed_flag_taxa_rows(
    assessed_output = assessed,
    all_feature_map = all_feature_map,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region
  )
}

.restore_unassessed_flag_taxa_rows <- function(
  assessed_output,
  all_feature_map,
  expected_environment,
  expected_habitat,
  expected_region
) {
  unassessed <- .build_unassessed_flag_taxa_rows(
    all_feature_map = all_feature_map,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region
  )

  rows <- Filter(Negate(is.null), list(assessed_output, unassessed))
  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- .bind_rows_fill(rows)
  out <- out[order(out$.flag_taxa_original_row), , drop = FALSE]
  out$.flag_taxa_original_row <- NULL
  rownames(out) <- NULL
  .finalise_flag_taxa_output(out, expected_region = expected_region)
}

.build_unassessed_flag_taxa_rows <- function(
  all_feature_map,
  expected_environment,
  expected_habitat,
  expected_region
) {
  if (nrow(all_feature_map) == 0) {
    return(NULL)
  }

  unassessed <- all_feature_map[!all_feature_map$useful_taxonomy, , drop = FALSE]
  if (nrow(unassessed) == 0) {
    return(NULL)
  }

  data.frame(
    .flag_taxa_original_row = as.integer(unassessed$original_row),
    feature_id = as.character(unassessed$feature_id),
    query_id = NA_character_,
    taxon_name = NA_character_,
    taxon_rank = NA_character_,
    lineage = as.character(unassessed$lineage),
    worms_environment = NA_character_,
    expected_environment = rep(expected_environment, nrow(unassessed)),
    expected_habitat = rep(.nullable_character(expected_habitat), nrow(unassessed)),
    expected_region = rep(.nullable_character(expected_region), nrow(unassessed)),
    expected_environment_status = "insufficient_taxonomic_resolution",
    expected_habitat_status = "insufficient_taxonomic_resolution",
    expected_region_status = if (is.null(expected_region)) "not_assessed" else "unknown",
    ecological_status = "insufficient_taxonomic_resolution",
    occurrence_interpretation = "not_assessed",
    recommended_action = "flag_for_review",
    rationale = "No useful taxonomic assignment was available, so this feature was not assessed by the LLM.",
    references = NA_character_,
    evidence_basis = "conservative_reasoning_only",
    llm_tool_used = FALSE,
    llm_tool_name = NA_character_,
    llm_tool_query = NA_character_,
    llm_tool_evidence_summary = NA_character_,
    llm_tool_references = NA_character_,
    evidence_sources = "none",
    evidence_summary = "No matched local evidence because no useful taxonomic assignment was available.",
    llm_prompt_selected = FALSE,
    llm_selected = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.validate_feature_level_taxon_identity <- function(result, feature_query_map) {
  if ("taxon_name" %in% colnames(result)) {
    result_name <- result$taxon_name
  } else if ("query_name" %in% colnames(result)) {
    result_name <- result$query_name
  } else {
    result_name <- NULL
  }

  if ("taxon_rank" %in% colnames(result)) {
    result_rank <- result$taxon_rank
  } else if ("query_rank" %in% colnames(result)) {
    result_rank <- result$query_rank
  } else {
    result_rank <- NULL
  }

  if (!is.null(result_name) && !is.null(result_rank)) {
    result_key <- .flag_taxa_result_query_key(result_name, result_rank)
    expected_key <- .flag_taxa_result_query_key(
      feature_query_map$query_name,
      feature_query_map$query_rank
    )
    if (!identical(result_key, expected_key)) {
      stop(
        "`result` taxon_name/taxon_rank values must match the feature-level ",
        "query taxa in the same row order.",
        call. = FALSE
      )
    }
  }

  invisible(result)
}

.validate_flag_taxa_result_query_ids <- function(result_query_id, expected_query_id) {
  result_query_id <- as.character(result_query_id)
  expected_query_id <- as.character(expected_query_id)

  duplicate_ids <- unique(result_query_id[duplicated(result_query_id)])
  if (length(duplicate_ids) > 0) {
    stop(
      "`result$query_id` contains duplicate values: ",
      paste(duplicate_ids, collapse = ", "),
      call. = FALSE
    )
  }

  missing_ids <- setdiff(expected_query_id, result_query_id)
  unexpected_ids <- setdiff(result_query_id, expected_query_id)
  if (length(missing_ids) > 0 || length(unexpected_ids) > 0) {
    pieces <- character()
    if (length(missing_ids) > 0) {
      pieces <- c(pieces, paste0("missing: ", paste(missing_ids, collapse = ", ")))
    }
    if (length(unexpected_ids) > 0) {
      pieces <- c(pieces, paste0("unexpected: ", paste(unexpected_ids, collapse = ", ")))
    }
    stop(
      "`result$query_id` must match the expected unique query IDs exactly (",
      paste(pieces, collapse = "; "),
      ").",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.enrich_flag_taxa_final_output <- function(
  result,
  feature_query_map,
  unique_query_table,
  expected_environment,
  expected_habitat,
  expected_region,
  taxon_evidence = NULL,
  tool_evidence = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  max_evidence_rows_per_query = 5,
  max_evidence_summary_chars = 160
) {
  if (nrow(result) != nrow(feature_query_map)) {
    stop(
      "Internal error: feature-level result and feature query map differ in row count.",
      call. = FALSE
    )
  }

  evidence_provenance <- .build_flag_taxa_evidence_provenance(
    unique_query_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  feature_key <- .flag_taxa_query_key(feature_query_map)
  evidence_index <- match(feature_key, evidence_provenance$query_key)
  evidence_sources <- rep("none", length(feature_key))
  evidence_summary <- rep("No matched local evidence.", length(feature_key))
  worms_environment <- rep(NA_character_, length(feature_key))
  matched_evidence <- !is.na(evidence_index)
  evidence_sources[matched_evidence] <- evidence_provenance$evidence_sources[evidence_index[matched_evidence]]
  evidence_summary[matched_evidence] <- evidence_provenance$evidence_summary[evidence_index[matched_evidence]]
  worms_environment[matched_evidence] <- evidence_provenance$worms_environment[evidence_index[matched_evidence]]

  out <- data.frame(
    .flag_taxa_original_row = as.integer(feature_query_map$original_row),
    feature_id = as.character(feature_query_map$feature_id),
    query_id = as.character(feature_query_map$query_id),
    taxon_name = as.character(feature_query_map$query_name),
    taxon_rank = as.character(feature_query_map$query_rank),
    lineage = as.character(feature_query_map$lineage),
    worms_environment = worms_environment,
    expected_environment = rep(expected_environment, nrow(feature_query_map)),
    expected_habitat = rep(.nullable_character(expected_habitat), nrow(feature_query_map)),
    expected_region = rep(.nullable_character(expected_region), nrow(feature_query_map)),
    expected_environment_status = as.character(result$expected_environment_status),
    expected_habitat_status = as.character(result$expected_habitat_status),
    expected_region_status = as.character(result$expected_region_status),
    ecological_status = as.character(result$ecological_status),
    occurrence_interpretation = as.character(result$occurrence_interpretation),
    recommended_action = as.character(result$recommended_action),
    rationale = as.character(result$rationale),
    references = as.character(result$references),
    evidence_basis = as.character(result$evidence_basis),
    llm_tool_used = as.logical(result$llm_tool_used),
    llm_tool_name = as.character(result$llm_tool_name),
    llm_tool_query = as.character(result$llm_tool_query),
    llm_tool_evidence_summary = as.character(result$llm_tool_evidence_summary),
    llm_tool_references = as.character(result$llm_tool_references),
    evidence_sources = evidence_sources,
    evidence_summary = evidence_summary,
    llm_prompt_selected = if ("query_id" %in% colnames(result)) {
      as.character(feature_query_map$query_id) %in% as.character(result$query_id)
    } else {
      rep(TRUE, nrow(feature_query_map))
    },
    llm_selected = if ("query_id" %in% colnames(result)) {
      as.character(feature_query_map$query_id) %in% as.character(result$query_id)
    } else {
      rep(TRUE, nrow(feature_query_map))
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  out
}

.normalise_flag_taxa_tool_columns <- function(result) {
  if (!("evidence_basis" %in% colnames(result))) {
    result$evidence_basis <- "conservative_reasoning_only"
  }
  if (!("llm_tool_used" %in% colnames(result))) {
    result$llm_tool_used <- FALSE
  }
  if (!("llm_tool_name" %in% colnames(result))) {
    result$llm_tool_name <- NA_character_
  }
  if (!("llm_tool_query" %in% colnames(result))) {
    result$llm_tool_query <- NA_character_
  }
  if (!("llm_tool_evidence_summary" %in% colnames(result))) {
    result$llm_tool_evidence_summary <- NA_character_
  }
  if (!("llm_tool_references" %in% colnames(result))) {
    result$llm_tool_references <- NA_character_
  }

  result
}

.normalise_flag_taxa_interpretation_columns <- function(result) {
  if (!inherits(result, "data.frame")) {
    return(result)
  }
  if (!("ecological_status" %in% colnames(result))) {
    result$ecological_status <- .derive_flag_taxa_ecological_status(result)
  }
  if (!("occurrence_interpretation" %in% colnames(result))) {
    result$occurrence_interpretation <- .derive_flag_taxa_occurrence_interpretation(result)
  }
  result
}

.derive_flag_taxa_ecological_status <- function(result) {
  n <- nrow(result)
  env <- .optional_result_column(result, "expected_environment_status", n)
  habitat <- .optional_result_column(result, "expected_habitat_status", n)
  action <- .optional_result_column(result, "recommended_action", n)

  out <- rep("unknown", n)
  out[env == "insufficient_taxonomic_resolution" | habitat == "insufficient_taxonomic_resolution"] <-
    "insufficient_taxonomic_resolution"
  out[env %in% c("compatible", "mixed_within_rank", "possible_likely") &
    habitat %in% c("compatible", "possible_likely")] <- "compatible"
  out[env %in% c("compatible", "mixed_within_rank", "possible_likely") &
    habitat %in% c("no_known_habitat_evidence", "unknown")] <- "plausible"
  out[habitat == "incompatible"] <- "incompatible_resident"
  out[env == "incompatible" | action %in% c("exclude")] <- "incompatible_resident"
  out
}

.derive_flag_taxa_occurrence_interpretation <- function(result) {
  n <- nrow(result)
  env <- .optional_result_column(result, "expected_environment_status", n)
  habitat <- .optional_result_column(result, "expected_habitat_status", n)
  action <- .optional_result_column(result, "recommended_action", n)

  out <- rep("uncertain", n)
  out[env == "insufficient_taxonomic_resolution" | habitat == "insufficient_taxonomic_resolution"] <-
    "not_assessed"
  out[env %in% c("compatible", "mixed_within_rank", "possible_likely") &
    habitat %in% c("compatible", "possible_likely")] <- "expected_resident"
  out[env %in% c("compatible", "mixed_within_rank", "possible_likely") &
    habitat %in% c("no_known_habitat_evidence", "unknown")] <- "plausible_resident"
  out[habitat == "incompatible"] <- "possible_transient_or_allochthonous"
  out[action %in% c("exclude")] <- "likely_contaminant_or_misassignment"
  out
}

.optional_result_column <- function(result, column_name, n) {
  if (column_name %in% colnames(result)) {
    return(as.character(result[[column_name]]))
  }
  rep(NA_character_, n)
}

.normalise_flag_taxa_categorical_columns <- function(result, expected_region = NULL) {
  if (!inherits(result, "data.frame")) {
    return(result)
  }
  allowed <- .flag_taxa_allowed_values(
    expected_region = expected_region,
    for_validation = FALSE
  )
  for (column_name in .flag_taxa_final_categorical_columns()) {
    if (!(column_name %in% colnames(result))) {
      next
    }
    values <- result[[column_name]]
    if (is.factor(values)) {
      values <- as.character(values)
    } else if (is.numeric(values) || is.integer(values)) {
      values <- .map_flag_taxa_enum_codes(
        values,
        column_name = column_name,
        allowed_values = allowed[[column_name]]
      )
    } else {
      values <- as.character(values)
      coded <- !is.na(values) & grepl("^[0-9]+$", values)
      if (any(coded)) {
        values[coded] <- .map_flag_taxa_enum_codes(
          values[coded],
          column_name = column_name,
          allowed_values = allowed[[column_name]]
        )
      }
    }
    if (identical(column_name, "recommended_action")) {
      values[!is.na(values) & values == "review"] <- "flag_for_review"
    }
    result[[column_name]] <- as.character(values)
  }

  result
}

.map_flag_taxa_enum_codes <- function(values, column_name, allowed_values) {
  if (is.null(allowed_values)) {
    stop(
      "Internal flag_taxa error: numeric enum codes were detected in `",
      column_name,
      "`, but no explicit enum map is available.",
      call. = FALSE
    )
  }
  codes <- suppressWarnings(as.integer(values))
  invalid <- is.na(codes) |
    codes < 1L |
    codes > length(allowed_values) |
    as.numeric(codes) != suppressWarnings(as.numeric(values))
  if (any(invalid)) {
    stop(
      "Internal flag_taxa error: numeric enum codes in `",
      column_name,
      "` could not be safely mapped to labels.",
      call. = FALSE
    )
  }

  allowed_values[codes]
}

.flag_taxa_final_categorical_columns <- function() {
  c(
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "evidence_basis"
  )
}

.finalise_flag_taxa_output <- function(result, expected_region = NULL) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(result)
  }
  result <- .normalise_flag_taxa_interpretation_columns(result)
  result <- .normalise_flag_taxa_categorical_columns(
    result,
    expected_region = expected_region
  )
  for (column_name in .flag_taxa_final_categorical_columns()) {
    if (!(column_name %in% colnames(result))) {
      next
    }
    if (!is.character(result[[column_name]])) {
      stop(
        "Internal flag_taxa error: `",
        column_name,
        "` must be returned as character labels, not enum codes.",
        call. = FALSE
      )
    }
    coded <- !is.na(result[[column_name]]) & grepl("^[0-9]+$", result[[column_name]])
    if (any(coded)) {
      stop(
        "Internal flag_taxa error: `",
        column_name,
        "` contains bare numeric enum codes after finalisation.",
        call. = FALSE
      )
    }
  }

  .order_flag_taxa_final_columns(result)
}

.flag_taxa_final_output_columns <- function() {
  c(
    "feature_id",
    "taxon_name",
    "taxon_rank",
    "lineage",
    "worms_environment",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "evidence_sources",
    "evidence_summary",
    "llm_prompt_selected",
    "llm_tool_used",
    "llm_tool_name",
    "llm_tool_query",
    "llm_tool_evidence_summary",
    "llm_tool_references",
    "expected_environment",
    "expected_habitat",
    "expected_region",
    "llm_selected"
  )
}

.order_flag_taxa_final_columns <- function(result) {
  leading <- intersect(.flag_taxa_final_output_columns(), colnames(result))
  trailing <- setdiff(colnames(result), leading)
  result[, c(leading, trailing), drop = FALSE]
}

.validate_flag_taxa_tool_traceability <- function(result) {
  if (!is.logical(result$llm_tool_used) || any(is.na(result$llm_tool_used))) {
    stop("`llm_tool_used` must contain only TRUE or FALSE values.", call. = FALSE)
  }

  tool_basis <- result$evidence_basis %in% c("tool_evidence", "local_and_tool_evidence")
  if (any(tool_basis & !result$llm_tool_used)) {
    stop(
      "`llm_tool_used` must be TRUE when `evidence_basis` is ",
      "`tool_evidence` or `local_and_tool_evidence`.",
      call. = FALSE
    )
  }

  used_tool <- result$llm_tool_used
  if (any(used_tool & !.is_non_empty_trace_value(result$llm_tool_name))) {
    stop(
      "`llm_tool_name` must not be missing or empty when `llm_tool_used` is TRUE.",
      call. = FALSE
    )
  }
  if (any(used_tool & !.is_non_empty_trace_value(result$llm_tool_evidence_summary))) {
    stop(
      "`llm_tool_evidence_summary` must not be missing or empty when ",
      "`llm_tool_used` is TRUE.",
      call. = FALSE
    )
  }
  if (any(used_tool & !.is_non_empty_trace_value(result$references))) {
    stop(
      "`references` must not be missing or empty when `llm_tool_used` is TRUE.",
      call. = FALSE
    )
  }
  if (any(used_tool & !.is_non_empty_trace_value(result$llm_tool_references))) {
    stop(
      "`llm_tool_references` must not be missing or empty when `llm_tool_used` is TRUE. ",
      "Use 'Tool used but no verifiable references were returned.' if the tool returned no usable references.",
      call. = FALSE
    )
  }
  vague_reference <- used_tool & .is_vague_flag_taxa_tool_reference(result$llm_tool_references)
  if (any(vague_reference)) {
    stop(
      "`llm_tool_references` must contain verifiable references, DOIs, URLs, ",
      "or the explicit statement 'Tool used but no verifiable references were returned.' ",
      "Vague placeholders such as 'scientific_literature (tool evidence summary)' are not sufficient.",
      call. = FALSE
    )
  }

  invisible(result)
}

.is_vague_flag_taxa_tool_reference <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value[is.na(value)] <- ""
  allowed_no_reference <- value == "tool used but no verifiable references were returned."
  vague <- value %in% c(
    "scientific_literature (tool evidence summary)",
    "scientific literature (tool evidence summary)",
    "tool evidence summary",
    "scientific_literature",
    "scientific literature"
  )
  vague & !allowed_no_reference
}

.is_non_empty_trace_value <- function(value) {
  value <- as.character(value)
  .is_non_empty_value(value) & !(tolower(trimws(value)) %in% c("na", "null", "none"))
}

.validate_flag_taxa_exclusion_evidence <- function(result, allow_llm_tools) {
  exclude <- result$recommended_action == "exclude"
  conservative_exclude <- exclude & result$evidence_basis == "conservative_reasoning_only"
  if (any(conservative_exclude)) {
    stop(
      "`recommended_action = \"exclude\"` cannot be used with ",
      "`evidence_basis = \"conservative_reasoning_only\"`.",
      call. = FALSE
    )
  }

  habitat_absence_exclude <- exclude &
    result$expected_habitat_status %in% c(
      "no_known_habitat_evidence",
      "possible_likely",
      "possible_unlikely"
    )
  if (any(habitat_absence_exclude)) {
    stop(
      "`recommended_action = \"exclude\"` cannot be used with ",
      "`expected_habitat_status = \"no_known_habitat_evidence\", ",
      "\"possible_likely\", or \"possible_unlikely\"`. These are uncertainty ",
      "or plausibility statuses, not explicit incompatibility.",
      call. = FALSE
    )
  }

  no_evidence_tool_exclude <- exclude &
    result$llm_tool_used &
    .flag_taxa_tool_summary_no_evidence(result$llm_tool_evidence_summary) &
    result$evidence_basis != "local_evidence"
  if (any(no_evidence_tool_exclude)) {
    stop(
      "`recommended_action = \"exclude\"` requires explicit positive ",
      "incompatibility evidence. Tool summaries indicating no evidence found ",
      "must be treated as uncertainty/review, not exclusion.",
      call. = FALSE
    )
  }

  no_evidence_tool_habitat_incompatible <-
    result$expected_habitat_status == "incompatible" &
    result$llm_tool_used &
    .flag_taxa_tool_summary_no_evidence(result$llm_tool_evidence_summary) &
    result$evidence_basis != "local_evidence"
  if (any(no_evidence_tool_habitat_incompatible)) {
    stop(
      "`expected_habitat_status = \"incompatible\"` requires explicit ",
      "positive habitat incompatibility evidence. Tool summaries indicating ",
      "no evidence found must be treated as uncertainty, not incompatibility.",
      call. = FALSE
    )
  }

  if (isTRUE(allow_llm_tools)) {
    supported_exclude <- result$evidence_basis == "local_evidence" |
      (result$evidence_basis %in% c("tool_evidence", "local_and_tool_evidence") &
        result$llm_tool_used)
    if (any(exclude & !supported_exclude)) {
      stop(
        "When `allow_llm_tools = TRUE`, `recommended_action = \"exclude\"` ",
        "requires local evidence or explicit tool evidence.",
        call. = FALSE
      )
    }
  }

  invisible(result)
}

.flag_taxa_tool_summary_no_evidence <- function(summary) {
  text <- tolower(as.character(summary))
  out <- grepl(
    "no_evidence_found|no explicit evidence|no evidence found|no useful evidence|no records found|no known habitat|no known regional|no distribution evidence",
    text
  )
  out[is.na(out)] <- FALSE
  out
}

.repair_flag_taxa_invalid_llm_output <- function(
  result,
  repair_invalid_llm = TRUE,
  expected_region = NULL,
  verbose = FALSE,
  start_time = NULL
) {
  if (!inherits(result, "data.frame")) {
    return(result)
  }
  repair_invalid_llm <- .validate_logical_scalar(repair_invalid_llm, "repair_invalid_llm")
  verbose <- .validate_logical_scalar(verbose, "verbose")
  result <- .normalise_flag_taxa_tool_columns(result)
  result <- .normalise_flag_taxa_interpretation_columns(result)
  result <- .normalise_flag_taxa_categorical_columns(
    result,
    expected_region = expected_region
  )
  missing_tool_references <- .flag_taxa_missing_tool_reference_rows(result)
  if (any(missing_tool_references) && isTRUE(repair_invalid_llm)) {
    repaired_count <- sum(missing_tool_references)
    result$llm_tool_references[missing_tool_references] <-
      .flag_taxa_missing_tool_reference_value()
    .flag_taxa_progress(
      verbose,
      start_time,
      "filled missing llm_tool_references for ",
      repaired_count,
      " tool-used row",
      if (identical(repaired_count, 1L)) "" else "s",
      "."
    )
  }
  unsupported <- .flag_taxa_repairable_overcall_rows(result)
  if (!any(unsupported) || !isTRUE(repair_invalid_llm)) {
    return(result)
  }

  if ("recommended_action" %in% colnames(result)) {
    result$recommended_action[unsupported] <- "flag_for_review"
  }
  if ("expected_habitat_status" %in% colnames(result)) {
    repair_habitat <- unsupported &
      result$expected_habitat_status %in% c("incompatible", "no_known_habitat_evidence", "possible_likely", "possible_unlikely")
    result$expected_habitat_status[repair_habitat] <- "no_known_habitat_evidence"
  }
  if ("rationale" %in% colnames(result)) {
    result$rationale[unsupported] <- paste(
      as.character(result$rationale[unsupported]),
      "biohelper repair: unsupported exclusion or habitat incompatibility was downgraded to flag_for_review because absence of evidence is not positive incompatibility evidence."
    )
  }
  if ("evidence_basis" %in% colnames(result)) {
    conservative <- unsupported &
      result$evidence_basis == "conservative_reasoning_only"
    result$evidence_basis[conservative] <- "conservative_reasoning_only"
  }

  result
}

.flag_taxa_missing_tool_reference_rows <- function(result) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(logical())
  }
  result <- .normalise_flag_taxa_tool_columns(result)
  used_tool <- as.logical(result$llm_tool_used)
  used_tool[is.na(used_tool)] <- FALSE
  out <- used_tool & !.is_non_empty_trace_value(result$llm_tool_references)
  out[is.na(out)] <- FALSE
  out
}

.flag_taxa_missing_tool_reference_value <- function() {
  "Tool used but no verifiable references were returned."
}

.flag_taxa_required_tool_missing_summary <- function() {
  "Tool evidence was required but no usable tool evidence was returned; result should be treated conservatively."
}

.flag_taxa_required_tool_missing_reference <- function() {
  "Tool evidence was required but no verifiable references were returned."
}

.flag_taxa_no_explicit_tool_evidence_summary <- function() {
  "Tool search found no explicit evidence relevant to the expected context; this is treated as missing evidence, not incompatibility."
}

.flag_taxa_default_tool_name <- function(prompt_tools = NULL) {
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  if (length(prompt_tools) > 0) {
    return(.format_prompt_tool_names(prompt_tools))
  }
  "registered_tool"
}

.flag_taxa_required_tool_query_text <- function(
  result,
  row_index,
  query_taxonomy = NULL,
  expected_environment = NULL,
  expected_habitat = NULL,
  expected_region = NULL
) {
  taxon_name <- if ("taxon_name" %in% colnames(result)) {
    as.character(result$taxon_name[[row_index]])
  } else if ("query_name" %in% colnames(result)) {
    as.character(result$query_name[[row_index]])
  } else {
    NA_character_
  }
  taxon_rank <- if ("taxon_rank" %in% colnames(result)) {
    as.character(result$taxon_rank[[row_index]])
  } else if ("query_rank" %in% colnames(result)) {
    as.character(result$query_rank[[row_index]])
  } else {
    NA_character_
  }
  query_id <- if ("query_id" %in% colnames(result)) {
    as.character(result$query_id[[row_index]])
  } else {
    NA_character_
  }
  if (!is.null(query_taxonomy) && .is_non_empty_value(query_id) && "query_id" %in% colnames(query_taxonomy)) {
    query_match <- match(query_id, as.character(query_taxonomy$query_id))
    if (!is.na(query_match)) {
      if ("query_name" %in% colnames(query_taxonomy)) {
        taxon_name <- as.character(query_taxonomy$query_name[[query_match]])
      }
      if ("query_rank" %in% colnames(query_taxonomy)) {
        taxon_rank <- as.character(query_taxonomy$query_rank[[query_match]])
      }
    }
  }
  pieces <- c(taxon_name, taxon_rank, expected_environment, expected_habitat, expected_region)
  pieces <- pieces[.is_non_empty_value(pieces)]
  if (length(pieces) == 0) {
    return("required tool evidence query")
  }
  paste(pieces, collapse = " ")
}

.repair_flag_taxa_required_tool_attempts <- function(
  result,
  tool_requirement = "optional",
  allow_llm_tools = FALSE,
  repair_invalid_llm = TRUE,
  query_taxonomy = NULL,
  prompt_tools = NULL,
  expected_environment = NULL,
  expected_habitat = NULL,
  expected_region = NULL,
  verbose = FALSE,
  start_time = NULL
) {
  if (!inherits(result, "data.frame") || nrow(result) == 0) {
    return(result)
  }
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  allow_llm_tools <- .validate_logical_scalar(allow_llm_tools, "allow_llm_tools")
  repair_invalid_llm <- .validate_logical_scalar(repair_invalid_llm, "repair_invalid_llm")
  verbose <- .validate_logical_scalar(verbose, "verbose")
  if (!isTRUE(allow_llm_tools) || !identical(tool_requirement, "required_for_llm")) {
    return(result)
  }

  result <- .normalise_flag_taxa_tool_columns(result)
  result <- .normalise_flag_taxa_interpretation_columns(result)
  used_tool <- as.logical(result$llm_tool_used)
  used_tool[is.na(used_tool)] <- FALSE
  missing_attempt <- !used_tool
  if (!any(missing_attempt)) {
    return(result)
  }
  if (!isTRUE(repair_invalid_llm)) {
    stop(
      "`tool_requirement = \"required_for_llm\"` requires every LLM-selected row ",
      "to report a tool attempt. Missing tool attempts were found for ",
      sum(missing_attempt),
      " row(s).",
      call. = FALSE
    )
  }

  tool_name <- .flag_taxa_default_tool_name(prompt_tools)
  row_indexes <- which(missing_attempt)
  for (row_index in row_indexes) {
    result$llm_tool_used[[row_index]] <- TRUE
    if ("llm_tool_status" %in% colnames(result)) {
      result$llm_tool_status[[row_index]] <- "failed"
    }
    if (!.is_non_empty_trace_value(result$llm_tool_name[[row_index]])) {
      result$llm_tool_name[[row_index]] <- tool_name
    }
    if (!.is_non_empty_trace_value(result$llm_tool_query[[row_index]])) {
      result$llm_tool_query[[row_index]] <- .flag_taxa_required_tool_query_text(
        result = result,
        row_index = row_index,
        query_taxonomy = query_taxonomy,
        expected_environment = expected_environment,
        expected_habitat = expected_habitat,
        expected_region = expected_region
      )
    }
    if (!.is_non_empty_trace_value(result$llm_tool_evidence_summary[[row_index]])) {
      result$llm_tool_evidence_summary[[row_index]] <- .flag_taxa_required_tool_missing_summary()
    }
    if (!.is_non_empty_trace_value(result$llm_tool_references[[row_index]])) {
      result$llm_tool_references[[row_index]] <- .flag_taxa_required_tool_missing_reference()
    }
    if ("references" %in% colnames(result) && !.is_non_empty_trace_value(result$references[[row_index]])) {
      result$references[[row_index]] <- .flag_taxa_required_tool_missing_reference()
    }
  }

  if ("recommended_action" %in% colnames(result)) {
    downgraded <- missing_attempt & result$recommended_action == "exclude"
    if (any(downgraded)) {
      result$recommended_action[downgraded] <- "flag_for_review"
    }
  }
  if ("rationale" %in% colnames(result)) {
    result$rationale[missing_attempt] <- paste(
      as.character(result$rationale[missing_attempt]),
      "biohelper repair: tool evidence was required for this LLM-selected row, but no tool attempt was reported; the row was marked as a failed/no-evidence tool attempt and interpreted conservatively."
    )
  }

  .flag_taxa_progress(
    verbose,
    start_time,
    "filled required tool-evidence metadata for ",
    sum(missing_attempt),
    " LLM-selected row",
    if (identical(sum(missing_attempt), 1L)) "" else "s",
    "."
  )

  result
}

.flag_taxa_repairable_overcall_rows <- function(result) {
  if (!all(c("recommended_action", "expected_habitat_status", "evidence_basis") %in% colnames(result))) {
    return(rep(FALSE, nrow(result)))
  }
  exclude <- result$recommended_action == "exclude"
  conservative_exclude <- exclude & result$evidence_basis == "conservative_reasoning_only"
  uncertain_habitat_exclude <- exclude &
    result$expected_habitat_status %in% c(
      "no_known_habitat_evidence",
      "possible_likely",
      "possible_unlikely"
    )
  no_evidence_tool_exclude <- exclude &
    result$llm_tool_used &
    .flag_taxa_tool_summary_no_evidence(result$llm_tool_evidence_summary) &
    result$evidence_basis != "local_evidence"
  no_evidence_tool_habitat_incompatible <-
    result$expected_habitat_status == "incompatible" &
    result$llm_tool_used &
    .flag_taxa_tool_summary_no_evidence(result$llm_tool_evidence_summary) &
    result$evidence_basis != "local_evidence"

  out <- conservative_exclude |
    uncertain_habitat_exclude |
    no_evidence_tool_exclude |
    no_evidence_tool_habitat_incompatible
  out[is.na(out)] <- FALSE
  out
}

.augment_flag_taxa_result_with_tool_evidence <- function(
  result,
  query_taxonomy,
  taxon_evidence,
  tool_evidence,
  tax_ranks,
  max_evidence_rows_per_query,
  max_evidence_summary_chars = 160
) {
  if (is.null(tool_evidence) || !inherits(result, "data.frame")) {
    return(result)
  }
  tool_evidence <- .validate_flag_taxa_tool_evidence(tool_evidence)
  if (nrow(tool_evidence) == 0) {
    return(result)
  }

  result <- .normalise_flag_taxa_tool_columns(result)
  result <- .normalise_flag_taxa_interpretation_columns(result)
  if (!("query_id" %in% colnames(result))) {
    result$query_id <- query_taxonomy$query_id
    added_query_id <- TRUE
  } else {
    added_query_id <- FALSE
  }

  local_sources <- .build_flag_taxa_evidence_provenance(
    unique_query_table = query_taxonomy,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  local_lookup <- stats::setNames(
    local_sources$evidence_sources,
    query_taxonomy$query_id
  )
  tool_lookup <- split(tool_evidence, as.character(tool_evidence$query_id))

  for (row_index in seq_len(nrow(result))) {
    query_id <- as.character(result$query_id[[row_index]])
    if (!(query_id %in% names(tool_lookup))) {
      next
    }
    tool_rows <- tool_lookup[[query_id]]
    explicit_rows <- tool_rows[tool_rows$tool_found_explicit_evidence | tool_rows$tool_used, , drop = FALSE]
    if (nrow(explicit_rows) == 0) {
      next
    }
    first_row <- explicit_rows[1, , drop = FALSE]
    result$llm_tool_used[[row_index]] <- TRUE
    if (!.is_non_empty_trace_value(result$llm_tool_name[[row_index]])) {
      result$llm_tool_name[[row_index]] <- .first_non_empty_value(first_row$tool_name)
    }
    if (!.is_non_empty_trace_value(result$llm_tool_query[[row_index]])) {
      result$llm_tool_query[[row_index]] <- .first_non_empty_value(first_row$tool_query)
    }
    if (!.is_non_empty_trace_value(result$llm_tool_evidence_summary[[row_index]])) {
      result$llm_tool_evidence_summary[[row_index]] <- .first_non_empty_value(first_row$tool_evidence_summary)
    }
    if (!.is_non_empty_trace_value(result$llm_tool_references[[row_index]])) {
      result$llm_tool_references[[row_index]] <- .first_non_empty_value(first_row$tool_references)
    }
    if (!.is_non_empty_trace_value(result$references[[row_index]])) {
      result$references[[row_index]] <- .first_non_empty_value(first_row$tool_references)
    }

    local_value <- if (query_id %in% names(local_lookup)) {
      local_lookup[[query_id]]
    } else {
      NA_character_
    }
    has_local <- .is_non_empty_value(local_value) &&
      !identical(local_value, "none")
    result$evidence_basis[[row_index]] <- if (isTRUE(has_local)) {
      "local_and_tool_evidence"
    } else {
      "tool_evidence"
    }
  }

  if (isTRUE(added_query_id)) {
    result$query_id <- NULL
  }
  result
}

.build_flag_taxa_evidence_provenance <- function(
  unique_query_table,
  taxon_evidence,
  tax_ranks,
  max_evidence_rows_per_query,
  max_evidence_summary_chars = 160
) {
  out <- data.frame(
    query_key = .flag_taxa_query_key(unique_query_table),
    evidence_sources = rep("none", nrow(unique_query_table)),
    evidence_summary = rep("No matched local evidence.", nrow(unique_query_table)),
    worms_environment = rep(NA_character_, nrow(unique_query_table)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (is.null(taxon_evidence) || nrow(taxon_evidence) == 0) {
    return(out)
  }

  relevant_evidence <- collect_relevant_taxon_evidence(
    tax_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  if (is.null(relevant_evidence) || nrow(relevant_evidence) == 0) {
    return(out)
  }

  for (query_index in seq_len(nrow(unique_query_table))) {
    group <- relevant_evidence[
      relevant_evidence$taxonomy_row == query_index,
      ,
      drop = FALSE
    ]
    if (nrow(group) == 0) {
      next
    }
    out$evidence_sources[[query_index]] <- .collapse_distinct_values(
      group$source,
      max_values = 12
    )
    out$evidence_summary[[query_index]] <- .format_flag_taxa_evidence_summary(group)
    out$worms_environment[[query_index]] <- .flag_taxa_worms_environment(group)
  }

  out
}

.flag_taxa_worms_environment <- function(evidence_rows) {
  if (is.null(evidence_rows) || nrow(evidence_rows) == 0) {
    return(NA_character_)
  }
  source_values <- tolower(as.character(.optional_evidence_column(evidence_rows, "source")))
  worms_rows <- startsWith(source_values, "worms")
  if (!any(worms_rows)) {
    return(NA_character_)
  }
  values <- .flag_taxa_evidence_values(evidence_rows[worms_rows, , drop = FALSE], "environment")
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(unique(values), collapse = "; ")
}

.format_flag_taxa_evidence_summary <- function(evidence_rows) {
  source_values <- unique(as.character(evidence_rows$source)[.is_non_empty_value(evidence_rows$source)])
  if (length(source_values) == 0) {
    return("No matched local evidence.")
  }

  pieces <- vapply(source_values, function(source_value) {
    group <- evidence_rows[as.character(evidence_rows$source) == source_value, , drop = FALSE]
    detail <- character()
    evidence_type <- .collapse_distinct_values(.optional_evidence_column(group, "evidence_type"), max_values = 3)
    environment <- .collapse_distinct_values(.optional_evidence_column(group, "environment"), max_values = 4)
    habitat <- .collapse_distinct_values(.optional_evidence_column(group, "habitat"), max_values = 4)
    region <- .collapse_distinct_values(.optional_evidence_column(group, "region"), max_values = 4)
    if (.is_non_empty_value(evidence_type)) {
      detail <- c(detail, paste0("evidence_type=", evidence_type))
    }
    if (.is_non_empty_value(environment)) {
      detail <- c(detail, paste0("environment=", environment))
    }
    if (.is_non_empty_value(habitat)) {
      detail <- c(detail, paste0("habitat=", habitat))
    }
    if (.is_non_empty_value(region)) {
      detail <- c(detail, paste0("region=", region))
    }
    if (any(grepl("^Summarised evidence:", group$evidence_summary))) {
      detail <- c(detail, "summarised evidence")
    }
    if (length(detail) == 0) {
      detail <- "matched local evidence"
    }

    paste0(source_value, ": ", paste(detail, collapse = ", "))
  }, character(1))

  paste(pieces, collapse = "; ")
}

.nullable_character <- function(value) {
  if (is.null(value)) {
    return(NA_character_)
  }

  as.character(value)
}

.flag_taxa_query_key <- function(query_table) {
  paste(
    .flag_taxa_result_query_key(query_table$query_name, query_table$query_rank),
    .taxon_match_key(query_table$lineage),
    sep = "\r"
  )
}

.flag_taxa_result_query_key <- function(query_name, query_rank) {
  paste(
    .taxon_match_key(query_name),
    tolower(trimws(as.character(query_rank))),
    sep = "\r"
  )
}

.flag_taxa_feature_id <- function(tax_table) {
  column_lookup <- stats::setNames(colnames(tax_table), tolower(colnames(tax_table)))
  for (column_name in c("feature_id", "taxon_id", "asv", "otu")) {
    lookup_name <- tolower(column_name)
    if (lookup_name %in% names(column_lookup)) {
      matched_column <- column_lookup[[lookup_name]]
      return(as.character(tax_table[[matched_column]]))
    }
  }

  as.character(seq_len(nrow(tax_table)))
}

.flag_taxa_reserved_input_columns <- function() {
  tolower(c(
    .flag_taxa_output_columns(),
    .flag_taxa_structured_output_columns(include_feature_id = FALSE)
  ))
}

#' Validate user-provided taxon evidence
#'
#' @param taxon_evidence Optional evidence table supplied to `flag_taxa()`.
#'
#' @return `NULL`, or a data frame with standardised column order.
#' @noRd
validate_taxon_evidence <- function(taxon_evidence) {
  if (is.null(taxon_evidence)) {
    return(NULL)
  }

  if (!inherits(taxon_evidence, "data.frame")) {
    stop("`taxon_evidence` must be a data.frame.", call. = FALSE)
  }

  required_columns <- .taxon_evidence_required_columns()
  missing_columns <- setdiff(required_columns, colnames(taxon_evidence))
  if (length(missing_columns) > 0) {
    stop(
      "`taxon_evidence` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(taxon_evidence) > 0) {
    empty_columns <- required_columns[
      vapply(taxon_evidence[required_columns], .is_entirely_empty_column, logical(1))
    ]
    if (
      "reference" %in% empty_columns &&
        all(as.character(taxon_evidence$evidence_type) %in% c(
          "not_queried",
          "taxonomy_lookup_failed"
        ))
    ) {
      empty_columns <- setdiff(empty_columns, "reference")
    }
    if (length(empty_columns) > 0) {
      stop(
        "Required `taxon_evidence` columns must not be entirely empty: ",
        paste(empty_columns, collapse = ", "),
        call. = FALSE
      )
    }
  }

  ordered_columns <- c(
    required_columns,
    intersect(.taxon_evidence_optional_columns(), colnames(taxon_evidence))
  )
  ordered_columns <- c(ordered_columns, setdiff(colnames(taxon_evidence), ordered_columns))

  taxon_evidence[, ordered_columns, drop = FALSE]
}

.validate_flag_taxa_evidence_sources <- function(evidence_sources) {
  if (is.null(evidence_sources)) {
    return(NULL)
  }

  evidence_sources <- .validate_character_vector(
    evidence_sources,
    "evidence_sources"
  )
  supported_sources <- "worms"
  unsupported_sources <- setdiff(evidence_sources, supported_sources)
  if (length(unsupported_sources) > 0) {
    stop(
      "`evidence_sources` contains unsupported value(s): ",
      paste(unsupported_sources, collapse = ", "),
      ". Currently supported value: worms.",
      call. = FALSE
    )
  }

  evidence_sources
}

.collect_flag_taxa_requested_evidence <- function(
  x,
  evidence_sources,
  taxon_evidence,
  ccz_evidence_path,
  worms_cache_path,
  worms_sleep,
  worms_batch_size,
  worms_max_tries,
  worms_retry_sleep,
  continue_on_worms_error,
  worms_verbose
) {
  if (is.null(evidence_sources) && is.null(ccz_evidence_path)) {
    return(taxon_evidence)
  }

  fetched_evidence <- NULL
  if ("worms" %in% evidence_sources) {
    taxa_to_query <- .flag_taxa_extract_taxa_for_evidence(x)
    fetched_evidence <- .flag_taxa_fetch_worms_evidence(
      taxa = taxa_to_query,
      by = "name",
      cache_path = worms_cache_path,
      sleep = worms_sleep,
      batch_size = worms_batch_size,
      max_tries = worms_max_tries,
      retry_sleep = worms_retry_sleep,
      continue_on_error = continue_on_worms_error,
      verbose = worms_verbose
    )
  }

  ccz_evidence <- NULL
  if (!is.null(ccz_evidence_path)) {
    ccz_evidence <- .flag_taxa_read_worms_ccz_evidence(ccz_evidence_path)
  }

  .combine_taxon_evidence_tables(taxon_evidence, fetched_evidence, ccz_evidence)
}

.flag_taxa_extract_taxa_for_evidence <- function(x) {
  extract_taxa_for_evidence(x)
}

.flag_taxa_fetch_worms_evidence <- function(
  taxa,
  by,
  cache_path,
  sleep,
  batch_size,
  max_tries,
  retry_sleep,
  continue_on_error,
  verbose
) {
  fetch_worms_evidence(
    taxa = taxa,
    by = by,
    cache_path = cache_path,
    sleep = sleep,
    batch_size = batch_size,
    max_tries = max_tries,
    retry_sleep = retry_sleep,
    continue_on_error = continue_on_error,
    verbose = verbose
  )
}

.flag_taxa_read_worms_ccz_evidence <- function(path) {
  read_worms_ccz_evidence(path)
}

call_flag_taxa_ellmer <- function(
  prompt,
  chat,
  include_feature_id = TRUE,
  expected_region = NULL,
  max_tries = 3,
  retry_sleep = 5
) {
  .require_ellmer_for_flag_taxa()
  if (is.null(chat) || !is.function(chat$chat_structured)) {
    stop(
      "`chat` must be an ellmer chat object with a `chat_structured()` method.",
      call. = FALSE
    )
  }
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  max_tries <- .validate_positive_integer_scalar(max_tries, "max_tries")
  retry_sleep <- .validate_non_negative_numeric_scalar(retry_sleep, "retry_sleep")

  output_type <- .flag_taxa_ellmer_output_type(
    include_feature_id = include_feature_id,
    expected_region = expected_region
  )
  for (attempt in seq_len(max_tries)) {
    result <- tryCatch(
      list(value = chat$chat_structured(prompt, type = output_type), error = NULL),
      error = function(error) {
        list(value = NULL, error = error)
      }
    )

    if (is.null(result$error)) {
      return(.as_flag_taxa_result_table(result$value))
    }

    if (
      !.is_retryable_ellmer_error(result$error) ||
        attempt >= max_tries
    ) {
      .abort_ellmer_structured_call(
        error = result$error,
        chat = chat,
        attempts = attempt
      )
    }

    if (retry_sleep > 0) {
      Sys.sleep(retry_sleep)
    }
  }

  stop("ellmer structured output call failed unexpectedly.", call. = FALSE)
}

.is_retryable_ellmer_error <- function(error) {
  status_code <- .ellmer_error_status_code(error)
  if (!is.na(status_code) && status_code %in% c(429L, 500L, 502L, 503L, 504L, 520L)) {
    return(TRUE)
  }

  message <- conditionMessage(error)
  classes <- paste(class(error), collapse = " ")
  has_retryable_http_status <- grepl(
    "\\b(429|500|502|503|504|520)\\b",
    message
  )
  has_network_error <- grepl(
    "timeout|timed out|curl|network|connection|failed to connect|could not resolve|ssl|tls|temporar",
    paste(message, classes),
    ignore.case = TRUE
  )

  isTRUE(has_retryable_http_status || has_network_error)
}

.ellmer_error_status_code <- function(error) {
  current <- error
  for (i in seq_len(5)) {
    status <- .status_code_from_condition(current)
    if (!is.na(status)) {
      return(status)
    }

    parent <- tryCatch(current$parent, error = function(error) NULL)
    if (!inherits(parent, "condition")) {
      break
    }
    current <- parent
  }

  NA_integer_
}

.status_code_from_condition <- function(error) {
  for (field in c("status_code", "status", "http_status")) {
    status <- .as_status_code(tryCatch(error[[field]], error = function(error) NULL))
    if (!is.na(status)) {
      return(status)
    }
  }

  response <- tryCatch(error$response, error = function(error) NULL)
  if (!is.null(response)) {
    for (field in c("status_code", "status")) {
      status <- .as_status_code(tryCatch(response[[field]], error = function(error) NULL))
      if (!is.na(status)) {
        return(status)
      }
    }
  }

  NA_integer_
}

.as_status_code <- function(value) {
  if (is.null(value) || length(value) != 1) {
    return(NA_integer_)
  }

  status <- suppressWarnings(as.integer(value))
  if (is.na(status)) {
    return(NA_integer_)
  }

  status
}

.abort_ellmer_structured_call <- function(error, chat, attempts) {
  backend <- .ellmer_chat_backend_label(chat)
  backend_text <- if (.is_non_empty_value(backend)) {
    paste0(" (", backend, ")")
  } else {
    ""
  }
  guidance <- .ellmer_status_guidance(error)

  stop(
    "ellmer structured output call failed after ",
    attempts,
    " attempt(s)",
    backend_text,
    ". Original error: ",
    conditionMessage(error),
    guidance,
    "\nRetry later, check OpenAI status, or test `chat$chat(\"hello\")`.",
    call. = FALSE
  )
}

.ellmer_status_guidance <- function(error) {
  status_code <- .ellmer_error_status_code(error)
  if (is.na(status_code)) {
    status_code <- .ellmer_status_code_from_message(conditionMessage(error))
  }

  if (identical(status_code, 413L)) {
    return(paste0(
      "\nFinal LLM judgement prompt was too large for this model route. ",
      "HTTP 413 indicates request payload/model context limit, not provider quota. ",
      "Try a larger-context model/provider. If you must use this route, reduce ",
      "`max_evidence_rows_per_query`, reduce `max_tool_taxa`, review fewer ",
      "taxa, disable tool evidence, or use `prompt_only = TRUE` for manual ",
      "prompt handling."
    ))
  }

  if (identical(status_code, 429L)) {
    return(paste0(
      "\nHTTP 429 is likely an LLM provider request quota/rate-limit. It may ",
      "occur before or during tool use. It does not necessarily mean Scite ",
      "or another registered tool reached quota, and it does not necessarily ",
      "indicate prompt-size problems."
    ))
  }

  if (.ellmer_error_indicates_truncated_output(error)) {
    return(paste0(
      "\nFinal LLM judgement output appears to be truncated or incomplete. ",
      "This usually means the prompt/output was too large for the model or ",
      "output limit. Try reducing `max_evidence_rows_per_query`, reducing ",
      "`max_tool_taxa`, reviewing fewer taxa, or using a larger-context / ",
      "larger-output model."
    ))
  }

  ""
}

.ellmer_error_indicates_truncated_output <- function(error) {
  text <- tolower(paste(
    conditionMessage(error),
    paste(class(error), collapse = " ")
  ))
  grepl(
    "parse error|premature eof|unexpected end of input|truncated|incomplete json|unterminated|unexpected end",
    text
  )
}

.ellmer_status_code_from_message <- function(message) {
  matched <- regmatches(
    message,
    regexpr("\\b(413|429|500|502|503|504|520)\\b", message)
  )
  if (length(matched) == 0 || !nzchar(matched)) {
    return(NA_integer_)
  }

  as.integer(matched)
}

.ellmer_chat_backend_label <- function(chat) {
  provider <- .first_non_empty_value(unlist(lapply(
    c("provider", "provider_name", "api_provider"),
    .chat_field,
    chat = chat
  ), use.names = FALSE))
  model <- .first_non_empty_value(unlist(lapply(
    c("model", "model_name", "deployment", "deployment_name"),
    .chat_field,
    chat = chat
  ), use.names = FALSE))

  pieces <- character()
  if (.is_non_empty_value(provider)) {
    pieces <- c(pieces, paste0("provider=", provider))
  }
  if (.is_non_empty_value(model)) {
    pieces <- c(pieces, paste0("model=", model))
  }

  paste(pieces, collapse = ", ")
}

.chat_field <- function(field, chat) {
  value <- tryCatch(chat[[field]], error = function(error) NULL)
  if (is.function(value)) {
    return(NULL)
  }
  value
}

.require_ellmer_for_flag_taxa <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop(
      "The ellmer package is required for LLM-backed flag_taxa(). ",
      "Install it with install.packages('ellmer') or use prompt_only = TRUE.",
      call. = FALSE
    )
  }
}

.flag_taxa_ellmer_output_type <- function(
  include_feature_id = TRUE,
  expected_region = NULL
) {
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  allowed <- .flag_taxa_allowed_values(expected_region = expected_region)
  row_fields <- list(
    query_id = ellmer::type_string("Copy of query_id for this unique query taxon."),
    taxon_name = ellmer::type_string("Copy of query_name for this input row."),
    taxon_rank = ellmer::type_string("Copy of query_rank for this input row."),
    expected_environment_status = ellmer::type_enum(
      allowed$expected_environment_status,
      "Assessment of compatibility with the expected environment."
    ),
    expected_habitat_status = ellmer::type_enum(
      allowed$expected_habitat_status,
      "Assessment of compatibility with the expected habitat."
    ),
    expected_region_status = ellmer::type_enum(
      allowed$expected_region_status,
      "Assessment of evidence for the expected region."
    ),
    ecological_status = ellmer::type_enum(
      allowed$ecological_status,
      "Ecological resident compatibility interpretation for this query taxon."
    ),
    occurrence_interpretation = ellmer::type_enum(
      allowed$occurrence_interpretation,
      "Interpretation of whether the sequence represents resident, allochthonous, contaminant, misassignment, or uncertain occurrence."
    ),
    recommended_action = ellmer::type_enum(
      allowed$recommended_action,
      "Recommended action for this taxon: retain, flag_for_review, flag_possible_exclusion, or exclude."
    ),
    rationale = ellmer::type_string(
      "Concise rationale. Unknown or missing evidence and no WoRMS match must not be treated as incompatibility."
    ),
    references = ellmer::type_string(
      "Reference labels from local evidence or registered tools only. Do not invent AphiaIDs, DOIs, URLs, article titles, author names, or citations."
    ),
    evidence_basis = ellmer::type_enum(
      allowed$evidence_basis,
      "Type of evidence supporting this taxon-level assessment."
    ),
    llm_tool_used = ellmer::type_boolean(
      "Whether any registered LLM tool was used for this taxon."
    ),
    llm_tool_name = ellmer::type_string(
      "Compact name of the registered tool used, such as scite, or NA if no tool was used."
    ),
    llm_tool_query = ellmer::type_string(
      "Compact summary of the tool query or search terms used, or NA if no tool was used."
    ),
    llm_tool_evidence_summary = ellmer::type_string(
      "Compact summary of tool-derived evidence used, or NA if no tool was used."
    ),
    llm_tool_references = ellmer::type_string(
      "Verifiable tool-derived references, titles, DOIs, URLs, or 'Tool used but no verifiable references were returned.' if no usable tool references were returned. Do not invent references."
    )
  )
  if (isTRUE(include_feature_id)) {
    row_fields <- c(
      list(feature_id = ellmer::type_string("Feature or ASV identifier copied from the input row.")),
      row_fields
    )
  }

  row_type <- do.call(
    ellmer::type_object,
    c(
      list(.description = "One taxon-level flag_taxa result row."),
      row_fields
    )
  )
  ellmer::type_array(
    row_type,
    description = "One result row for each unique query taxon, in the same order."
  )
}

.as_flag_taxa_result_table <- function(result) {
  if (inherits(result, "data.frame")) {
    return(as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (!is.list(result)) {
    stop(
      "The ellmer result could not be converted to a data.frame.",
      call. = FALSE
    )
  }
  if (length(result) == 0) {
    return(data.frame())
  }

  if (!is.null(names(result)) && all(vapply(result, length, integer(1)) == 1)) {
    return(as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (all(vapply(result, is.list, logical(1)))) {
    rows <- lapply(result, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })
    out <- .bind_rows_fill(rows)
    rownames(out) <- NULL
    return(out)
  }

  stop(
    "The ellmer result could not be converted to a data.frame.",
    call. = FALSE
  )
}

.combine_taxon_evidence_tables <- function(...) {
  tables <- list(...)
  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0) {
    return(NULL)
  }

  tables <- lapply(tables, validate_taxon_evidence)
  all_columns <- unique(unlist(lapply(tables, colnames), use.names = FALSE))
  tables <- lapply(tables, function(table) {
    missing_columns <- setdiff(all_columns, colnames(table))
    for (column_name in missing_columns) {
      table[[column_name]] <- NA_character_
    }
    table[, all_columns, drop = FALSE]
  })

  validate_taxon_evidence(.bind_rows_fill(tables))
}

.bind_rows_fill <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- dplyr::bind_rows(rows)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Extract a taxonomy table from supported inputs
#'
#' @param x A data frame, matrix, phyloseq taxonomy table, or phyloseq object.
#'
#' @return A data frame containing taxonomy.
#' @noRd
extract_tax_table <- function(x) {
  if (methods::is(x, "phyloseq")) {
    tax_table <- phyloseq::tax_table(x, errorIfNULL = FALSE)
    if (is.null(tax_table)) {
      stop("`x` does not contain a taxonomy table.", call. = FALSE)
    }
    tax_table <- as.data.frame(
      methods::as(tax_table, "matrix"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (inherits(x, "taxonomyTable")) {
    tax_table <- as.data.frame(
      methods::as(x, "matrix"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (inherits(x, "data.frame")) {
    tax_table <- as.data.frame(
      x,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (is.matrix(x)) {
    tax_table <- as.data.frame(
      x,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    stop(
      "`x` must be a data.frame taxonomy table or a phyloseq/speedyseq object.",
      call. = FALSE
    )
  }

  if (nrow(tax_table) < 1) {
    stop("`x` must contain at least one taxonomic row.", call. = FALSE)
  }
  if (ncol(tax_table) < 1) {
    stop("`x` must contain at least one taxonomy column.", call. = FALSE)
  }

  .add_taxon_id_column(tax_table)
}

#' Choose the most specific available taxon from a taxonomy row
#'
#' @param tax_row A named taxonomy row.
#' @param tax_ranks Character vector of ranks ordered from most specific to
#'   broadest.
#'
#' @return A list with `query_rank` and `query_name`.
#' @noRd
choose_query_taxon <- function(
  tax_row,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum", "kingdom")
) {
  tax_row <- .coerce_tax_row(tax_row)
  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")

  for (rank in tax_ranks) {
    taxon_name <- .tax_rank_value(tax_row, rank)
    if (.is_informative_taxon(taxon_name)) {
      return(list(
        query_rank = tolower(rank),
        query_name = .clean_taxon_name(taxon_name)
      ))
    }
  }

  list(
    query_rank = NA_character_,
    query_name = NA_character_
  )
}

#' Build the flag_taxa prompt
#'
#' @param tax_table Data frame taxonomy table.
#' @inheritParams flag_taxa
#'
#' @return A character scalar prompt.
#' @noRd
build_flag_taxa_prompt <- function(
  tax_table,
  expected_environment,
  expected_habitat = NULL,
  expected_region = NULL,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  evidence_sources = NULL,
  taxon_evidence = NULL,
  tool_evidence = NULL,
  query_table = NULL,
  preliminary_judgement = NULL,
  max_evidence_rows_per_query = 5,
  max_evidence_summary_chars = 160,
  max_tool_evidence_summary_chars = 300,
  max_lineage_chars = 300,
  allow_llm_tools = FALSE,
  prompt_tools = NULL,
  tool_requirement = "required_for_llm"
) {
  if (!inherits(tax_table, "data.frame")) {
    stop("`tax_table` must be a data.frame.", call. = FALSE)
  }

  tax_ranks <- .validate_character_vector(tax_ranks, "tax_ranks")
  evidence_sources <- .validate_flag_taxa_evidence_sources(evidence_sources)
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  tool_evidence <- .validate_flag_taxa_tool_evidence(tool_evidence, query_table = NULL)
  max_evidence_rows_per_query <- .validate_positive_integer_scalar(
    max_evidence_rows_per_query,
    "max_evidence_rows_per_query"
  )
  max_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_evidence_summary_chars,
    "max_evidence_summary_chars"
  )
  max_tool_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_tool_evidence_summary_chars,
    "max_tool_evidence_summary_chars"
  )
  max_lineage_chars <- .validate_positive_integer_scalar(
    max_lineage_chars,
    "max_lineage_chars"
  )
  allow_llm_tools <- .validate_logical_scalar(allow_llm_tools, "allow_llm_tools")
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  if (is.null(query_table)) {
    feature_query_map <- .prepare_flag_taxa_feature_query_map(
      tax_table = tax_table,
      tax_ranks = tax_ranks
    )
    query_table <- .unique_flag_taxa_query_table(feature_query_map)
  } else {
    query_table <- .validate_flag_taxa_query_table(query_table)
  }

  relevant_taxon_evidence <- collect_relevant_taxon_evidence(
    tax_table = query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )

  allowed_values <- .flag_taxa_allowed_values(expected_region = expected_region)
  expected_habitat_text <- .format_optional_prompt_value(expected_habitat)
  expected_region_text <- .format_optional_prompt_value(expected_region)
  expected_region_instruction <- .format_expected_region_instruction(expected_region)
  evidence_text <- .format_flag_taxa_evidence_sources(evidence_sources)
  json_fields <- .flag_taxa_json_output_fields(include_feature_id = FALSE)
  user_evidence_section <- .format_taxon_evidence_prompt_section(
    relevant_taxon_evidence = relevant_taxon_evidence,
    taxon_evidence_provided = !is.null(taxon_evidence)
  )
  preliminary_section <- .format_flag_taxa_preliminary_prompt_section(
    query_table = query_table,
    preliminary_judgement = preliminary_judgement
  )
  tool_evidence_section <- .format_flag_taxa_tool_evidence_prompt_section(
    tool_evidence = tool_evidence,
    query_table = query_table,
    max_tool_evidence_summary_chars = max_tool_evidence_summary_chars
  )
  tool_instruction_section <- .format_flag_taxa_tool_instruction_section(
    allow_llm_tools = allow_llm_tools,
    prompt_tools = prompt_tools,
    tool_requirement = tool_requirement
  )

  prompt_query_table <- query_table[, c("query_id", "query_name", "query_rank", "lineage"), drop = FALSE]
  prompt_query_table$lineage <- .truncate_prompt_value(
    prompt_query_table$lineage,
    max_chars = max_lineage_chars
  )
  query_id_list <- paste(as.character(prompt_query_table$query_id), collapse = ", ")
  if (nrow(prompt_query_table) == 0) {
    query_id_list <- "none; no taxa were selected for LLM/manual judgement under the current judgement_mode"
  }
  selection_instruction <- if (nrow(prompt_query_table) == 0) {
    "No query taxa were selected for LLM/manual judgement under the current judgement_mode. The unique query table below contains only headers."
  } else {
    "Assess only the unique query taxa listed in this prompt; deterministic rows not listed here should remain unchanged."
  }

  paste(
    "Review metabarcoding taxonomy for ecological/geographic compatibility.",
    "",
    "Task:",
    "Assess each unique query taxon against the expected context. The input table is unique taxa, not duplicated feature/ASV rows.",
    selection_instruction,
    "Use query_rank/query_name as the group to assess and lineage for context. Copy query_id exactly, query_name to taxon_name, and query_rank to taxon_rank.",
    paste0("Required query_id values for this prompt: ", query_id_list),
    "Return exactly one JSON object for every listed query_id. Do not omit listed query_id values and do not add extra query_id values.",
    "",
    "Expected context:",
    paste0("- expected_environment: ", expected_environment),
    paste0("- expected_habitat: ", expected_habitat_text),
    paste0("- expected_region: ", expected_region_text),
    paste0("- preferred evidence sources: ", evidence_text),
    preliminary_section,
    user_evidence_section,
    tool_evidence_section,
    tool_instruction_section,
    "",
    "Create a downloadable JSON file named `flag_taxa_llm_result.json`.",
    "Do not just print explanatory text. The file content must be a valid JSON array of objects.",
    "If your interface cannot create a downloadable file, output only the raw JSON array so the user can save it as `flag_taxa_llm_result.json`.",
    "The JSON array must contain one object per query taxon in the same order.",
    "Do not return markdown, a markdown table, a plain text table, code fences, comments, or text outside the JSON file/raw JSON array.",
    "Use JSON null, true, and false; do not use R values such as NA, TRUE, or FALSE.",
    "Avoid nested quotes in string fields.",
    "Do not include feature_id; results are joined back to features after taxon-level assessment. Each object must contain:",
    paste(paste0("- ", json_fields), collapse = "\n"),
    "",
    "Allowed values:",
    .format_allowed_value_prompt_block(
      "expected_environment_status",
      allowed_values$expected_environment_status
    ),
    "",
    .format_allowed_value_prompt_block(
      "expected_habitat_status",
      allowed_values$expected_habitat_status
    ),
    "",
    .format_allowed_value_prompt_block(
      "expected_region_status",
      allowed_values$expected_region_status
    ),
    "",
    .format_allowed_value_prompt_block(
      "ecological_status",
      allowed_values$ecological_status
    ),
    "",
    .format_allowed_value_prompt_block(
      "occurrence_interpretation",
      allowed_values$occurrence_interpretation
    ),
    "",
    .format_allowed_value_prompt_block(
      "recommended_action",
      allowed_values$recommended_action
    ),
    "",
    .format_allowed_value_prompt_block(
      "evidence_basis",
      allowed_values$evidence_basis
    ),
    "",
    "Tool traceability: evidence_basis describes local/tool/conservative support. If a registered tool is used, set llm_tool_used = TRUE and fill llm_tool_name, llm_tool_query, llm_tool_evidence_summary, and llm_tool_references; otherwise use FALSE/NA. If a tool was used but returned no usable references, set llm_tool_references to 'Tool used but no verifiable references were returned.' Do not use vague placeholders such as 'scientific_literature (tool evidence summary)' as the only tool reference.",
    "",
    "Conservative biological rules:",
    "- Do not classify taxa as incompatible based on majority ecology.",
    "- A taxon is environmentally incompatible only when no known member of the queried taxonomic group is compatible with the expected environment.",
    "- If a broad taxonomic rank contains both compatible and incompatible organisms, use mixed_within_rank.",
    "- If taxonomy is too coarse to judge, use insufficient_taxonomic_resolution.",
    "- If evidence is missing, use unknown or no_distribution_evidence. Do not invent evidence.",
    "- Distinguish living habitat incompatibility from transient or allochthonous DNA.",
    "- Use biological judgement based on known ecology, not only occurrence records.",
    "- Missing records alone are uncertainty, not incompatibility; known ecology or defining biology can support unlikely_resident or incompatible_resident.",
    "- Consider whether the sequence could be transported or allochthonous DNA rather than a resident organism.",
    "- Consider whether the sequence is more consistent with contamination or taxonomic misassignment.",
    "- For photosynthetic taxa in aphotic deep-sea benthic habitats, resident status is incompatible or unlikely, but allochthonous DNA may be plausible.",
    "- For taxa with marine evidence and an expected marine habitat subtype, do not call habitat incompatible just because exact habitat records are missing.",
    "- For broad ranks, avoid hard exclude unless the whole group's defining biology is incompatible or the evidence is very strong and rank-appropriate.",
    "- Use no_known_habitat_evidence when available evidence does not document the taxon in the expected habitat but does not demonstrate biological incompatibility.",
    "- Use exclude only for obvious incompatibilities supported by explicit positive evidence.",
    "- Include concise references that justify categorisations.",
    "- Missing evidence is not evidence of incompatibility; use unknown and review when evidence is insufficient.",
    "- The preliminary evidence-first judgement is the baseline. Do not override review to exclude unless explicit positive incompatibility evidence exists.",
    "- Absence of evidence from tools is not a reason to make the call stricter than the preliminary judgement.",
    "- For broad taxa, prefer cautious statuses such as mixed_within_rank, insufficient_taxonomic_resolution, no_known_habitat_evidence, possible_likely, possible_unlikely, unknown, and review.",
    "- Failure to find explicit evidence for a taxon in a habitat or region is not evidence that the taxon is incompatible with that habitat or region.",
    "- Never infer incompatibility from lack of records alone. 'No evidence found' means uncertainty, not exclusion.",
    "- Before assigning recommended_action = \"exclude\", identify explicit evidence that demonstrates incompatibility. If you cannot identify explicit incompatibility evidence, use flag_for_review or flag_possible_exclusion.",
    "- Do not use expected_habitat_status = \"incompatible\" or recommended_action = \"exclude\" only because no deep-sea records, no CCZ records, no Scite evidence, or no tool evidence were found.",
    "- If you assign expected_habitat_status = \"incompatible\", ecological_status = \"unlikely_resident\" or \"incompatible_resident\", occurrence_interpretation indicating allochthony, contamination, or misassignment, or recommended_action = \"flag_possible_exclusion\" or \"exclude\", your rationale must explain the specific biological/ecological reason. Do not merely repeat status values.",
    "- Absence of records alone is not incompatibility. If the issue is only missing habitat or region records, use no_known_habitat_evidence/no_distribution_evidence and flag_for_review rather than incompatible or exclude.",
    "- Known ecology can be positive evidence. For example, photosynthetic taxa in aphotic deep-sea benthic habitats may be incompatible as residents, although their DNA may be possible allochthonous material.",
    "- If WoRMS or another database includes the expected environment among mixed flags, and habitat/region records are missing, use mixed_within_rank or compatible for environment as appropriate, no_known_habitat_evidence or unknown for habitat, no_distribution_evidence or unknown for region, and review.",
    "- Never use recommended_action = \"exclude\" with evidence_basis = \"conservative_reasoning_only\".",
    "- Do not exclude a taxon based only on general model knowledge.",
    "",
    "If expected_habitat is not specified, use unknown for expected_habitat_status unless evidence supports a more specific interpretation.",
    expected_region_instruction,
    "Write concise rationales and include references for each row. Do not add free-text commentary outside the JSON.",
    "",
    "Unique query taxa, tab-separated:",
    .format_prompt_table(prompt_query_table),
    sep = "\n"
  )
}

.validate_flag_taxa_query_table <- function(query_table) {
  if (!inherits(query_table, "data.frame")) {
    stop("`query_table` must be a data.frame.", call. = FALSE)
  }
  required_columns <- c("query_name", "query_rank", "lineage")
  missing_columns <- setdiff(required_columns, colnames(query_table))
  if (length(missing_columns) > 0) {
    stop(
      "`query_table` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  out <- query_table[, unique(c("query_id", required_columns)), drop = FALSE]
  if (!("query_id" %in% colnames(out))) {
    out <- cbind(
      data.frame(
        query_id = .flag_taxa_query_ids(nrow(out)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      out,
      stringsAsFactors = FALSE
    )
  }
  out
}

.format_allowed_value_prompt_block <- function(name, values) {
  paste(
    c(
      paste0(name, ":"),
      paste0("- ", values)
    ),
    collapse = "\n"
  )
}

.format_expected_region_instruction <- function(expected_region) {
  if (is.null(expected_region)) {
    return("If expected_region is not specified, use not_assessed for expected_region_status.")
  }

  paste(
    "Because expected_region is provided, do not use not_assessed.",
    "If evidence is insufficient, use unknown or no_distribution_evidence as appropriate."
  )
}

.format_flag_taxa_tool_instruction_section <- function(
  allow_llm_tools,
  prompt_tools = NULL,
  tool_requirement = "required_for_llm"
) {
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  if (isTRUE(allow_llm_tools)) {
    reference_policy <- paste(
      "Do not invent AphiaIDs, DOIs, URLs, article titles, author names, or citations.",
      "Use only references provided in the local evidence table or returned by tools.",
      "If using general biological knowledge, state it in the rationale without fabricating a reference."
    )
    if (identical(tool_requirement, "required_for_llm")) {
      named_tool_instruction <- if (length(prompt_tools) > 0) {
        tool_label <- .format_prompt_tool_names(prompt_tools)
        paste0(
          "Use the ",
          tool_label,
          " tool",
          if (length(prompt_tools) == 1) "" else "s",
          " for every query_id in this prompt."
        )
      } else {
        "Use the registered literature/search tool available through the chat interface for every query_id in this prompt."
      }
      return(paste(
        "",
        "Required registered tool use",
        "Registered tools may be available through the chat object, for example Scite or another literature/search tool.",
        "Tool use is mandatory for every query_id in this prompt. Use the registered tool(s), such as Scite if available, before returning the final JSON. For each query_id, report the tool query, tool evidence summary, and tool references. If no evidence is found, say so explicitly. Do not invent references, AphiaIDs, DOIs, URLs, or citations.",
        named_tool_instruction,
        "For each query_id, fill llm_tool_used, llm_tool_name, llm_tool_query, llm_tool_evidence_summary, and llm_tool_references.",
        "If the tool returns no explicit evidence, set llm_tool_used = TRUE, report the query, set llm_tool_evidence_summary to a no-evidence statement, and set llm_tool_references to 'Tool used but no verifiable references were returned.'",
        "A tool search that finds no explicit evidence is missing evidence, not incompatibility evidence. Do not make the final interpretation stricter because a tool found no records.",
        reference_policy,
        "Use the full compact lineage to build precise search queries and disambiguate taxa. Search terms should include the query taxon name, relevant higher taxonomy, and expected context such as marine, deep sea, Clarion-Clipperton Zone, or the user-provided expected region/habitat/environment.",
        "Summarise any tool-derived evidence in rationale, references, and the tool-use fields.",
        "Tool use must happen before the final structured result is returned. The final response must still be valid JSON matching the requested schema.",
        sep = "\n"
      ))
    }
    named_tool_instruction <- if (length(prompt_tools) > 0) {
      tool_label <- .format_prompt_tool_names(prompt_tools)
      paste(
        paste0("Tool use requested: use the ", tool_label, " tool", if (length(prompt_tools) == 1) "" else "s", " for each query taxon listed in this prompt unless local evidence already clearly resolves the taxon."),
        "For each tool search, report the exact search query, whether explicit evidence was found, and the references used.",
        "When a listed query taxon remains selected for manual/LLM judgement under judgement_mode, treat it as eligible for the named tool search; deterministic retain rows are not included here unless judgement_mode = \"llm_all\".",
        "If Scite is one of the named tools, use it to find explicit literature evidence and report verifiable article titles, short citations, DOIs, URLs, or state that no verifiable references were returned.",
        sep = "\n"
      )
    } else {
      ""
    }
    return(paste(
      "",
      "Optional registered tool use",
      "Registered tools may be available through the chat object, for example Scite or another literature/search tool. Use local evidence first. Use registered tools selectively, especially when the decision would otherwise be flag_for_review, flag_possible_exclusion, or exclude.",
      named_tool_instruction,
      "Do not use tools for taxa where local evidence clearly supports retain.",
      "Use registered tools only for taxa selected by judgement_mode: missing-evidence taxa for llm_missing_evidence, flag_for_review and flag_possible_exclusion taxa for llm_flagged, flag_possible_exclusion taxa for llm_possible_exclusion, or all taxa for llm_all.",
      "Use registered tools for selected taxa that would otherwise be flag_for_review or flag_possible_exclusion because local evidence is missing, weak, ambiguous, broad, contradictory, or indicates a possible resident-community concern.",
      "Before assigning recommended_action = \"exclude\", use registered tools to check for explicit supporting evidence unless the provided local evidence already clearly demonstrates incompatibility.",
      "Do not exclude a taxon based only on general model knowledge. If no explicit local or tool-derived evidence supports exclusion, use flag_for_review or flag_possible_exclusion instead.",
      "A tool search that finds no explicit evidence is missing evidence, not incompatibility evidence. No deep-sea records, no CCZ records, no Scite evidence, or no known habitat records should usually produce unknown, no_known_habitat_evidence, no_distribution_evidence, and flag_for_review rather than exclude.",
      "If tool evidence does not resolve the case, keep flag_for_review or flag_possible_exclusion as appropriate.",
      "Do not use tool evidence to overrule broad-rank uncertainty unless the evidence is appropriate to the query rank. For broad ranks, use mixed_within_rank or insufficient_taxonomic_resolution when appropriate.",
      "Use the full compact lineage to build precise search queries and disambiguate taxa. Search terms should include the query taxon name, relevant higher taxonomy, and expected context such as marine, deep sea, Clarion-Clipperton Zone, or the user-provided expected region/habitat/environment.",
      "Use tools mainly for flagged and possible-exclusion boundary cases. Do not spend tool calls confirming easy retains with clear local evidence. Use tools to try to resolve cases that would otherwise remain flag_for_review or flag_possible_exclusion, and to verify any potential exclude decision. If tool evidence does not clearly resolve the case, keep flag_for_review or flag_possible_exclusion as appropriate.",
      "If registered tools are unavailable or return no useful evidence, report llm_tool_used = FALSE or TRUE with a no_evidence_found summary as appropriate, and use unknown/flag_for_review conservatively unless local evidence clearly supports another action.",
      reference_policy,
      "Summarise any tool-derived evidence in rationale, references, and the tool-use fields.",
      "Tool use must happen before the final structured result is returned. The final response must still be valid JSON matching the requested schema.",
      sep = "\n"
    ))
  }

  paste(
    "",
    "External tools and web search",
    "Do not assume access to external tools or web search. Use only the evidence provided in this prompt and conservative biological reasoning. Missing evidence is not evidence of incompatibility.",
    sep = "\n"
  )
}

.format_prompt_tool_names <- function(prompt_tools) {
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  if (length(prompt_tools) == 0) {
    return("")
  }
  if (length(prompt_tools) == 1) {
    return(prompt_tools[[1]])
  }
  paste(prompt_tools, collapse = ", ")
}

.prepare_flag_taxa_taxonomy <- function(tax_table, tax_ranks) {
  reserved_columns <- .flag_taxa_output_columns()
  collisions <- colnames(tax_table)[tolower(colnames(tax_table)) %in% reserved_columns]
  if (length(collisions) > 0) {
    stop(
      "`x` already contains columns reserved by `flag_taxa()`: ",
      paste(collisions, collapse = ", "),
      call. = FALSE
    )
  }

  matched_ranks <- tolower(tax_ranks) %in% tolower(colnames(tax_table))
  if (!any(matched_ranks)) {
    stop(
      "`x` must contain at least one requested taxonomic rank column. ",
      "Expected one of: ",
      paste(tax_ranks, collapse = ", "),
      call. = FALSE
    )
  }

  query_taxa <- lapply(seq_len(nrow(tax_table)), function(row_index) {
    choose_query_taxon(tax_table[row_index, , drop = FALSE], tax_ranks = tax_ranks)
  })

  tax_table$query_rank <- vapply(
    query_taxa,
    function(query_taxon) query_taxon$query_rank,
    character(1)
  )
  tax_table$query_name <- vapply(
    query_taxa,
    function(query_taxon) query_taxon$query_name,
    character(1)
  )
  tax_table$lineage <- vapply(
    seq_len(nrow(tax_table)),
    function(row_index) .build_lineage(tax_table[row_index, , drop = FALSE], tax_ranks),
    character(1)
  )

  tax_table
}

collect_relevant_taxon_evidence <- function(
  tax_table,
  taxon_evidence,
  tax_ranks = c("species", "genus", "family", "order", "class", "phylum"),
  max_evidence_rows_per_query = 5,
  max_evidence_summary_chars = 160
) {
  taxon_evidence <- validate_taxon_evidence(taxon_evidence)
  if (is.null(taxon_evidence)) {
    return(NULL)
  }

  if (!inherits(tax_table, "data.frame")) {
    stop("`tax_table` must be a data.frame.", call. = FALSE)
  }
  max_evidence_rows_per_query <- .validate_positive_integer_scalar(
    max_evidence_rows_per_query,
    "max_evidence_rows_per_query"
  )
  max_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_evidence_summary_chars,
    "max_evidence_summary_chars"
  )

  if (!all(c("query_rank", "query_name") %in% colnames(tax_table))) {
    tax_table <- .unique_flag_taxa_query_table(
      .prepare_flag_taxa_feature_query_map(tax_table, tax_ranks)
    )
  } else if (!("lineage" %in% colnames(tax_table))) {
    tax_table$lineage <- NA_character_
  }

  if (nrow(taxon_evidence) == 0 || nrow(tax_table) == 0) {
    return(.empty_relevant_taxon_evidence_table(tax_table))
  }

  rows <- lapply(seq_len(nrow(tax_table)), function(row_index) {
    query_name <- tax_table$query_name[[row_index]]
    if (!.is_informative_taxon(query_name)) {
      return(NULL)
    }

    matched <- .taxon_evidence_match_indices(
      tax_table = tax_table,
      taxon_evidence = taxon_evidence,
      row_index = row_index
    )
    if (length(matched) == 0) {
      return(NULL)
    }

    evidence_rows <- taxon_evidence[matched, , drop = FALSE]
    if (nrow(evidence_rows) > max_evidence_rows_per_query) {
      evidence_rows <- .summarise_relevant_taxon_evidence(
        evidence_rows = evidence_rows,
        query_name = query_name,
        query_rank = tax_table$query_rank[[row_index]]
      )
    }
    context <- .taxon_evidence_match_context(tax_table, row_index)
    context <- context[rep(1, nrow(evidence_rows)), , drop = FALSE]
    cbind(context, evidence_rows)
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(.empty_relevant_taxon_evidence_table(tax_table))
  }

  out <- .bind_rows_fill(rows)
  rownames(out) <- NULL
  .select_relevant_taxon_evidence_columns(
    out,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
}

.taxon_evidence_match_indices <- function(tax_table, taxon_evidence, row_index) {
  query_name <- tax_table$query_name[[row_index]]
  query_rank <- tolower(as.character(tax_table$query_rank[[row_index]]))
  query_key <- .taxon_match_key(query_name)
  matched <- integer()

  name_keys <- .taxon_match_key(taxon_evidence$taxon_name)
  matched <- c(matched, which(name_keys == query_key))

  if ("accepted_name" %in% colnames(taxon_evidence)) {
    accepted_keys <- .taxon_match_key(taxon_evidence$accepted_name)
    matched <- c(matched, which(accepted_keys == query_key))
  }

  lineage_columns <- c(
    phylum = "phylum",
    class = "class",
    order = "order",
    family = "family",
    genus = "genus"
  )
  lineage_column <- NULL
  if (query_rank %in% names(lineage_columns)) {
    lineage_column <- lineage_columns[[query_rank]]
  }
  if (!is.null(lineage_column) && lineage_column %in% colnames(taxon_evidence)) {
    lineage_keys <- .taxon_match_key(taxon_evidence[[lineage_column]])
    matched <- c(matched, which(lineage_keys == query_key))
  }

  query_ids <- .taxon_evidence_query_ids(tax_table[row_index, , drop = FALSE])
  if (length(query_ids) > 0) {
    evidence_ids <- .taxon_evidence_row_ids(taxon_evidence)
    matched <- c(matched, which(evidence_ids %in% query_ids))
  }

  unique(matched)
}

.taxon_match_key <- function(taxon_name) {
  tolower(trimws(as.character(taxon_name)))
}

.taxon_evidence_query_ids <- function(tax_row) {
  id_columns <- intersect(
    c(
      "source_taxon_id",
      "source_record_id",
      "AphiaID",
      "aphia_id",
      "aphiaID",
      "taxonID",
      "taxon_id"
    ),
    colnames(tax_row)
  )
  if (length(id_columns) == 0) {
    return(character())
  }

  ids <- unlist(tax_row[id_columns], use.names = FALSE)
  ids <- .taxon_id_match_key(ids)
  unique(ids[.is_non_empty_value(ids)])
}

.taxon_evidence_row_ids <- function(taxon_evidence) {
  id_columns <- intersect(c("source_taxon_id", "source_record_id"), colnames(taxon_evidence))
  if (length(id_columns) == 0) {
    return(rep(NA_character_, nrow(taxon_evidence)))
  }

  ids <- .row_first_non_empty(taxon_evidence, id_columns)
  .taxon_id_match_key(ids)
}

.taxon_id_match_key <- function(value) {
  tolower(trimws(as.character(value)))
}

.summarise_relevant_taxon_evidence <- function(evidence_rows, query_name, query_rank) {
  .summarise_relevant_taxon_evidence_group(
    group = evidence_rows,
    query_name = query_name,
    query_rank = query_rank
  )
}

.summarise_relevant_taxon_evidence_group <- function(group, query_name, query_rank) {
  out <- as.data.frame(
    stats::setNames(
      replicate(
        length(.taxon_evidence_prompt_columns()),
        NA_character_,
        simplify = FALSE
      ),
      .taxon_evidence_prompt_columns()
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  examples <- .collapse_distinct_values(.optional_evidence_column(group, "taxon_name"), max_values = 5)
  environments <- .collapse_distinct_values(.optional_evidence_column(group, "environment"), max_values = 8)
  habitats <- .collapse_distinct_values(.optional_evidence_column(group, "habitat"), max_values = 5)
  regions <- .collapse_distinct_values(.optional_evidence_column(group, "region"), max_values = 5)
  sources <- .collapse_distinct_values(.optional_evidence_column(group, "source"), max_values = 5)
  evidence_types <- .collapse_distinct_values(.optional_evidence_column(group, "evidence_type"), max_values = 5)
  references <- .collapse_distinct_values(.optional_evidence_column(group, "reference"), max_values = 3)

  out$taxon_name <- paste0("Summarised matches for ", tolower(query_rank), ": ", query_name)
  out$taxon_rank <- .collapse_distinct_values(group$taxon_rank, max_values = 5)
  out$source <- sources
  out$evidence_type <- evidence_types
  out$environment <- environments
  out$habitat <- habitats
  out$region <- regions
  out$evidence_summary <- paste0(
    "Summarised evidence: ",
    nrow(group),
    " matching evidence row(s) for ",
    tolower(query_rank),
    " '",
    query_name,
    "'. Examples: ",
    .format_summary_value(examples),
    ". Distinct environment(s): ",
    .format_summary_value(environments),
    ". Distinct habitat(s): ",
    .format_summary_value(habitats),
    ". Distinct region(s): ",
    .format_summary_value(regions),
    ". Source/evidence type(s): ",
    .format_summary_value(sources),
    " / ",
    .format_summary_value(evidence_types),
    ". This is summarised evidence; not all matching rows are shown."
  )
  out$reference <- references

  out
}

.optional_evidence_column <- function(evidence_rows, column_name) {
  if (column_name %in% colnames(evidence_rows)) {
    return(evidence_rows[[column_name]])
  }

  rep(NA_character_, nrow(evidence_rows))
}

.collapse_distinct_values <- function(value, max_values = 5) {
  value <- as.character(value)
  value <- unlist(strsplit(value, "\\s*;\\s*"), use.names = FALSE)
  value <- trimws(value)
  value <- unique(value[.is_non_empty_value(value)])
  if (length(value) == 0) {
    return(NA_character_)
  }
  if (length(value) <= max_values) {
    return(paste(value, collapse = "; "))
  }

  paste0(
    paste(utils::head(value, max_values), collapse = "; "),
    "; ... (",
    length(value),
    " total)"
  )
}

.first_non_empty_value <- function(value) {
  value <- as.character(value)
  value <- value[.is_non_empty_value(value)]
  if (length(value) == 0) {
    return(NA_character_)
  }

  value[[1]]
}

.format_summary_value <- function(value) {
  if (!.is_non_empty_value(value)) {
    return("not available")
  }

  value
}

.taxon_evidence_match_context <- function(tax_table, row_index) {
  out <- data.frame(
    taxonomy_row = row_index,
    query_id = if ("query_id" %in% colnames(tax_table)) {
      as.character(tax_table$query_id[[row_index]])
    } else {
      NA_character_
    },
    query_rank = tax_table$query_rank[[row_index]],
    query_name = tax_table$query_name[[row_index]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if ("feature_id" %in% colnames(tax_table)) {
    out$feature_id <- as.character(tax_table$feature_id[[row_index]])
  } else if ("taxon_id" %in% colnames(tax_table)) {
    out$taxon_id <- as.character(tax_table$taxon_id[[row_index]])
  }

  out
}

.empty_relevant_taxon_evidence_table <- function(tax_table) {
  context_columns <- c("taxonomy_row", "query_id", "query_rank", "query_name")
  if ("feature_id" %in% colnames(tax_table)) {
    context_columns <- c(context_columns, "feature_id")
  } else if ("taxon_id" %in% colnames(tax_table)) {
    context_columns <- c(context_columns, "taxon_id")
  }
  out <- as.data.frame(
    stats::setNames(
      replicate(length(context_columns), character(), simplify = FALSE),
      context_columns
    ),
    stringsAsFactors = FALSE
  )
  for (column_name in .taxon_evidence_prompt_columns()) {
    out[[column_name]] <- character()
  }

  out
}

.select_relevant_taxon_evidence_columns <- function(
  taxon_evidence,
  max_evidence_summary_chars = 160
) {
  for (column_name in .taxon_evidence_prompt_columns()) {
    if (!(column_name %in% colnames(taxon_evidence))) {
      taxon_evidence[[column_name]] <- NA_character_
    }
  }

  context_columns <- intersect(
    c("taxonomy_row", "feature_id", "taxon_id", "query_id", "query_rank", "query_name"),
    colnames(taxon_evidence)
  )
  out <- taxon_evidence[
    ,
    c(context_columns, .taxon_evidence_prompt_columns()),
    drop = FALSE
  ]
  .compact_taxon_evidence_prompt_values(
    out,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
}

.flag_taxa_output_columns <- function() {
  c(
    "query_rank",
    "query_name",
    "lineage",
    "expected_environment",
    "expected_environment_status",
    "expected_habitat",
    "expected_habitat_status",
    "expected_region",
    "expected_region_status",
    "recommended_action",
    "rationale",
    "references"
  )
}

.flag_taxa_json_output_fields <- function(include_feature_id) {
  .flag_taxa_structured_output_columns(include_feature_id)
}

.flag_taxa_legacy_output_columns <- function(include_feature_id) {
  fields <- .flag_taxa_output_columns()
  if (isTRUE(include_feature_id)) {
    fields <- c("feature_id", fields)
  }

  fields
}

.flag_taxa_structured_output_columns <- function(
  include_feature_id,
  include_query_id = TRUE
) {
  fields <- c(
    "query_id",
    "taxon_name",
    "taxon_rank",
    "expected_environment_status",
    "expected_habitat_status",
    "expected_region_status",
    "ecological_status",
    "occurrence_interpretation",
    "recommended_action",
    "rationale",
    "references",
    "evidence_basis",
    "llm_tool_used",
    "llm_tool_name",
    "llm_tool_query",
    "llm_tool_evidence_summary",
    "llm_tool_references"
  )
  if (!isTRUE(include_query_id)) {
    fields <- setdiff(fields, "query_id")
  }
  if (isTRUE(include_feature_id)) {
    fields <- c("feature_id", fields)
  }

  fields
}

.flag_taxa_taxon_result_columns <- function(result) {
  with_query_id <- .flag_taxa_structured_output_columns(
    include_feature_id = FALSE,
    include_query_id = TRUE
  )
  without_query_id <- .flag_taxa_structured_output_columns(
    include_feature_id = FALSE,
    include_query_id = FALSE
  )
  if (all(with_query_id %in% colnames(result))) {
    return(with_query_id)
  }
  without_query_id
}

.taxon_evidence_required_columns <- function() {
  c(
    "taxon_name",
    "taxon_rank",
    "source",
    "evidence_type",
    "evidence_summary",
    "reference"
  )
}

.taxon_evidence_optional_columns <- function() {
  c(
    "accepted_name",
    "accepted_rank",
    "source_taxon_id",
    "source_record_id",
    "environment",
    "habitat",
    "region",
    "locality",
    "decimal_latitude",
    "decimal_longitude",
    "basis_of_record",
    "occurrence_count",
    "reference_url",
    "doi",
    "checked_at",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus"
  )
}

.taxon_evidence_prompt_columns <- function() {
  c(
    "taxon_name",
    "taxon_rank",
    "source",
    "evidence_type",
    "environment",
    "habitat",
    "region",
    "evidence_summary",
    "reference"
  )
}

.compact_taxon_evidence_prompt_values <- function(
  taxon_evidence,
  max_evidence_summary_chars = 160
) {
  max_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_evidence_summary_chars,
    "max_evidence_summary_chars"
  )
  if ("evidence_summary" %in% colnames(taxon_evidence)) {
    taxon_evidence$evidence_summary <- .truncate_prompt_value(
      taxon_evidence$evidence_summary,
      max_chars = max_evidence_summary_chars
    )
  }
  if ("reference" %in% colnames(taxon_evidence)) {
    taxon_evidence$reference <- .truncate_prompt_value(
      taxon_evidence$reference,
      max_chars = 120
    )
  }

  taxon_evidence
}

.truncate_prompt_value <- function(value, max_chars) {
  value <- as.character(value)
  too_long <- !is.na(value) & nchar(value, type = "chars") > max_chars
  value[too_long] <- paste0(
    substr(value[too_long], 1, max_chars),
    "... [truncated]"
  )
  value
}

.taxon_evidence_lineage_columns <- function() {
  c("kingdom", "phylum", "class", "order", "family", "genus")
}

.is_entirely_empty_column <- function(column) {
  values <- as.character(column)
  all(is.na(values) | !nzchar(trimws(values)))
}

.format_taxon_evidence_prompt_section <- function(
  relevant_taxon_evidence,
  taxon_evidence_provided
) {
  if (!isTRUE(taxon_evidence_provided)) {
    return("")
  }

  evidence_body <- if (
    is.null(relevant_taxon_evidence) ||
      nrow(relevant_taxon_evidence) == 0
  ) {
    "No rows in `taxon_evidence` matched the query taxa in the input taxonomy table by taxon_name, accepted_name, taxon/source ID, or supported lineage columns."
  } else {
    paste(
      "Compact relevant evidence table, tab-separated. Large lineage matches may be summarised rather than shown row-by-row:",
      .format_prompt_table(relevant_taxon_evidence),
      sep = "\n"
    )
  }

  paste(
    "",
    "User-provided taxon evidence",
    "Only evidence rows relevant to the input query taxa are included below.",
    "The prompt includes compact evidence fields only; full references, URLs, identifiers, coordinates, and DOIs remain in the evidence object but are not pasted here by default.",
    "Use this evidence before performing online searches.",
    "Do not search online for taxa where the provided evidence is sufficient to assign the requested statuses.",
    "Use online sources only when the provided evidence is incomplete, missing, or conflicting.",
    "Treat user-provided taxon evidence, including WoRMS evidence, as supporting evidence rather than an automatic decision.",
    "Treat context-specific regional or habitat evidence, such as source = worms_ccz, as distinct from generic WoRMS evidence; do not replace it with less-specific rows.",
    "Do not treat generic WoRMS evidence as more authoritative than regional or deep-sea checklist evidence when assessing matching regions or habitats.",
    "A no-match, lookup failure, unknown, empty environment, or missing WoRMS evidence row is missing evidence, not evidence of incompatibility.",
    evidence_body,
    sep = "\n"
  )
}

.format_flag_taxa_preliminary_prompt_section <- function(
  query_table,
  preliminary_judgement
) {
  if (is.null(preliminary_judgement) || nrow(query_table) == 0) {
    return("")
  }
  query_ids <- as.character(query_table$query_id)
  rows <- preliminary_judgement[
    as.character(preliminary_judgement$query_id) %in% query_ids,
    ,
    drop = FALSE
  ]
  if (nrow(rows) == 0) {
    return("")
  }
  rows <- rows[match(query_ids, as.character(rows$query_id)), , drop = FALSE]
  prompt_rows <- rows[
    ,
    c(
      "query_id",
      "query_name",
      "query_rank",
      "prelim_environment_status",
      "prelim_habitat_status",
      "prelim_region_status",
      "prelim_ecological_status",
      "prelim_occurrence_interpretation",
      "prelim_recommended_action",
      "prelim_evidence_basis",
      "prelim_rationale"
    ),
    drop = FALSE
  ]
  prompt_rows$prelim_rationale <- .truncate_prompt_value(
    prompt_rows$prelim_rationale,
    max_chars = 260
  )

  paste(
    "",
    "Preliminary evidence-first judgement",
    "The table below is a deterministic local-evidence baseline. Treat it as the starting point for judgement.",
    "Do not make a preliminary review or retain stricter unless explicit positive incompatibility evidence is available.",
    "Absence of tool evidence, no records, no known habitat evidence, or no regional record should not make the final call stricter by itself.",
    "Preliminary judgement table, tab-separated:",
    .format_prompt_table(prompt_rows),
    sep = "\n"
  )
}

.validate_flag_taxa_tool_evidence <- function(tool_evidence, query_table = NULL) {
  if (is.null(tool_evidence)) {
    return(NULL)
  }
  if (!inherits(tool_evidence, "data.frame")) {
    stop("`tool_evidence` must be a data.frame.", call. = FALSE)
  }
  required_columns <- .flag_taxa_tool_evidence_columns()
  missing_columns <- setdiff(required_columns, colnames(tool_evidence))
  if (length(missing_columns) > 0) {
    stop(
      "`tool_evidence` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  tool_evidence <- as.data.frame(
    tool_evidence,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  tool_evidence$tool_used <- as.logical(tool_evidence$tool_used)
  tool_evidence$tool_found_explicit_evidence <- as.logical(tool_evidence$tool_found_explicit_evidence)
  if (any(is.na(tool_evidence$tool_used))) {
    stop("`tool_evidence$tool_used` must contain only TRUE or FALSE values.", call. = FALSE)
  }
  if (any(is.na(tool_evidence$tool_found_explicit_evidence))) {
    stop("`tool_evidence$tool_found_explicit_evidence` must contain only TRUE or FALSE values.", call. = FALSE)
  }
  if (!is.null(query_table) && nrow(tool_evidence) > 0) {
    .validate_flag_taxa_result_query_ids(
      result_query_id = tool_evidence$query_id,
      expected_query_id = query_table$query_id
    )
  }
  tool_evidence[, required_columns, drop = FALSE]
}

.flag_taxa_tool_evidence_columns <- function() {
  c(
    "query_id",
    "query_name",
    "query_rank",
    "tool_used",
    "tool_name",
    "tool_query",
    "tool_evidence_summary",
    "tool_references",
    "tool_found_explicit_evidence"
  )
}

.format_flag_taxa_tool_evidence_prompt_section <- function(
  tool_evidence,
  query_table,
  max_tool_evidence_summary_chars = 300
) {
  max_tool_evidence_summary_chars <- .validate_positive_integer_scalar(
    max_tool_evidence_summary_chars,
    "max_tool_evidence_summary_chars"
  )
  tool_evidence <- .filter_flag_taxa_tool_evidence_to_query(
    tool_evidence = tool_evidence,
    query_table = query_table
  )
  if (is.null(tool_evidence) || nrow(tool_evidence) == 0) {
    return("")
  }

  prompt_columns <- c(
    "query_id",
    "query_name",
    "query_rank",
    "tool_used",
    "tool_name",
    "tool_query",
    "tool_evidence_summary",
    "tool_references",
    "tool_found_explicit_evidence"
  )
  tool_evidence <- tool_evidence[, prompt_columns, drop = FALSE]
  tool_evidence$tool_evidence_summary <- .truncate_prompt_value(
    tool_evidence$tool_evidence_summary,
    max_chars = max_tool_evidence_summary_chars
  )
  tool_evidence$tool_references <- .truncate_prompt_value(
    tool_evidence$tool_references,
    max_chars = min(160L, max_tool_evidence_summary_chars)
  )

  paste(
    "",
    "Tool-derived evidence",
    "Tool-derived evidence may have been gathered in a prior tool evidence pass. Use it when available.",
    "A tool-derived no_evidence_found or no explicit evidence summary is missing evidence, not incompatibility evidence, and must not by itself justify expected_habitat_status = incompatible or recommended_action = exclude.",
    "During the final structured-output step, you may use registered tools only if necessary, but the primary tool evidence should come from the provided tool-derived evidence section.",
    "Compact tool evidence table, tab-separated:",
    .format_prompt_table(tool_evidence),
    sep = "\n"
  )
}

.filter_flag_taxa_tool_evidence_to_query <- function(tool_evidence, query_table) {
  tool_evidence <- .validate_flag_taxa_tool_evidence(tool_evidence)
  if (is.null(tool_evidence) || nrow(tool_evidence) == 0) {
    return(tool_evidence)
  }
  if (is.null(query_table) || !("query_id" %in% colnames(query_table))) {
    return(tool_evidence)
  }
  out <- tool_evidence[as.character(tool_evidence$query_id) %in% as.character(query_table$query_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.should_run_flag_taxa_tool_pass <- function(
  allow_llm_tools,
  tool_use_policy,
  prompt_only,
  run_tools_in_prompt_only
) {
  isTRUE(allow_llm_tools) &&
    !identical(tool_use_policy, "never") &&
    (!isTRUE(prompt_only) || isTRUE(run_tools_in_prompt_only))
}

.run_flag_taxa_tool_evidence_pass <- function(
  unique_query_table,
  preliminary_judgement = NULL,
  judgement_mode = "llm_all",
  taxon_evidence,
  expected_environment,
  expected_habitat,
  expected_region,
  tax_ranks,
  max_evidence_rows_per_query,
  allow_llm_tools,
  tool_use_policy,
  max_tool_taxa,
  tool_batch_size,
  chat,
  llm_max_tries,
  llm_retry_sleep,
  verbose,
  start_time
) {
  candidates <- .select_flag_taxa_tool_candidates(
    unique_query_table = unique_query_table,
    preliminary_judgement = preliminary_judgement,
    judgement_mode = judgement_mode,
    taxon_evidence = taxon_evidence,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    tool_use_policy = tool_use_policy
  )
  candidate_count <- nrow(candidates)
  if (!is.infinite(max_tool_taxa) && candidate_count > max_tool_taxa) {
    candidates <- candidates[seq_len(max_tool_taxa), , drop = FALSE]
  }
  skipped <- candidate_count - nrow(candidates)
  .flag_taxa_progress(
    verbose,
    start_time,
    "selected ",
    nrow(candidates),
    " / ",
    nrow(unique_query_table),
    " unique taxa for tool evidence using policy '",
    tool_use_policy,
    "'",
    if (skipped > 0) paste0("; skipped ", skipped, " due to max_tool_taxa") else "",
    "."
  )
  if (nrow(candidates) == 0) {
    return(.empty_flag_taxa_tool_evidence())
  }

  batches <- split(
    candidates,
    ceiling(seq_len(nrow(candidates)) / tool_batch_size)
  )
  rows <- vector("list", length(batches))
  for (batch_index in seq_along(batches)) {
    batch <- batches[[batch_index]]
    .flag_taxa_progress(
      verbose,
      start_time,
      "running tool evidence batch ",
      batch_index,
      " / ",
      length(batches),
      " (",
      nrow(batch),
      " taxa)."
    )
    rows[[batch_index]] <- tryCatch(
      call_flag_taxa_tool_evidence(
        query_table = batch,
        expected_environment = expected_environment,
        expected_habitat = expected_habitat,
        expected_region = expected_region,
        chat = chat,
        max_tries = llm_max_tries,
        retry_sleep = llm_retry_sleep
      ),
      error = function(error) {
        stop(
          "Tool evidence pass failed for batch ",
          batch_index,
          " / ",
          length(batches),
          ": ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    rows[[batch_index]] <- .validate_flag_taxa_tool_evidence(
      rows[[batch_index]],
      query_table = batch
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed tool evidence batch ",
      batch_index,
      " / ",
      length(batches),
      "."
    )
  }

  .validate_flag_taxa_tool_evidence(.bind_rows_fill(rows))
}

.select_flag_taxa_tool_candidates <- function(
  unique_query_table,
  preliminary_judgement = NULL,
  judgement_mode = "llm_all",
  taxon_evidence,
  expected_environment,
  expected_habitat,
  expected_region,
  tax_ranks,
  max_evidence_rows_per_query,
  tool_use_policy
) {
  if (identical(tool_use_policy, "never")) {
    return(unique_query_table[0, , drop = FALSE])
  }
  if (identical(judgement_mode, "evidence_only")) {
    return(unique_query_table[0, , drop = FALSE])
  }
  if (!is.null(preliminary_judgement) && nrow(preliminary_judgement) > 0) {
    selected <- .flag_taxa_judgement_mode_selection(
      preliminary_judgement = preliminary_judgement,
      judgement_mode = judgement_mode
    )
    selected_ids <- as.character(preliminary_judgement$query_id[selected])
    out <- unique_query_table[
      as.character(unique_query_table$query_id) %in% selected_ids,
      ,
      drop = FALSE
    ]
    rownames(out) <- NULL
    return(out)
  }

  if (judgement_mode %in% c("llm_missing_evidence", "llm_flagged", "llm_possible_exclusion", "llm_all")) {
    return(unique_query_table)
  }

  stop("Unsupported judgement_mode.", call. = FALSE)
}

.evidence_summary_contains_expected_context <- function(
  evidence_summary,
  expected_environment,
  expected_habitat,
  expected_region
) {
  evidence_summary <- tolower(as.character(evidence_summary))
  contains <- rep(TRUE, length(evidence_summary))
  if (!is.null(expected_environment) && .is_non_empty_value(expected_environment)) {
    contains <- contains & grepl(tolower(expected_environment), evidence_summary, fixed = TRUE)
  }
  if (!is.null(expected_habitat) && .is_non_empty_value(expected_habitat)) {
    contains <- contains & grepl(tolower(expected_habitat), evidence_summary, fixed = TRUE)
  }
  if (!is.null(expected_region) && .is_non_empty_value(expected_region)) {
    contains <- contains & grepl(tolower(expected_region), evidence_summary, fixed = TRUE)
  }
  contains[is.na(contains)] <- FALSE
  contains
}

.evidence_summary_suggests_possible_incompatibility <- function(
  evidence_summary,
  expected_environment,
  expected_habitat,
  expected_region
) {
  text <- tolower(as.character(evidence_summary))
  out <- grepl("\\bincompatible\\b|restricted_elsewhere|known_elsewhere_only", text)
  if (identical(tolower(expected_environment), "marine")) {
    out <- out | (grepl("freshwater|terrestrial", text) & !grepl("marine", text))
  }
  out[is.na(out)] <- FALSE
  out
}

.empty_flag_taxa_tool_evidence <- function() {
  out <- as.data.frame(
    stats::setNames(
      replicate(length(.flag_taxa_tool_evidence_columns()), character(), simplify = FALSE),
      .flag_taxa_tool_evidence_columns()
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$tool_used <- logical()
  out$tool_found_explicit_evidence <- logical()
  out
}

.flag_taxa_required_tool_failure_evidence <- function(
  query_table,
  error_message,
  prompt_tools = NULL,
  expected_environment = NULL,
  expected_habitat = NULL,
  expected_region = NULL
) {
  query_table <- .validate_flag_taxa_query_table(query_table)
  if (nrow(query_table) == 0) {
    return(.empty_flag_taxa_tool_evidence())
  }
  tool_name <- .flag_taxa_default_tool_name(prompt_tools)
  tool_query <- vapply(seq_len(nrow(query_table)), function(row_index) {
    pieces <- c(
      query_table$query_name[[row_index]],
      query_table$query_rank[[row_index]],
      expected_environment,
      expected_habitat,
      expected_region
    )
    pieces <- pieces[.is_non_empty_value(pieces)]
    paste(pieces, collapse = " ")
  }, character(1))
  summary <- paste(
    "tool_failed:",
    "Tool evidence was required but the tool evidence pass failed:",
    error_message,
    .flag_taxa_required_tool_missing_summary()
  )
  out <- data.frame(
    query_id = as.character(query_table$query_id),
    query_name = as.character(query_table$query_name),
    query_rank = as.character(query_table$query_rank),
    tool_used = TRUE,
    tool_name = tool_name,
    tool_query = tool_query,
    tool_evidence_summary = summary,
    tool_references = .flag_taxa_required_tool_missing_reference(),
    tool_found_explicit_evidence = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  .validate_flag_taxa_tool_evidence(out, query_table = query_table)
}

call_flag_taxa_tool_evidence <- function(
  query_table,
  expected_environment,
  expected_habitat = NULL,
  expected_region = NULL,
  chat,
  max_tries = 3,
  retry_sleep = 5
) {
  .require_ellmer_for_flag_taxa()
  if (is.null(chat) || !is.function(chat$chat_structured)) {
    stop(
      "`chat` must be an ellmer chat object with a `chat_structured()` method.",
      call. = FALSE
    )
  }
  prompt <- .build_flag_taxa_tool_evidence_prompt(
    query_table = query_table,
    expected_environment = expected_environment,
    expected_habitat = expected_habitat,
    expected_region = expected_region
  )
  output_type <- .flag_taxa_tool_evidence_output_type()
  for (attempt in seq_len(max_tries)) {
    result <- tryCatch(
      list(value = chat$chat_structured(prompt, type = output_type), error = NULL),
      error = function(error) list(value = NULL, error = error)
    )
    if (is.null(result$error)) {
      return(.as_flag_taxa_result_table(result$value))
    }
    if (!.is_retryable_ellmer_error(result$error) || attempt >= max_tries) {
      .abort_ellmer_structured_call(
        error = result$error,
        chat = chat,
        attempts = attempt
      )
    }
    if (retry_sleep > 0) {
      Sys.sleep(retry_sleep)
    }
  }
  stop("ellmer tool evidence call failed unexpectedly.", call. = FALSE)
}

.build_flag_taxa_tool_evidence_prompt <- function(
  query_table,
  expected_environment,
  expected_habitat,
  expected_region
) {
  paste(
    "Use registered tools, if available, to find explicit evidence for these taxonomy queries.",
    "Search only for evidence about environment, habitat, or geographic occurrence relevant to the expected context.",
    "Use the compact lineage to disambiguate taxa and build precise search terms.",
    "Return valid JSON only. If no explicit evidence is found, set tool_found_explicit_evidence = FALSE and use a conservative summary beginning with 'no_evidence_found:'.",
    "No explicit tool evidence is missing evidence, not incompatibility evidence. Do not claim absence of records proves absence from the habitat or region.",
    "Do not invent AphiaIDs, DOIs, URLs, article titles, author names, or citations. Use only references returned by registered tools.",
    "",
    "Expected context:",
    paste0("- expected_environment: ", expected_environment),
    paste0("- expected_habitat: ", .format_optional_prompt_value(expected_habitat)),
    paste0("- expected_region: ", .format_optional_prompt_value(expected_region)),
    "",
    "Return one object per query_id with fields:",
    paste(paste0("- ", .flag_taxa_tool_evidence_columns()), collapse = "\n"),
    "",
    "Query taxa, tab-separated:",
    .format_prompt_table(query_table[, c("query_id", "query_name", "query_rank", "lineage"), drop = FALSE]),
    sep = "\n"
  )
}

.flag_taxa_tool_evidence_output_type <- function() {
  row_fields <- list(
    query_id = ellmer::type_string("Copy of query_id for this query taxon."),
    query_name = ellmer::type_string("Copy of query_name."),
    query_rank = ellmer::type_string("Copy of query_rank."),
    tool_used = ellmer::type_boolean("Whether a registered tool was used."),
    tool_name = ellmer::type_string("Compact tool/source name, such as scite, or NA."),
    tool_query = ellmer::type_string("Compact summary of tool query terms, or NA."),
    tool_evidence_summary = ellmer::type_string("Compact summary of explicit tool-derived evidence, or no-evidence statement."),
    tool_references = ellmer::type_string("References, URLs, or citation labels from tool-derived evidence, or NA."),
    tool_found_explicit_evidence = ellmer::type_boolean("Whether explicit relevant evidence was found.")
  )
  row_type <- do.call(
    ellmer::type_object,
    c(list(.description = "One tool evidence row for a query taxon."), row_fields)
  )
  ellmer::type_array(
    row_type,
    description = "One tool evidence row for each query taxon in the batch."
  )
}

.flag_taxa_allowed_values <- function(expected_region = NULL, for_validation = FALSE) {
  expected_region <- .validate_optional_scalar_character(
    expected_region,
    "expected_region"
  )
  if (!is.logical(for_validation) || length(for_validation) != 1 || is.na(for_validation)) {
    stop("`for_validation` must be a logical scalar.", call. = FALSE)
  }

  expected_region_status <- c(
    "known_in_region",
    "known_near_region",
    "known_elsewhere_only",
    "restricted_elsewhere",
    "no_distribution_evidence",
    "unknown"
  )
  if (isTRUE(for_validation)) {
    expected_region_status <- c(expected_region_status, "not_assessed")
  } else if (is.null(expected_region)) {
    expected_region_status <- "not_assessed"
  }

  list(
    expected_environment_status = c(
      "compatible",
      "incompatible",
      "mixed_within_rank",
      "unknown",
      "insufficient_taxonomic_resolution"
    ),
    expected_habitat_status = c(
      "compatible",
      "incompatible",
      "no_known_habitat_evidence",
      "possible_likely",
      "possible_unlikely",
      "transient_or_allochthonous_possible",
      "mixed_within_rank",
      "unknown",
      "insufficient_taxonomic_resolution"
    ),
    expected_region_status = expected_region_status,
    ecological_status = c(
      "compatible",
      "plausible",
      "unlikely_resident",
      "incompatible_resident",
      "unknown",
      "insufficient_taxonomic_resolution"
    ),
    occurrence_interpretation = c(
      "expected_resident",
      "plausible_resident",
      "possible_transient_or_allochthonous",
      "possible_contaminant_or_misassignment",
      "likely_contaminant_or_misassignment",
      "uncertain",
      "not_assessed"
    ),
    recommended_action = c(
      "retain",
      "flag_for_review",
      "flag_possible_exclusion",
      "exclude"
    ),
    evidence_basis = c(
      "local_evidence",
      "tool_evidence",
      "local_and_tool_evidence",
      "conservative_reasoning_only"
    )
  )
}

.validate_flag_taxa_enum_column <- function(result, column_name, allowed_values) {
  values <- as.character(result[[column_name]])
  invalid_values <- unique(values[is.na(values) | !(values %in% allowed_values)])

  if (length(invalid_values) > 0) {
    invalid_values[is.na(invalid_values)] <- "<NA>"
    stop(
      "`",
      column_name,
      "` contains invalid values: ",
      paste(invalid_values, collapse = ", "),
      ". Allowed values are: ",
      paste(allowed_values, collapse = ", "),
      call. = FALSE
    )
  }
}

.validate_flag_taxa_region_assessment <- function(result, expected_region) {
  region_status <- as.character(result$expected_region_status)
  not_assessed <- region_status == "not_assessed"

  if (is.null(expected_region) && any(!not_assessed)) {
    invalid_values <- unique(region_status[!not_assessed])
    stop(
      "`expected_region_status` must be `not_assessed` when ",
      "`expected_region` is NULL. Found: ",
      paste(invalid_values, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(expected_region) && any(not_assessed)) {
    stop(
      "`expected_region_status` must not be `not_assessed` when ",
      "`expected_region` is provided.",
      call. = FALSE
    )
  }
}

.add_taxon_id_column <- function(tax_table) {
  row_names <- rownames(tax_table)
  default_row_names <- as.character(seq_len(nrow(tax_table)))
  if (
    !is.null(row_names) &&
      !identical(row_names, default_row_names) &&
      !("taxon_id" %in% colnames(tax_table))
  ) {
    tax_table <- cbind(
      taxon_id = row_names,
      tax_table,
      stringsAsFactors = FALSE
    )
  }

  rownames(tax_table) <- NULL
  tax_table
}

.coerce_tax_row <- function(tax_row) {
  if (inherits(tax_row, "data.frame")) {
    if (nrow(tax_row) != 1) {
      stop("`tax_row` must contain exactly one row.", call. = FALSE)
    }
    out <- vapply(tax_row, function(value) as.character(value[[1]]), character(1))
  } else if (is.list(tax_row)) {
    out <- vapply(tax_row, function(value) as.character(value[[1]]), character(1))
  } else {
    out <- as.character(tax_row)
    names(out) <- names(tax_row)
  }

  if (is.null(names(out)) || any(names(out) == "")) {
    stop("`tax_row` must be a named taxonomy row.", call. = FALSE)
  }

  out
}

.tax_rank_value <- function(tax_row, rank) {
  rank_index <- which(tolower(names(tax_row)) == tolower(rank))
  if (length(rank_index) == 0) {
    return(NA_character_)
  }

  as.character(tax_row[[rank_index[[1]]]])
}

.build_lineage <- function(tax_row, tax_ranks) {
  tax_row <- .coerce_tax_row(tax_row)
  lineage_ranks <- rev(tax_ranks)
  lineage <- character()

  for (rank in lineage_ranks) {
    taxon_name <- .tax_rank_value(tax_row, rank)
    if (.is_informative_taxon(taxon_name)) {
      lineage <- c(
        lineage,
        paste0(tolower(rank), ": ", .clean_taxon_name(taxon_name))
      )
    }
  }

  if (length(lineage) == 0) {
    return(NA_character_)
  }

  paste(lineage, collapse = "; ")
}

.build_compact_lineage <- function(tax_row, tax_ranks) {
  tax_row <- .coerce_tax_row(tax_row)
  lineage_ranks <- .flag_taxa_lineage_ranks(tax_ranks)
  lineage <- character(length(lineage_ranks))

  for (rank_index in seq_along(lineage_ranks)) {
    rank <- lineage_ranks[[rank_index]]
    taxon_name <- .tax_rank_value(tax_row, rank)
    taxon_name <- if (.is_informative_taxon(taxon_name)) {
      .clean_taxon_name(taxon_name)
    } else {
      "NA"
    }
    lineage[[rank_index]] <- paste0(
      .format_tax_rank_label(rank),
      "=",
      taxon_name
    )
  }

  paste(lineage, collapse = "; ")
}

.flag_taxa_lineage_ranks <- function(tax_ranks) {
  tax_ranks <- tolower(trimws(as.character(tax_ranks)))
  default_lineage_ranks <- c("phylum", "class", "order", "family", "genus", "species")
  lineage_ranks <- default_lineage_ranks[default_lineage_ranks %in% tax_ranks]
  if (length(lineage_ranks) > 0) {
    return(lineage_ranks)
  }

  rev(tax_ranks)
}

.format_tax_rank_label <- function(rank) {
  rank <- tolower(trimws(as.character(rank[[1]])))
  labels <- c(
    domain = "Domain",
    superkingdom = "Superkingdom",
    kingdom = "Kingdom",
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "Species"
  )
  if (rank %in% names(labels)) {
    return(labels[[rank]])
  }

  paste0(toupper(substr(rank, 1, 1)), substring(rank, 2))
}

.is_informative_taxon <- function(taxon_name) {
  if (length(taxon_name) == 0 || is.na(taxon_name[[1]])) {
    return(FALSE)
  }

  cleaned_name <- tolower(.clean_taxon_name(taxon_name[[1]]))
  unknown_values <- c(
    "",
    "na",
    "nan",
    "null",
    "none",
    "unknown",
    "uncultured",
    "unclassified",
    "unassigned",
    "undetermined",
    "unidentified",
    "environmental sample",
    "metagenome",
    "incertae sedis",
    "not assigned"
  )

  !(
    cleaned_name %in% unknown_values ||
      grepl("^[[:alpha:]]__\\s*$", cleaned_name) ||
      grepl("^(asv|otu)[_-]?[0-9]+$", cleaned_name) ||
      grepl(
        "\\b(uncultured|unidentified|unclassified|unknown|environmental sample|metagenome|incertae sedis)\\b",
        cleaned_name
      )
  )
}

.clean_taxon_name <- function(taxon_name) {
  taxon_name <- trimws(as.character(taxon_name[[1]]))
  taxon_name <- gsub("[\r\n\t]+", " ", taxon_name)
  taxon_name <- sub("^[[:alpha:]]__", "", taxon_name)
  trimws(taxon_name)
}

.format_prompt_table <- function(tax_table) {
  formatted <- as.data.frame(
    tax_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  formatted[] <- lapply(formatted, function(column) {
    column <- as.character(column)
    column[is.na(column)] <- ""
    gsub("[\r\n\t]+", " ", column)
  })
  if (nrow(formatted) == 0) {
    return(paste(colnames(formatted), collapse = "\t"))
  }

  rows <- apply(formatted, 1, paste, collapse = "\t")
  paste(c(paste(colnames(formatted), collapse = "\t"), rows), collapse = "\n")
}

.format_optional_prompt_value <- function(value) {
  if (is.null(value)) {
    return("not specified")
  }

  value
}

.format_flag_taxa_evidence_sources <- function(evidence_sources) {
  if (is.null(evidence_sources) || length(evidence_sources) == 0) {
    return("none requested")
  }

  paste(evidence_sources, collapse = ", ")
}

.check_flag_taxa_prompt_size <- function(
  prompt,
  max_prompt_chars,
  chunk_index = NULL,
  chunk_taxa = NULL
) {
  if (is.null(max_prompt_chars)) {
    return(invisible(prompt))
  }

  prompt_chars <- nchar(prompt, type = "chars")
  if (prompt_chars > max_prompt_chars) {
    chunk_text <- if (!is.null(chunk_index)) {
      paste0(" for chunk ", chunk_index)
    } else {
      ""
    }
    taxa_text <- if (!is.null(chunk_taxa)) {
      paste0(" containing ", chunk_taxa, " taxa")
    } else {
      ""
    }
    stop(
      "`flag_taxa()` generated a prompt",
      chunk_text,
      taxa_text,
      " with ",
      prompt_chars,
      " characters, which exceeds `max_prompt_chars = ",
      max_prompt_chars,
      "`. Reduce `max_evidence_rows_per_query`, reduce `max_tool_taxa`, ",
      "provide fewer taxa, set `judgement_mode` to review fewer taxa, or use ",
      "a larger-context model.",
      call. = FALSE
    )
  }

  invisible(prompt)
}

.build_flag_taxa_prompt_chunks <- function(
  tax_table,
  unique_query_table,
  preliminary_judgement,
  expected_environment,
  expected_habitat,
  expected_region,
  tax_ranks,
  evidence_sources,
  taxon_evidence,
  tool_evidence,
  max_evidence_rows_per_query,
  max_evidence_summary_chars,
  max_tool_evidence_summary_chars,
  max_lineage_chars,
  allow_llm_tools,
  tool_requirement = "required_for_llm",
  llm_chunk_size,
  max_prompt_chars,
  auto_reduce_llm_chunk_size = TRUE,
  prompt_tools = NULL
) {
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  query_chunks <- .split_flag_taxa_query_table(unique_query_table, llm_chunk_size)
  if (length(query_chunks) == 0) {
    query_chunks <- list(chunk_001 = unique_query_table)
  }

  build_one_prompt <- function(query_table) {
    build_flag_taxa_prompt(
      tax_table = tax_table,
      expected_environment = expected_environment,
      expected_habitat = expected_habitat,
      expected_region = expected_region,
      tax_ranks = tax_ranks,
      evidence_sources = evidence_sources,
      taxon_evidence = taxon_evidence,
      tool_evidence = tool_evidence,
      query_table = query_table,
      preliminary_judgement = preliminary_judgement,
      max_evidence_rows_per_query = max_evidence_rows_per_query,
      max_evidence_summary_chars = max_evidence_summary_chars,
      max_tool_evidence_summary_chars = max_tool_evidence_summary_chars,
      max_lineage_chars = max_lineage_chars,
      allow_llm_tools = allow_llm_tools,
      prompt_tools = prompt_tools,
      tool_requirement = tool_requirement
    )
  }

  if (isTRUE(auto_reduce_llm_chunk_size) && !is.null(max_prompt_chars)) {
    query_chunks <- .auto_reduce_flag_taxa_query_chunks(
      query_chunks = query_chunks,
      build_prompt = build_one_prompt,
      max_prompt_chars = max_prompt_chars
    )
  }

  prompts <- lapply(seq_along(query_chunks), function(chunk_index) {
    query_table <- query_chunks[[chunk_index]]
    prompt <- build_one_prompt(query_table)
    .check_flag_taxa_prompt_size(
      prompt = prompt,
      max_prompt_chars = max_prompt_chars,
      chunk_index = chunk_index,
      chunk_taxa = nrow(query_table)
    )
    prompt
  })
  names(prompts) <- names(query_chunks)
  attr(prompts, "query_chunks") <- query_chunks
  attr(prompts, "prompt_chars") <- vapply(prompts, nchar, integer(1), type = "chars")
  prompts
}

.auto_reduce_flag_taxa_query_chunks <- function(
  query_chunks,
  build_prompt,
  max_prompt_chars
) {
  queue <- unname(query_chunks)
  out <- list()

  while (length(queue) > 0) {
    chunk <- queue[[1]]
    queue <- queue[-1]
    prompt <- build_prompt(chunk)
    prompt_chars <- nchar(prompt, type = "chars")
    if (prompt_chars <= max_prompt_chars || nrow(chunk) <= 1) {
      out[[length(out) + 1L]] <- chunk
      next
    }

    split_at <- floor(nrow(chunk) / 2)
    first <- chunk[seq_len(split_at), , drop = FALSE]
    second <- chunk[seq.int(split_at + 1L, nrow(chunk)), , drop = FALSE]
    rownames(first) <- NULL
    rownames(second) <- NULL
    queue <- c(list(first, second), queue)
  }

  names(out) <- sprintf("chunk_%03d", seq_along(out))
  out
}

.report_flag_taxa_prompt_chunks <- function(
  prompt_chunks,
  query_chunks,
  verbose,
  start_time
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  if (is.null(query_chunks)) {
    query_chunks <- attr(prompt_chunks, "query_chunks", exact = TRUE)
  }
  prompt_chars <- vapply(prompt_chunks, nchar, integer(1), type = "chars")
  for (chunk_index in seq_along(prompt_chunks)) {
    query_count <- if (!is.null(query_chunks) && length(query_chunks) >= chunk_index) {
      nrow(query_chunks[[chunk_index]])
    } else {
      NA_integer_
    }
    .flag_taxa_progress(
      TRUE,
      start_time,
      "prompt chunk ",
      chunk_index,
      " / ",
      length(prompt_chunks),
      " contains ",
      query_count,
      " taxa and ",
      prompt_chars[[chunk_index]],
      " chars."
    )
  }
  invisible(NULL)
}

.write_flag_taxa_prompt_file <- function(prompt_chunks, prompt_path = NULL) {
  prompt_path <- .validate_optional_path(prompt_path, "prompt_path")
  prompt_path <- .resolve_flag_taxa_prompt_path(prompt_path)
  parent <- dirname(prompt_path)
  if (!identical(parent, ".") && !dir.exists(parent)) {
    stop(
      "The parent directory for `prompt_path` does not exist: ",
      parent,
      call. = FALSE
    )
  }

  prompt <- .combine_flag_taxa_prompt_chunks(prompt_chunks)
  writeLines(as.character(prompt), con = prompt_path, useBytes = TRUE)
  prompt_path
}

.resolve_flag_taxa_prompt_path <- function(prompt_path = NULL) {
  default_name <- "flag_taxa_prompt.txt"
  if (is.null(prompt_path)) {
    return(file.path(getwd(), default_name))
  }
  prompt_path <- path.expand(as.character(prompt_path))
  if (dir.exists(prompt_path)) {
    return(file.path(prompt_path, default_name))
  }
  extension <- tolower(tools::file_ext(prompt_path))
  if (identical(extension, "txt")) {
    return(prompt_path)
  }
  paste0(prompt_path, ".txt")
}

.combine_flag_taxa_prompt_chunks <- function(prompt_chunks) {
  if (length(prompt_chunks) == 0) {
    return("")
  }
  if (length(prompt_chunks) == 1) {
    return(as.character(prompt_chunks[[1]]))
  }
  pieces <- lapply(seq_along(prompt_chunks), function(chunk_index) {
    paste(
      paste0("===== CHUNK ", chunk_index, " / ", length(prompt_chunks), " ====="),
      as.character(prompt_chunks[[chunk_index]]),
      sep = "\n"
    )
  })
  paste(unlist(pieces, use.names = FALSE), collapse = "\n\n")
}

.warn_flag_taxa_large_prompt_chunks <- function(
  prompt_chunks,
  query_chunks,
  warn_prompt_chars,
  llm_chunk_size
) {
  if (is.null(warn_prompt_chars) || length(prompt_chunks) == 0) {
    return(invisible(NULL))
  }
  prompt_chars <- vapply(prompt_chunks, nchar, integer(1), type = "chars")
  large <- which(prompt_chars > warn_prompt_chars)
  if (length(large) == 0) {
    return(invisible(NULL))
  }

  if (length(prompt_chunks) == 1) {
    warning(
      "flag_taxa(): final judgement prompt is ",
      prompt_chars[[1]],
      " chars in one chunk. This may exceed model input/output limits. ",
      "Consider reducing `max_evidence_rows_per_query`, reducing tool evidence, ",
      "reviewing fewer taxa, or using a larger-context / larger-output model.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  for (chunk_index in large) {
    query_count <- if (!is.null(query_chunks) && length(query_chunks) >= chunk_index) {
      nrow(query_chunks[[chunk_index]])
    } else {
      NA_integer_
    }
    warning(
      "flag_taxa(): final judgement prompt chunk ",
      chunk_index,
      " / ",
      length(prompt_chunks),
      " contains ",
      query_count,
      " taxa and ",
      prompt_chars[[chunk_index]],
      " chars. This may exceed model input/output limits. Consider reducing ",
      "`max_evidence_rows_per_query`, tool evidence, or using a larger-context / ",
      "larger-output model.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

.split_flag_taxa_query_table <- function(unique_query_table, llm_chunk_size) {
  if (nrow(unique_query_table) == 0) {
    return(list())
  }
  if (is.infinite(llm_chunk_size) || nrow(unique_query_table) <= llm_chunk_size) {
    return(stats::setNames(list(unique_query_table), "chunk_001"))
  }

  chunk_index <- ceiling(seq_len(nrow(unique_query_table)) / llm_chunk_size)
  chunks <- split(unique_query_table, chunk_index)
  names(chunks) <- sprintf("chunk_%03d", seq_along(chunks))
  lapply(chunks, function(chunk) {
    rownames(chunk) <- NULL
    chunk
  })
}

.call_flag_taxa_ellmer_chunks <- function(
  prompt_chunks,
  query_chunks,
  chat,
  expected_region,
  max_tries,
  retry_sleep,
  allow_llm_tools,
  tool_requirement = "optional",
  prompt_tools = NULL,
  repair_invalid_llm = TRUE,
  verbose,
  start_time
) {
  tool_requirement <- .validate_flag_taxa_tool_requirement(tool_requirement)
  prompt_tools <- .validate_optional_prompt_tools(prompt_tools)
  rows <- vector("list", length(prompt_chunks))
  for (chunk_index in seq_along(prompt_chunks)) {
    prompt_chars <- nchar(prompt_chunks[[chunk_index]], type = "chars")
    .flag_taxa_progress(
      verbose,
      start_time,
      "running LLM judgement chunk ",
      chunk_index,
      " / ",
      length(prompt_chunks),
      " (",
      nrow(query_chunks[[chunk_index]]),
      " taxa, ",
      prompt_chars,
      " chars)."
    )
    chunk_result <- call_flag_taxa_ellmer(
      prompt = prompt_chunks[[chunk_index]],
      chat = chat,
      include_feature_id = FALSE,
      expected_region = expected_region,
      max_tries = max_tries,
      retry_sleep = retry_sleep
    )
    chunk_result <- .repair_flag_taxa_invalid_llm_output(
      result = chunk_result,
      repair_invalid_llm = repair_invalid_llm,
      expected_region = expected_region,
      verbose = verbose,
      start_time = start_time
    )
    rows[[chunk_index]] <- validate_flag_taxa_taxon_output(
      result = chunk_result,
      query_taxonomy = query_chunks[[chunk_index]],
      expected_region = expected_region,
      allow_llm_tools = allow_llm_tools
    )
    .flag_taxa_progress(
      verbose,
      start_time,
      "completed LLM judgement chunk ",
      chunk_index,
      " / ",
      length(prompt_chunks),
      "."
    )
  }

  .bind_rows_fill(rows)
}

.format_flag_taxa_llm_start_message <- function(unique_taxa, chunks, llm_chunk_size) {
  if (chunks <= 1) {
    return(paste0(
      "starting LLM structured-output call for ",
      unique_taxa,
      " unique taxa in 1 chunk."
    ))
  }

  paste0(
    "starting LLM structured-output calls for ",
    unique_taxa,
    " unique taxa in ",
    chunks,
    " chunks of up to ",
    llm_chunk_size,
    " taxa."
  )
}

.flag_taxa_progress <- function(verbose, start_time, ..., total = FALSE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  label <- if (isTRUE(total)) " Total time: " else " Elapsed: "
  message(
    "flag_taxa(): ",
    paste0(..., collapse = ""),
    label,
    .format_elapsed_time(start_time),
    "."
  )
  invisible(NULL)
}

.format_elapsed_time <- function(start_time) {
  seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  sprintf("%.1f sec", seconds)
}

.flag_taxa_evidence_match_stats <- function(
  unique_query_table,
  taxon_evidence,
  tax_ranks,
  max_evidence_rows_per_query,
  max_evidence_summary_chars = 160
) {
  relevant <- collect_relevant_taxon_evidence(
    tax_table = unique_query_table,
    taxon_evidence = taxon_evidence,
    tax_ranks = tax_ranks,
    max_evidence_rows_per_query = max_evidence_rows_per_query,
    max_evidence_summary_chars = max_evidence_summary_chars
  )
  if (is.null(relevant) || nrow(relevant) == 0) {
    return(list(evidence_rows = 0L, query_taxa_with_evidence = 0L))
  }

  if ("query_id" %in% colnames(relevant)) {
    covered <- length(unique(as.character(relevant$query_id)))
  } else {
    covered <- length(unique(.flag_taxa_result_query_key(
      relevant$query_name,
      relevant$query_rank
    )))
  }
  list(evidence_rows = nrow(relevant), query_taxa_with_evidence = covered)
}

.validate_required_scalar_character <- function(value, name) {
  if (
    missing(value) ||
      !is.character(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !nzchar(trimws(value))
  ) {
    stop("`", name, "` must be a non-empty character scalar.", call. = FALSE)
  }

  trimws(value)
}

.validate_optional_scalar_character <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }

  .validate_required_scalar_character(value, name)
}

.validate_optional_path <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }

  .validate_required_scalar_character(value, name)
}

.validate_character_vector <- function(value, name) {
  if (
    !is.character(value) ||
      length(value) < 1 ||
      any(is.na(value)) ||
      any(!nzchar(trimws(value)))
  ) {
    stop("`", name, "` must be a non-empty character vector.", call. = FALSE)
  }

  unique(tolower(trimws(value)))
}

.validate_optional_prompt_tools <- function(value) {
  if (is.null(value)) {
    return(character())
  }
  if (is.character(value) && length(value) == 0) {
    return(character())
  }
  if (
    !is.character(value) ||
      length(value) < 1 ||
      any(is.na(value)) ||
      any(!nzchar(trimws(value)))
  ) {
    stop("`prompt_tools` must be NULL or a non-empty character vector.", call. = FALSE)
  }

  unique(trimws(value))
}

.validate_flag_taxa_tool_requirement <- function(tool_requirement) {
  match.arg(tool_requirement, choices = c("required_for_llm", "optional"))
}

.validate_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }

  value
}

.validate_optional_positive_integer_scalar <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }

  .validate_positive_integer_scalar(value, name)
}

.validate_positive_count_or_inf <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be a positive number or Inf.", call. = FALSE)
  }
  if (is.infinite(value)) {
    return(Inf)
  }
  if (value < 1 || value != floor(value)) {
    stop("`", name, "` must be a positive integer or Inf.", call. = FALSE)
  }
  as.integer(value)
}

.validate_non_negative_count_or_inf <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be a non-negative number or Inf.", call. = FALSE)
  }
  if (is.infinite(value)) {
    return(Inf)
  }
  if (value < 0 || value != floor(value)) {
    stop("`", name, "` must be a non-negative integer or Inf.", call. = FALSE)
  }
  as.integer(value)
}

.validate_flag_taxa_tool_use_policy <- function(
  tool_use_policy,
  allow_llm_tools,
  use_default = FALSE
) {
  if (isTRUE(use_default)) {
    return(if (isTRUE(allow_llm_tools)) "review_or_exclude" else "never")
  }
  policy <- match.arg(
    tool_use_policy,
    choices = c("never", "review_or_exclude", "exclude_only", "always")
  )
  if (!isTRUE(allow_llm_tools)) {
    return("never")
  }
  policy
}

.validate_flag_taxa_judgement_mode <- function(
  judgement_mode,
  use_default,
  prompt_only,
  chat,
  supplied_llm_result
) {
  choices <- c(
    "evidence_only",
    "llm_missing_evidence",
    "llm_flagged",
    "llm_possible_exclusion",
    "llm_all",
    "llm_review"
  )
  if (isTRUE(use_default)) {
    if (isTRUE(prompt_only)) {
      return("evidence_only")
    }
    if (!is.null(supplied_llm_result)) {
      return("llm_all")
    }
    if (!is.null(chat)) {
      return("llm_all")
    }
    return("evidence_only")
  }

  judgement_mode <- match.arg(judgement_mode, choices = choices)
  if (identical(judgement_mode, "llm_review")) {
    return("llm_flagged")
  }
  judgement_mode
}
