# biohelper 0.0.24.9000 (development)

## Taxonomy normalisation

- Initial multi-ID name lookup warnings are handled internally. Ambiguity
  warnings describe conflicts remaining after lineage matching and optional
  BLAST, including unassigned cases with no distinguishing broader name.
  Unrelated lookup warnings and errors remain visible.
- BLAST ambiguity checking now defaults to `fallback = "lca"`; use an explicit
  `fallback = "first"` for the previous fallback behaviour. Name-only ambiguity
  defaults are unchanged.
- Progress omits repeated per-round name lookup messages. Local BLAST heartbeats
  show elapsed minutes every five minutes by default (`progress_seconds = 300`).
  Both help pages and the README include complete local BLAST examples.
- `taxo_merge()` plans all NCBI name matches before running one pooled BLAST
  batch for the uncached distinct sequences selected across all input tables.
  FASTA reading and hit ancestry lookup are shared. `refresh = TRUE` also
  deduplicates across inputs, and feature subsets may span different input IDs.
- Both normalisation entry points group repeated ordered taxonomy profiles,
  reuse name/lineage matching, and keep separate evidence for distinct DNA.
  BLAST hit screening and merging avoid repeated full-table scans.
- `blast_options` adds `verbose = TRUE` and `progress_seconds = 300`, showing
  cache/search counts, active settings, elapsed-time updates and fallback stages.
  Local searches use processx for progress and cleanup on interruption; numerical
  screening settings, task and one-thread default are unchanged.
- `taxo_merge()` now accepts `ambiguity`, `sequences`, `db`, `blastapp_path`,
  and `blast_options`, forwarding them to each input's normalisation before
  merging. Local BLAST checking with an LCA fallback can be enabled directly
  in the merge call. Reports for non-default ambiguity strategies are retained
  in the `normalisation` attribute without adding taxonomy columns. The existing
  `spnc` argument is now forwarded as supplied. Calls with no assigned ranks
  return an empty table instead of failing during merging.
- `taxo_normalisation()` searches ordered names from right to left across any
  number of taxonomic columns. Input labels do not impose NCBI ranks.
- Conflicting names are matched against full ordered NCBI lineages. Ambiguity
  options are `first`, `lca`, `unassigned`, and `blast_check`.
- Optional sequence checking uses `blast_options`; local BLAST is the default.
  `db` and `blastapp_path` follow `blastn_taxo_assignment()` conventions. Missing
  sequence IDs automatically use the selected fallback.
- Default output remains a taxonomy table or phyloseq object. Reports are
  available through attributes or the explicit `return_details = TRUE` list.
  Alignment/provenance columns are preserved separately from taxonomy names.
- `spnc = TRUE` explicitly joins the final two taxonomic columns, independently
  of their labels.

## BLAST/LCA assignment

- `blastn_taxo_assignment()` adds `resume = TRUE` and `overwrite = FALSE`, with
  independent validated checkpoints for both searches, both LCA tables and
  consensus. Completed searches survive downstream failures and threshold changes.
- Superseded results are archived by default. An explicit fresh run uses
  `resume = FALSE, overwrite = TRUE`. Older files can be structurally adopted
  with explicit provenance uncertainties.
- Taxonomy updates run once per invocation and use a staged HTTPS download and
  validated compact biohelper cache. Existing shared taxizedb databases are not
  replaced. Optional `taxonomy_db` accepts existing taxizedb/taxonomizr SQLite
  databases. `update = FALSE` never downloads taxonomy.
- Representative hits are selected deterministically by E-value, bitscore,
  identity, coverage and subject ID. Consensus keeps statistics from one hit.
- LCA conflicts are now checked broad-to-narrow. Additional TaxIDs attached to a
  reference are no longer discarded, and unresolved or identity-masked
  representatives do not suppress known lineage evidence. These bug fixes can
  change assignments; the existing numeric identity/coverage thresholds are retained.
- The final genuine CSV is `blastn_taxo_assignment.csv`; the old misspelled
  filename is recognised during resumption. Every FASTA query is retained.
- Five fields follow `species`: `pident`, `qcovs`, `evalue`, `assignment_method`,
  and `support_method`. Select rank columns by name rather than assuming nine
  total columns. The wrapper now returns its documented dataframe, rather than
  the previous accidental `NULL`.
- `lcaPident()` retains its original column aliases and adds consistent numeric
  alignment fields and a representative subject identifier. Both assignment
  functions record database/software/parameter provenance and validate inputs.

## Other corrections

- `taxo_bar_plot()` preserves supplied `n_rank1` values and uses `n_rank2` to
  control lower-rank taxa per group.
- `taxo_merge()` returns lowercase taxonomy rank columns by default.

## Validation

- The complete test suite passed 1,949 assertions before publication preparation.
- Real local BLAST checks used three supplied ASVs and a small authentic NCBI
  reference set. Separate R sessions verified checkpoint reuse and LCA-only
  recalculation after threshold changes. These are software checks, not exhaustive
  biological validation or species-identification guarantees.
