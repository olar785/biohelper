<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

Helper functions for various bioinformatics, data visualisation and analysis operations, including taxonomic normalisation, ASV decontamination, barplot of taxonomic composition, etc.

The package can be installed with: 

if (!require("remotes")) install.packages("remotes");
remotes::install_github("olar785/biohelper")

The development version containing the taxonomy and resumable BLAST updates is
available from the `development` branch:

```r
remotes::install_github("olar785/biohelper", ref = "development", upgrade = "never")
```

Restart R after installation and check `packageVersion("biohelper")`. The current
development version is `0.0.24.9000`. See [NEWS.md](NEWS.md) for changes and output
compatibility notes. This development version is provided for testing.

## Restartable BLAST taxonomy assignment

`blastn_taxo_assignment()` independently checkpoints each megablast/blastn search,
each LCA assignment and final merging. Rerunning an interrupted call reuses valid
completed stages. Changing LCA thresholds reuses raw searches; changing one
method's search options invalidates only that method and its downstream results.

```r
arguments <- list(
  blastapp_path = "/tools/blastn", method = "both", queries = "ASVs.fasta",
  db = "/databases/nt",                 # formatted BLAST sequence database
  taxonomy_db = "accessionTaxa.sql",    # separate local taxonomy SQLite file
  output_path = "blast_results", nthreads = 4
)
result <- do.call(biohelper::blastn_taxo_assignment, arguments)
# After interruption: same call; completed stages print SKIPPED.
result <- do.call(biohelper::blastn_taxo_assignment, arguments)
# Recalculate taxonomy without repeating searches.
arguments$minSim <- 98
revised <- do.call(biohelper::blastn_taxo_assignment, arguments)
```

`resume = TRUE, overwrite = FALSE` is the default. Superseded files and their
checkpoint provenance are preserved under `.biohelper/archive` after a validated
replacement succeeds. Hard links or renaming avoid copying large files solely
for checkpointing. Archives can grow and are not automatically purged.
`resume = TRUE, overwrite = TRUE` still reuses valid checkpoints, but permits
discarding superseded results. An explicit fresh run uses
`resume = FALSE, overwrite = TRUE`; `FALSE/FALSE` refuses to replace requested
existing results. Output directories have a process lock; dead same-host locks
are recovered, while cross-host/unverifiable locks require confirming the previous
process stopped. Use a different output directory for concurrent analyses.

The final genuine CSV is **`blastn_taxo_assignment.csv`**, containing every query
ASV, sorted by identifier, with no row-name column and NA for missing values:

```text
asv,domain,kingdom,phylum,class,order,family,genus,species,pident,qcovs,evalue,assignment_method,support_method
```

The five new fields follow `species`. Statistics are numeric, but
`phyloseq::tax_table()` can coerce them to character alongside taxonomy names.
Representative statistics come from **one hit**, ranked by lowest E-value,
highest bitscore, highest identity, highest coverage, then subject ID, with a
canonical final tie-breaker. They are alignment measurements, not LCA confidence.
Processed tables retain the selected subject, bitscore and TaxIDs, plus the
legacy identity/coverage aliases. `taxo_normalisation()` preserves the new
metadata and excludes it from its ordered-name search.

Compatible classifications retain the more resolved method and its statistics.
Identical or equally resolved compatible results use the representative-hit
ranking to choose their source; complete method ties prefer megablast.
`assignment_method` is megablast/blastn for an unmodified chosen classification.
Actual disagreements truncate the first conflicting rank and all finer ranks,
setting `assignment_method = "consensus"`. `support_method` identifies the single
method providing statistics, even for a truncated consensus; it is NA when no
statistics exist. ASVs without qualifying hits have missing taxonomy/statistics.

Small RDS records under `.biohelper/stages`, per-search logs and
`blastn_taxo_assignment_provenance.rds` record file hashes, parameters, versions,
reference identity, stage completion, hit/query counts, unresolved TaxIDs and
assignment summaries. Sequence database identity uses available file names,
sizes and modification times rather than hashing an entire large database.
After replacing a reference without changing its recorded identity, use a new
database prefix or an explicit fresh run.

Legacy raw files without metadata require structural/query-ID validation and
are adopted with an explicit warning that their original parameters/reference
versions are unknown. Empty raw files without completion metadata are rerun;
completed empty stages with metadata are reusable. Old processed files without
new statistics are recalculated from reusable raw hits. The misspelled
`blastn_taxo_assingment.csv` is detected and preserved. A final-only legacy
analysis can be migrated without searches, but unavailable statistics/provenance
remain NA and are marked unverified. Name-based rank selection remains compatible;
code assuming exactly nine output columns must account for five appended fields.
The wrapper now returns its documented data.frame; previously it inadvertently
returned NULL.

## Local taxonomy reliability and LCA corrections

`update = FALSE` never downloads. `taxonomy_db` accepts existing taxizedb or
taxonomizr SQLite names/nodes databases. NULL uses biohelper's managed cache if
present, otherwise the existing taxizedb `NCBI.sql`. Missing/inaccessible/corrupt
databases fail early; unresolved TaxIDs remain missing. Lookups batch unique IDs
and shared ancestors through one read-only SQLite connection per LCA call. They
do not query taxonomy online. Concurrent database changes cause a retryable error.

`update = TRUE` performs **one** staged update before assignment, not one per
search method. It uses the [official HTTPS NCBI taxonomy archive](https://support.nlm.nih.gov/kbArticle/?pn=KA-05543),
builds compact names/nodes tables and validates SQLite/name/parent integrity
before publication. Failed transfers/builds preserve the current cache; successful
replacements archive it. Installed taxizedb code and its shared database are
untouched. The compact cache omits taxizedb's large hierarchy table and uses
biohelper's own reader. To update, omit an explicitly selected external database
or choose a new database path; external databases are not replaced in place.

Thresholds remain species 97%, genus 95%, family 87%, order 83%, class 81%,
phylum 79%, kingdom 71%; minCov=80 filters qcovs. minSim is the species threshold,
not an overall hit-identity filter. `pident="before"` masks each hit before LCA;
`"after"` masks LCA ranks using the representative hit; `"no"` disables identity
masking. Standalone lcaPident defaults to no; the wrapper defaults to before.
Literal legacy search defaults remain unchanged: BLAST interprets perc_identity
values as percentages, not fractions.

Assignments can change where the old fine-to-broad loop concealed a broader
conflict, a merged hit's additional TaxIDs were discarded, or an unresolved or identity-masked
representative suppressed other known lineage evidence. Conflicts are now
checked broad-to-narrow, all reference TaxIDs participate, and known evidence is
retained. Deterministic representative selection can also change after-mode
masking when E-value ties have different identities. Other threshold methodology
is preserved. `ftbl_file` and `taxonly` remain accepted for compatibility: the
former remains unused, and the latter does not remove required statistics.
See `?blastn_taxo_assignment` and `?lcaPident` for complete details.
