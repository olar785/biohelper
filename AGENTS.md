# AGENTS.md

## Project

This repository contains the R package `biohelper`, a helper package for molecular ecology, metabarcoding, taxonomy handling, phyloseq/speedyseq workflows, visualisation, and ecological interpretation.

## General R package conventions

- Follow standard R package structure.
- Use roxygen2 documentation for exported functions.
- Add or update testthat tests for new user-facing behaviour.
- Use explicit namespaces for non-base R functions.
- Avoid adding heavy dependencies unless justified.
- Validate inputs early with informative errors.
- Do not silently modify user objects.

## Phyloseq/speedyseq conventions

- Support phyloseq/speedyseq objects where relevant.
- Do not modify phyloseq/speedyseq objects unless the function explicitly promises to do so.
- For taxonomic summaries and melting workflows, prefer speedyseq equivalents where available.
- Keep examples compact and biologically meaningful.

## Testing and documentation

Before considering work complete:

- Run `devtools::document()` after changing roxygen comments.
- Run `devtools::test()` after changing code.
- Run `devtools::check()` before finalising larger changes.
- Mock external API calls in tests.
- Tests should not require live API keys.
