test_that("read_worms_ccz_evidence reads a small CSV fixture", {
  ccz <- data.frame(
    AphiaID = c("1001", "1002"),
    ScientificName = c("Aaptos kanuux", "Suberitidae"),
    AphiaID_accepted = c("1001", "1002"),
    ScientificName_accepted = c("Aaptos kanuux", "Suberitidae"),
    taxonRank = c("Species", "Family"),
    Kingdom = c("Animalia", "Animalia"),
    Phylum = c("Porifera", "Porifera"),
    Class = c("Demospongiae", "Demospongiae"),
    Order = c("Suberitida", "Suberitida"),
    Family = c("Suberitidae", "Suberitidae"),
    Genus = c("Aaptos", NA_character_),
    Marine = c(1, 1),
    Brackish = c(0, 1),
    Fresh = c(0, 0),
    Terrestrial = c(0, 0),
    taxonomicStatus = c("accepted", "accepted"),
    LSID = c("urn:lsid:marinespecies.org:taxname:1001", "urn:lsid:marinespecies.org:taxname:1002"),
    Citation = c(
      "WoRMS CCZ checklist https://example.org/ccz/aaptos",
      "WoRMS CCZ checklist"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  path <- tempfile(fileext = ".csv")
  utils::write.csv(ccz, path, row.names = FALSE, na = "")

  evidence <- read_worms_ccz_evidence(path, checked_at = "2026-06-20")

  expect_true(all(biohelper:::.taxon_evidence_required_columns() %in% colnames(evidence)))
  expect_equal(nrow(evidence), 2)
  expect_equal(evidence$taxon_name, c("Aaptos kanuux", "Suberitidae"))
  expect_equal(evidence$taxon_rank, c("Species", "Family"))
  expect_true(all(evidence$source == "worms_ccz"))
  expect_true(all(evidence$evidence_type == "regional_deepsea_checklist"))
  expect_true(all(evidence$habitat == "deep sea"))
  expect_true(all(evidence$region == "Clarion-Clipperton Zone"))
  expect_equal(evidence$environment, c("marine", "marine; brackish"))
  expect_equal(evidence$reference, c("WoRMS CCZ checklist https://example.org/ccz/aaptos", "WoRMS CCZ checklist"))
  expect_equal(evidence$reference_url[[1]], "https://example.org/ccz/aaptos")
  expect_equal(
    evidence$reference_url[[2]],
    "https://www.marinespecies.org/deepsea/CCZ/aphia.php?p=taxdetails&id=1002"
  )
  expect_equal(evidence$kingdom, c("Animalia", "Animalia"))
  expect_equal(evidence$phylum, c("Porifera", "Porifera"))
  expect_equal(evidence$class, c("Demospongiae", "Demospongiae"))
  expect_equal(evidence$order, c("Suberitida", "Suberitida"))
  expect_equal(evidence$family, c("Suberitidae", "Suberitidae"))
  expect_equal(evidence$genus, c("Aaptos", NA_character_))
  expect_equal(evidence$checked_at, c("2026-06-20", "2026-06-20"))
  expect_no_error(biohelper:::validate_taxon_evidence(evidence))
})
