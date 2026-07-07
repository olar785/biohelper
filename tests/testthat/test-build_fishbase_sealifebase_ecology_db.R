mock_fb_species <- function() {
  data.frame(
    SpecCode = c(1, 2),
    Genus = c("Thunnus", "Abyssogena"),
    Species = c("albacares", "southwardae"),
    FBname = c("Yellowfin tuna", NA_character_),
    Saltwater = c(1, 1),
    DemersPelag = c("pelagic-oceanic", NA_character_),
    DepthRangeShallow = c(1, 2800),
    DepthRangeDeep = c(1602, 6400),
    DepthRangeComShallow = c(1, NA),
    DepthRangeComDeep = c(250, NA),
    Comments = c("species note", NA_character_),
    GenCode = c(10, 20),
    FamCode = c(100, 200),
    Subfamily = c("Scombrinae", NA_character_),
    stringsAsFactors = FALSE
  )
}

mock_fb_taxonomy <- function() {
  list(
    genera = data.frame(
      GenCode = c(10, 20),
      GenName = c("Thunnus", "Abyssogena"),
      GenComName = c("tunas", "vent clams"),
      FamCode = c(100, 200),
      Subfamily = c("Scombrinae", NA_character_),
      Habitat = c("open ocean", "vent-associated"),
      WaterSalinity = c("marine", "marine"),
      Distribution = c("global warm seas", "deep sea vents"),
      Diagnosis = c("tuna genus diagnosis", "clam genus diagnosis"),
      Remark = c("tuna genus remark", "clam genus remark"),
      Etymology = c("Thunnus etymology", "Abyssogena etymology"),
      stringsAsFactors = FALSE
    ),
    families = data.frame(
      FamCode = c(100, 200),
      Family = c("Scombridae", "Vesicomyidae"),
      CommonName = c("mackerels", "clam family"),
      Order = c("Scombriformes", "Venerida"),
      Ordnum = c(1000, 2000),
      Class = c("Actinopterygii", "Bivalvia"),
      ClassNum = c(1, 2),
      Habitat = c("open water", "chemosynthetic habitats"),
      BodyShapeI = c("fusiform", "bivalved"),
      WaterSalinity = c("marine", "marine"),
      Remark = c("family remark", "clam family remark"),
      Etymology = c("Scombridae etymology", "Vesicomyidae etymology"),
      stringsAsFactors = FALSE
    ),
    orders = data.frame(
      Ordnum = c(1000, 2000),
      Order = c("Scombriformes", "Venerida"),
      CommonName = c("mackerel-like fishes", "venerids"),
      ClassNum = c(1, 2),
      Class = c("Actinopterygii", "Bivalvia"),
      WaterSalinity = c("marine", "marine"),
      ClassificationRemark = c("order classification remark", NA_character_),
      Etymology = c("Scombriformes etymology", "Venerida etymology"),
      stringsAsFactors = FALSE
    ),
    classes = data.frame(
      ClassNum = c(1, 2),
      Class = c("Actinopterygii", "Bivalvia"),
      CommonName = c("ray-finned fishes", "bivalves"),
      SuperClass = c("Osteichthyes", NA_character_),
      Subclass = c("Neopterygii", NA_character_),
      WaterSalinity = c("marine; freshwater", "marine"),
      Remarks = c("class remarks", "bivalve class remarks"),
      Etymology = c("Actinopterygii etymology", "Bivalvia etymology"),
      stringsAsFactors = FALSE
    )
  )
}

mock_slb_species <- function() {
  data.frame(
    SpecCode = c(1, 3),
    Genus = c("Psychropotes", "Calocalanus"),
    Species = c("longicauda", "pavo"),
    Saltwater = c(1, 1),
    DemersPelag = c("benthic", "pelagic"),
    DepthRangeDeep = c(5173, 200),
    Comments = c("sea cucumber note", "copepod note"),
    GenCode = c(30, 40),
    FamCode = c(300, 400),
    stringsAsFactors = FALSE
  )
}

mock_slb_taxonomy <- function() {
  list(
    genera = data.frame(
      GenCode = c(30, 40),
      GenName = c("Psychropotes", "Calocalanus"),
      CommonName = c("sea pigs", "copepods"),
      Famcode = c(300, 400),
      Subfamily = c(NA_character_, "Calocalaninae"),
      Habitat = c("abyssal benthos", "plankton"),
      WaterSalinity = c("marine", "marine"),
      Distribution = c("deep Pacific", "global ocean"),
      Diagnosis = c("sea pig diagnosis", "copepod diagnosis"),
      Comment = c("sea pig comment", "copepod comment"),
      stringsAsFactors = FALSE
    ),
    families = data.frame(
      FamCode = c(300, 400),
      Family = c("Psychropotidae", "Paracalanidae"),
      CommonName = c("psychropotids", "paracalanids"),
      Order = c("Holothuriida", "Calanoida"),
      Ordnum = c(3000, 4000),
      Class = c("Holothuroidea", "Hexanauplia"),
      ClassNum = c(30, 40),
      Phylum = c("Echinodermata", "Arthropoda"),
      Distribution = c("abyssal", "marine plankton"),
      Remark = c("family remarks", "copepod family remarks"),
      stringsAsFactors = FALSE
    ),
    orders = data.frame(
      Ordnum = c(3000, 4000),
      Order = c("Holothuriida", "Calanoida"),
      CommonName = c("sea cucumbers", "calanoid copepods"),
      Class = c("Holothuroidea", "Hexanauplia"),
      ClassNum = c(30, 40),
      Phylum = c("Echinodermata", "Arthropoda"),
      PhylumNum = c(3, 4),
      Remark = c("order remarks", "calanoid remarks"),
      Etymology = c("Holothuriida etymology", "Calanoida etymology"),
      stringsAsFactors = FALSE
    ),
    classes = data.frame(
      ClassNum = c(30, 40),
      Class = c("Holothuroidea", "Hexanauplia"),
      CommonName = c("sea cucumbers", "copepod class"),
      Phylum = c("Echinodermata", "Arthropoda"),
      PhylumNum = c(3, 4),
      WaterSalinity = c("marine", "marine"),
      Remarks = c("class remarks", "copepod class remarks"),
      Etymology = c("Holothuroidea etymology", "Hexanauplia etymology"),
      stringsAsFactors = FALSE
    ),
    phylums = data.frame(
      PhylumId = c(3, 4),
      Kingdom = c("Metazoa", "Metazoa"),
      Phylum = c("Echinodermata", "Arthropoda"),
      CommonName = c("echinoderms", "arthropods"),
      Etymology = c("Echinodermata etymology", "Arthropoda etymology"),
      stringsAsFactors = FALSE
    ),
    phylumclass = data.frame(
      Phylum = c("Echinodermata", "Arthropoda"),
      Class = c("Holothuroidea", "Hexanauplia"),
      stringsAsFactors = FALSE
    )
  )
}

test_that("clean_fishbase_species_table keeps species comments and renamed depth fields", {
  species <- mock_fb_species()
  species$Land <- NULL
  species$FBname <- NULL

  out <- biohelper:::clean_fishbase_species_table(species, source = "fishbase")

  expect_s3_class(out, "tbl_df")
  expect_equal(out$common_name, rep(NA_character_, 2))
  expect_equal(out$ecology_environment[[1]], "marine")
  expect_equal(out$ecology_position, c("pelagic", "unknown"))
  expect_equal(out$ecology_habitat_broad, c("pelagic", "unknown"))
  expect_equal(out$ecology_depth_zone[[1]], "bathyal_possible")
  expect_equal(out$ecology_depth_zone[[2]], "hadal_or_abyssal_possible")
  expect_equal(out$depth_min[[1]], 1)
  expect_equal(out$depth_max[[1]], 1602)
  expect_equal(out$common_depth_min[[1]], 1)
  expect_equal(out$common_depth_max[[1]], 250)
  expect_equal(out$species_comments[[1]], "species note")
  expect_false("comments" %in% colnames(out))
})

test_that("clean_fishbase_species_table handles FBname and casing in DemersPelag", {
  species <- data.frame(
    SpecCode = 10,
    Genus = "Gallus",
    Species = "gallus",
    FBname = "Chicken",
    Land = 1,
    DemersPelag = "Sessile",
    DepthRangeDeep = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- biohelper:::clean_fishbase_species_table(species, source = "sealifebase")

  expect_equal(out$common_name, "Chicken")
  expect_equal(out$taxon_name, "Gallus gallus")
  expect_equal(out$ecology_environment, "terrestrial")
  expect_equal(out$ecology_position, "benthic")
  expect_equal(out$ecology_habitat_broad, "benthic")
})

test_that("clean_fishbase_ecology_table collapses flags and keeps ecology comments", {
  ecology <- data.frame(
    SpecCode = c(1, 1, 2, 999),
    Neritic = c(-1, 0, NA, 0),
    Oceanic = c(0, -1, 0, 0),
    Epipelagic = c(-1, 0, 0, 0),
    Mesopelagic = c(0, -1, 0, 0),
    Benthic = c(-1, 0, 0, 0),
    Demersal = c(0, -1, 0, 0),
    Pelagic = c(0, 0, -1, 0),
    Sessile = c(-1, 0, 0, 0),
    Mobile = c(0, -1, 0, 0),
    Endofauna = c(-1, 0, 0, 0),
    Megabenthos = c(0, -1, 0, 0),
    SoftBottom = c(-1, 0, 0, 0),
    Mud = c(0, -1, 0, 0),
    Seamounts = c(-1, 0, 0, 0),
    DeepWaterCorals = c(0, -1, 0, 0),
    AddRems = c("first note", "second note", NA_character_, "unmatched note"),
    stringsAsFactors = FALSE
  )

  out <- biohelper:::clean_fishbase_ecology_table(ecology, source = "fishbase")
  row1 <- out[out$source == "fishbase" & out$spec_code == "1", , drop = FALSE]

  expect_equal(nrow(out), 3)
  expect_equal(row1$ecology_zone, "neritic; oceanic")
  expect_equal(row1$ecology_water_column_zone, "epipelagic; mesopelagic")
  expect_false("ecology_habitat_broad" %in% colnames(out))
  expect_equal(row1$ecology_mobility, "mobile; sessile")
  expect_equal(row1$ecology_size_class, "endofauna; megabenthos")
  expect_equal(row1$ecology_substrate, "mud; soft_bottom")
  expect_equal(row1$ecology_special_habitat, "deep_water_corals; seamounts")
  expect_equal(row1$ecology_comments, "first note; second note")
  expect_false("add_rems" %in% colnames(out))
  expect_false(any(c("benthic", "demersal", "sessile", "soft_bottom", "seamounts") %in% colnames(out)))
})

test_that("native FishBase and SeaLifeBase taxonomy joins correctly", {
  fb_tax <- mock_fb_taxonomy()
  slb_tax <- mock_slb_taxonomy()

  fb <- biohelper:::clean_fishbase_taxonomy_table(
    species = mock_fb_species(),
    genera = fb_tax$genera,
    families = fb_tax$families,
    orders = fb_tax$orders,
    classes = fb_tax$classes
  )
  slb <- biohelper:::clean_sealifebase_taxonomy_table(
    species = mock_slb_species(),
    genera = slb_tax$genera,
    families = slb_tax$families,
    orders = slb_tax$orders,
    classes = slb_tax$classes,
    phylums = slb_tax$phylums,
    phylumclass = slb_tax$phylumclass
  )

  expect_equal(fb$kingdom, rep("Metazoa", 2))
  expect_equal(fb$phylum, rep("Chordata", 2))
  expect_equal(fb$taxonomy_source, rep("fishbase_native_inferred_phylum", 2))
  expect_equal(fb$family[[1]], "Scombridae")
  expect_equal(fb$class[[1]], "Actinopterygii")
  expect_equal(fb$superclass[[1]], "Osteichthyes")

  expect_equal(slb$kingdom, rep("Metazoa", 2))
  expect_equal(slb$phylum, c("Echinodermata", "Arthropoda"))
  expect_equal(slb$taxonomy_source, rep("sealifebase_native", 2))
  expect_equal(slb$order[[2]], "Calanoida")
  expect_equal(slb$family[[1]], "Psychropotidae")

  fb_meta <- attr(fb, "rank_metadata", exact = TRUE)
  slb_meta <- attr(slb, "rank_metadata", exact = TRUE)
  fb_genus <- fb_meta[fb_meta$taxon_rank == "genus" & fb_meta$taxon_name == "Thunnus", , drop = FALSE]
  fb_family <- fb_meta[fb_meta$taxon_rank == "family" & fb_meta$taxon_name == "Scombridae", , drop = FALSE]
  fb_order <- fb_meta[fb_meta$taxon_rank == "order" & fb_meta$taxon_name == "Scombriformes", , drop = FALSE]
  slb_class <- slb_meta[slb_meta$taxon_rank == "class" & slb_meta$taxon_name == "Holothuroidea", , drop = FALSE]
  slb_phylum <- slb_meta[slb_meta$taxon_rank == "phylum" & slb_meta$taxon_name == "Echinodermata", , drop = FALSE]

  expect_equal(fb_genus$common_name, "tunas")
  expect_equal(fb_genus$rank_habitat_note, "open ocean")
  expect_equal(fb_genus$rank_etymology, "Thunnus etymology")
  expect_equal(fb_family$rank_body_shape, "fusiform")
  expect_equal(fb_family$rank_remarks, "family remark")
  expect_equal(fb_family$rank_etymology, "Scombridae etymology")
  expect_equal(fb_order$rank_remarks, "order classification remark")
  expect_equal(fb_order$rank_etymology, "Scombriformes etymology")
  expect_equal(slb_class$rank_water_salinity, "marine")
  expect_equal(slb_class$rank_remarks, "class remarks")
  expect_equal(slb_class$rank_etymology, "Holothuroidea etymology")
  expect_equal(slb_phylum$common_name, "echinoderms")
  expect_equal(slb_phylum$rank_etymology, "Echinodermata etymology")
})

test_that("ecology_habitat_broad is derived only from species DemersPelag", {
  species <- data.frame(
    SpecCode = c(1, 2, 3),
    Genus = c("A", "B", "C"),
    Species = c("one", "two", "three"),
    Saltwater = 1,
    DemersPelag = c("pelagic", "demersal", NA_character_),
    stringsAsFactors = FALSE
  )
  ecology <- data.frame(
    SpecCode = c(1, 2, 3),
    Benthic = c(-1, 0, 0),
    Pelagic = c(0, -1, -1),
    stringsAsFactors = FALSE
  )

  out <- biohelper:::.join_fishbase_species_ecology(
    species = species,
    ecology = ecology,
    taxonomy = biohelper:::.empty_fishbase_taxonomy_table("fishbase"),
    source = "fishbase"
  )

  expect_equal(out$ecology_habitat_broad, c("pelagic", "demersal", "unknown"))
})

test_that("ecology_depth_zone follows depth_max thresholds", {
  expect_equal(
    biohelper:::.fishbase_depth_zone(c(NA, 199, 200, 999, 1000, 2999, 3000, 5999, 6000)),
    c(
      "unknown",
      "shallow_only",
      "mesophotic_or_upper_bathyal_possible",
      "mesophotic_or_upper_bathyal_possible",
      "bathyal_possible",
      "bathyal_possible",
      "abyssal_possible",
      "abyssal_possible",
      "hadal_or_abyssal_possible"
    )
  )
})

test_that("ecology evidence strength uses documented habitat support rule", {
  expect_equal(biohelper:::.fishbase_evidence_strength(0, NA_real_, NA_character_), "unknown")
  expect_equal(biohelper:::.fishbase_evidence_strength(10, 0.9, "pelagic"), "high")
  expect_equal(biohelper:::.fishbase_evidence_strength(5, 0.6, "pelagic; demersal"), "moderate")
  expect_equal(biohelper:::.fishbase_evidence_strength(10, 0.5, "pelagic; demersal"), "mixed")
  expect_equal(biohelper:::.fishbase_evidence_strength(4, 1, "benthic"), "low")
})

test_that("fishbase ecology aggregation de-duplicates terms and keeps one dominant value", {
  collapsed <- biohelper:::.collapse_fishbase_text_values(c(
    "benthic; host_associated; pelagic; benthic",
    "host_associated; unknown; pelagic"
  ))
  dominant <- biohelper:::.dominant_fishbase_value(c(
    "pelagic; benthic",
    "benthic",
    "unknown"
  ))
  tied <- biohelper:::.dominant_fishbase_value(c("pelagic", "benthic"))
  missing <- biohelper:::.dominant_fishbase_value(c(NA_character_, ""))

  expect_equal(collapsed, "benthic; host_associated; pelagic; unknown")
  expect_equal(dominant$value, "benthic")
  expect_equal(dominant$prop, 0.5)
  expect_false(grepl(";", dominant$value, fixed = TRUE))
  expect_equal(tied$value, "benthic")
  expect_equal(tied$prop, 0.5)
  expect_true(is.na(missing$value))
  expect_true(is.na(missing$prop))
})

test_that("fishbase ecology combiner keeps source/spec_code keys distinct and comments separate", {
  fb_tax <- mock_fb_taxonomy()
  slb_tax <- mock_slb_taxonomy()
  fb_ecology <- data.frame(
    SpecCode = c(1, 2, 2, 999),
    Pelagic = c(-1, 0, 0, -1),
    Benthic = c(0, -1, 0, 0),
    Demersal = c(0, 0, -1, 0),
    SoftBottom = c(0, -1, 0, 0),
    HydrothermalVents = c(0, -1, 0, 0),
    AddRems = c("pelagic ecology", "benthic ecology", "demersal ecology", "unmatched ecology"),
    stringsAsFactors = FALSE
  )
  slb_ecology <- data.frame(
    SpecCode = c(1, 1),
    Benthic = c(-1, 0),
    Pelagic = c(0, -1),
    Mud = c(-1, 0),
    Mobile = c(-1, 0),
    Meiobenthos = c(0, -1),
    AddRems = c("deep-sea note", NA_character_),
    stringsAsFactors = FALSE
  )

  out <- biohelper:::.combine_fishbase_sealifebase_ecology_tables(
    fishbase_species = mock_fb_species(),
    fishbase_ecology = fb_ecology,
    fishbase_genera = fb_tax$genera,
    fishbase_families = fb_tax$families,
    fishbase_orders = fb_tax$orders,
    fishbase_classes = fb_tax$classes,
    sealifebase_species = mock_slb_species(),
    sealifebase_ecology = slb_ecology,
    sealifebase_genera = slb_tax$genera,
    sealifebase_families = slb_tax$families,
    sealifebase_orders = slb_tax$orders,
    sealifebase_classes = slb_tax$classes,
    sealifebase_phylums = slb_tax$phylums,
    sealifebase_phylumclass = slb_tax$phylumclass
  )
  species <- out$species
  report <- out$diagnostics

  expect_equal(nrow(species), nrow(mock_fb_species()) + nrow(mock_slb_species()))
  expect_equal(report$expected_rows, 4)
  expect_equal(report$final_rows, 4)
  expect_equal(report$fishbase_ecology_unmatched_spec_code, 1)
  expect_equal(report$sealifebase_ecology_unmatched_spec_code, 0)
  expect_equal(report$fishbase_taxonomy_unmatched_spec_code, 0)
  expect_equal(report$sealifebase_taxonomy_unmatched_spec_code, 0)

  expected_columns <- c(
    "source", "spec_code", "taxon_rank", "taxon_name", "kingdom", "phylum",
    "superclass", "class", "subclass", "order", "family", "subfamily",
    "genus", "species", "common_name", "taxonomy_source",
    "ecology_environment", "ecology_position", "ecology_habitat_broad",
    "ecology_depth_zone", "depth_min", "depth_max", "common_depth_min",
    "common_depth_max", "ecology_zone", "ecology_water_column_zone",
    "ecology_mobility", "ecology_size_class", "ecology_substrate",
    "ecology_special_habitat", "raw_demers_pelag", "species_comments",
    "ecology_comments"
  )
  expect_equal(colnames(species), expected_columns)

  fb_one <- species[species$source == "fishbase" & species$spec_code == "1", , drop = FALSE]
  slb_one <- species[species$source == "sealifebase" & species$spec_code == "1", , drop = FALSE]
  expect_equal(fb_one$taxon_name, "Thunnus albacares")
  expect_equal(fb_one$family, "Scombridae")
  expect_equal(fb_one$ecology_habitat_broad, "pelagic")
  expect_equal(fb_one$species_comments, "species note")
  expect_equal(fb_one$ecology_comments, "pelagic ecology")
  expect_equal(slb_one$taxon_name, "Psychropotes longicauda")
  expect_equal(slb_one$phylum, "Echinodermata")
  expect_equal(slb_one$ecology_habitat_broad, "benthic")
  expect_equal(slb_one$ecology_mobility, "mobile")
  expect_equal(slb_one$ecology_size_class, "meiobenthos")
  expect_equal(slb_one$ecology_comments, "deep-sea note")

  higher <- build_fishbase_sealifebase_higher_rank_ecology_db(species)
  slb_phylum <- higher[
    higher$source == "sealifebase" &
      higher$taxon_rank == "phylum" &
      higher$taxon_name == "Echinodermata",
    ,
    drop = FALSE
  ]
  fb_phylum <- higher[
    higher$source == "fishbase" &
      higher$taxon_rank == "phylum" &
      higher$taxon_name == "Chordata",
    ,
    drop = FALSE
  ]
  expect_equal(slb_phylum$common_name, "echinoderms")
  expect_false(all(is.na(higher$common_name)))
  expect_true(is.na(fb_phylum$common_name))
})

test_that("higher-rank aggregation collapses ecology and preserves rank metadata", {
  species <- data.frame(
    source = rep("fishbase", 10),
    spec_code = as.character(seq_len(10)),
    taxon_rank = "species",
    taxon_name = paste("Testus", paste0("sp", seq_len(10))),
    kingdom = "Metazoa",
    phylum = "Chordata",
    superclass = NA_character_,
    class = "Actinopterygii",
    subclass = NA_character_,
    order = "Testiformes",
    family = "Testidae",
    subfamily = NA_character_,
    genus = "Testus",
    species = paste0("sp", seq_len(10)),
    common_name = NA_character_,
    taxonomy_source = "fishbase_native_inferred_phylum",
    ecology_environment = c(rep("marine", 9), "marine; brackish"),
    ecology_position = rep("benthic", 10),
    ecology_habitat_broad = rep("benthic", 10),
    ecology_depth_zone = c(rep("bathyal_possible", 8), "abyssal_possible", "unknown"),
    depth_min = rep(100, 10),
    depth_max = seq(1000, 1900, by = 100),
    common_depth_min = NA_real_,
    common_depth_max = NA_real_,
    ecology_zone = c(rep("oceanic", 8), "neritic", NA_character_),
    ecology_water_column_zone = c(rep("mesopelagic", 5), rep(NA_character_, 5)),
    ecology_mobility = NA_character_,
    ecology_size_class = NA_character_,
    ecology_substrate = NA_character_,
    ecology_special_habitat = NA_character_,
    raw_demers_pelag = NA_character_,
    species_comments = NA_character_,
    ecology_comments = NA_character_,
    stringsAsFactors = FALSE
  )
  species <- tibble::as_tibble(species)
  attr(species, "rank_metadata") <- tibble::tibble(
    source = "fishbase",
    taxon_rank = "family",
    taxon_name = "Testidae",
    common_name = "test family",
    rank_water_salinity = NA_character_,
    rank_body_shape = NA_character_,
    rank_habitat_note = "metadata habitat note",
    rank_distribution = NA_character_,
    rank_diagnosis = NA_character_,
    rank_remarks = "metadata remarks"
  )

  higher <- build_fishbase_sealifebase_higher_rank_ecology_db(species)
  family <- higher[higher$taxon_rank == "family" & higher$taxon_name == "Testidae", , drop = FALSE]

  expect_setequal(unique(higher$taxon_rank), c("kingdom", "phylum", "class", "order", "family", "genus"))
  expect_equal(family$ecology_environment, "brackish; marine")
  expect_equal(family$ecology_habitat_broad, "benthic")
  expect_equal(family$n_species, 10)
  expect_equal(family$n_species_with_habitat_broad, 10)
  expect_equal(family$dominant_habitat_broad, "benthic")
  expect_equal(family$dominant_habitat_broad_prop, 1)
  expect_equal(family$ecology_evidence_strength, "high")
  expect_equal(family$common_name, "test family")
  expect_equal(family$rank_habitat_note, "metadata habitat note")
  expect_true("rank_remarks" %in% colnames(higher))

  combined <- biohelper:::.combine_species_and_higher_rank_fishbase_db(species, higher)
  expect_setequal(
    unique(combined$taxon_rank),
    c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  )
  expect_equal(
    setdiff(
      unique(combined$taxon_rank),
      c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    ),
    character()
  )
})

test_that("fishbase combined ecology DB filters non-standard taxon_rank rows", {
  species <- data.frame(
    source = "fishbase",
    spec_code = "1",
    taxon_rank = "species",
    taxon_name = "Testus alpha",
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Testiformes",
    family = "Testidae",
    genus = "Testus",
    species = "alpha",
    stringsAsFactors = FALSE
  )
  species_extra <- species
  species_extra$taxon_rank <- "subclass"
  species_extra$taxon_name <- "Actinopteri"
  higher <- data.frame(
    source = "fishbase",
    taxon_rank = c("family", "subfamily"),
    taxon_name = c("Testidae", "Testinae"),
    kingdom = "Metazoa",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Testiformes",
    family = c("Testidae", "Testidae"),
    genus = NA_character_,
    stringsAsFactors = FALSE
  )

  combined <- biohelper:::.combine_species_and_higher_rank_fishbase_db(
    dplyr::bind_rows(species, species_extra),
    higher
  )

  expect_equal(
    setdiff(
      unique(combined$taxon_rank),
      c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    ),
    character()
  )
  expect_false(any(combined$taxon_rank %in% c("superclass", "subclass", "subfamily", "superfamily", "tribe")))
  expect_true(any(combined$taxon_rank == "family"))
  expect_true(any(combined$taxon_rank == "species"))
})
