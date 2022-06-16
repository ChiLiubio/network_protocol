
# use pacman package to load or install packages
tmp <- "pacman"
if(!require(tmp, character.only = TRUE)) {
	install.packages(tmp, dependencies = TRUE)
	require(tmp, character.only = TRUE)
}

p_load(curatedMetagenomicData, dplyr, mia)

alcoholStudy <-
    filter(sampleMetadata, age >= 18) |>
    filter(!is.na(alcohol)) |>
    filter(body_site == "stool") |>
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance", counts = TRUE, rownames = "short")

View(alcoholStudy@assays@data$relative_abundance)
# convert TreeSummarizedExperiment object to phyloseq object
tmp <- makePhyloseqFromTreeSummarizedExperiment(alcoholStudy, abund_values = "relative_abundance")

# load file2meco package
p_load(file2meco)
stool_raw <- phyloseq2meco(tmp)
save(stool_raw, file = "stool_raw.RData")

# filter data
p_load(microeco, magrittr)

load("stool_raw.RData")
stool_met <- clone(stool_raw)

# generate species name in tax_table
tmp <- stool_met$tax_table
rownames(tmp) %<>% gsub("\\[|\\]", "", .)
tmp$species <- paste0("s__", rownames(tmp))
stool_met$tax_table <- tmp

# rename taxonomic names
colnames(stool_met$tax_table) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
tmp <- ifelse(stool_met$tax_table$Phylum == "p__Euryarchaeota", "k__Archaea", "k__Bacteria")
stool_met$tax_table <- data.frame(Kingdom = tmp, stool_met$tax_table)

# filter the species with low abundance
stool_met$otu_table %<>% .[apply(., 1, sum)/sum(.) > 0.0005, ]
stool_met$tidy_dataset()

# filter two species with low occurrence frequency
stool_met$otu_table %<>% .[rownames(.) != "Prevotella sp. CAG:1124", ]
stool_met$otu_table %<>% .[rownames(.) != "Bacteroides sp. CAG:927", ]
stool_met$tidy_dataset()

# use samples have 0 count in less than 70 species
tmp <- apply(stool_met$otu_table, 2, function(x){sum(x == 0)})
tmp <- tmp[tmp < 70]
stool_met$otu_table %<>% .[, colnames(.) %in% names(tmp)]
stool_met$tidy_dataset()
# otu_table have 144 rows and 198 columns

# use species have 0 count in less than 120 samples
tmp <- apply(stool_met$otu_table, 1, function(x){sum(x == 0)})
tmp <- tmp[tmp < 120]
stool_met$otu_table %<>% .[rownames(.) %in% names(tmp), ]
stool_met$tidy_dataset()


save(stool_met, file = "stool_met.RData")
rm(stool_raw)

