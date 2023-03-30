options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/bhklab_orcestra/snakemake/PSet_gCSI2019/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/bhklab_orcestra/snakemake/PSet_gCSI2019/processed"

cell_all <- read.csv(file.path(download_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))
drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))


curationCell <- cell_all[apply(!is.na(cell_all[, c("gCSI.cellid", "GNE.cellid")]), 1, any), ]
curationTissue <- cell_all[apply(!is.na(cell_all[, c("gCSI.cellid", "GNE.cellid")]), 1, any), ]
curationCell <- curationCell[, c("unique.cellid", "gCSI.cellid", "GNE.cellid")]
curationTissue <- curationTissue[, c("unique.tissueid", "gCSI.tissueid", "GNE.tissueid")]

rownames(curationTissue) <- curationCell[, "unique.cellid"]
rownames(curationCell) <- curationCell[, "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[, "gCSI.drugid"])), ]
curationDrug <- curationDrug[, c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[, "unique.drugid"]

saveRDS(curationCell, file.path(processed_dir, 'curationCell.rds'))
saveRDS(curationTissue, file.path(processed_dir, 'curationTissue.rds'))
saveRDS(curationDrug, file.path(processed_dir, 'curationDrug.rds'))