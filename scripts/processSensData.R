# library(PharmacoGx)
# library(readr)
# library(tximport)
# library(rhdf5)
# library(gdata)
# library(readxl)
# library(openxlsx)
# library(CoreGx)
library(data.table)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/processed"

gCSI_GR_AOC <- read.csv(file.path(download_dir, "gCSI_GRvalues_v1.3.tsv"), sep = "\t")

## Using data.table for efficency and consistency with other psets. If you
## arent familiar with syntax (it implements inplace modifications), please
## read up on it.
gCSI_GR_AOC <- data.table(gCSI_GR_AOC)
gCSI_GR_AOC[, cellid := CellLineName]
gCSI_GR_AOC[, drugid := DrugName]

gCSI_GR_AOC[, expid := paste(cellid, drugid, ExperimentNumber, sep = "_")]
max_conc <- max(gCSI_GR_AOC[, .N, expid][, N])
uids <- unique(gCSI_GR_AOC$expid)

## Create matrices for the data
doses_final <- matrix(NA_real_, nrow = length(uids), ncol = max_conc, dimnames = list(uids, paste0("dose", seq_len(max_conc))))
viability_final <- matrix(NA_real_, nrow = length(uids), ncol = max_conc, dimnames = list(uids, paste0("dose", seq_len(max_conc))))

setorder(gCSI_GR_AOC, expid, log10Concentration)
gCSI_GR_AOC_list <- split(gCSI_GR_AOC, by = "expid")

for (exp in names(gCSI_GR_AOC_list)) {
  xx <- gCSI_GR_AOC_list[[exp]]
  concentrations <- 10^(xx$log10Concentration) # remove log10
  concentrations <- concentrations / 1.0e-6 # convert molar to uM
  viability <- xx$relative_cell_count * 100

  doses_final[exp, 1:length(concentrations)] <- concentrations
  viability_final[exp, 1:length(viability)] <- viability
}


raw.sensitivity <- array(c(as.matrix(as.numeric(doses_final)), as.matrix(as.numeric(viability_final))),
  c(length(uids), max_conc, 2),
  dimnames = list(
    rownames(doses_final),
    colnames(doses_final),
    c("Dose", "Viability")
  )
)



sensitivityInfo_2018 <- as.data.frame(gCSI_GR_AOC[, c("expid", "cellid", "drugid", "ExperimentNumber", "TrtDuration", "doublingtime")])
sensitivityInfo_2018 <- unique(sensitivityInfo_2018)
rownames(sensitivityInfo_2018) <- sensitivityInfo_2018$expid

## Read in published data.

gCSI_GR_AOC_Pub <- read.csv(file.path(download_dir, "gCSI_GRmetrics_v1.3.tsv"), sep = "\t")



uids2 <- unique(sprintf(
  "%s_%s_%s", gCSI_GR_AOC_Pub$CellLineName,
  gCSI_GR_AOC_Pub$DrugName,
  gCSI_GR_AOC_Pub$ExperimentNumber
))

rownames(gCSI_GR_AOC_Pub) <- uids2

sensitivity.info <- sensitivityInfo_2018
sensitivity.published <- gCSI_GR_AOC_Pub

save(raw.sensitivity, sensitivity.info, sensitivity.published, file = file.path(processed_dir, "sens.data.RData"))

raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity) / 1000))




dir.create(file.path(processed_dir, "slices"))

for (i in seq_along(raw.sensitivity.x)) {
  slce <- raw.sensitivity[raw.sensitivity.x[[i]], , ]
  saveRDS(slce, file = file.path(processed_dir, "slices", paste0("gcsi2017_raw_sens_", i, ".rds")))
}

zip(zipfile = file.path(processed_dir, "raw_sense_slices.zip"), files = list.files(file.path(processed_dir, "slices"), full.names = TRUE))
unlink(file.path(processed_dir, "slices"), recursive = TRUE)
