library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
# library(gdata)
library(readxl)
library(openxlsx)
library(Biobase)
library(reshape2)
library(abind)
library(data.table)
library(parallel)
library(CoreGx)
library(SummarizedExperiment)

verbose <- FALSE
nthread <- 1

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")
output_dir <- args[[1]]
filename <- args[[2]]
tools <- args[[3]]
transcriptome <- args[[4]]

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/processed"
# output_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/"
# tools <- "Kallisto-0.46.1"
# transcriptome <- "Gencode_v33"

tools <- gsub("-", "_", tools)
# untar(file.path(download_dir, paste0(tools, ".tar.gz")), exdir = download_dir)
tool_path <- expand.grid(a = tools, b = transcriptome)
tool_path <- paste0(tool_path$a, "_", tool_path$b)

# untar(tarfile = file.path(download_dir, paste0(tools, ".tar.gz")), exdir = processed_dir)

rnaseq_results <- readRDS(file.path(processed_dir, 'rnaseq_results.rds'))
curationCell <- readRDS(file.path(processed_dir, 'curationCell.rds'))
curationTissue <- readRDS(file.path(processed_dir, 'curationTissue.rds'))
curationDrug <- readRDS(file.path(processed_dir, 'curationDrug.rds'))

# myDirPrefix <- "/pfs/"
# args = commandArgs(trailingOnly=TRUE)
# rnaseq_select <- args
# print(rnaseq_select)
# ORCESTRA_ID = tail(rnaseq_select, n=1)


# tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
# tools <- rnaseq_select[tools]
# tools <- gsub("-", "_", tools)
# transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
# transcriptome <- rnaseq_select[transcriptome]
# tool_path = expand.grid(a = tools,b = transcriptome)
# tool_path = paste0(tool_path$a, "_",tool_path$b)

print(tool_path)

matchToIDTable <- function(ids, tbl, column, returnColumn = "unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)", Hmisc::escapeRegex(trimws(x)), "((///)|$)"), tbl[, column])
    if (length(myx) > 1) {
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if (length(myx) == 0) {
      return(NA_character_)
    }
    return(tbl[myx, returnColumn])
  })
}

load(file.path(processed_dir, "profiles.RData"))
load(file.path(processed_dir, "sens.data.RData"))
load(file.path(download_dir, "gCSI_molData.RData"))

cell_all <- read.csv(file.path(download_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))
drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))

# curationCell <- cell_all[apply(!is.na(cell_all[, c("gCSI.cellid", "GNE.cellid")]), 1, any), ]
# curationTissue <- cell_all[apply(!is.na(cell_all[, c("gCSI.cellid", "GNE.cellid")]), 1, any), ]
# curationCell <- curationCell[, c("unique.cellid", "gCSI.cellid", "GNE.cellid")]
# curationTissue <- curationTissue[, c("unique.tissueid", "gCSI.tissueid", "GNE.tissueid")]
# 
# rownames(curationTissue) <- curationCell[, "unique.cellid"]
# rownames(curationCell) <- curationCell[, "unique.cellid"]
# 
# curationDrug <- drug_all[which(!is.na(drug_all[, "gCSI.drugid"])), ]
# curationDrug <- curationDrug[, c("unique.drugid", "gCSI.drugid")]
# rownames(curationDrug) <- curationDrug[, "unique.drugid"]

stopifnot(setequal(rownames(sensitivity.info), rownames(sensitivity.published)))
stopifnot(setequal(rownames(raw.sensitivity), rownames(sensitivity.published)))
stopifnot(setequal(rownames(res), rownames(sensitivity.published)))


recomputed_2018 <- res

sensitivityProfiles_2018 <- data.frame("aac_recomputed" = NA, "ic50_recomputed" = NA, "HS" = NA, "E_inf" = NA, "EC50" = NA)
sensitivityProfiles_2018[nrow(sensitivity.info), ] <- NA
sensitivityProfiles_2018[, "aac_recomputed"] <- as.numeric(recomputed_2018[, "aac_recomputed"])
sensitivityProfiles_2018[, "ic50_recomputed"] <- as.numeric((recomputed_2018[, "ic50_recomputed"]))
sensitivityProfiles_2018[, "HS"] <- as.numeric((recomputed_2018[, "HS"]))
sensitivityProfiles_2018[, "E_inf"] <- as.numeric((recomputed_2018[, "E_inf"]))
sensitivityProfiles_2018[, "EC50"] <- as.numeric((recomputed_2018[, "EC50"]))

rownames(sensitivityProfiles_2018) <- rownames(recomputed_2018)

sensitivityProfiles_2018 <- sensitivityProfiles_2018[rownames(sensitivity.info), ]
sensitivity.published <- sensitivity.published[rownames(sensitivity.info), ]
raw.sensitivity <- raw.sensitivity[rownames(sensitivity.info), , ]

sensitivityProfiles_2018$GR_AOC_published <- sensitivity.published[, c("GR_AOC")]
sensitivityProfiles_2018$meanviability_published <- sensitivity.published[, c("meanviability")]
sensitivityProfiles_2018$GRmax_published <- sensitivity.published[, c("GRmax")]
sensitivityProfiles_2018$Emax_published <- sensitivity.published[, c("Emax")]
sensitivityProfiles_2018$GRinf_published <- sensitivity.published[, c("GRinf")]
sensitivityProfiles_2018$GEC50_published <- sensitivity.published[, c("GEC50")]
sensitivityProfiles_2018$GR50_published <- sensitivity.published[, c("GR50")]


## Only doing this for data added to the pset.

reps <- matchToIDTable(mut$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
mut$cellid <- reps
# mut$tissueid <- curationTissue[mut$cellid, "unique.tissueid"]

reps <- matchToIDTable(cnv$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
cnv$cellid <- reps
# cnv$tissueid <- curationTissue[cnv$cellid, "unique.tissueid"]

message("mapping sensitivity cells")
mapInfo <- data.frame(
  gCSI.cellid = unique(sensitivity.info$cellid),
  unique.cellid = matchToIDTable(unique(sensitivity.info$cellid), curationCell, "gCSI.cellid", "unique.cellid")
)
stopifnot(!anyNA(mapInfo[, 2]))
sensitivity.info$cellid <- mapInfo[match(sensitivity.info$cellid, mapInfo[, 1]), 2]

message("mapping sensitivity drugs")

mapInfo <- data.frame(
  gCSI.drugid = unique(sensitivity.info$drugid),
  unique.drugid = matchToIDTable(unique(sensitivity.info$drugid), curationDrug, "gCSI.drugid", "unique.drugid")
)
stopifnot(!anyNA(mapInfo[, 2]))
sensitivity.info$drugid <- mapInfo[match(sensitivity.info$drugid, mapInfo[, 1]), 2]

# summarizeRnaSeq <- function(dir,
#                             features_annotation,
#                             samples_annotation,
#                             method) {
#   library(Biobase)
#   library(readr)
#   library(tximport)
# 
#   load(features_annotation)
# 
#   tx2gene <- as.data.frame(cbind("transcript" = tx2gene$transcripts, "gene" = tx2gene$genes))
# 
#   files <- list.files(dir, recursive = TRUE, full.names = T)
#   if (method == "kallisto") {
#     resFiles <- grep("abundance.h5", files)
#   } else {
#     resFiles <- grep("quant.sf", files)
#   }
#   resFiles <- files[resFiles]
#   length(resFiles)
#   names(resFiles) <- basename(dirname(resFiles))
# 
#   if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
#     txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
#   } else {
#     txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)
#   }
# 
#   head(txi$counts[, 1:5])
#   dim(txi$counts)
# 
#   xx <- txi$abundance
#   gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
#   fData(gene.exp) <- features_gene[featureNames(gene.exp), ]
#   pData(gene.exp) <- samples_annotation[sampleNames(gene.exp), ]
#   annotation(gene.exp) <- "rnaseq"
# 
#   xx <- txi$counts
#   gene.count <- Biobase::ExpressionSet(log2(xx + 1))
#   fData(gene.count) <- features_gene[featureNames(gene.count), ]
#   pData(gene.count) <- samples_annotation[sampleNames(gene.count), ]
#   annotation(gene.count) <- "rnaseq"
# 
#   txii <- tximport(resFiles, type = method, txOut = T)
# 
#   if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
#     # remove non-coding transcripts in ensembl
#     rownames(txii$abundance) <- gsub("\\..*", "", rownames(txii$abundance))
#     txii$abundance[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
#     missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
#     txii$abundance <- txii$abundance[-which(rownames(txii$abundance) %in% missing_transcript), ]
#   }
# 
#   xx <- txii$abundance
#   transcript.exp <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 0.001))
#   if (features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")) {
#     featureNames(transcript.exp) <- gsub("\\|.*", "", featureNames(transcript.exp))
#     fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
#   } else {
#     fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
#   }
#   pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp), ]
#   annotation(transcript.exp) <- "isoforms"
# 
# 
#   if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
#     # remove non-coding transcripts in ensembl
#     rownames(txii$counts) <- gsub("\\..*", "", rownames(txii$counts))
#     txii$counts <- txii$counts[-which(rownames(txii$counts) %in% missing_transcript), ]
#   }
#   xx <- txii$counts
#   transcript.count <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 1))
#   if (features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")) {
#     featureNames(transcript.count) <- gsub("\\|.*", "", featureNames(transcript.count))
#     fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
#   } else {
#     fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
#   }
#   pData(transcript.count) <- samples_annotation[sampleNames(transcript.count), ]
#   annotation(transcript.count) <- "isoforms"
# 
#   return(list(
#     "rnaseq" = gene.exp,
#     "rnaseq.counts" = gene.count,
#     "isoforms" = transcript.exp,
#     "isoforms.counts" = transcript.count
#   ))
# }
# 
# 
# rnaseq.sampleinfo <- read.csv(file = file.path(download_dir, "gCSI_rnaseq_meta.csv"), stringsAsFactors = FALSE, row.names = 1)
# rnaseq.sampleinfo[, "cellid"] <- matchToIDTable(ids = rnaseq.sampleinfo[, "Cell_line"], tbl = curationCell, column = "GNE.cellid", returnColumn = "unique.cellid")
# 
# 
# for (r in 1:length(tool_path)) {
#   print(tool_path[r])
#   if (length(grep(pattern = "Kallisto", x = tool_path[r])) > 0) {
#     tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
#     # tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
#     rnatool <- "kallisto"
#   } else {
#     tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
#     # tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
#     rnatool <- "salmon"
#   }
# 
# 
#   if (length(grep(pattern = "lift37", x = tool_path[r])) > 0) {
#     annot <- file.path(download_dir, "Gencode.v33lift37.annotation.RData")
#   } else if (length(grep(pattern = "v33", x = tool_path[r])) > 0) {
#     annot <- file.path(download_dir, "Gencode.v33.annotation.RData")
#   } else {
#     annot <- file.path(download_dir, "Ensembl.v99.annotation.RData")
#   }
#   print(annot)
# 
#   rnaseq <- summarizeRnaSeq(
#     dir = file.path(processed_dir, tool, tool_path[r]),
#     features_annotation = annot,
#     samples_annotation = rnaseq.sampleinfo,
#     method = rnatool
#   )
# 
# 
#   reps <- matchToIDTable(rnaseq$rnaseq$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
#   stopifnot(!anyNA(reps))
#   rnaseq$rnaseq$cellid <- reps
# 
#   reps <- matchToIDTable(rnaseq$rnaseq.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
#   stopifnot(!anyNA(reps))
#   rnaseq$rnaseq.counts$cellid <- reps
# 
#   reps <- matchToIDTable(rnaseq$isoforms$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
#   stopifnot(!anyNA(reps))
#   rnaseq$isoforms$cellid <- reps
# 
#   reps <- matchToIDTable(rnaseq$isoforms.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
#   stopifnot(!anyNA(reps))
#   rnaseq$isoforms.counts$cellid <- reps
# 
#   rnaseq_results <- c(rnaseq_results, c(
#     rnaseq <- setNames(rnaseq, paste0(tool, ".", names(rnaseq)))
#   ))
# }

rnaseq_cellid_all <- pData(rnaseq_results[[1]])[, "cellid"]
reps <- matchToIDTable(rownames(cellInfo), curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rownames(cellInfo) <- reps

# add missing cell lines to cell info
cellInfo <- as.data.frame(cellInfo)

cellinall <- CoreGx::.unionList(rnaseq_cellid_all, sensitivity.info$cellid, cnv$cellid, mut$cellid)

newCells <- setdiff(cellinall, rownames(cellInfo))

newrows <- cellInfo[newCells, ]
rownames(newrows) <- newCells
cellInfo <- rbind(cellInfo, newrows)

cellInfo$tissueid <- curationTissue[rownames(cellInfo), "unique.tissueid"]
cellInfo$cellid <- rownames(cellInfo)
curationTissue <- curationTissue[rownames(cellInfo), ]
curationCell <- curationCell[rownames(cellInfo), ]

print("rows of cellinfo")
print(nrow(cellInfo))

# removed_exp_ids <- sub("__[^_]+$", "", rownames(sensitivityProfiles_2018))

# tt <- removed_exp_ids
# for(i in 1:length(tt)) {
#   xx <- which(tt == tt[i])
#   if(length(xx) > 1) {
#     for(j in 1:length(xx)) {
#       tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
#     }
#   }
# }
# tt <- sub("__", "_", tt)

# rownames(sensitivityProfiles_2018) <- tt

# removed_exp_ids <- sub("__[^_]+$", "", rownames(sensitivityRaw_2018))

# tt <- removed_exp_ids
# for(i in 1:length(tt)) {
#   xx <- which(tt == tt[i])
#   if(length(xx) > 1) {
#     for(j in 1:length(xx)) {
#       tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
#     }
#   }
# }
# tt <- sub("__", "_", tt)

# rownames(sensitivityRaw_2018) <- tt


# removed_exp_ids <- sub("__[^_]+$", "", rownames(sensitivityInfo_2018))

# tt <- removed_exp_ids
# for(i in 1:length(tt)) {
#   xx <- which(tt == tt[i])
#   if(length(xx) > 1) {
#     for(j in 1:length(xx)) {
#       tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
#     }
#   }
# }
# tt <- sub("__", "_", tt)

# rownames(sensitivityInfo_2018) <- tt


drugsPresent <- sort(unique(sensitivity.info$drugid))

druginfo <- curationDrug[drugsPresent, ]

drug_all <- read.csv(file.path(download_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[, "gCSI.drugid"])), ]
drug_all <- drug_all[, c("unique.drugid", "gCSI.drugid", "smiles", "inchikey", "cid", "FDA")]
rownames(drug_all) <- drug_all[, "unique.drugid"]

drug_all <- drug_all[rownames(druginfo), ]
druginfo[, c("smiles", "inchikey", "cid", "FDA")] <- drug_all[, c("smiles", "inchikey", "cid", "FDA")]
druginfo$gCSI.drugid <- NULL
colnames(druginfo)[which(names(druginfo) == "unique.drugid")] <- "drugid"
curationDrug <- curationDrug[rownames(druginfo), ]

z <- list()

z <- c(z, c(
  rnaseq_results,
  "cnv" = cnv,
  "mutation" = mut
))

.converteSetToSE <- function(eSets) {
  SEfinal <- lapply(
    eSets,
    function(eSet) {
      # Change rownames from probes to EnsemblGeneId for rna data type
      if (grepl("^rna$", Biobase::annotation(eSet))) {
        rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
      }

      # Build summarized experiment from eSet
      SE <- SummarizedExperiment::SummarizedExperiment(
        ## TODO:: Do we want to pass an environment for better memory efficiency?
        assays = S4Vectors::SimpleList(as.list(Biobase::assayData(eSet))),
        # Switch rearrange columns so that IDs are first, probes second
        rowData = S4Vectors::DataFrame(Biobase::fData(eSet),
          rownames = rownames(Biobase::fData(eSet))
        ),
        colData = S4Vectors::DataFrame(Biobase::pData(eSet),
          rownames = rownames(Biobase::pData(eSet))
        ),
        metadata = list(
          "experimentData" = eSet@experimentData,
          "annotation" = Biobase::annotation(eSet),
          "protocolData" = Biobase::protocolData(eSet)
        )
      )
      ## TODO:: Determine if this can be done in the SE constructor?
      # Extract names from expression set
      SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
      mDataType <- Biobase::annotation(eSet)
      eSets[[mDataType]] <- SE
    }
  )
  # setNames(pSet@molecularProfiles, names(eSets))
  return(SEfinal)
}

z <- .converteSetToSE(z)

sensitivity.info <- as.data.frame(sensitivity.info)

standardizeRawDataConcRange <- function(sens.info, sens.raw) {
  unq.drugs <- unique(sens.info$drugid)

  conc.m <- data.table(melt(sens.raw[, , 1], as.is = TRUE))
  conc.m[, drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[, .(l = min(value, na.rm = T), r = max(value, na.rm = T)), c("drugid", "Var1")]
  conc.ranges[, Var1 := NULL]
  conc.ranges <- conc.ranges[, unique(.SD), drugid]
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[
    ,
    {
      sq <- sort(unique(c(l, r)))
      l <- sq[seq(1, length(sq) - 1)]
      r <- sq[seq(2, length(sq))]
      .(l = l, r = r)
    },
    drugid
  ]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N > 0)
    unlist(sapply(seq(1, N), function(ii) {
      return(sapply(seq(ii, N), function(jj) {
        return(seq(ii, jj))
      }))
    }), recursive = FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl) {
    if (length(indicies) == 1) {
      return(TRUE)
    }
    sq <- seq(indicies[1], indicies[length(indicies)] - 1)
    all(lr.tbl[["l"]][sq + 1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[, .N, drugid][, N], returnConsInts)

  names(per.drug.range.indicies) <- conc.ranges.disj[, unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]


  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug) {
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug) {
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) {
      return(c(lr.tbl[x[1], l], lr.tbl[x[length(x)], r]))
    }))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify = FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)

  conc.ranges <- conc.m[, .(l = min(value, na.rm = T), r = max(value, na.rm = T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)
  # tic()
  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug) {
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid == drug, .(l, r)], 1, function(rng) {
      conc.m[drugid == drug][conc.ranges[drugid == drug][l <= rng["l"]][r >= rng["r"], Var1], on = "Var1"][value >= rng["l"]][value <= rng["r"], .N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid == drug][which(num.points.in.range == max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })
  # toc()
  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug) {
    rng <- unlist(chosen.drug.ranges[[drug]][, .(l, r)])
    exp.out.range <- conc.ranges[drugid == drug][l > rng["l"] | r < rng["r"], Var1]
    return(exp.out.range)
  }))

  sens.raw[removed.experiments, , ] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]

  for (drug in unq.drugs) {
    rng <- unlist(chosen.drug.ranges[[drug]][, .(l, r)])
    myx <- conc.ranges.kept[drugid == drug, Var1]
    doses <- sens.raw[myx, , "Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, , "Dose"][which(which.remove, arr.ind = TRUE)] <- NA_real_
    sens.raw[myx, , "Viability"][which(which.remove, arr.ind = TRUE)] <- NA_real_

    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid == drug, "chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid == drug, "chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments, "rm.by.conc.range"] <- TRUE

  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}

# add cellosaurus disease type to cell-info

disease <- cell_all$Cellosaurus.Disease.Type[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Cellosaurus.Disease.Type <- disease

# add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Cellosaurus.Accession.id <- assession

# add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$PharmacoDB.id <- pdb

# add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$unique.tissueid.fromstudies <- study_tissue

# add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$CellLine.Type <- cell_type

# add metastatic info to cell_info
metastatic <- cell_all$Metastatic[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Metastatic <- metastatic

# standardize <- standardizeRawDataConcRange(sens.info = sensitivity.info, sens.raw = raw.sensitivity)


gCSI_2018 <- PharmacoSet(
  molecularProfiles = z,
  name = "gCSI",
  cell = cellInfo,
  drug = druginfo,
  sensitivityInfo = sensitivity.info,
  sensitivityRaw = raw.sensitivity,
  sensitivityProfiles = sensitivityProfiles_2018,
  curationCell = curationCell,
  curationDrug = curationDrug,
  curationTissue = curationTissue,
  datasetType = "sensitivity"
)


# filter noisy curves from PSet (modified function to take into account standardized conc range)

filterNoisyCurves2 <- function(pSet, epsilon = 25, positive.cutoff.percent = .80, mean.viablity = 200, nthread = 1) {
  acceptable <- mclapply(rownames(PharmacoGx::sensitivityInfo(pSet)), function(xp) {
    # for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp, , ], 2, as.numeric), stringsAsFactors = FALSE)
    if (!all(is.na(drug.responses))) {
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[, "delta"] <- .computeDelta(drug.responses$Viability)

      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)

      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)

      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
        (delta.sum < epsilon) &
        (max.cum.sum < (2 * epsilon)) &
        (mean(drug.responses$Viability) < mean.viablity)) {
        return(xp)
      }
    }
  }, mc.cores = nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy" = noisy, "ok" = acceptable))
}

.computeDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if (trunc) {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx) - 1]), 0))
  } else {
    return(c(xx[2:length(xx)] - xx[1:length(xx) - 1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if (trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2, simplify = TRUE))
  tt <- tt[which(((tt[, 2] - tt[, 1]) >= 2) == TRUE), ]
  cum.sum <- unlist(lapply(1:nrow(tt), function(x) {
    xx[tt[x, 2]] - xx[tt[x, 1]]
  }))
  return(max(cum.sum))
}

noisy_out <- filterNoisyCurves2(gCSI_2018)
print("filter done")
gCSI_2018@sensitivity$profiles[noisy_out$noisy, ] <- NA

saveRDS(gCSI_2018, file = paste0(output_dir, filename))

# unlink(file.path(processed_dir, tool), recursive = TRUE)

# dataset <- "gCSI2"
# #output ORCESTRA_ID and Pachyderm commit id
# write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
# write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
# pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
# write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
