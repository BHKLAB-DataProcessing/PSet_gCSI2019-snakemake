library(Biobase)
library(readr)
library(tximport)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")
tools <- args[[2]]
transcriptome <- args[[3]]

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/bhklab_orcestra/snakemake/PSet_gCSI2019/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/bhklab_orcestra/snakemake/PSet_gCSI2019/processed"
# tools <- "Kallisto-0.46.1"
# transcriptome <- "Gencode_v33"

tools <- gsub("-", "_", tools)
untar(file.path(download_dir, paste0(tools, ".tar.gz")), exdir = download_dir)
tool_path <- expand.grid(a = tools, b = transcriptome)
tool_path <- paste0(tool_path$a, "_", tool_path$b)

# untar(tarfile = file.path(download_dir, paste0(tools, ".tar.gz")), exdir = processed_dir)
curationCell <- readRDS(file.path(processed_dir, "curationCell.rds"))
rnaseq_results <- list()

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

summarizeRnaSeq <- function(dir,
                            features_annotation,
                            samples_annotation,
                            method) {
  load(features_annotation)

  tx2gene <- as.data.frame(cbind("transcript" = tx2gene$transcripts, "gene" = tx2gene$genes))

  files <- list.files(dir, recursive = TRUE, full.names = T)
  if (method == "kallisto") {
    resFiles <- grep("abundance.h5", files)
  } else {
    resFiles <- grep("quant.sf", files)
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))

  if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
    txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else {
    txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)
  }

  head(txi$counts[, 1:5])
  dim(txi$counts)

  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp), ]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp), ]
  annotation(gene.exp) <- "rnaseq"

  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count), ]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count), ]
  annotation(gene.count) <- "rnaseq"

  txii <- tximport(resFiles, type = method, txOut = T)

  if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
    # remove non-coding transcripts in ensembl
    rownames(txii$abundance) <- gsub("\\..*", "", rownames(txii$abundance))
    txii$abundance[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
    missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
    txii$abundance <- txii$abundance[-which(rownames(txii$abundance) %in% missing_transcript), ]
  }

  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 0.001))
  if (features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")) {
    featureNames(transcript.exp) <- gsub("\\|.*", "", featureNames(transcript.exp))
    fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
  } else {
    fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp), ]
  annotation(transcript.exp) <- "isoforms"


  if (features_annotation == file.path(download_dir, "Ensembl.v99.annotation.RData")) {
    # remove non-coding transcripts in ensembl
    rownames(txii$counts) <- gsub("\\..*", "", rownames(txii$counts))
    txii$counts <- txii$counts[-which(rownames(txii$counts) %in% missing_transcript), ]
  }
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 1))
  if (features_annotation == file.path(download_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(download_dir, "Gencode.v33lift37.annotation.RData")) {
    featureNames(transcript.count) <- gsub("\\|.*", "", featureNames(transcript.count))
    fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
  } else {
    fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count), ]
  annotation(transcript.count) <- "isoforms"

  return(list(
    "rnaseq" = gene.exp,
    "rnaseq.counts" = gene.count,
    "isoforms" = transcript.exp,
    "isoforms.counts" = transcript.count
  ))
}

rnaseq.sampleinfo <- read.csv(file = file.path(download_dir, "gCSI_rnaseq_meta.csv"), stringsAsFactors = FALSE, row.names = 1)
rnaseq.sampleinfo[, "cellid"] <- matchToIDTable(ids = rnaseq.sampleinfo[, "Cell_line"], tbl = curationCell, column = "GNE.cellid", returnColumn = "unique.cellid")


for (r in 1:length(tool_path)) {
  print(tool_path[r])
  if (length(grep(pattern = "Kallisto", x = tool_path[r])) > 0) {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool <- "kallisto"
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool <- "salmon"
  }


  if (length(grep(pattern = "lift37", x = tool_path[r])) > 0) {
    annot <- file.path(download_dir, "Gencode.v33lift37.annotation.RData")
  } else if (length(grep(pattern = "v33", x = tool_path[r])) > 0) {
    annot <- file.path(download_dir, "Gencode.v33.annotation.RData")
  } else {
    annot <- file.path(download_dir, "Ensembl.v99.annotation.RData")
  }
  print(annot)

  rnaseq <- summarizeRnaSeq(
    dir = file.path(download_dir, tool, tool_path[r]),
    features_annotation = annot,
    samples_annotation = rnaseq.sampleinfo,
    method = rnatool
  )


  reps <- matchToIDTable(rnaseq$rnaseq$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$rnaseq$cellid <- reps

  reps <- matchToIDTable(rnaseq$rnaseq.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$rnaseq.counts$cellid <- reps

  reps <- matchToIDTable(rnaseq$isoforms$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$isoforms$cellid <- reps

  reps <- matchToIDTable(rnaseq$isoforms.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$isoforms.counts$cellid <- reps

  rnaseq_results <- c(rnaseq_results, c(
    rnaseq <- setNames(rnaseq, paste0(tool, ".", names(rnaseq)))
  ))
}

unlink(file.path(download_dir, tools), recursive = TRUE)

saveRDS(rnaseq_results, file = file.path(processed_dir, "rnaseq_results.rds"))
