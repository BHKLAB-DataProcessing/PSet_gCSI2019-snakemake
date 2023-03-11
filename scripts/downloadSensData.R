library(downloader)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_gCSI2019-snakemake/download"

download('http://research-pub.gene.com/gCSI_GRvalues2019/gCSI_GRdata_v1.3.tsv.tar.gz', destfile = file.path(download_dir, "gCSI_GRdata_v1.3.tsv.tar.gz"))
untar(tarfile=file.path(download_dir, "gCSI_GRdata_v1.3.tsv.tar.gz"), exdir=download_dir)
file.copy(file.path(download_dir, 'OUTPUT', 'gCSI_GRmetrics_v1.3.tsv'), file.path(download_dir, 'gCSI_GRmetrics_v1.3.tsv'))
file.copy(file.path(download_dir, 'OUTPUT', 'gCSI_GRvalues_v1.3.tsv'), file.path(download_dir, 'gCSI_GRvalues_v1.3.tsv'))
unlink(file.path(download_dir, 'OUTPUT'), recursive=TRUE)
file.remove(file.path(download_dir, "gCSI_GRdata_v1.3.tsv.tar.gz"))