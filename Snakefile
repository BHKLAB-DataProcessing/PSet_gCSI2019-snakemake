from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]
rna_tool = config["rna_tool"]
rna_ref = config["rna_ref"]

basePath = "https://orcestradata.blob.core.windows.net/gcsi/gCSI/2018"

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_gCSI2019:
    input:
        S3.remote(prefix + "download/gCSI_molData.RData"),
        S3.remote(prefix + "processed/profiles.RData"),
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz'),
        S3.remote(prefix + "processed/sens.data.RData"),
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + 'download/' + rna_ref_file),
        S3.remote(prefix + "download/gCSI_rnaseq_meta.csv")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb = 10000,
        disk_mb = 8000
    shell:
        """
        Rscript scripts/getgCSI2019.R {prefix} {filename} {rna_tool} {rna_ref}
        """

rule recalculate_and_assemble:
    input:
        S3.remote(prefix + "processed/raw_sense_slices.zip")
    output:
        S3.remote(prefix + "processed/profiles.RData")
    shell:
        """
        Rscript scripts/recalculateAndAssemble.R {prefix}
        """

rule process_sens_data:
    input:
        S3.remote(prefix + "download/gCSI_GRmetrics_v1.3.tsv"),
        S3.remote(prefix + "download/gCSI_GRvalues_v1.3.tsv")
    output:
        S3.remote(prefix + "processed/raw_sense_slices.zip"),
        S3.remote(prefix + "processed/sens.data.RData")
    shell:
        """
        Rscript scripts/processSensData.R {prefix}
        """

rule download_annotation:
    output:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + 'download/' + rna_ref_file),
        S3.remote(prefix + "download/gCSI_rnaseq_meta.csv")
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/{rna_ref_file}' \
            -O {prefix}download/{rna_ref_file}
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/gCSI_rnaseq_meta.csv' \
            -O {prefix}download/gCSI_rnaseq_meta.csv
        """

rule download_sens_data:
    output:
        S3.remote(prefix + "download/gCSI_GRmetrics_v1.3.tsv"),
        S3.remote(prefix + "download/gCSI_GRvalues_v1.3.tsv")
    shell:
        """
        Rscript scripts/downloadSensData.R {prefix}
        """

rule download_data:
    output:
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz'),
        S3.remote(prefix + "download/gCSI_molData.RData")
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        wget '{basePath}/Mutation/gCSI_molData.RData' -O {prefix}download/gCSI_molData.RData
        """
