library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)
library(Biobase)
library(reshape2)
library(abind)
library(data.table)
library(parallel)

verbose=FALSE
nthread=1
            
options(stringsAsFactors = FALSE)
z <- list()

dir.prefix <- "pfs"

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(trimws(x)),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
        if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}

# load(file.path(dir.prefix, "downloadrnagCSI/gCSI_2017_molecprofile.RData"))
load(file.path(dir.prefix, "gcsi2018ProfilesAssemble/profiles.RData"))
# load(file.path(dir.prefix, "getgCSI2017/gcsidrugpost.RData"))
load(file.path(dir.prefix, "gcsi2018SensData/sens.data.RData"))
load(file.path(dir.prefix, "gcsi2017RawSensitivity/gCSI_molData.RData"))
load(file.path(dir.prefix, "downloadrnagCSI/gCSI_2017_molecprofile.RData"))


cell_all <- read.csv(file.path(dir.prefix, "downAnnotations/cell_annotation_all.csv"), na.strings=c("", " ", "NA"))
drug_all <- read.csv(file.path(dir.prefix, "downAnnotations/drugs_with_ids.csv"), na.strings=c("", " ", "NA"))


curationCell <- cell_all[apply(!is.na(cell_all[,c("gCSI.cellid", "GNE.cellid")]),1,any),]
curationTissue <- cell_all[apply(!is.na(cell_all[,c("gCSI.cellid", "GNE.cellid")]),1,any),]
curationCell <- curationCell[ , c("unique.cellid", "gCSI.cellid", "GNE.cellid")]
curationTissue <- curationTissue[ , c("unique.tissueid", "gCSI.tissueid", "GNE.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

stopifnot(setequal(rownames(sensitivity.info), rownames(sensitivity.published)))
stopifnot(setequal(rownames(raw.sensitivity), rownames(sensitivity.published)))
stopifnot(setequal(rownames(res), rownames(sensitivity.published)))


recomputed_2018 <- res

sensitivityProfiles_2018 <- data.frame("aac_recomputed" = NA, "ic50_recomputed"=NA, "HS"=NA, "E_inf"=NA, "EC50"=NA)
sensitivityProfiles_2018[nrow(sensitivity.info),] <- NA
sensitivityProfiles_2018[,"aac_recomputed"] <- as.numeric(recomputed_2018[,"aac_recomputed"])
sensitivityProfiles_2018[,"ic50_recomputed"] <- as.numeric((recomputed_2018[,"ic50_recomputed"]))
sensitivityProfiles_2018[,"HS"] <- as.numeric((recomputed_2018[,"HS"]))
sensitivityProfiles_2018[,"E_inf"] <- as.numeric((recomputed_2018[,"E_inf"]))
sensitivityProfiles_2018[,"EC50"] <- as.numeric((recomputed_2018[,"EC50"]))

rownames(sensitivityProfiles_2018) <- rownames(recomputed_2018)

sensitivityProfiles_2018 <- sensitivityProfiles_2018[rownames(sensitivity.info),]
sensitivity.published <- sensitivity.published[rownames(sensitivity.info),]
raw.sensitivity <- raw.sensitivity[rownames(sensitivity.info),,]

sensitivityProfiles_2018$GR_AOC_published <- sensitivity.published[,c("GR_AOC")]
sensitivityProfiles_2018$meanviability_published <- sensitivity.published[,c("meanviability")]
sensitivityProfiles_2018$GRmax_published <- sensitivity.published[,c("GRmax")]
sensitivityProfiles_2018$Emax_published <- sensitivity.published[,c("Emax")]
sensitivityProfiles_2018$GRinf_published <- sensitivity.published[,c("GRinf")]
sensitivityProfiles_2018$GEC50_published <- sensitivity.published[,c("GEC50")]
sensitivityProfiles_2018$GR50_published <- sensitivity.published[,c("GR50")]


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
mapInfo <- data.frame(gCSI.cellid = unique(sensitivity.info$cellid), 
                      unique.cellid = matchToIDTable(unique(sensitivity.info$cellid), curationCell, "gCSI.cellid", "unique.cellid"))
stopifnot(!anyNA(mapInfo[,2]))
sensitivity.info$cellid <- mapInfo[match(sensitivity.info$cellid, mapInfo[,1]),2]

message("mapping sensitivity drugs")

mapInfo <- data.frame(gCSI.drugid = unique(sensitivity.info$drugid), 
                      unique.drugid = matchToIDTable(unique(sensitivity.info$drugid), curationDrug, "gCSI.drugid", "unique.drugid"))
stopifnot(!anyNA(mapInfo[,2]))
sensitivity.info$drugid <- mapInfo[match(sensitivity.info$drugid, mapInfo[,1]),2]


rnaseq <- gCSI@molecularProfiles$rnaseq

rnaseq$cellid <- rnaseq$Cell_line

reps <- matchToIDTable(rnaseq$cellid, curationCell, "GNE.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rnaseq$cellid <- reps
# rnaseq$tissueid <- curationTissue[rnaseq$cellid, "unique.tissueid"]



rnaseq.counts <- gCSI@molecularProfiles$rnaseq.counts

rnaseq.counts$cellid <- rnaseq.counts$Cell_line

reps <- matchToIDTable(rnaseq.counts$cellid, curationCell, "GNE.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rnaseq.counts$cellid <- reps


isoforms <- gCSI@molecularProfiles$isoforms

isoforms$cellid <- isoforms$Cell_line

reps <- matchToIDTable(isoforms$cellid, curationCell, "GNE.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
isoforms$cellid <- reps
# rnaseq$tissueid <- curationTissue[rnaseq$cellid, "unique.tissueid"]


isoforms.counts <- gCSI@molecularProfiles$isoforms.counts

isoforms.counts$cellid <- isoforms.counts$Cell_line

reps <- matchToIDTable(isoforms.counts$cellid, curationCell, "GNE.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
isoforms.counts$cellid <- reps



reps <- matchToIDTable(rownames(cellInfo), curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rownames(cellInfo) <- reps

cellInfo <- as.data.frame(cellInfo)

cellinall <- unionList(rnaseq$cellid, sensitivity.info$cellid, cnv$cellid, mut$cellid)

newCells <- setdiff(cellinall, rownames(cellInfo))

newrows <- cellInfo[newCells,]
rownames(newrows) <- newCells
cellInfo <- rbind(cellInfo, newrows)

cellInfo$tissueid <- curationTissue[rownames(cellInfo), "unique.tissueid"]

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

druginfo <- curationDrug[drugsPresent,]

drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "gCSI.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(druginfo),]
druginfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]

z <- list()

z <- c(z,c(
  "rnaseq"=rnaseq,
  "rnaseq.counts" = rnaseq.counts,
  "isoforms" = isoforms,
  "isoforms.counts" = isoforms.counts,
  "cnv"=cnv,
  "mutation" = mut)
)

sensitivity.info <- as.data.frame(sensitivity.info)

standardizeRawDataConcRange <- function(sens.info, sens.raw){
	unq.drugs <- unique(sens.info$drugid)

	conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
	conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
	conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
	conc.ranges[,Var1 := NULL]
	conc.ranges <- conc.ranges[,unique(.SD), drugid]	
	# conc.ranges[,N := .N, drugid]
	conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
				                       l = sq[seq(1,length(sq)-1)];
				                       r = sq[seq(2,length(sq))];
				                       .(l=l,r=r)}, drugid]
    ## Function below returns all consecutive ranges of ints between 1 and N
    returnConsInts <- function(N) {
        stopifnot(N>0)
        unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
    }
    rangeNoHoles <- function(indicies, lr.tbl){
        if(length(indicies) == 1) return(TRUE)
        sq <- seq(indicies[1], indicies[length(indicies)]-1)
        all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
    }
    per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][,N], returnConsInts)

    names(per.drug.range.indicies) <- conc.ranges.disj[,unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
    

    # Check if there are any holes in the chosen range combination
    per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug){

        lr.tbl <- conc.ranges.disj[drugid == drug]
        per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]

        })
    per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug){

        lr.tbl <- conc.ranges.disj[drugid == drug]
        res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
        colnames(res) <- c("l", "r")
        res <- data.frame(res)
        res <- cbind(drugid = drug, res)
        }, simplify=FALSE)
    per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
    
    conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
    setkey(conc.m, Var1)
    conc.m <- na.omit(conc.m)
    setkey(conc.m, drugid, Var1, value)
    setkey(conc.ranges, drugid, l, r)
    # tic()
    ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
    ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
    chosen.drug.ranges <- lapply(unq.drugs, function(drug){
        num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
            conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
            # conc.m[drugid==drug][, Var1]
            })
        max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
        max.ranges[which.max(log10(r) - log10(l)), ]
    })
    # toc()
    names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
    removed.experiments <- unlist(lapply(unq.drugs, function(drug){
        rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
        exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
        return(exp.out.range)
        }))

    sens.raw[removed.experiments,,] <- NA_real_
    conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]

    for(drug in unq.drugs){
        rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
        myx <- conc.ranges.kept[drugid==drug,Var1]
        doses <- sens.raw[myx, ,"Dose"]
        which.remove <- (doses < rng["l"] | doses > rng["r"])
        sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
        sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_

        ## Annotate sens info with chosen range
        sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
        sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
    }
    sens.info$rm.by.conc.range <- FALSE
    sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE

    return(list("sens.info" = sens.info, sens.raw = sens.raw))
}
		 
		 
		 
standardize <- standardizeRawDataConcRange(sens.info = sensitivity.info, sens.raw = raw.sensitivity)
		 

gCSI_2018 <- PharmacoSet(molecularProfiles=z,
                         name="gCSI",
                         cell=cellInfo,
                         drug=druginfo,
                         sensitivityInfo=standardize$sens.info,
                         sensitivityRaw=standardize$sens.raw,
                         sensitivityProfiles=sensitivityProfiles_2018,
                         curationCell=curationCell,
                         curationDrug=curationDrug,
                         curationTissue=curationTissue,
                         datasetType="sensitivity")

                 
 #filter noisy curves from PSet (modified function to take into account standardized conc range)
		 

filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
acceptable <- mclapply(rownames(PharmacoGx::sensitivityInfo(pSet)), function(xp) {
  #for(xp in rownames(sensitivityInfo(pSet))){
  drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
  if (!all(is.na(drug.responses))){
    
  
  drug.responses <- drug.responses[complete.cases(drug.responses), ]
  doses.no <- nrow(drug.responses)
  drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
  
  delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
  
  max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
  
  if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
      (delta.sum < epsilon) &
      (max.cum.sum < (2 * epsilon)) &
      (mean(drug.responses$Viability) < mean.viablity)) {
    return (xp)
  }
  }
  
}, mc.cores=nthread)
acceptable <- unlist(acceptable)
noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}
		 
noisy_out <- filterNoisyCurves2(gCSI_2018)
print("filter done")
gCSI_2018@sensitivity$profiles[noisy_out$noisy, ] <- NA                
                          
saveRDS(gCSI_2018, file="/pfs/out/gCSI_2018.rds")

# #rna-seq (processed - Zhaleh)

# load("/Users/anthonymammoliti/Desktop/gCSI_kallisto.RData")

# gCSI@molecularProfiles$rnaseq@phenoData@data$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$rnaseq$cellid[1] <- "888-mel"

# gCSI@molecularProfiles$rnaseq$cellid[which(!gCSI@molecularProfiles$rnaseq$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")
# gCSI@molecularProfiles$rnaseq@phenoData@data$cellid[which(!gCSI@molecularProfiles$rnaseq@phenoData@data$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")

# gCSI@molecularProfiles$rnaseq$cellid <- cellosaurus$cellosaurus.name[match(gCSI@molecularProfiles$rnaseq$cellid, cellosaurus$lab_annotation)]

# gCSI_cellname <- gCSI@molecularProfiles$rnaseq@phenoData@data$Cell_line[which(!gCSI@molecularProfiles$rnaseq$cellid %in% curationCell$unique.cellid)]

# curationCell[570:788,] <- NA

# curationCell$gCSI.cellid[570:788] <- gCSI_cellname
# curationCell$unique.cellid[570:788] <- gCSI@molecularProfiles$rnaseq$cellid[which(!gCSI@molecularProfiles$rnaseq$cellid %in% curationCell$unique.cellid)]

# #rna-seq isoforms (processed - Zhaleh)
# gCSI@molecularProfiles$isoforms@phenoData@data$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$isoforms$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$isoforms$cellid[which(!gCSI@molecularProfiles$isoforms$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")
# gCSI@molecularProfiles$isoforms@phenoData@data$cellid[which(!gCSI@molecularProfiles$isoforms@phenoData@data$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")

# gCSI@molecularProfiles$isoforms$cellid <- cellosaurus$cellosaurus.name[match(gCSI@molecularProfiles$isoforms$cellid, cellosaurus$lab_annotation)]

# #rna-seq counts (processed - Zhaleh)
# gCSI@molecularProfiles$rnaseq.counts@phenoData@data$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$rnaseq.counts$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$rnaseq.counts$cellid[which(!gCSI@molecularProfiles$rnaseq.countscellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")
# gCSI@molecularProfiles$rnaseq.counts@phenoData@data$cellid[which(!gCSI@molecularProfiles$rnaseq.counts@phenoData@data$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")

# gCSI@molecularProfiles$rnaseq.counts$cellid <- cellosaurus$cellosaurus.name[match(gCSI@molecularProfiles$rnaseq.counts$cellid, cellosaurus$lab_annotation)]


# #rna-seq isoform counts (processed - Zhaleh)

# gCSI@molecularProfiles$isoforms.counts@phenoData@data$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$isoforms.counts$cellid[1] <- "888-mel"
# gCSI@molecularProfiles$isoforms.counts$cellid[which(!gCSI@molecularProfiles$isoforms.counts$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")
# gCSI@molecularProfiles$isoforms.counts@phenoData@data$cellid[which(!gCSI@molecularProfiles$isoforms.counts@phenoData@data$cellid %in% cellosaurus$lab_annotation)] <- c("DOV13","NCI-H322","Hs 695T","NCI-H322","OVCAR-3","OVCA433","Rh30","SR")

# gCSI@molecularProfiles$isoforms.counts$cellid <- cellosaurus$cellosaurus.name[match(gCSI@molecularProfiles$isoforms.counts$cellid, cellosaurus$lab_annotation)]

# #cnv

# gCSI@molecularProfiles$cnv$cellid[231] <- "SNU-16"
# gCSI@molecularProfiles$cnv$cellid[282] <- "U266B1"
# gCSI@molecularProfiles$cnv$cellid <- cellosaurus$cellosaurus.name[match(gCSI@molecularProfiles$cnv$cellid, cellosaurus$lab_annotation)]



# curationCell[789,] <- NA
# curationCell$gCSI.cellid[789] <- "HEC-1"
# curationCell$unique.cellid[789] <- "HEC-1"

# curationCell <- curationCell[c(-675,-727),] 

# rownames(curationCell) <- curationCell$unique.cellid

# curationDrug <- data.frame("unique.drugid"=unique(sensitivityInfo_2018$drugid), "gCSI.drugid"=unique(sensitivityInfo_2018$DrugName))
# rownames(curationDrug) <- curationDrug$unique.drugid

# curationTissue <- gCSI_2018@curation$tissue

# curationTissue_missing <- curationCell$gCSI.cellid[which(!curationCell$gCSI.cellid %in% rownames(curationTissue))]
# curationTissue <- curationTissue[which(rownames(curationTissue) %in% curationCell$gCSI.cellid),]
# curationTissue_uniqueid <- curationCell$unique.cellid[match(rownames(curationTissue),curationCell$gCSI.cellid)]
# rownames(curationTissue) <- curationTissue_uniqueid

# cell_all <- read.csv("~/Desktop/getgcsi_2018-master/annotation/cell_annotation_all.csv", sep=",", comment.char="#", na.strings=c("", " ", "NA"))


# matched_missing <- matchToIDTable(curationTissue_missing, cell_all, "gCSI.cellid")
# matched_missing[matched_missing=="character(0)"] <- NA
# still_missing <- which(is.na(matched_missing))
# still_missing <- c('HCC70','A3-KAW','A4-Fuk','COLO-849','DMS-53','DMS-79','HCC-1359','MCF10A','MT3','OAW-28','SW1417')
# matched_missing[which(is.na(matched_missing))] <- still_missing
# matched_missing <- as.character(matched_missing)
# curationTissue[673:787,] <- NA
# curationTissue$gCSI.tissueid[673:787] <- cell_all$unique.tissueid[which(matched_missing %in% cell_all$unique.cellid)]
# curationTissue$unique.tissueid[673:787] <- cell_all$unique.tissueid[which(matched_missing %in% cell_all$unique.cellid)]
# rownames(curationTissue)[673:787] <- curationCell$unique.cellid[match(curationTissue_missing,curationCell$gCSI.cellid)]

# cell <- gCSI_2017@cell
# cell[411:787,] <- NA
# cell$unique.id[411:787] <- curationCell$unique.cellid[which(!curationCell$gCSI.cellid %in% rownames(cell))]
# cell$CellLineName[411:787] <- curationCell$gCSI.cellid[which(!curationCell$unique.cellid %in% rownames(cell))]
# rownames(cell)[411:787] <- cell$unique.id[411:787]

# ####PUBLISHED METRICS#####

# gCSI_GR_AOC_Pub <- read.csv("/Users/anthonymammoliti/Desktop/getgcsi_2018-master/2018/sensitivity/gCSI_GRmetrics_v1.2.tsv", sep="\t")


# #match to cellosaurus

# cellosaurus <- read.xlsx("~/Desktop/annotations/cellosaurus_names.xlsx", na.strings=c("", " ", "NA"))


# gCSI_GR_AOC_Pub$CellLineName[grep("1-10", gCSI_GR_AOC_Pub$CellLineName)] <- "NB-TU-1-10"

# cello_temp <- cellosaurus
# cello_temp$lab_annotation <-  tolower(gsub(badchars, "",cello_temp$lab_annotation))

# gCSI_GR_AOC_Temp <-  gCSI_GR_AOC_Pub
# gCSI_GR_AOC_Temp$CellLineName <- tolower(gsub(badchars, "",gCSI_GR_AOC_Temp$CellLineName))

# cello_matched <- cellosaurus$cellosaurus.name[match(gCSI_GR_AOC_Temp$CellLineName, cello_temp$lab_annotation)]
# gCSI_GR_AOC_Pub$CellLineName <- cello_matched


# #match to pubchem

# pubchem <- read.xlsx("~/Desktop/annotations/drugsWithids_pub.xlsx", na.strings=c("", " ", "NA"))
# pubchem$Selected.name[which(is.na(pubchem$Selected.name))] <- pubchem$unique.drugid[which(is.na(pubchem$Selected.name))]

# gCSI_GR_AOC_Pub$DrugName[grep("Vinblastine", gCSI_GR_AOC_Pub$DrugName)] <- "Vincaleukoblastine"
# gCSI_GR_AOC_Pub$DrugName[grep("TGX-221", gCSI_GR_AOC_Pub$DrugName)] <- "TGX221"
# gCSI_GR_AOC_Pub$DrugName[grep("GDC-0941", gCSI_GR_AOC_Pub$DrugName)] <- "Pictilisib"
# gCSI_GR_AOC_Pub$DrugName[grep("CAL-101", gCSI_GR_AOC_Pub$DrugName)] <- "Idelalisib"
# gCSI_GR_AOC_Pub$DrugName[grep("AZD-7762", gCSI_GR_AOC_Pub$DrugName)] <- "AZD7762"
# gCSI_GR_AOC_Pub$DrugName[grep("PLX-4720", gCSI_GR_AOC_Pub$DrugName)] <- "PLX4720"
# gCSI_GR_AOC_Pub$DrugName[grep("5-FU", gCSI_GR_AOC_Pub$DrugName)] <- "5-Fluorouracil"
# gCSI_GR_AOC_Pub$DrugName[grep("AZD-8055", gCSI_GR_AOC_Pub$DrugName)] <- "AZD8055"
# gCSI_GR_AOC_Pub$DrugName[grep("JQ1", gCSI_GR_AOC_Pub$DrugName)] <- "JQ1 compound"
# gCSI_GR_AOC_Pub$DrugName[grep("CHIR-99021", gCSI_GR_AOC_Pub$DrugName)] <- "Chir-99021"
# gCSI_GR_AOC_Pub$DrugName[grep("NU 7441", gCSI_GR_AOC_Pub$DrugName)] <- "NU-7441"
# gCSI_GR_AOC_Pub$DrugName[grep("17-AAG", gCSI_GR_AOC_Pub$DrugName)] <- "Tanespimycin"
# gCSI_GR_AOC_Pub$DrugName[grep("MLN_2480", gCSI_GR_AOC_Pub$DrugName)] <- "MLN2480"

# pubchem_matched <- pubchem$Selected.name[match(gCSI_GR_AOC_Pub$DrugName, pubchem$Selected.name)]
# gCSI_GR_AOC_Pub$DrugName <- pubchem_matched

# uids <- unique(sprintf("%s__%s__%s",gCSI_GR_AOC_Pub$CellLineName,
#                        gCSI_GR_AOC_Pub$DrugName,
#                        gCSI_GR_AOC_Pub$ExperimentNumber))

# rownames(gCSI_GR_AOC_Pub) <- uids

# sensitivityProfiles_2018[,"GR_AOC_published"] <- as.numeric(gCSI_GR_AOC_Pub$GR_AOC[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])
# sensitivityProfiles_2018[,"meanviability_published"] <- as.numeric(gCSI_GR_AOC_Pub$meanviability[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])
# sensitivityProfiles_2018[,"GRmax_published"] <- as.numeric(gCSI_GR_AOC_Pub$GRmax[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])
# sensitivityProfiles_2018[,"GRinf_published"] <- as.numeric(gCSI_GR_AOC_Pub$GRinf[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])
# sensitivityProfiles_2018[,"GEC50_published"] <- as.numeric(gCSI_GR_AOC_Pub$GEC50[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])
# sensitivityProfiles_2018[,"GR50_published"] <- as.numeric(gCSI_GR_AOC_Pub$GR50[match(rownames(sensitivityProfiles_2018), rownames(gCSI_GR_AOC_Pub))])


