#  =============================================================================
#  TCGA-Assembler version 2
#
#  Copyright (C) <2017>  <Yitan Zhu>
#  This file is part of TCGA-Assembler.
#
#  TCGA-Assembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TCGA-Assembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TCGA-Assembler.  If not, see <http://www.gnu.org/licenses/>.
#  =============================================================================


#  =============================================================================
#  TCGA-Assembler Version 2 Module A
#  =============================================================================


#  =============================================================================
#  variable prefix example
#  =============================================================================
#  s : string
#  n : number
#  v : vector
#  d : data.frame
#  m : matrix
#  l : list
#  lv : list of vector
#  ld : list of data.frame


#  =============================================================================
#  internal functions, NOT used directly by user
#  =============================================================================

#' Get fields names of GDC entities
#'
#' @param arch String of archive type: "legacy" or "".
#' @param endp String of endpoint: "files".
#' @return Character vector of all available field names.
#' @examples
#' fieldsList <- FieldsList("legacy")
#' fieldsList <- FieldsList("")
FieldsList <- function(arch = "legacy",
                       endp = "files") {
    library(rjson)
    #url <- paste("https://gdc-api.nci.nih.gov/",
    url <- paste("https://api.gdc.cancer.gov/",
                 arch,
                 ifelse(arch == "", "", "/"),
                 endp,
                 "/_mapping",
                 sep = "")
    opt <- "--silent --show-error"
    arg <- paste(opt, url)
    jsn <- paste(system2("curl", arg, stdout = T), collapse = "")
    return(fromJSON(jsn)$fields)
}

#' Define the default fields used in metadata file
#'
#' @return Character vector of selected fields names for metadata.
FieldsMeta <- function() {
    return(c(# "access",
             "archive.file_name",
             # "cases.case_id",
             "cases.project.project_id",
             "cases.samples.portions.analytes.aliquots.submitter_id",  # 1st
             # "cases.samples.portions.analytes.submitter_id",
             "cases.samples.portions.submitter_id",  # 2nd
             "cases.samples.sample_type",
             "cases.samples.sample_type_id",  #  sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
             "data_category",
             "data_type",
             "experimental_strategy",
             "file_id",
             "file_name",
             "file_size",
             # "md5sum",
             "platform",
             "updated_datetime"))
}

#' Get count of GDC entities in /archive/endpoint
#'
#' @param arch String of archive type: "legacy" or "".
#' @param endp String of endpoint: "files".
#' @return Count of entities in specified archive & endpoint.
#' @examples
#' entityCount <- EntityCount("legacy")
#' entityCount <- EntityCount("")
EntityCount <- function(arch = "legacy",
                        endp = "files") {
    library(rjson)
    #url <- paste("https://gdc-api.nci.nih.gov/",
    url <- paste("https://api.gdc.cancer.gov/",
                 arch,
                 ifelse(arch == "", "", "/"),
                 endp,
                 sep = "")
    opt <- "--silent --show-error"
    arg <- paste(opt, url)
    jsn <- paste(system2("curl", arg, stdout = T), collapse = "")
    # stopifnot(length(fromJSON(jsn)$warnings) == 0)
    return(fromJSON(jsn)$data$pagination$total)
}

#' Get a string of current time (YYYYMMDDhhmmss)
#'
#' @return String of current time (YYYYMMDDhhmmss).
#' @examples
#' TimeNow()
TimeNow <- function() {
    return(gsub("[- :]", "", as.character(Sys.time())))
}

#' Get index of the entity with newest archive file version
#'
#' @param archiveNames Character vector of "archive.file_name". All of
#' these entities have same "file_name".
#' @return Index of the newest one.
#' @examples
#' v <- c("mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.114.1.0.tar.gz",
#'        "mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.2.1.0.tar.gz",
#'        "mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.10.1.0.tar.gz")
#' n <- ArchiveNewest(v)
ArchiveNewest <- function(archiveNames) {
    m <- sapply(strsplit(archiveNames, split = "\\."),
                function(x){rev(x)[5 : 3]})
    mode(m) <- "numeric"  # OR class(m) <- "numeric"
    vOrder <- do.call(order, lapply(seq(nrow(m)), function(x){m[x, ]}))
    return(vOrder[length(vOrder)])
}

#' Get vector of bool indicate the one with newest archive file version or not
#' in one group
#'
#' @param archiveNames Character vector of "archive.file_name". All of
#' these entities have same "file_name".
#' @return Bool vector to indicate the newest one in one group
ArchiveNewestInGroup <- function(archiveNames) {
    b <- seq(length(archiveNames)) == ArchiveNewest(archiveNames)
    return(ifelse(b, T, F))
}

#' Get vector of bool indicate the one with newest archive file version or not
#' in every group
#'
#' @param archiveNames Character vector of "archive.file_name", could be
#' split into groups, each group with a unique "file_name".
#' @param splitFactor Character vector of "file_name" (sorted).
#' @return Bool vector to indicate the newest one in every group.
ArchiveNewestInGroups <- function(archiveNames,
                                  splitFactor) {
    stopifnot(all(splitFactor == splitFactor[order(splitFactor)])) # splitFactor should be sorted before split
    return(unlist(lapply(split(archiveNames, splitFactor),
                         ArchiveNewestInGroup)))
}

#' Choose columns from a file
#'
#' @param fileName String of filename.
#' @param fileId String of fileId with same length of fileName.
#' @param colNames Character vector of colname, specify the columns choosed.
#' @param sortBy Name one column which acts as the rownames.
#' @param skipLines Number of lines skipped in \code{read.csv}.
#' @param naStrings String indicating the \code{NA} in \code{read.csv}.
#' @return A \code{data.frame} with specified columns from one file.
ColumnsFromFile <- function(fileName,
                            fileId,
                            colNames,
                            sortBy = "",
                            skipLines = 0,
                            naStrings) {
    d <- read.csv(fileName,
                  sep = "\t",
                  row.names = NULL,
                  as.is = T,
                  skip = skipLines,
                  na.strings = naStrings)
    if (sortBy != "") {
        if (length(d[, sortBy]) != length(unique(d[, sortBy]))) { # probe duplicated, check same value for each probe group
            stopifnot(all(tapply(Reduce(paste, d[, colNames]), d[, sortBy],
                                 function(x){length(unique(x)) == 1})))
            d <- d[!duplicated(d[, sortBy]), ]
        }
        stopifnot(length(d[, sortBy]) == length(unique(d[, sortBy])))
        rownames(d) <- d[, sortBy]
        d <- d[order(d[, sortBy]), colNames]
    } else {
        d <- cbind("fileId" = rep(fileId, dim(d)[1]), d[, colNames])
    }
    return(d)
}

#' Choose columns from each file of filenames
#'
#' @param fileName String of filename, named with file_id.
#' @param colNames Character vector of colname, specify the columns choosed.
#' @param sortBy Name one column which acts as the rownames.
#' @param skipLines Number of lines skipped in \code{read.csv}.
#' @param naStrings String indicating the \code{NA} in \code{read.csv}.
#' @return List of \code{data.frame} with specified columns from each file.
ColumnsFromFiles <- function(fileNameById,
                             colNames,
                             sortBy = "",
                             skipLines = 0,
                             naStrings = "NA") {
    l <- lapply(names(fileNameById),
                function(x){
                    ColumnsFromFile(fileNameById[x],
                                    x,
                                    colNames,
                                    sortBy,
                                    skipLines,
                                    naStrings)
                })
    names(l) <- names(fileNameById)
    return(l)
}

#' Strip characters at left or right end of each string in a vector
#'
#' @param vUnstripped Character vector of unstripped strings.
#' @param stripNum Number of characters need to be stripped.
#' @param stripEnd Sting indicating the terminal: "right" or "left".
#' @return Character vector of stripped strings.
StripEnd <- function(vUnstripped,
                     stripNum,
                     stripEnd = "right") {
    if (stripNum == 0) {
        vStripped <- vUnstripped
    } else  {
        if (stripEnd %in% c("r", "right")) {
            vStripped <- sapply(strsplit(vUnstripped, split = ""),
                                function(x){
                                    paste(x[-((length(x) - stripNum + 1) : length(x))],
                                          collapse = "")
                                })
        } else if (stripEnd %in% c("l", "left")) {
            vStripped <- sapply(strsplit(vUnstripped, split = ""),
                                function(x){
                                    paste(x[-(1 : stripNum)], collapse = "")
                                })
        } else {
            print("stripEnd should be one of 'right', 'r', 'l' or 'left'")
        }
    }
    return(vStripped)
}

#' Make a vector of values named with probes
#'
#' @param dProbeValue A \code{data.frame}, usually read from filename.
#' @param colProbe String of column name of probe, used as name.
#' @param colValue String of column name of value.
#' @param stripNum Number of characters need to be stripped from "probe".
#' @param stripEnd Sting indicating the terminal: "right" or "left".
#' @return Named vector.
ProbeValue <- function(dProbeValue,
                       colProbe,
                       colValue,
                       stripNum = 0,
                       stripEnd = "right") {
    v <- dProbeValue[, colValue]
    names(v) <- StripEnd(dProbeValue[, colProbe], stripNum, stripEnd)
    return(v)
}

#' Cut the "*\\.1\\.*" columns and rows with \code{NA} only
#'
#' @param metaData A \code{data.frame} read from metadata file.
#' @param colPattern String of regular expression pattern to filter colname.
#' @return A \code{data.frame} after cutting.
MetaCut <- function(metaData,
                    colPattern = "\\.1\\.") {
    colIdx <- grep(colPattern, colnames(metaData))
    rowIdx <- sapply(seq(dim(metaData)[1]),
                     function(x){
                         ifelse(all(is.na(metaData[x, colIdx])), T, F)
                     })
    metaCut <- metaData[rowIdx, seq(dim(metaData)[2])[-colIdx]]
    return(metaCut)
}

#' Define the filters with specified assay platform
#'
#' @param sAssay String of assay platform.
#' @return List of filter.
#' @example
#' filter <- Filter("methylation_450")
Filter <- function(sAssay) {
    filter <- list()
    vAssay <- c(# Copy number segmentation
                "cna_cnv.hg18",
                "cna_cnv.hg19",
                "cna_nocnv.hg18",
                "cna_nocnv.hg19",
                # Exon junction quantification
                "exonJunction_RNAseq",
                "exonJunction_TotalRNAseq",
                # Exon quantification
                "exon_RNAseq",
                "exon_TotalRNAseq",
                # Gene expression quantification
                "gene_Array",
                "gene.normalized_RNAseq",
                "gene_RNAseq",
                "gene.normalized_TotalRNAseq",
                "gene_TotalRNAseq",
                # Isoform expression quantification
                "isoform.normalized_RNAseq",
                "isoform_RNAseq",
                "isoform.normalized_TotalRNAseq",
                "isoform_TotalRNAseq",
                # Methylation beta value
                "methylation_27",
                "methylation_450",
                # miRNA gene quantification
                "mir_GA.hg18",
                "mir_GA.hg19",
                "mir_GA.hg19.mirbase20",
                # miRNA gene quantification
                "mir_HiSeq.hg18",
                "mir_HiSeq.hg19",
                "mir_HiSeq.hg19.mirbase20",
                # miRNA isoform quantification
                "mirIsoform_GA.hg18",
                "mirIsoform_GA.hg19",
                "mirIsoform_GA.hg19.mirbase20",
                # miRNA isoform quantification
                "mirIsoform_HiSeq.hg18",
                "mirIsoform_HiSeq.hg19",
                "mirIsoform_HiSeq.hg19.mirbase20",
                # Protein expression quantification
                "protein_RPPA",
                # Simple somatic mutation
                "somaticMutation_DNAseq",
                # CPTAC
                "glycoproteome_iTRAQ",
                "phosphoproteome_iTRAQ",
                "proteome_iTRAQ")
    if (!sAssay %in% vAssay) {  # assayPlatform = sAssay
        print(paste(c("assayPlatform should be one of:", vAssay), collapse = " "))
    } else if (sAssay == "cna_cnv.hg18") {
        filter$data_category         <- c("Copy number variation")
        filter$data_type             <- c("Copy number segmentation")
        filter$experimental_strategy <- "Genotyping array"
        filter$platform              <- "Affymetrix SNP Array 6.0"
        filter$file_name             <- "\\.hg18\\.seg\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "cna_cnv.hg19") {
        filter$data_category         <- c("Copy number variation")
        filter$data_type             <- c("Copy number segmentation")
        filter$experimental_strategy <- "Genotyping array"
        filter$platform              <- "Affymetrix SNP Array 6.0"
        filter$file_name             <- "\\.hg19\\.seg\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "cna_nocnv.hg18") {
        filter$data_category         <- c("Copy number variation")
        filter$data_type             <- c("Copy number segmentation")
        filter$experimental_strategy <- "Genotyping array"
        filter$platform              <- "Affymetrix SNP Array 6.0"
        filter$file_name             <- "\\.nocnv_hg18\\.seg\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "cna_nocnv.hg19") {
        filter$data_category         <- c("Copy number variation")
        filter$data_type             <- c("Copy number segmentation")
        filter$experimental_strategy <- "Genotyping array"
        filter$platform              <- "Affymetrix SNP Array 6.0"
        filter$file_name             <- "\\.nocnv_hg19\\.seg\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "exon_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Exon quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.bt\\.exon_quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "exon_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Exon quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.bt\\.exon_quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "exonJunction_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Exon junction quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.junction_quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "exonJunction_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Exon junction quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.junction_quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "gene_Array") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Gene expression quantification"
        filter$experimental_strategy <- "Gene expression array"
        filter$platform              <- "AgilentG4502A_07_3"
        filter$file_name             <- "\\.txt_lmean\\.out\\.logratio\\.gene\\.tcga_level3\\.data\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "gene.normalized_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Gene expression quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.genes\\.normalized_results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "gene_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Gene expression quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.genes\\.results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "gene.normalized_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Gene expression quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.genes\\.normalized_results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "gene_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Gene expression quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.genes\\.results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "isoform.normalized_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Isoform expression quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.isoforms\\.normalized_results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "isoform_RNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Isoform expression quantification"
        filter$experimental_strategy <- "RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.isoforms\\.results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "isoform.normalized_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Isoform expression quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.isoforms\\.normalized_results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "isoform_TotalRNAseq") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "Isoform expression quantification"
        filter$experimental_strategy <- "Total RNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "\\.rsem\\.isoforms\\.results"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_GA.hg18") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_GA.hg19") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_GA.hg19.mirbase20") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_HiSeq.hg18") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_HiSeq.hg19") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mir_HiSeq.hg19.mirbase20") {
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA gene quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.mirna\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_GA.hg18") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_GA.hg19") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_GA.hg19.mirbase20") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina GA"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_HiSeq.hg18") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_HiSeq.hg19") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "mirIsoform_HiSeq.hg19.mirbase20") {  # rows different
        filter$data_category         <- c("Gene expression")
        filter$data_type             <- "miRNA isoform quantification"
        filter$experimental_strategy <- "miRNA-Seq"
        filter$platform              <- "Illumina HiSeq"
        filter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.isoform\\.quantification\\.txt"
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "methylation_27") {
        filter$data_category         <- "DNA methylation"
        filter$data_type             <- "Methylation beta value"
        filter$experimental_strategy <- "Methylation array"
        filter$platform              <- "Illumina Human Methylation 27"
        filter$file_name             <- "jhu-usc\\.edu_.*\\.HumanMethylation27\\."
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "methylation_450") {
        filter$data_category         <- "DNA methylation"
        filter$data_type             <- "Methylation beta value"
        filter$experimental_strategy <- "Methylation array"
        filter$platform              <- "Illumina Human Methylation 450"
        filter$file_name             <- "jhu-usc\\.edu_.*\\.HumanMethylation450\\."
        filter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (sAssay == "protein_RPPA") {  # differnt "archive.file_name"
        filter$data_category         <- "Protein expression"
        filter$data_type             <- "Protein expression quantification"
        filter$experimental_strategy <- "Protein expression array"
        filter$platform              <- "MDA_RPPA_Core"
        filter$file_name             <-
            "mdanderson\\.org_.*\\.MDA_RPPA_Core\\.protein_expression\\.Level_3\\."
        filter$submitter_id          <- "cases.0.samples.0.portions.0.submitter_id"
    } else if (sAssay == "somaticMutation_DNAseq") {
        filter$data_category         <- "Simple nucleotide variation"
        filter$data_type             <- "Simple somatic mutation"
        filter$experimental_strategy <- "DNA-Seq"
        filter$platform              <- c("Illumina GA", "Illumina HiSeq", "Mixed platforms")
        filter$file_name             <- "\\.somatic\\.maf"
        filter$submitter_id          <- c("cases.0.samples.0.portions.0.submitter_id",  # both exist
                                          "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id")
    } else if (sAssay == "glycoproteome_iTRAQ") {  # from CTPAC not GDC
        filter$data_category         <- NA
        filter$data_type             <- NA
        filter$experimental_strategy <- NA
        filter$platform              <- "glycoproteome_iTRAQ"
        filter$file_name             <- "_Glycoproteome\\.glycosite\\.itraq\\.tsv"
        filter$submitter_id          <- NA
    } else if (sAssay == "phosphoproteome_iTRAQ") {  # from CTPAC not GDC
        filter$data_category         <- NA
        filter$data_type             <- NA
        filter$experimental_strategy <- NA
        filter$platform              <- "phosphoproteome_iTRAQ"
        filter$file_name             <- "_Phosphoproteome\\.phosphosite\\.itraq\\.tsv"
        filter$submitter_id          <- NA
    } else if (sAssay == "proteome_iTRAQ") {  # from CTPAC not GDC
        filter$data_category         <- NA
        filter$data_type             <- NA
        filter$experimental_strategy <- NA
        filter$platform              <- "proteome_iTRAQ"
        filter$file_name             <- "_Proteome.(itraq|spectral_counts|CDAP\\.r2\\.spectral_counts)\\.tsv"
        filter$submitter_id          <- NA
    }
    return(filter)
}

#' Get metadata of biospecimen & clinical data files with defined filter
#'
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
#' @example
#' metadata <- MetaDataClin(tmpDir = ".")
MetaDataClin <- function(tmpDir = ".",
                         arch = "legacy",
                         fieldsMeta = "",
                         entityCount = (-1),
                         endp = "files") {
    library(rjson)
    if (entityCount == (-1)) {
        entityCount <- EntityCount(arch, endp)
    }
    if (fieldsMeta == "") {
        fieldsMeta <- FieldsMeta()
    } else {
        fieldsList <- FieldsList(arch, endp)
        print("fields names: checking ...")
        stopifnot(all(fieldsMeta %in% fieldsList))
        print("fields names: checking done!")
    }
    print("metadata file: preparing ...")
    out <- paste(tmpDir, "/tmp_metadata_", endp, ".tsv", sep = "")
    #url <- paste('"https://gdc-api.nci.nih.gov/',
    url <- paste("https://api.gdc.cancer.gov/",
                 arch,
                 ifelse(arch == '', '', '/'), endp,
                 sep = '')
    opt <- paste("-o ", out,
                 " --silent --show-error --request POST",
                 " --header Content-Type:application/json --data @",
                 tmpDir, "/tmp_metadata.json",
                 sep = "")
    arg <- paste(opt, url)
    filter2 <- list()
    filter2$access <- list(op = "=",
                           content = list(field = "access", value = "open"))
    filter2$data_format <- list(op = "in",
                                content = list(field = "data_format",
                                               value = "Biotab"))
    filterAll <- list(op = "and",
                      content = list(filter2$access,
                                     filter2$data_format))
    payload <- list(filters = filterAll,
                    format = "TSV",
                    sort = "file_id",
                    from = 0,
                    size = entityCount,
                    fields = paste(fieldsMeta, collapse = ","))
    cat(toJSON(payload), file = paste(tmpDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", arg, stdout = T)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    metaData <- tryCatch(read.csv(out,
                                  sep = "\t",
                                  row.names = NULL,
                                  as.is = T,
                                  na.strings = ""),
                         error = function(e){print(e$message); return(NULL)})
    if (!is.null(metaData)) {  # > 0 lines available in input
        if (nrow(metaData) == 0) {
            metaData <- NULL
        } else {
            write.table(metaData,
                        file = paste(tmpDir,
                                     "/tmp_metadata_clinical__",
                                     arch,
                                     "__",
                                     endp,
                                     ".tsv",
                                     sep = ""),
                        quote = F,
                        sep = "\t",
                        col.names = NA,
                        row.names = T)
        }
    }
    print("metadata file: preparing done!")
    return(metaData)
}

#' Get metadata of somatic mutation data files with defined filter
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
MetaDataSoma <- function(vCancer = "BRCA",
                         sAssay,
                         sampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
                         tmpDir = ".",
                         arch = "legacy",
                         fieldsMeta = "",
                         entityCount = (-1),
                         endp = "files") {
    library(rjson)
    if (entityCount == (-1)) {
        entityCount <- EntityCount(arch, endp)
    }
    if (fieldsMeta == "") {
        fieldsMeta <- FieldsMeta()
    } else {
        fieldsList <- FieldsList(arch, endp)
        print("fields names: checking ...")
        stopifnot(all(fieldsMeta %in% fieldsList))
        print("fields names: checking done!")
    }
    filter <- Filter(sAssay)
    print("metadata file: preparing ...")
    out <- paste(tmpDir, "/tmp_metadata_", endp, ".tsv", sep = "")
    #url <- paste('"https://gdc-api.nci.nih.gov/',
    url <- paste("https://api.gdc.cancer.gov/",
                 arch,
                 ifelse(arch == '', '', '/'), endp,
                 sep = '')
    opt <- paste("-o ", out, " --silent --show-error --request POST ",
                 "--header Content-Type:application/json --data @",
                 tmpDir, "/tmp_metadata.json", sep = "")
    arg <- paste(opt, url)
    filter2 <- list()
    filter2$access <- list(op = "=" ,
                           content = list(field = "access",
                                          value = "open"))
    filter2$data_format <- list(op = "in",
                                content = list(field = "data_format",
                                               value = "MAF"))
    filter2$project_id <- list(op = "in",
                               content = list(field = "cases.project.project_id",
                                              value = paste("TCGA", vCancer,
                                                            sep = "-")))
    filter2$data_category <- list(op = "in",
                                  content = list(field = "data_category",
                                                 value = filter$data_category))
    filter2$data_type <- list(op = "in",
                              content = list(field = "data_type",
                                             value = filter$data_type))
    filter2$experimental_strategy <- list(op = "in",
                                          content = list(field = "experimental_strategy",
                                                         value = filter$experimental_strategy))
    filter2$platform <- list(op = "in",
                             content = list(field = "platform",
                                            value = filter$platform))
    filterAll <- list(op = "and",
                      content = list(filter2$access,
                                     filter2$data_format,
                                     filter2$project_id,
                                     filter2$data_category,
                                     filter2$data_type,
                                     filter2$experimental_strategy,
                                     filter2$platform))
    payload <- list(filters = filterAll,
                    format = "TSV",
                    sort = "file_id",
                    from = 0,
                    size = entityCount,
                    fields = paste(fieldsMeta, collapse = ","))
    cat(toJSON(payload), file = paste(tmpDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", arg, stdout = T)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    metaData <- tryCatch(read.csv(out,
                                  sep = "\t",
                                  row.names = NULL,
                                  as.is = T,
                                  na.strings = ""),
                         error = function(e){print(e$message); return(NULL)})
    if (!is.null(metaData)) {  # > 0 lines available in input
        colNames <- c("archive.file_name",
                      "cases.0.project.project_id",
                      "data_category",
                      "data_type",
                      "experimental_strategy",
                      "file_id",
                      "file_name",
                      "file_size",
                      "platform",
                      "updated_datetime")
        metaData <- metaData[, colNames]
        rownames(metaData) <- metaData[, "file_id"]
        write.table(metaData,
                    file = paste(tmpDir,
                                 "/tmp_metadata_byType__",
                                 arch,
                                 "__",
                                 endp,
                                 "__",
                                 paste(vCancer, collapse = "_"),
                                 "__",
                                 sAssay,
                                 ".tsv",
                                 sep = ""),
                    quote = F,
                    sep = "\t",
                    col.names = NA,
                    row.names = T)
    }
    print("metadata file: preparing done!")
    return(metaData)
}


#' Get metadata of spefified assay platform files with defined filter
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
#' @example
#' metaData <- function(vCancer = "BRCA",
#' 										 sAssay = "gene_RNAseq",
#'										 sampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
#' 										 tmpDir = ".",
#' 										 arch = "legacy",
#' 										 fieldsMeta = "",
#' 										 entityCount = (-1),
#' 										 endp = "files")
MetaData <- function(vCancer = "BRCA",
                     sAssay,
                     sampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
                     tmpDir = ".",
                     arch = "legacy",
                     fieldsMeta = "",
                     entityCount = (-1),
                     endp = "files") {
    library(rjson)
    if (entityCount == (-1)) {
        entityCount <- EntityCount(arch, endp)
    }
    if (fieldsMeta == "") {
        fieldsMeta <- FieldsMeta()
    } else {
        fieldsList <- FieldsList(arch, endp)
        print("fields names: checking ...")
        stopifnot(all(fieldsMeta %in% fieldsList))
        print("fields names: checking done!")
    }
    filter <- Filter(sAssay)
    print("metadata file: preparing ...")
    out <- paste(tmpDir, "/tmp_metadata_", endp, ".tsv", sep = "")
    #url <- paste('"https://gdc-api.nci.nih.gov/',
    url <- paste("https://api.gdc.cancer.gov/",
                 arch,
                 ifelse(arch == '', '', '/'), endp,
                 sep = '')
    opt <- paste("-o ", out,
                 " --silent --show-error --request POST",
                 " --header Content-Type:application/json --data @",
                 tmpDir, "/tmp_metadata.json",
                 sep = "")
    arg <- paste(opt, url)
    filter2 <- list()
    filter2$access <- list(op = "=" ,
                           content = list(field = "access",
                                          value = "open"))
    filter2$data_format <- list(op = "in",
                                content = list(field = "data_format",
                                               value = "TXT"))
    filter2$project_id <- list(op = "in",
                               content = list(field = "cases.project.project_id",
                                              value = paste("TCGA", vCancer,
                                                            sep = "-")))
    filter2$data_category <- list(op = "in",
                                  content = list(field = "data_category",
                                                 value = filter$data_category))
    filter2$data_type <- list(op = "in",
                              content = list(field = "data_type",
                                             value = filter$data_type))
    filter2$experimental_strategy <- list(op = "in",
                                          content = list(field = "experimental_strategy",
                                                         value = filter$experimental_strategy))
    filter2$platform <- list(op = "in",
                             content = list(field = "platform",
                                            value = filter$platform))
    names(sampleTypeId) <- NULL  # avoid names added into tmp_metadata.json
    fieldSampleTypeId <- "cases.samples.sample_type_id"
    filterSampleTypeId <- list(op = "in",
                               content = list(field = fieldSampleTypeId,
                                              value = sampleTypeId))
    filterAll <- list(op = "and",
                      content = list(filter2$access,
                                     filter2$data_format,
                                     filter2$project_id,
                                     filter2$data_category,
                                     filter2$data_type,
                                     filter2$experimental_strategy,
                                     filter2$platform,
                                     filterSampleTypeId))
    payload <- list(filters = filterAll,
                    format = "TSV",
                    sort = "file_id",
                    from = 0,
                    size = entityCount,
                    fields = paste(fieldsMeta, collapse = ","))
    cat(toJSON(payload), file = paste(tmpDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", arg, stdout = T)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    metaData <- tryCatch(read.csv(out,
                                  sep = "\t",
                                  row.names = NULL,
                                  as.is = T,
                                  na.strings = "",
                                  colClasses = "character"),
                         error = function(e){print(e$message); return(NULL)})
    if (!is.null(metaData)) {  # > 0 lines available in input
        if (length(grep("\\.1\\.", colnames(metaData))) > 0) {
            metaData <- MetaCut(metaData, colPattern = "\\.1\\.")
        }
        metaCut <- metaData[grep(filter$file_name, metaData$file_name),
                            sort(colnames(metaData))]  # filter with filename
        if (nrow(metaCut) == 0) {
            metaData <- NULL
        } else {
            metaCut <- metaCut[order(metaCut$file_name), ]  # sort before split !!
            if (any(table(metaCut[, filter$submitter_id]) != 1)) { # duplicated files with same barcode
                boolRow <- ArchiveNewestInGroups(metaCut$archive.file_name,
                                                 metaCut$file_name)
                # not filter$submitter_id because of one barcode to multi files
                metaCut <- metaCut[boolRow, ]
            }
            stopifnot(all(table(metaCut[, filter$submitter_id]) == 1)) # duplicated files with same barcode
            rownames(metaCut) <- metaCut[, filter$submitter_id]
            metaData <- metaCut[order(rownames(metaCut)), order(colnames(metaCut))] # sort by barcode
            stopifnot(all(rownames(metaData) == metaData[filter$submitter_id]))
            write.table(metaData,
                        file = paste(tmpDir,
                                     "/tmp_metadata_byType__",
                                     arch,
                                     "__",
                                     endp,
                                     "__",
                                     paste(vCancer, collapse = "_"),
                                     "__",
                                     sAssay,
                                     ".tsv",
                                     sep = ""),
                        quote = F, sep = "\t", col.names = NA, row.names = T)
        }
    }
    print("metadata file: preparing done!")
    return(metaData)
}

#' Download files by \code{file_id}
#'
#' @param fileId2Bar Character vector of barcode with file_id as the name.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @return Character vector of downloaded filename with file_id as the name.
#' @example
#' fileNameById <- FileNameById(metaData$file_id, tmpDir = ".")
FileNameById <- function(fileId2Bar,
                         tmpDir = ".",
                         arch = "legacy") {
    if (length(fileId2Bar) == 0) {  # no files available for this filter
        fileNameById <- NULL
    } else {
        stopifnot(length(fileId2Bar) > 0)  # no files available for this filter
        library(rjson)
        tmpTar <- paste(tmpDir, "/gdc_download_", TimeNow(), ".tar.gz", sep = "")
        #url <- paste("https://gdc-api.nci.nih.gov/",
        url <- paste("https://api.gdc.cancer.gov/",
                     arch,
                     ifelse(arch == "", "", "/"),
                     "data",
                     sep = "")
        opt <- paste("-o ", tmpTar,
                     " --silent --show-error --request POST ",
                     "--header Content-Type:application/json --data @",
                     tmpDir, "/tmp_id.json",
                     sep = "")
        arg <- paste(opt, url)
        cat(toJSON(list(ids = names(fileId2Bar))),
            file = paste(tmpDir, "/tmp_id.json", sep = ""))
        print("*.tar.gz file: downloading & unzipping ...")
        err <- "error"
        while (err != 0) {
            stdOut <- system2("curl", arg, stdout = T)
            if (!is.null(attr(stdOut, "status"))) {
                print("error (download): check the proxy")
            }
            stopifnot(is.null(attr(stdOut, "status")))
            tmpUntar <- strsplit(tmpTar, split = "\\.")[[1]][1]
            err <- tryCatch(untar(tmpTar, exdir = tmpUntar, tar = "internal"),
                            error = function(e){return(e$message)})
            if (length(fileId2Bar) == 1) {  # only 1 unzipped file in *.tar.gz
                file.rename(from = tmpTar,
                            to = paste(tmpUntar, "/", fileId2Bar[1], ".tsv", sep = ""))
                manifest <- data.frame(filename = paste(fileId2Bar[1],
                                                        ".tsv",
                                                        sep = ""))
                write.table(manifest,
                            file = paste(tmpUntar, "/MANIFEST.txt", sep = ""),
                            quote = F,
                            sep = "\t",
                            col.names = T,
                            row.names = F)
                err <- 0
            }
        }
        manifest <- read.csv(paste(tmpUntar, "MANIFEST.txt", sep = "/"),
                             sep = "\t",
                             row.names = NULL,
                             as.is = T,
                             na.strings = "")
        fileNameById <- paste(tmpUntar, manifest$filename, sep = "/")
        names(fileNameById) <- fileId2Bar[manifest$id]
    }
    if ("data" %in% dir()) {file.remove("data")}
    print("*.tar.gz file: downloading & unzipping done!")
    return(fileNameById[order(names(fileNameById))])  # sort by file_id
}


#  =============================================================================
#  merge functions, NOT used directly by user
#  =============================================================================

#' Merge copy number variations ("cna") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeCopy <- function(fileNameById) {
    colNames <- c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
    ld <- ColumnsFromFiles(fileNameById,
                           colNames,
                           sortBy = "",
                           skipLines = 0,
                           naStrings = "NA")
    dMerged <- Reduce(rbind, ld)
    dMerged <- as.data.frame(as.matrix(dMerged), stringsAsFactors = F)
    return(dMerged)
}

#' Merge methyloation ("methy") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMethy <- function(fileNameById) {
    vValue <- c("Beta_value")
    vProbe <- c("Composite.Element.REF",
                "Gene_Symbol",
                "Chromosome",
                "Genomic_Coordinate")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 1,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c("CpG",
                           "Gene_Symbol",
                           "Chromosome",
                           "Genomic_Coordinate",
                           names(fileNameById))
    return(dMerged)
}

#' Merge microRNA ("mir") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMir <- function(fileNameById) {
    vValue <- c("read_count", "reads_per_million_miRNA_mapped")
    vProbe <- c("miRNA_ID")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById) * 2)
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("miRNA_ID",
                       rep(names(fileNameById), each = 2)),
                     c("miRNA_ID",
                       rep(c("read_count",
                             "reads_per_million_miRNA_mapped"),
                           length(fileNameById))),
                     dMerged)
    return(dMerged)
}

#' Merge microRNA isoform ("mirIsoform") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMirIso <- function(fileNameById) {
    vValue <- c("read_count", "reads_per_million_miRNA_mapped")
    vProbe <- c("isoform_coords", "miRNA_ID", "miRNA_region", "cross.mapped")
    ld <- ColumnsFromFiles(fileNameById,
                           c(vProbe, vValue),
                           sortBy = vProbe[1],
                           skipLines = 0,
                           naStrings = "NA")
    lv <- lapply(ld, function(x){paste(x[, "isoform_coords"],
                                       x[, "miRNA_ID"],
                                       x[, "miRNA_region"],
                                       x[, "cross.mapped"],
                                       sep = "|")})
    vMirIsoId <- sort(unique(unlist(lv)))
    for (n in seq(length(ld))) {
        rownames(ld[[n]]) <- lv[[n]]
    }
    m <- matrix(nrow = length(vMirIsoId),
                ncol = length(fileNameById) * 2)
    for (n in seq(length(fileNameById))) {
        m[, 2 * n - 1] <- ld[[n]][vMirIsoId, "read_count"]
        m[, 2 * n  ] <- ld[[n]][vMirIsoId, "reads_per_million_miRNA_mapped"]
    }
    dMerged <- cbind(Reduce(rbind, strsplit(vMirIsoId, split = "\\|")), m)
    dMerged <- rbind(c(vProbe, rep(names(fileNameById), each = 2)),
                     c(vProbe,
                       rep(c("read_count",
                             "reads_per_million_miRNA_mapped"),
                           length(ld))),
                     dMerged)
    return(dMerged)
}

#' Merge gene array ("gene_Array") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneArray <- function(fileNameById) {
    vValue <- c("log2.lowess.normalized..cy5.cy3..collapsed.by.gene.symbol")
    vProbe <- c("Composite.Element.REF")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 1,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("gene_id", rep(names(fileNameById), each = 1)),
                     dMerged)
    return(dMerged)
}

#' Merge gene ("gene.normalized_RNAseq") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneRnaSeqNorm <- function(fileNameById) {
    vValue <- c("normalized_count")
    vProbe <- c("gene_id")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("gene_id", rep(names(fileNameById), each = 1)), dMerged)
    return(dMerged)
}

#' Merge unnormalized gene ("gene_RNAseq") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneRnaSeqUnnorm <- function(fileNameById) {
    vValue <- c("raw_count", "scaled_estimate")
    vProbe <- c("gene_id")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById) * 2)
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("gene_id", rep(names(fileNameById), each = 2)),
                     c("gene_id", rep(c("raw_count", "scaled_estimate"),
                                      length(fileNameById))),
                     dMerged)
    return(dMerged)
}

#' Merge gene ("isoform.normalized_RNAseq") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneIsoRnaSeqNorm <- function(fileNameById) {
    vValue <- c("normalized_count")
    vProbe <- c("isoform_id")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("isoform_id", rep(names(fileNameById), each = 1)),
                     dMerged)
    return(dMerged)
}

#' Merge unnormalized gene ("isoform_RNAseq") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneIsoRnaSeqUnnorm <- function(fileNameById) {
    vValue <- c("raw_count", "scaled_estimate")
    vProbe <- c("isoform_id")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById) * 2)
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, 2 * n - 1] <- d[vName, vValue[1]]
            m[, 2 * n  ] <- d[vName, vValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("isoform_id", rep(names(fileNameById), each = 2)),
                     c("isoform_id", rep(c("raw_count", "scaled_estimate"),
                                         length(fileNameById))),
                     dMerged)
    return(dMerged)
}

#' Merge exon ("exon") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeExon <- function(fileNameById) {
    vValue <- c("RPKM")
    vProbe <- c("exon")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe,
                               vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("exon", rep(names(fileNameById), each = 1)),
                     dMerged)
    return(dMerged)
}

#' Merge exon ("exonJunction") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeExonJunction <- function(fileNameById) {  # duplicated rows
    vValue <- c("raw_counts")
    vProbe <- c("junction")
    for (n in seq(length(fileNameById))) {
        d <- ColumnsFromFile(fileNameById[n],
                             names(fileNameById)[n],
                             c(vProbe, vValue),
                             sortBy = vProbe[1],
                             skipLines = 0,
                             naStrings = "NA")
        if (n == 1) {
            vName <- sort(rownames(d))
            stopifnot(length(vName) == length(unique(vName)))
            dProbe <- d[vName, vProbe]
            m <- matrix(nrow = nrow(d), ncol = length(fileNameById))
            m[, n] <- d[vName, vValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == vName))
            m[, n] <- d[vName, vValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    dMerged <- rbind(c("exonJunction", rep(names(fileNameById), each = 1)),
                     dMerged)
    return(dMerged)
}

#' Merge protein ("protein_RPPA") files, distributed by \code{Merge}
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeProtein <- function(fileNameById) {  # duplicated rows
    vPr <- c("ABL1|c-Abl",
             "ACACA ACACB|ACC_pS79",
             "ACACA|ACC1",
             "ACVRL1|ACVRL1",
             "ADAR|ADAR1",
             "AKT1 AKT2 AKT3|Akt",
             "AKT1 AKT2 AKT3|Akt_pS473",
             "AKT1 AKT2 AKT3|Akt_pT308",
             "AKT1S1|PRAS40_pT246",
             "ANXA1|Annexin-1",
             "ANXA7|Annexin_VII",
             "ARAF|A-Raf",
             "ARAF|A-Raf_pS299",
             "ARID1A|ARID1A",
             "AR|AR",
             "ASNS|ASNS",
             "ATM|ATM",
             "ATP5A1|Oxphos-complex-V_subunitb",
             "AXL|Axl",
             "BAD|Bad_pS112",
             "BAK1|Bak",
             "BAP1|Bap1-c-4",
             "BAX|Bax",
             "BCL2A1|Bcl2A1",
             "BCL2L11|Bim",
             "BCL2L1|Bcl-xL",
             "BCL2|Bcl-2",
             "BECN1|Beclin",
             "BID|Bid",
             "BIRC2|cIAP",
             "BRAF|B-Raf",
             "BRAF|B-Raf_pS445",
             "BRCA2|BRCA2",
             "BRD4|BRD4",
             "CA9|CA9",
             "CASP3|Caspase-3",
             "CASP7|Caspase-7_cleavedD198",
             "CASP8|Caspase-8",
             "CASP9|Caspase-9",
             "CAV1|Caveolin-1",
             "CCNB1|Cyclin_B1",
             "CCND1|Cyclin_D1",
             "CCNE1|Cyclin_E1",
             "CCNE2|Cyclin_E2",
             "CD274|CD274",
             "CD274|PD-L1",
             "CDH1|E-Cadherin",
             "CDH2|N-Cadherin",
             "CDH3|P-Cadherin",
             "CDK1|CDK1",
             "CDK1|CDK1_pY15",
             "CDKN1A|p21",
             "CDKN1B|p27",
             "CDKN1B|p27_pT157",
             "CDKN1B|p27_pT198",
             "CDKN2A|p16_INK4a",
             "CHEK1|Chk1",
             "CHEK1|Chk1_pS296",
             "CHEK1|Chk1_pS345",
             "CHEK2|Chk2",
             "CHEK2|Chk2_pT68",
             "CHGA|Chromogranin-A-N-term",
             "CLDN7|Claudin-7",
             "COG3|COG3",
             "COL6A1|Collagen_VI",
             "COPS5|JAB1",
             "CTLA4|CTLA4",
             "CTNNA1|alpha-Catenin",
             "CTNNB1|beta-Catenin",
             "DIABLO|Smac",
             "DIRAS3|DIRAS3",
             "DPP4|CD26",
             "DVL3|Dvl3",
             "E2F1|E2F1",
             "EEF2K|eEF2K",
             "EEF2|eEF2",
             "EGFR|EGFR",
             "EGFR|EGFR_pY1068",
             "EGFR|EGFR_pY1173",
             "EIF4EBP1|4E-BP1",
             "EIF4EBP1|4E-BP1_pS65",
             "EIF4EBP1|4E-BP1_pT37_T46",
             "EIF4EBP1|4E-BP1_pT70",
             "EIF4E|eIF4E",
             "EIF4G1|eIF4G",
             "ENY2|ENY2",
             "EPPK1|EPPK1",
             "ERBB2|HER2",
             "ERBB2|HER2_pY1248",
             "ERBB3|HER3",
             "ERBB3|HER3_pY1289",
             "ERCC1|ERCC1",
             "ERCC5|ERCC5",
             "ERRFI1|MIG-6",
             "ESR1|ER-alpha",
             "ESR1|ER-alpha_pS118",
             "ETS1|ETS-1",
             "EZH2|EZH2",
             "FASN|FASN",
             "FN1|Fibronectin",
             "FOXM1|FoxM1",
             "FOXO3|FOXO3a",
             "FOXO3|FOXO3a_pS318_S321",
             "G6PD|G6PD",
             "GAB2|GAB2",
             "GAPDH|GAPDH",
             "GATA3|GATA3",
             "GATA6|GATA6",
             "GSK3A GSK3B|GSK3-alpha-beta",
             "GSK3A GSK3B|GSK3-alpha-beta_pS21_S9",
             "GSK3A GSK3B|GSK3_pS9",
             "GUSP4|DUSP4",
             "GYG1|GYG-Glycogenin1",
             "GYS1|GYS",
             "GYS1|GYS_pS641",
             "HIF1A|HIF-1_alpha",
             "HSPA1A|HSP70",
             "IGFBP2|IGFBP2",
             "IGFR1|IGF1R_pY1135_Y1136",
             "INPP4B|INPP4B",
             "IRF1|IRF-1",
             "IRS1|IRS1",
             "ITGA2|CD49b",
             "JAK2|Jak2",
             "JUN|c-Jun_pS73",
             "KAT2A|GCN5L2",
             "KDR|VEGFR2",
             "KEAP1|KEAP1",
             "KIT|c-Kit",
             "KRT5|CK5",
             "LCK|Lck",
             "LCN2|LCN2a",
             "LDHA|LDHA",
             "LDHB|LDHB",
             "MACC1|MACC1",
             "MAP2K1|MEK1",
             "MAP2K1|MEK1_pS217_S221",
             "MAPK1 MAPK3|MAPK_pT202_Y204",
             "MAPK14|p38_MAPK",
             "MAPK14|p38_pT180_Y182",
             "MAPK1|ERK2",
             "MAPK8|JNK_pT183_pY185",
             "MAPK9|JNK2",
             "MET|c-Met",
             "MET|c-Met_pY1235",
             "MRE11A|Mre11",
             "MS4A1|CD20",
             "MSH2|MSH2",
             "MSH6|MSH6",
             "MTCO2|Mitochondria",
             "MTOR|mTOR",
             "MTOR|mTOR_pS2448",
             "MYC|c-Myc",
             "MYH11|MYH11",
             "MYH9|Myosin-IIa",
             "MYH9|Myosin-IIa_pS1943",
             "NAPSA|Napsin-A",
             "NDRG1|NDRG1_pT346",
             "NF2|NF2",
             "NFE2L2|Nrf2",
             "NFKB1|NF-kB-p65_pS536",
             "NKX2-1|TTF1",
             "NOTCH1|Notch1",
             "NRAS|N-Ras",
             "NRG1|Heregulin",
             "PARK7|DJ-1",
             "PARP1|PARP-Ab-3",
             "PARP1|PARP1",
             "PARP1|PARP_cleaved",
             "PCNA|PCNA",
             "PDCD1|PDCD1",
             "PDCD4|PDCD4",
             "PDK1|PDK1",
             "PDK1|PDK1_pS241",
             "PEA15|PEA15",
             "PEA15|PEA15_pS116",
             "PECAM1|CD31",
             "PGR|PR",
             "PIK3CA|PI3K-p110-alpha",
             "PIK3R1|PI3K-p85",
             "PKM2|PKM2",
             "PRDX1|PRDX1",
             "PREX1|PREX1",
             "PRKAA1|AMPK_alpha",
             "PRKAA1|AMPK_pT172",
             "PRKCA|PKC-alpha",
             "PRKCA|PKC-alpha_pS657",
             "PRKCB|PKC-pan_BetaII_pS660",
             "PRKCD|PKC-delta_pS664",
             "PTEN|PTEN",
             "PTPN11|SHP-2_pY542",
             "PXN|Paxillin",
             "PYGB|PYGB",
             "PYGB|PYGB-AB2",
             "PYGL|PYGL",
             "PYGM|PYGM",
             "RAB11A RAB11B|Rab11",
             "RAB25|Rab25",
             "RAD50|Rad50",
             "RAD51|Rad51",
             "RAF1|C-Raf",
             "RAF1|C-Raf_pS338",
             "RB1|Rb",
             "RB1|Rb_pS807_S811",
             "RBM15|RBM15",
             "RET|Ret_pY905",
             "RICTOR|Rictor",
             "RICTOR|Rictor_pT1135",
             "RPS6KA1|p90RSK",
             "RPS6KA1|p90RSK_pT359_S363",
             "RPS6KB1|p70S6K",
             "RPS6KB1|p70S6K_pT389",
             "RPS6|S6",
             "RPS6|S6_pS235_S236",
             "RPS6|S6_pS240_S244",
             "RPTOR|Raptor",
             "SCD|SCD",
             "SDHB|Complex-II_subunit30",
             "SERPINE1|PAI-1",
             "SETD2|SETD2",
             "SHC1|Shc_pY317",
             "SLC1A5|SLC1A5",
             "SMAD1|Smad1",
             "SMAD3|Smad3",
             "SMAD4|Smad4",
             "SNAI1|Snail",
             "SQSTM1|p62-LCK-ligand",
             "SRC|Src",
             "SRC|Src_pY416",
             "SRC|Src_pY527",
             "SRSF1|SF2",
             "STAT3|STAT3_pY705",
             "STAT5A|STAT5-alpha",
             "STK11|LKB1",
             "STMN1|Stathmin",
             "SYK|Syk",
             "SYP|Synaptophysin",
             "TFRC|TFRC",
             "TGM2|Transglutaminase",
             "TIGAR|TIGAR",
             "TP53BP1|53BP1",
             "TP53|p53",
             "TP63|p63",
             "TSC1|TSC1",
             "TSC2|Tuberin",
             "TSC2|Tuberin_pT1462",
             "TUBA1B|Acetyl-a-Tubulin-Lys40",
             "TYMS|Thymidilate-Synthase",
             "WWTR1|TAZ",
             "XBP1|XBP1",
             "XRCC1|XRCC1",
             "XRCC5|Ku80",
             "YAP1|YAP",
             "YAP1|YAP_pS127",
             "YBX1|YB-1",
             "YBX1|YB-1_pS102",
             "YWHAB|14-3-3_beta",
             "YWHAE|14-3-3_epsilon",
             "YWHAZ|14-3-3_zeta" )
    names(vPr) <- sapply(strsplit(vPr, split = "\\|"), function(x){x[2]})
    colNames <- c("Composite.Element.REF", "Protein.Expression")
    ld  <- ColumnsFromFiles(fileNameById,
                            colNames,
                            sortBy = "Composite.Element.REF",
                            skipLines = 1,
                            naStrings = "NA")
    lvAbValue <- lapply(ld,
                        function(x){
                            ProbeValue(x,
                                       colProbe = "Composite.Element.REF",
                                       colValue = "Protein.Expression",
                                       stripNum = 0)})
    vAb <- sort(unique(unlist(lapply(lvAbValue, names))))
    vPr4Ab <- vPr[StripEnd(vUnstripped = vAb, stripNum = 4,
                           stripEnd = "right")]
    vPr2Ab <- paste(sapply(strsplit(vPr4Ab, split = "\\|"),
                           function(x){x[1]}),
                    vAb, sep = "|")
    names(vPr2Ab) <- vAb
    m <- matrix(nrow = length(vPr2Ab), ncol = length(fileNameById))
    rownames(m) <- vPr2Ab
    colnames(m) <- names(fileNameById)
    for (s in names(fileNameById)) {
        m[, s] <- lvAbValue[[s]][vAb]
    }
    dMerged <- cbind(protein = vPr2Ab, m)
    return(dMerged)
}

#' Main function of \code{Merge}, distribute filenames by assay platform
#'
#' @param fileNameById Character vector of filename (named with file_id).
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @return A \code{data.frame} of merged table.
Merge <- function(fileNameById,
                  sAssay) {
    print("merging files: merging unzipped data files ...")
    if (is.null(fileNameById)) {
        print("merging files: no file satisfies the filter!")
        dMerged <- NULL
        return(dMerged)
    } else if (sAssay %in% c("cna_cnv.hg18",
                             "cna_cnv.hg19",
                             "cna_nocnv.hg18",
                             "cna_nocnv.hg19")) {
        dMerged <- MergeCopy(fileNameById)
    } else if (sAssay %in% c("exonJunction_RNAseq")) {
        dMerged <- MergeExonJunction(fileNameById)
    } else if (sAssay %in% c("exon_RNAseq")) {
        dMerged <- MergeExon(fileNameById)
    } else if (sAssay %in% c("gene_Array")) {
        dMerged <- MergeGeneArray(fileNameById)
    } else if (sAssay %in% c("gene.normalized_RNAseq")) {
        dMerged <- MergeGeneRnaSeqNorm(fileNameById)
    } else if (sAssay %in% c("gene_RNAseq")) {
        dMerged <- MergeGeneRnaSeqUnnorm(fileNameById)
    } else if (sAssay %in% c("isoform.normalized_RNAseq")) {
        dMerged <- MergeGeneIsoRnaSeqNorm(fileNameById)
    } else if (sAssay %in% c("isoform_RNAseq")) {
        dMerged <- MergeGeneIsoRnaSeqUnnorm(fileNameById)
    } else if (sAssay %in% c("methylation_27",
                             "methylation_450")) {
        dMerged <- MergeMethy(fileNameById)
    } else if (sAssay %in% c("mir_GA.hg18",
                             "mir_GA.hg19",
                             "mir_GA.hg19.mirbase20",
                             "mir_HiSeq.hg18",
                             "mir_HiSeq.hg19",
                             "mir_HiSeq.hg19.mirbase20")) {
        dMerged <- MergeMir(fileNameById)
    } else if (sAssay %in% c("mirIsoform_GA.hg18",
                             "mirIsoform_GA.hg19",
                             "mirIsoform_GA.hg19.mirbase20",
                             "mirIsoform_HiSeq.hg18",
                             "mirIsoform_HiSeq.hg19",
                             "mirIsoform_HiSeq.hg19.mirbase20")) {
        dMerged <- MergeMirIso(fileNameById)
    } else if (sAssay %in% c("protein_RPPA")) {
        dMerged <- MergeProtein(fileNameById)
    }
    print("merging files: merging unzipped data files done!")
    return(dMerged)
}


#  =============================================================================
#  adapting functions, internal <> interface, NOT used directly by user
#  =============================================================================

#' Check the user specified parameters
#'
#' @param vCancer String of cancer type.
#' @param vAssay Character vector of assay platform.
#' @param sampleTypeName Character vector of name for sample_type_id.
#' @param assayGroup String of assay platform goup.
#' @return List of checked parameters.
CheckParam <- function(vCancer,
                       vAssay,
                       sampleTypeName,
                       assayGroup) {
    vCancerAll <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                    "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                    "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
                    "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                    "THCA", "THYM", "UCEC", "UCS", "UVM")
    vSampleTypeIdAll <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    names(vSampleTypeIdAll) <-
        c("TP",   # 01, 'Primary Tumor'
          "TR",   # 02, 'Recurrent Tumor'
          "TB",   # 03, 'Primary Blood Derived Cancer - Peripheral Blood'
          "TRBM", # 04, 'Recurrent Blood Derived Cancer - Bone Marrow'
          "TAP",  # 05, 'Additional - New Primary'
          "TM",   # 06, 'Metastatic'
          "TAM",  # 07, 'Additional Metastatic'
          "THOC", # 08, 'Human Tumor Original Cells'
          "TBM",  # 09, 'Primary Blood Derived Cancer - Bone Marrow'
          "NB",   # 10, 'Blood Derived Normal'
          "NT",   # 11, 'Solid Tissue Normal'
          "NBC",  # 12, 'Buccal Cell Normal'
          "NEBV", # 13, 'EBV Immortalized Normal'
          "NBM",  # 14, 'Bone Marrow Normal'
          "CELLC",# 20, 'Control Analyte'
          "TRB",  # 40, 'Recurrent Blood Derived Cancer - Peripheral Blood'
          "CELL", # 50, 'Cell Lines'
          "XP",   # 60, 'Primary Xenograft Tissue'
          "XCL")  # 61, 'Cell Line Derived Xenograft Tissue'
    lAssayGroup <- list(cna = c("cna_cnv.hg18",
                                "cna_cnv.hg19",
                                "cna_nocnv.hg18",
                                "cna_nocnv.hg19"),
                        gene = c("gene_Array",
                                 "gene.normalized_RNAseq",
                                 "gene_RNAseq",
                                 "isoform.normalized_RNAseq",
                                 "isoform_RNAseq",
                                 "exon_RNAseq",
                                 "exonJunction_RNAseq"),
                        methy = c("methylation_27",
                                  "methylation_450"),
                        mir = c("mir_GA.hg18",
                                "mir_GA.hg19",
                                "mir_GA.hg19.mirbase20",
                                "mir_HiSeq.hg18",
                                "mir_HiSeq.hg19",
                                "mir_HiSeq.hg19.mirbase20"),
                        mirIsoform = c("mirIsoform_GA.hg18",
                                       "mirIsoform_GA.hg19",
                                       "mirIsoform_GA.hg19.mirbase20",
                                       "mirIsoform_HiSeq.hg18",
                                       "mirIsoform_HiSeq.hg19",
                                       "mirIsoform_HiSeq.hg19.mirbase20"),
                        protein = c("protein_RPPA"),
                        somatic = c("somaticMutation_DNAseq"),
                        itraq = c("glycoproteome_iTRAQ",
                                  "phosphoproteome_iTRAQ",
                                  "proteome_iTRAQ"))
    vAssaySub <- lAssayGroup[[assayGroup]]
    if (is.null(vCancer)) {
        vCancer <- vCancerAll
    } else if (!all(vCancer %in% vCancerAll)) {
        print(c("cancerType should be 'NULL' (all) or one of: ", vCancerAll))
        stopifnot(!all(vCancer %in% vCancerAll))
    }
    if (is.null(vAssay)) {
        vAssay <- vAssaySub
    } else if (!all(vAssay %in% vAssaySub)) {
        print(c("assayPlatform should be 'NULL' (all) or one of: ", vAssaySub))
        stopifnot(!all(vAssay %in% vAssaySub))
    }
    if (is.null(sampleTypeName)) {
        sampleTypeId <- vSampleTypeIdAll
    } else if (!all(sampleTypeName %in% names(vSampleTypeIdAll))) {
        print(paste("tissueType should be 'NULL' (all) or one of:",
                    paste(vSampleTypeIdAll, collapse = ","),
                    "01) TP = 'Primary Tumor';",
                    "02) TR = 'Recurrent Tumor';",
                    "03) TB = 'Primary Blood Derived Cancer - Peripheral Blood';",
                    "04) TRBM = 'Recurrent Blood Derived Cancer - Bone Marrow';",
                    "05) TAP = 'Additional - New Primary';",
                    "06) TM = 'Metastatic';",
                    "07) TAM = 'Additional Metastatic';",
                    "08) THOC = 'Human Tumor Original Cells';",
                    "09) TBM = 'Primary Blood Derived Cancer - Bone Marrow';",
                    "10) NB = 'Blood Derived Normal';",
                    "11) NT = 'Solid Tissue Normal';",
                    "12) NBC = 'Buccal Cell Normal';",
                    "13) NEBV = 'EBV Immortalized Normal';",
                    "14) NBM = 'Bone Marrow Normal';",
                    "20) CELLC = 'Control Analyte';",
                    "40) TRB = 'Recurrent Blood Derived Cancer - Peripheral Blood';",
                    "50) CELL = 'Cell Lines';",
                    "60) XP = 'Primary Xenograft Tissue';",
                    "61) XCL = 'Cell Line Derived Xenograft Tissue'.",
                    sep = " "))
        stopifnot(all(sampleTypeName %in% names(vSampleTypeIdAll)))
    } else {
        sampleTypeId <- vSampleTypeIdAll[sampleTypeName]
    }
    return(list(vCancer = vCancer, vAssay = vAssay,
                sampleTypeId = sampleTypeId))
}

#' Pipe of metadata, download and merge (general function)
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param barCode Character vector of barcode, to specify patients.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
Pipe <- function(vCancer,
                 sAssay,
                 sampleTypeId,
                 barCode = NULL,
                 tmpDir = ".",
                 arch = "legacy",
                 fieldsMeta = "",
                 entityCount = (-1),
                 endp = "files") {
    vTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(vTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(vCancer = vCancer,
                      sAssay = sAssay,
                      sampleTypeId = sampleTypeId,
                      tmpDir = tmpDir,
                      arch = arch,
                      fieldsMeta = fieldsMeta,
                      entityCount = entityCount,
                      endp = endp)
    if (!is.null(dMeta)) {
        fileId2Bar <- rownames(dMeta)
        names(fileId2Bar) <- dMeta$file_id
        if (!is.null(barCode)) {
            fileId2Bar <- fileId2Bar[ifelse(substr(fileId2Bar, 1, 12) %in%
                                            substr(barCode, 1, 12), T, F)]
        }
        if (length(fileId2Bar) > 0) {
            fileNameById <- FileNameById(fileId2Bar,
                                         tmpDir = tmpDir,
                                         arch = arch)
            dPiped <- Merge(fileNameById, sAssay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        # 						paste(vCancer, collapse = "|"),
        # 						" & assayPlatform = ", sAssay, " & tissueType = ",
        # 						paste(vTissueType[sampleTypeId], collapse = "|"),
        # 						sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (with batch download)
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param barCode Character vector of barcode, to specify patients.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeBatch <- function(vCancer,
                      sAssay,
                      sampleTypeId,
                      barCode = NULL,
                      tmpDir = ".",
                      arch = "legacy",
                      fieldsMeta = "",
                      entityCount = (-1),
                      endp = "files") {
    vTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(vTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(vCancer = vCancer,
                      sAssay = sAssay,
                      sampleTypeId = sampleTypeId,
                      tmpDir = tmpDir,
                      arch = arch,
                      fieldsMeta = fieldsMeta,
                      entityCount = entityCount,
                      endp = endp)
    if (!is.null(dMeta)) {
        fileId2Bar <- rownames(dMeta)
        names(fileId2Bar) <- dMeta$file_id
        if (!is.null(barCode)) {
            fileId2Bar <- fileId2Bar[ifelse(substr(fileId2Bar, 1, 12) %in%
                                            substr(barCode, 1, 12), T, F)]
        }
        if (length(fileId2Bar) > 0) {
            fileNameById <- vector()  # batch download
            n1Batch <- 50
            for (n in seq(ceiling(length(fileId2Bar)/n1Batch)) - 1) {
                nStart <- n * n1Batch + 1
                nStop <- ifelse(n == (ceiling(length(fileId2Bar)/n1Batch) - 1),
                                length(fileId2Bar),
                                (n + 1) * n1Batch)
                vStep <- FileNameById(fileId2Bar[nStart : nStop],
                                      tmpDir = tmpDir,
                                      arch = arch)
                fileNameById <- c(fileNameById, vStep)
            }  # batch download
            dPiped <- Merge(fileNameById, sAssay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        # 						paste(vCancer, collapse = "|"),
        # 						" & assayPlatform = ", sAssay, " & tissueType = ",
        # 						paste(vTissueType[sampleTypeId], collapse = "|"),
        # 						sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for mir data with 705 probes)
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param barCode Character vector of barcode, to specify patients.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeMirLt23k <- function(vCancer,
                         sAssay,
                         sampleTypeId,
                         barCode = NULL,
                         tmpDir = ".",
                         arch = "legacy",
                         fieldsMeta = "",
                         entityCount = (-1),
                         endp = "files") {
    vTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(vTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(vCancer = vCancer,
                      sAssay = sAssay,
                      sampleTypeId = sampleTypeId,
                      tmpDir = tmpDir,
                      arch = arch,
                      fieldsMeta = fieldsMeta,
                      entityCount = entityCount,
                      endp = endp)
    if (!is.null(dMeta)) {
        fileId2Bar <- rownames(dMeta)
        names(fileId2Bar) <- dMeta$file_id
        vLt23k <- as.numeric(dMeta$file_size) < 23000  # file_size < 23000
        fileId2Bar <- fileId2Bar[vLt23k]
        if (!is.null(barCode)) {
            fileId2Bar <- fileId2Bar[ifelse(substr(fileId2Bar, 1, 12) %in%
                                            substr(barCode, 1, 12), T, F)]
        }
        if (length(fileId2Bar) > 0) {
            fileNameById <- FileNameById(fileId2Bar,
                                         tmpDir = tmpDir,
                                         arch = arch)
            dPiped <- Merge(fileNameById, sAssay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        # 						paste(vCancer, collapse = "|"),
        # 						" & assayPlatform = ", sAssay, " & tissueType = ",
        # 						paste(vTissueType[sampleTypeId], collapse = "|"),
        # 						sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for mir data with more than 705 probes)
#'
#' @param vCancer String of cancer type.
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param barCode Character vector of barcode, to specify patients.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeMirGt23k <- function(vCancer,
                         sAssay,
                         sampleTypeId,
                         barCode = NULL,
                         tmpDir = ".",
                         arch = "legacy",
                         fieldsMeta = "",
                         entityCount = (-1),
                         endp = "files") {
    vTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(vTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(vCancer = vCancer,
                      sAssay = sAssay,
                      sampleTypeId = sampleTypeId,
                      tmpDir = tmpDir,
                      arch = arch,
                      fieldsMeta = fieldsMeta,
                      entityCount = entityCount,
                      endp = endp)
    if (!is.null(dMeta)) {
        fileId2Bar <- rownames(dMeta)
        names(fileId2Bar) <- dMeta$file_id
        vGt23k <- as.numeric(dMeta$file_size) > 23000  # file_size > 23000
        fileId2Bar <- fileId2Bar[vGt23k]
        if (!is.null(barCode)) {
            fileId2Bar <- fileId2Bar[ifelse(substr(fileId2Bar, 1, 12) %in%
                                            substr(barCode, 1, 12), T, F)]
        }
        if (length(fileId2Bar) > 0) {
            fileNameById <- FileNameById(fileId2Bar,
                                         tmpDir = tmpDir,
                                         arch = arch)
            dPiped <- Merge(fileNameById, sAssay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        # 						paste(vCancer, collapse = "|"),
        # 						" & assayPlatform = ", sAssay, " & tissueType = ",
        # 						paste(vTissueType[sampleTypeId], collapse = "|"),
        # 						sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for somatic mutation data)
#'
#' @param vCancer String of cancer type
#' @param sAssay String of assay platform, used in \code{Filter}.
#' @param sampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param barCode Character vector of barcode, to specify patients.
#' @param tmpDir String of directory for temporary files.
#' @param arch String of archive type: "legacy" or "".
#' @param fieldsMeta Character vector of colnames in metadata.
#' @param entityCount Number of entity (row) in the metadata: "-1" means all.
#' @param endp String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeSomatic <- function(vCancer,
                        sAssay,
                        sampleTypeId,
                        barCode = NULL,
                        tmpDir = ".",
                        arch = "legacy",
                        fieldsMeta = "",
                        entityCount = (-1),
                        endp = "files") {
    vTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(vTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaDataSoma(vCancer = vCancer,
                          sAssay = sAssay,
                          sampleTypeId = sampleTypeId,
                          tmpDir = tmpDir,
                          arch = arch,
                          fieldsMeta = fieldsMeta,
                          entityCount = entityCount,
                          endp = endp)
    fileId2Bar <- dMeta$file_name
    names(fileId2Bar) <- dMeta$file_id
    fileNameById <- FileNameById(fileId2Bar,
                                 tmpDir = tmpDir,
                                 arch = arch)
    colNames <- c("hugo_symbol",
                  "entrez_gene_id",
                  "center",
                  "ncbi_build",
                  "chromosome",
                  "start_position",
                  "end_position",
                  "strand",
                  "variant_classification",
                  "variant_type",
                  "reference_allele",
                  "tumor_seq_allele1",
                  "tumor_seq_allele2",
                  "dbsnp_rs",
                  "dbsnp_val_status",
                  "tumor_sample_barcode",
                  "matched_norm_sample_barcode",
                  "match_norm_seq_allele1",
                  "match_norm_seq_allele2",
                  "tumor_validation_allele1",
                  "tumor_validation_allele2",
                  "match_norm_validation_allele1",
                  "match_norm_validation_allele2",
                  "verification_status",
                  "validation_status",
                  "mutation_status",
                  "sequencing_phase",
                  "sequence_source",
                  "validation_method",
                  "score",
                  "bam_file",
                  "sequencer",
                  "tumor_sample_uuid",
                  "matched_norm_sample_uuid")
    names(colNames) <- c("Hugo_Symbol",
                         "Entrez_Gene_Id",
                         "Center",
                         "NCBI_Build",
                         "Chromosome",
                         "Start_Position",
                         "End_Position",
                         "Strand",
                         "Variant_Classification",
                         "Variant_Type",
                         "Reference_Allele",
                         "Tumor_Seq_Allele1",
                         "Tumor_Seq_Allele2",
                         "dbSNP_RS",
                         "dbSNP_Val_Status",
                         "Tumor_Sample_Barcode",
                         "Matched_Norm_Sample_Barcode",
                         "Match_Norm_Seq_Allele1",
                         "Match_Norm_Seq_Allele2",
                         "Tumor_Validation_Allele1",
                         "Tumor_Validation_Allele2",
                         "Match_Norm_Validation_Allele1",
                         "Match_Norm_Validation_Allele2",
                         "Verification_Status",
                         "Validation_Status",
                         "Mutation_Status",
                         "Sequencing_Phase",
                         "Sequence_Source",
                         "Validation_Method",
                         "Score",
                         "BAM_File",
                         "Sequencer",
                         "Tumor_Sample_UUID",
                         "Matched_Norm_Sample_UUID")
    ldPiped <- list()
    for (sMaf in fileNameById) {
        dMaf <- read.csv(sMaf,
                         sep = "\t",
                         row.names = NULL,
                         as.is = T,
                         na.strings = "",
                         comment.char = "#")
        if (!is.null(dMaf)) {
            colnames(dMaf) <- tolower(colnames(dMaf))
            dPiped <- cbind(dMaf[, colNames])
            colnames(dPiped) <- names(colNames)
            if (!is.null(barCode)) {
                vbBar <- ifelse(substr(dPiped[, "Tumor_Sample_Barcode"], 1, 12) %in%
                                substr(barCode, 1, 12), T, F)
                if (any(vbBar)) {
                    dPiped <- dPiped[vbBar, , drop = F]
                } else {
                    dPiped <- NULL
                }
            }
            if (!is.null(dPiped) & length(sampleTypeId) < 14) {
                vbSampleTypeId <- ifelse(substr(dPiped[, "Tumor_Sample_Barcode"],
                                                14, 15) %in% sampleTypeId, T, F)
                if (any(vbSampleTypeId)) {
                    dPiped <- dPiped[vbSampleTypeId, , drop = F]
                } else {
                    dPiped <- NULL
                }
            }
            if (!is.null(dPiped)) {
                vbSymbol <- ifelse(dPiped[, "Hugo_Symbol"] %in% c("."), F, T)
                if (any(vbSymbol)) {
                    dPiped <- dPiped[vbSymbol, , drop = F]
                }
            }
        } else {
            # print(paste("metadata = NULL, when cancerType = ",
            # 						paste(vCancer, collapse = "|"),
            # 						" & assayPlatform = ", sAssay, " & tissueType = ",
            # 						paste(vTissueType[sampleTypeId], collapse = "|"),
            # 						sep = ""))
            dPiped <- NULL
        }
        ldPiped[[sMaf]] <- dPiped
    }
    return(ldPiped)
}


#  =============================================================================
#  interface functions, used directly by user
#  =============================================================================

#' DownloadmiRNASeqData: get miRNASeq data, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadmiRNASeqData <- function(cancerType = NULL,
                                 assayPlatform = NULL,
                                 tissueType = NULL,
                                 saveFolderName = ".",
                                 outputFileName = "",
                                 inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "mir")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        dLt23k <- PipeMirLt23k(vCancer = l$vCancer,
                               sAssay = sAssay,
                               sampleTypeId = l$sampleTypeId,
                               barCode = inputPatientIDs,
                               tmpDir = tmpDir,
                               arch = arch,
                               fieldsMeta = fieldsMeta,
                               entityCount = entityCount,
                               endp = endp)
        dGt23k <- PipeMirGt23k(vCancer = l$vCancer,
                               sAssay = sAssay,
                               sampleTypeId = l$sampleTypeId,
                               barCode = inputPatientIDs,
                               tmpDir = tmpDir,
                               arch = arch,
                               fieldsMeta = fieldsMeta,
                               entityCount = entityCount,
                               endp = endp)
        if (!is.null(dLt23k)) {
            sFileName705 <- paste(saveFolderName,
                                  "/",
                                  ifelse(outputFileName == "", "",
                                         paste(outputFileName, "__", sep = "")),
                                  paste(l$vCancer, collapse = "_"),
                                  "__",
                                  sAssay,
                                  "__705mir__",
                                  ifelse(is.null(tissueType), "tissueTypeAll",
                                         paste(tissueType, collapse = "_")),
                                  "__",
                                  TimeNow(),
                                  ".txt",
                                  sep = "")
            write.table(dLt23k,
                        file = sFileName705,
                        quote = F,
                        sep = "\t",
                        col.names = F,
                        row.names = F,
                        na = "")
            rm(dLt23k); gc()  # clear the memory
            vFileName[paste(sAssay, "_705", sep = "")] <- sFileName705
        }
        if (!is.null(dGt23k)) {
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName,
                                           "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dGt23k,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = F,
                        row.names = F,
                        na = "")
            rm(dGt23k); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadmiRisoformData: get miRisoform data, assayPlatform %in% c("mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadmiRisoformData <- function(cancerType = NULL,
                                   assayPlatform = NULL,
                                   tissueType = NULL,
                                   saveFolderName = ".",
                                   outputFileName = "",
                                   inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "mirIsoform")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        dPiped <- Pipe(vCancer = l$vCancer,
                       sAssay = sAssay,
                       sampleTypeId = l$sampleTypeId,
                       barCode = inputPatientIDs,
                       tmpDir = tmpDir,
                       arch = arch,
                       fieldsMeta = fieldsMeta,
                       entityCount = entityCount,
                       endp = endp)
        if (!is.null(dPiped)) {
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName, "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dPiped,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = F,
                        row.names = F,
                        na = "")
            rm(dPiped); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadmiRNASeqDataIncludeIsoform: get miRNASeq data, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20", "mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_GA.hg19.mirbase20", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19", "mirIsoform_HiSeq.hg19.mirbase20")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20", "mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_GA.hg19.mirbase20", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19", "mirIsoform_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadmiRNASeqDataIncludeIsoform <- function(cancerType = NULL,
                                               assayPlatform = NULL,
                                               tissueType = NULL,
                                               saveFolderName = ".",
                                               outputFileName = "",
                                               inputPatientIDs = NULL) {
    vMir <- c("mir_GA.hg18",
              "mir_GA.hg19",
              "mir_HiSeq.hg18",
              "mir_HiSeq.hg19")
    vMirIso <- c("mirIsoform_GA.hg18",
                 "mirIsoform_GA.hg19",
                 "mirIsoform_HiSeq.hg18",
                 "mirIsoform_HiSeq.hg19")
    assayPlatformMir <- assayPlatform[assayPlatform %in% vMir]
    assayPlatformIso <- assayPlatform[assayPlatform %in% vMirIso]
    vFileNameMir <- DownloadmiRNASeqData(cancerType,
                                         assayPlatformMir,
                                         tissueType,
                                         saveFolderName,
                                         outputFileName,
                                         inputPatientIDs)
    vFileNameMirIso <- DownloadmiRisoformData(cancerType,
                                              assayPlatformIso,
                                              tissueType,
                                              saveFolderName,
                                              outputFileName,
                                              inputPatientIDs)
    vFileName <- c(vFileNameMir, vFileNameMirIso)
    return(vFileName)
}

#' DownloadMethylationData: get methylation data, assayPlatform %in% c("methylation_27", "methylation_450")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("methylation_27", "methylation_450")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadMethylationData <- function(cancerType = NULL,
                                    assayPlatform = NULL,
                                    tissueType = NULL,
                                    saveFolderName = ".",
                                    outputFileName = "",
                                    inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "methy")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        dPiped <- PipeBatch(vCancer = l$vCancer,
                            sAssay = sAssay,
                            sampleTypeId = l$sampleTypeId,
                            barCode = inputPatientIDs,
                            tmpDir = tmpDir,
                            arch = arch,
                            fieldsMeta = fieldsMeta,
                            entityCount = entityCount,
                            endp = endp)
        if (!is.null(dPiped)) {
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName, "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dPiped,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = T,
                        row.names = F,
                        na = "")
            rm(dPiped); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadCNAData: get copy number data, assayPlatform %in% c("cna_cnv.hg18", "cna_cnv.hg19", "cna_nocnv.hg18", "cna_nocnv.hg19")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("cna_cnv.hg18", "cna_cnv.hg19", "cna_nocnv.hg18", "cna_nocnv.hg19")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadCNAData <- function(cancerType = NULL,
                            assayPlatform = NULL,
                            tissueType = NULL,
                            saveFolderName = ".",
                            outputFileName = "",
                            inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "cna")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        dPiped <- Pipe(vCancer = l$vCancer,
                       sAssay = sAssay,
                       sampleTypeId = l$sampleTypeId,
                       barCode = inputPatientIDs,
                       tmpDir = tmpDir,
                       arch = arch,
                       fieldsMeta = fieldsMeta,
                       entityCount = entityCount,
                       endp = endp)
        if (!is.null(dPiped)) {
            colnames(dPiped)[1] <- c("Sample")  # for Module_B
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName, "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dPiped,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = T,
                        row.names = F,
                        na = "")
            rm(dPiped); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadRNASeqData: get gene expression data, assayPlatform %in% c("gene_Array", "gene.normalized_RNAseq", "gene_RNAseq", "isoform.normalized_RNAseq", "isoform_RNAseq", "exon_RNAseq", "exonJunction_RNAseq")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("gene_Array", "gene.normalized_RNAseq", "gene_RNAseq", "isoform.normalized_RNAseq", "isoform_RNAseq", "exon_RNAseq", "exonJunction_RNAseq")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadRNASeqData <- function(cancerType = NULL,
                               assayPlatform = NULL,
                               tissueType = NULL,
                               saveFolderName = ".",
                               outputFileName = "",
                               inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "gene")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        if (sAssay %in% c("isoform.normalized_RNAseq",
                          "isoform_RNAseq",
                          "exon_RNAseq",
                          "exonJunction_RNAseq")) {
            dPiped <- PipeBatch(vCancer = l$vCancer,
                                sAssay = sAssay,
                                sampleTypeId = l$sampleTypeId,
                                barCode = inputPatientIDs,
                                tmpDir = tmpDir,
                                arch = arch,
                                fieldsMeta = fieldsMeta,
                                entityCount = entityCount,
                                endp = endp)
        } else {
            dPiped <- Pipe(vCancer = l$vCancer,
                           sAssay = sAssay,
                           sampleTypeId = l$sampleTypeId,
                           barCode = inputPatientIDs,
                           tmpDir = tmpDir,
                           arch = arch,
                           fieldsMeta = fieldsMeta,
                           entityCount = entityCount,
                           endp = endp)
        }
        if (!is.null(dPiped)) {
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName, "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dPiped,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = F,
                        row.names = F,
                        na = "")
            rm(dPiped); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadRPPAData: get protein expression data, assayPlatform %in% c("protein_RPPA")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("protein_RPPA")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
DownloadRPPAData <- function(cancerType = NULL,
                             assayPlatform = NULL,
                             tissueType = NULL,
                             saveFolderName = ".",
                             outputFileName = "",
                             inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "protein")
    vFileName <- rep(NA, length(l$vAssay))
    names(vFileName) <- l$vAssay
    for (sAssay in l$vAssay) {
        dPiped <- Pipe(vCancer = l$vCancer,
                       sAssay = sAssay,
                       sampleTypeId = l$sampleTypeId,
                       barCode = inputPatientIDs,
                       tmpDir = tmpDir,
                       arch = arch,
                       fieldsMeta = fieldsMeta,
                       entityCount = entityCount,
                       endp = endp)
        if (!is.null(dPiped)) {
            fileName <- paste(saveFolderName,
                              "/",
                              ifelse(outputFileName == "", "",
                                     paste(outputFileName, "__", sep = "")),
                              paste(l$vCancer, collapse = "_"),
                              "__",
                              sAssay,
                              "__",
                              ifelse(is.null(tissueType), "tissueTypeAll",
                                     paste(tissueType, collapse = "_")),
                              "__",
                              TimeNow(),
                              ".txt",
                              sep = "")
            write.table(dPiped,
                        file = fileName,
                        quote = F,
                        sep = "\t",
                        col.names = T,
                        row.names = F,
                        na = "")
            rm(dPiped); gc()  # clear the memory
            vFileName[sAssay] <- fileName
        }
    }
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadSomaticMutationData: get somatic mutation, assayPlatform %in% c("somaticMutation_DNAseq")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("somaticMutation_DNAseq")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
#' @examples
#' v <- DownloadSomaticMutationData(cancerType = "BRCA", assayPlatform = NULL, tissueType = NULL, saveFolderName = ".", outputFileName = "", inputPatientIDs = NULL)
DownloadSomaticMutationData <- function(cancerType = NULL,
                                        assayPlatform = NULL,
                                        tissueType = NULL,
                                        saveFolderName = ".",
                                        outputFileName = "",
                                        inputPatientIDs = NULL) {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "somatic")
    vFileName <- NULL
    for (sCancer in l$vCancer) {
        ldPiped <- PipeSomatic(vCancer = sCancer,
                               sAssay = l$vAssay,
                               sampleTypeId = l$sampleTypeId,
                               barCode = inputPatientIDs,
                               tmpDir = tmpDir,
                               arch = arch,
                               fieldsMeta = fieldsMeta,
                               entityCount = entityCount,
                               endp = endp)
        for (sPathname in names(ldPiped)) {
            dPiped <- ldPiped[[sPathname]]
            if (!is.null(dPiped)) {
                fileName <- paste(saveFolderName,
                                  "/",
                                  ifelse(outputFileName == "", "",
                                         paste(outputFileName, "__", sep = "")),
                                  paste(sCancer, collapse = "_"),
                                  "__",
                                  l$vAssay,
                                  "__",
                                  ifelse(is.null(tissueType), "tissueTypeAll",
                                         paste(tissueType, collapse = "_")),
                                  "__",
                                  TimeNow(),
                                  "__",
                                  rev(strsplit(sPathname, split = "/")[[1]])[1],
                                  ".txt",
                                  sep = "")
                write.table(dPiped,
                            file = fileName,
                            quote = F,
                            sep = "\t",
                            col.names = T,
                            row.names = F,
                            na = "")
                vFileName <- c(vFileName, fileName)
            }
            if (.Platform$OS.type == "windows") {
                dir.create(paste(saveFolderName, "/originalSomaticMutationFiles", sep = ""), recursive = T)
                bFileRename <-file.rename(from = sPathname,
                                          to = paste(saveFolderName, "/originalSomaticMutationFiles/",
                                                     rev(strsplit(sPathname, split = "/")[[1]])[1],
                                                     sep = ""))
                stopifnot(all(bFileRename))
            }
        }
    }
    rm(ldPiped); gc()  # clear the memory
    unlink(tmpDir, recursive = T)
    options(warn = 0)
    return(vFileName)
}

#' DownloadCPTACData: get CPTAC data, assayPlatform %in% c("glycoproteome_iTRAQ", "phosphoproteome_iTRAQ", "proteome_iTRAQ")
#'
#' @param canerType (i.e. vCancer, vector of cancer type), length(vCancer)> = 1. Now only c("BRCA", "OV", "COAD", "READ"), "BRCA"->"Breast", "OV"->"OV", c("COAD", "READ")->"Colorectal"
#' @param assayPlatform (i.e. vAssay, vector of type), length(vAssay)> = 1, assayPlatform %in% c("glycoproteome_iTRAQ", "phosphoproteome_iTRAQ", "proteome_iTRAQ")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...). Now only "TP"(01)
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. barCode, vector of barcode), find specified patients
#' @return vFileName (vector of path/filename), could be used by module B
#' @examples
#' v <- DownloadCPTACData(cancerType = NULL, assayPlatform = NULL, tissueType = NULL, saveFolderName = ".", outputFileName = "", inputPatientIDs = NULL)
DownloadCPTACData <- function(cancerType = NULL,
                              assayPlatform = NULL,
                              tissueType = NULL,
                              saveFolderName = ".",
                              outputFileName = "",
                              inputPatientIDs = NULL) {
    barCode <- inputPatientIDs
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    l <- CheckParam(vCancer = cancerType,
                    vAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    assayGroup = "itraq")
    vCancerCptac <- c("BRCA", "OV", "COAD", "READ")
    vCancer <- intersect(l$vCancer, vCancerCptac)
    if (is.null(vCancer) | length(vCancer) == 0) {
        print(c("cancerType should be 'NULL' (for all cancerType) or one of: ",
                vCancerCptac))
        stopifnot(is.null(vCancer) | length(vCancer) == 0)
    } else {
        urlPre <- "https://cptc-xfer.uis.georgetown.edu/publicData/Phase_II_Data/"
        print("CPTAC files  : downloading ...")
        vFileName <- NULL
        for (sCancer in vCancer) {
            if (sCancer %in% c("BRCA")) {
                vPathname <- c("TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Proteome.itraq.tsv",
                               ##"TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv",
                               # "TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r2/TCGA_Breast_BI_Proteome_CDAP.r2.peptides.tsv",
                               # "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.peptides.tsv",
                               # "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.phosphopeptide.itraq.tsv",
                               "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r4/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
                ##"TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
            } else if (sCancer %in% c("COAD", "READ")) {
                vPathname <- c(# "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.peptides.tsv",
                               # "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv",
                               #"TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome.spectral_counts.tsv",
                               "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv")
            } else if (sCancer %in% c("OV")) {
                vPathname <- c("TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Proteome.itraq.tsv",
                               ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.peptides.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.glycopeptide.itraq.tsv",
                               "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r4/TCGA_Ovarian_JHU_Glycoproteome.glycosite.itraq.tsv",
                               ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.glycosite.itraq.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.peptides.tsv",
                               "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Proteome.itraq.tsv",
                               ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.peptides.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.peptides.tsv",
                               # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.phosphopeptide.itraq.tsv",
                               "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r4/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")
                ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")
            }
            for (sAssay in l$vAssay) {
                filter <- Filter(sAssay)
                for (sPathname in vPathname[grep(filter$file_name,vPathname)]) {
                    fileName <- rev(strsplit(sPathname, split = "/")[[1]])[1]
                    url <- paste(urlPre, sPathname, sep = "")
                    out <- paste(tmpDir, "/", fileName, sep = "")
                    opt <- paste("--silent --show-error -o", out)
                    arg <- paste(opt, url)
                    stdOut <- system2("curl", arg, stdout = T)
                    if (!is.null(attr(stdOut, "status"))) {
                        print("error (download): check the proxy")
                    }
                    stopifnot(is.null(attr(stdOut, "status")))
                    # sPathName <- paste(saveFolderName, "/", fileName, sep = "")
                    # bFileRename <- file.rename(from = out, to = sPathName)
                    # stopifnot(all(bFileRename))
                    # fileName <- rev(strsplit(sPathName, split = "/")[[1]])[1]
                    d <- read.csv(out,
                                  skip = 0,
                                  sep = "\t",
                                  row.names = NULL,
                                  as.is = T,
                                  na.strings = "")
                    if (sAssay == "proteome_iTRAQ") {
                        if (sCancer %in% c("COAD", "READ")) {
                            vRowNot <- c("Total")
                            vColNot <- c("Total.Spectral.Counts",
                                         "Total.Unshared.Spectral.Counts")
                            vColDes <- c("Gene",
                                         "Description",
                                         "Organism",
                                         "Chromosome",
                                         "Locus")
                        } else {
                            vRowNot <- c("Mean", "Median", "StdDev")
                            vColNot <- NULL
                            vColDes <- c("Gene",
                                         "Description",
                                         "Organism",
                                         "Chromosome",
                                         "Locus")
                        }
                    } else if (sAssay == "phosphoproteome_iTRAQ") {
                        vRowNot <- NULL
                        vColNot <- NULL
                        vColDes <- c("Phosphosite","Peptide", "Gene", "Organism")
                    } else if (sAssay == "glycoproteome_iTRAQ") {
                        vRowNot <- NULL
                        vColNot <- NULL
                        vColDes <- c("Glycosite","Peptide", "Gene", "Organism")
                    }
                    d <- d[!(d[,vColDes[1]] %in% vRowNot), !(colnames(d) %in% vColNot)]
                    mInfo <- as.matrix(d[, vColDes])
                    mData <- as.matrix(d[, setdiff(colnames(d), vColDes)])
                    for (n in grep("^X[0-9A-Za-z][0-9A-Za-z]\\.", colnames(mData))) {
                        colnames(mData)[n] <- sub("^X", "", colnames(mData)[n])
                    }
                    for (n in grep("^[0-9A-Za-z][0-9A-Za-z]\\.", colnames(mData))) {
                        colnames(mData)[n] <- paste("TCGA-",
                                                    gsub("\\.", "-", colnames(mData)[n]),
                                                    sep = "")
                    }
                    for (n in grep("-[0-9][0-9]*-(Log|Unshared|Spectral)", colnames(mData))) {
                        substr(colnames(mData)[n], 17, 17) <- '.'
                    }
                    for (n in grep("^OVARIAN\\.CONTROL\\.", colnames(mData))) {
                        colnames(mData)[n] <- gsub("\\.", "-", colnames(mData)[n])
                        substr(colnames(mData)[n], 16, 16) <- '.'
                    }
                    if (!is.null(barCode)) {
                        vbCol <- ifelse(substr(colnames(mData), 1, 12) %in%
                                        substr(barCode, 1, 12), T, F)
                        if (!any(vbCol)) {
                            print("CPTAC files  : no inputPatientIDs found!")
                        }
                        mData <- mData[, vbCol, drop = F]
                    }
                    if (length(l$sampleTypeId) < 14) {
                        vbSampleTypeId <- ifelse(substr(colnames(mData), 14, 15) %in%
                                                 l$sampleTypeId, T, F)
                        if (!any(vbSampleTypeId)) {
                            print("CPTAC files  : no tissueType found!")
                        }
                        mData <- mData[, vbSampleTypeId, drop = F]
                    }
                    if (sCancer == "COAD") {
                        vPatient <- 
                            c("TCGA-A6-3807", "TCGA-A6-3808", "TCGA-A6-3810",
                              "TCGA-AA-3518", "TCGA-AA-3525", "TCGA-AA-3526",
                              "TCGA-AA-3529", "TCGA-AA-3531", "TCGA-AA-3534",
                              "TCGA-AA-3552", "TCGA-AA-3554", "TCGA-AA-3558",
                              "TCGA-AA-3561", "TCGA-AA-3664", "TCGA-AA-3666",
                              "TCGA-AA-3672", "TCGA-AA-3684", "TCGA-AA-3695",
                              "TCGA-AA-3710", "TCGA-AA-3715", "TCGA-AA-3818",
                              "TCGA-AA-3848", "TCGA-AA-3864", "TCGA-AA-3986",
                              "TCGA-AA-3989", "TCGA-AA-A004", "TCGA-AA-A00A",
                              "TCGA-AA-A00E", "TCGA-AA-A00F", "TCGA-AA-A00J",
                              "TCGA-AA-A00K", "TCGA-AA-A00N", "TCGA-AA-A00O",
                              "TCGA-AA-A00R", "TCGA-AA-A00U", "TCGA-AA-A010",
                              "TCGA-AA-A017", "TCGA-AA-A01C", "TCGA-AA-A01D",
                              "TCGA-AA-A01F", "TCGA-AA-A01I", "TCGA-AA-A01K",
                              "TCGA-AA-A01P", "TCGA-AA-A01R", "TCGA-AA-A01S",
                              "TCGA-AA-A01T", "TCGA-AA-A01V", "TCGA-AA-A01X",
                              "TCGA-AA-A01Z", "TCGA-AA-A022", "TCGA-AA-A024",
                              "TCGA-AA-A029", "TCGA-AA-A02E", "TCGA-AA-A02H",
                              "TCGA-AA-A02J", "TCGA-AA-A02O", "TCGA-AA-A02R",
                              "TCGA-AA-A02Y", "TCGA-AA-A03F", "TCGA-AA-A03J")
                        vb <- substr(colnames(mData), 1, 12) %in% vPatient
                        mData <- mData[, vb, drop = F]
                    }
                    if (sCancer == "READ") {
                        vPatient <- 
                            c("TCGA-AF-2691", "TCGA-AF-2692", "TCGA-AF-3400",
                              "TCGA-AF-3913", "TCGA-AG-3574", "TCGA-AG-3580",
                              "TCGA-AG-3584", "TCGA-AG-3593", "TCGA-AG-3594",
                              "TCGA-AG-4007", "TCGA-AG-A002", "TCGA-AG-A008",
                              "TCGA-AG-A00C", "TCGA-AG-A00H", "TCGA-AG-A00Y",
                              "TCGA-AG-A011", "TCGA-AG-A014", "TCGA-AG-A015",
                              "TCGA-AG-A016", "TCGA-AG-A01J", "TCGA-AG-A01L",
                              "TCGA-AG-A01N", "TCGA-AG-A01W", "TCGA-AG-A01Y",
                              "TCGA-AG-A020", "TCGA-AG-A026", "TCGA-AG-A02N",
                              "TCGA-AG-A02X", "TCGA-AG-A032", "TCGA-AG-A036")
                        vb <- substr(colnames(mData), 1, 12) %in% vPatient
                        mData <- mData[, vb, drop = F]
                    }
                    fileName1 <- paste(saveFolderName,
                                       "/",
                                       ifelse(outputFileName == "", "",
                                              paste(outputFileName, "__", sep = "")),
                                       paste(sCancer, collapse = "_"),
                                       "__",
                                       sAssay,
                                       "__",
                                       ifelse(is.null(tissueType), "tissueTypeAll",
                                              paste(tissueType, collapse = "_")),
                                       "__",
                                       # ifelse(sCancer == "OV",
                                       strsplit(fileName, split = "_")[[1]][3],
                                       #			""),
                                       "__",
                                       TimeNow(),
                                       ".txt",
                                       sep = "")
                    write.table(cbind(mInfo, mData),
                                file = fileName1,
                                quote = F,
                                sep = "\t",
                                col.names = T,
                                row.names = F,
                                na = "")
                    rm(mInfo,mData); gc()  # clear the memory
                    vFileName <- c(vFileName, fileName1)
                    # #
                    # vnUnshared <- grep("Unshared", colnames(mData))
                    # mDataShUn <- mData[, vnUnshared - 1]  # barcode.Log.Ratio
                    # mDataUnsh <- mData[, vnUnshared    ]  # barcode.Unshared.Log.Ratio
                    # colnames(mDataShUn) <- unlist(strsplit(colnames(mDataShUn),
                    # 																			 split = "-Log-Ratio"))
                    # colnames(mDataUnsh) <- unlist(strsplit(colnames(mDataUnsh),
                    # 																			 split = "-Unshared-Log-Ratio"))
                    # substr(colnames(mDataShUn), 17, 17) <- '.'
                    # substr(colnames(mDataUnsh), 17, 17) <- '.'
                    # mDataShUn <- mDataShUn[, order(colnames(mDataShUn))]
                    # mDataUnsh <- mDataUnsh[, order(colnames(mDataUnsh))]
                    # stopifnot(all(colnames(mDataUnsh) == colnames(mDataShUn)))
                    # if (!is.null(barCode)) {
                    # 	vbCol <- ifelse(substr(colnames(mDataShUn), 1, 12) %in%
                    # 									substr(barCode, 1, 12), T, F)
                    # 	if (!any(vbCol)) {
                    # 		print("CPTAC files  : no inputPatientIDs found!")
                    # 	}
                    # 	mDataShUn <- mDataShUn[, vbCol, drop = F]
                    # 	mDataUnsh <- mDataUnsh[, vbCol, drop = F]
                    # }
                    # if (length(l$sampleTypeId) < 14) {
                    # 	vbSampleTypeId <- ifelse(substr(colnames(mDataShUn), 14, 15) %in%
                    # 													 l$sampleTypeId, T, F)
                    # 	if (!any(vbSampleTypeId)) {
                    # 		print("CPTAC files  : no tissueType found!")
                    # 	}
                    # 	mDataShUn <- mDataShUn[, vbSampleTypeId, drop = F]
                    # 	mDataUnsh <- mDataUnsh[, vbSampleTypeId, drop = F]
                    # }
                    # sFileNameShUn <- paste(saveFolderName,
                    # 											 "/",
                    # 											 ifelse(outputFileName == "", "",
                    # 															paste(outputFileName, "__", sep = "")),
                    # 											 paste(sCancer, collapse = "_"),
                    # 											 "__",
                    # 											 l$vAssay,
                    # 											 "__",
                    # 											 strsplit(fileName, split = "\\.")[[1]][1],
                    # 											 "__LogRatio",
                    # 											 "__",
                    # 											 TimeNow(),
                    # 											 ".txt",
                    # 											 sep = "")
                    # sFileNameUnsh <- paste(saveFolderName,
                    # 											 "/",
                    # 											 ifelse(outputFileName == "", "",
                    # 															paste(outputFileName, "__", sep = "")),
                    # 											 paste(sCancer, collapse = "_"),
                    # 											 "__",
                    # 											 l$vAssay,
                    # 											 "__",
                    # 											 strsplit(fileName, split = "\\.")[[1]][1],
                    # 											 "__LogRatio_Unshared",
                    # 											 "__",
                    # 											 TimeNow(),
                    # 											 ".txt",
                    # 											 sep = "")
                    # write.table(cbind(mInfo, mDataShUn),
                    # 						file = sFileNameShUn,
                    # 						quote = F,
                    # 						sep = "\t",
                    # 						col.names = T,
                    # 						row.names = F,
                    # 						na = "")
                    # write.table(cbind(mInfo, mDataUnsh),
                    # 						file = sFileNameUnsh,
                    # 						quote = F,
                    # 						sep = "\t",
                    # 						col.names = T,
                    # 						row.names = F,
                    # 						na = "")
                    # vFileName <- c(vFileName, sFileNameShUn, sFileNameUnsh)
                    # #
            }}
        }
        unlink(tmpDir, recursive = T)
        options(warn = 0)
        print("CPTAC files  : downloading done!")
        return(vFileName)
    }
}

#' DownloadBiospecimenClinicalData: get biospecimen and clinical data
#'
#' @param canerType String indicating the specified cancer type
#' for which data should be downloaded.
#' Its value can be one of the cancer type abbreviations:
#' \code{c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
#' "GBM", HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM. Please refer to TCGA (https://tcga-data.nci.nih.gov/docs/publications/tcga/) for information about cancer type. The cancer type abbreviation table (Table 1) shows the full cancer type name.
#' @param saveFolderName (string of path to save the merged data): absolute or relative path
#' @param outputFileName (string of filename prefix)
#' @return vFileName (vector of path/filename), could be used by module B
DownloadBiospecimenClinicalData <- function(cancerType = NULL,
                                            saveFolderName =
                                                "./BiospecimenClinicalData",
                                            outputFileName = "") {
    arch <- "legacy"; fieldsMeta <- ""; entityCount <- (-1); endp <- "files"
    options(warn = -1)
    if (saveFolderName != ".") {dir.create(saveFolderName, recursive = T)}
    tmpDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(tmpDir)
    #	if (saveFolderName != ".") {
    #		if (!dir.exists(saveFolderName)) {
    #			dir.create(saveFolderName, recursive = T)
    #		}
    #	} # dir.exists # R >= 3.2
    #	tmpDir <- paste("tmp", TimeNow(), sep = "_")
    #	if (!dir.exists(tmpDir)) {dir.create(tmpDir)} # dir.exists # R >= 3.2
    d <- MetaDataClin(tmpDir = tmpDir,
                      arch = "legacy",
                      fieldsMeta = "",
                      entityCount = (-1),
                      endp = "files")
    vIdName <- d$file_name
    stopifnot(length(vIdName) == length(unique(vIdName))) # duplicated file_name
    names(vIdName) <- d$file_id
    vCancerAll <- c("ACC" , "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                    "ESCA", "GBM" , "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                    "LGG" , "LIHC", "LUAD", "LUSC", "MESO", "OV"  , "PAAD",
                    "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                    "THCA", "THYM", "UCEC", "UCS" , "UVM")
    if (is.null(cancerType)) {
        cancerType <- vCancerAll
    } else if (!all(cancerType %in% vCancerAll)) {
        print(c("cancerType should be 'NULL' (for all cancerType) or one of: ",
                vCancerAll))
        stopifnot(!all(cancerType %in% vCancerAll))
    } else {
        vbCancer <- ifelse(
                           toupper(
                                   sapply(
                                          strsplit(
                                                   sapply(
                                                          strsplit(vIdName,
                                                                   split = "\\."),
                                                          function(x){rev(x)[2]}),
                                                   split = "_"),
                                          function(y){rev(y)[1]})
                                   ) %in% cancerType,
                           T, F)
        vIdName <- vIdName[vbCancer]  # cancerType
    }
    fileNameById <- FileNameById(vIdName,
                                 tmpDir = tmpDir,
                                 arch = "legacy")
    vFileRename <- file.rename(from = fileNameById,
                               to = paste(saveFolderName,
                                          sapply(strsplit(fileNameById,
                                                          split = "/"),
                                                 function(x){rev(x)[1]}),
                                          sep = "/"))
    stopifnot(all(vFileRename))
    if (outputFileName != "") {
        vFileRename <- F
        vFileRename <- file.rename(from = paste(saveFolderName,
                                                vIdName,
                                                sep = "/"),
                                   to = paste(saveFolderName,
                                              paste(outputFileName,
                                                    vIdName,
                                                    sep = "__"),
                                              sep = "/"))
        stopifnot(all(vFileRename))
    }
    unlink(tmpDir, recursive = T)
    vFileName <- paste(saveFolderName, dir(saveFolderName), sep = "/")
    options(warn = 0)
    return(vFileName)
}


#  =============================================================================
#  Check whether this is the most updated version of TCGA-Assembler
#  =============================================================================

library(httr)
library(stringr)

VCwebContent <- try(content(GET("http://www.compgenome.org/TCGA-Assembler/"), as = "text"), silent = TRUE)
if (class(VCwebContent) == "try-error") {
    rm(VCwebContent)
} else {
    VCstartID <- str_locate_all(string = VCwebContent, pattern = "CheckVersionNumber1")
    if (dim(VCstartID[[1]])[1] == 0) {
        rm(VCstartID, VCwebContent)
    } else {
        VCstartID <- VCstartID[[1]][1, "end"]+1
        VCwebContent <- substr(VCwebContent, VCstartID, nchar(VCwebContent))
        VCstartID <- str_locate_all(string = VCwebContent, pattern = "\">")[[1]][1, "end"]+1
        VCendID <- str_locate_all(string = VCwebContent, pattern = "</span>")[[1]][1, "start"]-1
        VCnewestVersionNum <- substr(VCwebContent, VCstartID, VCendID)
        if (VCnewestVersionNum != "2.0.6") {
            writeLines("\n")
            writeLines("***************************************************************")
            writeLines("A new version of TCGA-Assembler is available!")
            writeLines(paste("Please download version ", VCnewestVersionNum,
                             " at www.compgenome.org/TCGA-Assembler/", sep = ""))
            writeLines("***************************************************************")
            writeLines("\n")
        }
        rm(VCstartID, VCwebContent, VCendID, VCnewestVersionNum)
    }
}

#  end
