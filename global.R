library(limma)
library(DT)
library(xlsx)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# GO DB
library(org.Hs.eg.db)
library(GO.db)
x <- Term(GOTERM)
# Swap names and values of x
goAll <- c(names(x))
names(goAll) <- x
names(goAll) <- firstup(names(goAll))
GO <- org.Hs.egGO2EG

# KEGG DB # Note: this action takes long and slows down the app to initialize
# definitions2pathways <- getKEGGPathwayNames(species.KEGG = "hsa", remove = TRUE)
# keggAll <- definitions2pathways$PathwayID
# names(keggAll) <- c(definitions2pathways$Description)
# genes2pathways <- getGeneKEGGLinks(species.KEGG = "hsa", convert = FALSE)
# ll <- split(genes2pathways, genes2pathways$PathwayID)
# ll <- lapply(ll, function(x) x[names(x) == "GeneID"])
# KEGG <- lapply(ll, function(x) as.character(unlist(x)))

keggAll <- NULL
KEGG <- NULL

# REACTOME DB
library("reactome.db")
file <- "../../Google Drive/Files/annotationdatabases/ReactomePathways.gmt"
REACTOME <- readLines(file)
REACTOME <- strsplit(REACTOME, split = "\t")
pathways <- lapply(REACTOME, function(x) x[1])
REACTOME <- lapply(REACTOME, function(x) x[-c(1, 2)])
reactomeAll <- unlist(pathways)
reactomeAll <- firstup(tolower(gsub("_"," ",gsub("REACTOME_", "", reactomeAll))))
names(REACTOME) <- reactomeAll

# MsigDB HALLMARK - Hs.H object
load("../../Google Drive/Files/annotationdatabases/human_H_v5p2.rdata")
Hs.H <- Hs.H
names(Hs.H) <- firstup(tolower(gsub("_", " ", gsub("HALLMARK_", "", names(Hs.H)))))
msigdbAll <- names(Hs.H)

# MsigDB HALLMARK - Hs.c5 object
load("../../Google Drive/Files/annotationdatabases/human_c5_v5p2.rdata")
Hs.c5 <- Hs.c5
msigAll_cancer <- names(Hs.c5)

