# prelim.R

# Load the libraries
library("Biobase")
library("GEOquery")
library("limma")
library("genefilter")
library("org.Hs.eg.db")
library("hgu133plus2.db")
library("illuminaHumanv2.db")
library("GSEABase")
library("Category")
library("SPIA")
library("KEGGREST")
library("KEGGprofile")
library("graph")
library("rgl")
library("pca3d")
require("annotate")
library("gplots")
library("stringr")
library("sva")
library("impute")

options(digits=3)

#################################################
# Some definitions
#################################################

# Source the function definition files
source("Rcode/funcdefs/makelog.R")
source("Rcode/funcdefs/hyperg.test.R")
source("Rcode/funcdefs/Fisher.test.R")

# Data directory
datapath <- "." # lab
keggxmlpath <- "./Metadata/keggxml" # lab

# Read the metadata files; get the list of study ids; also sample informations
gseinfo <- read.table(file="Metadata/gsegpl.txt", 
	header=T, sep="\t") # list of studies; gse, gpl
gplinfo <- read.table(file="Metadata/gplbioc.txt", 
	header=T, sep="\t") # list of platforms; gpl, biocdb
gsminfo <- read.table(file="Metadata/sample_meta.csv",
	header=T, sep="\t") # list of samples
study.ids <- as.character(gseinfo$gse)

# Define the study ids
ssids <- c("GSE13904", "GSE26378", "GSE26440", "GSE4607", "GSE8121", "GSE9692")
sepsisids <- c("GSE28750", "GSE32707", "GSE5772", "GSE6535", "GSE9960")

#####################################################################
# Load pathway information
#####################################################################
#source("Rcode/get_gbyp.R") 

