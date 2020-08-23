#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(rhdf5))
source("/usr/src/app/mutations.R")
load("/usr/src/app/constants.RData")
genomePath <- "/usr/src/app/genome.fa.gz" # path to the reference genome

# load functions and ranges please adjust!
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

outputPath <- args[length(args)] # where the output tensor is saved to

readVcfSave <- function(path) {
    vcf <- readVcf(path)
    vcf <- vcf[seqnames(vcf) %in% c(1:22, "X", "Y")] # filt(vcf) == "PASS" &
    elementLengths <- elementNROWS(alt(vcf))
    vcf <- vcf[elementLengths==1]
}

processVcf <- function(vcf) {
    print("Extracting trinculeotide context")
    tnc <- getTrinucleotideContext(vcf, genomePath)
    print("Extracting substitution types")
    sub <- getTrinucleotideSubs(vcf, tnc)
    print("Extracting transcription states")
    ts <- getStrandOrientation(vcf, TS)
    print("Extracting replication states")
    rs <- getStrandOrientation(vcf, RT)
    print("Extracting epigenetic states")
    ep <- getChromatinState(vcf, EPI)
    print("Extracting nucleosome states")
    nu <- getNucleosomeState(vcf, NUC)
    print("Extracting clustering states")
    cl <- getClustering(vcf)
    print("Creating tensor")
    t <- table(ts=ts, rs=rs, ep=ep, nu=nu, cl=cl, sub=sub)
    t[,,,,,SUB]
}

df2vcf <- function(df) {
  c <- df$chr
  s <- df$pos
  e <- as.integer(s+sapply(as.character(df$ref), nchar, simplify="array")-1)
  g <- GRanges(c, IRanges(s, e), "*")
  genome(g) <- "hg19"
  v <- VCF(rowRanges=g)
  ref(v) <- DNAStringSet(df$ref)
  alt(v) <- DNAStringSetList(lapply(df$alt, function(x) x))
  toNCBI(v)
}

print("Loading vcf(s) ...")
vcf <- lapply(args[1:length(args)-1], readVcfSave)

print("Processing ...")
snvTensor <- sapply(vcf, function(x) processVcf(x[isSNV(x)]), simplify="array")
indelTable <- sapply(vcf, function(x) getIndels(x[isIndel(x)]), simplify="array")
mnvTable <- sapply(vcf, function(x) getMNV(x), simplify="array")

print("Trying to save ...")
status <- h5createFile(outputPath)

print("Writing")
h5write(snvTensor, "SNVR", file=outputPath)
h5write(indelTable, "INDELS", file=outputPath)
h5write(mnvTable, "MNV", file=outputPath)
