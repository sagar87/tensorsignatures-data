library(rtracklayer)
library(Biostrings)

genomePath <- "genome.fa.gz"

# Set constants
ORI <- c("+", "-", "*")
NUC <- c("A", "C", "G", "T")
CHROMATINSTATES = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

PYR <- paste0(rep(rep(NUC, each=4), 6),
             "[", c(rep("C", 96/2), rep("T", 96/2)),
             ">", c(rep(c("A", "G", "T"), each=16), rep(c("A", "C", "G"), each=16)),
             "]", rep(NUC, 24))
TMP <- gsub("\\[|\\]|>", "", PYR)
PUR <- paste0(
  complement(DNAStringSet(substr(TMP, 4, 4))),
  "[", complement(DNAStringSet(substr(TMP, 2, 2))), ">", complement(DNAStringSet(substr(TMP, 3, 3))),
  "]", complement(DNAStringSet(substr(TMP, 1, 1))))
SUB <- c(PYR, PUR)

# Constants for other mutation types
NUCLEOTIDES <- c("A","C","G","T")
d <- as.character(t(outer(c("C","T","A","G"),NUCLEOTIDES,paste, sep="")))
e <- as.character(reverseComplement(DNAStringSet(d)))

DINUCLEOTIDES <- character()
while(length(d)>0){
  dd <- d[1]
  ee <- e[1]
  DINUCLEOTIDES <- c(DINUCLEOTIDES,dd)
  w <- c(1, which(d==ee))
  d <- d[-w]
  e <- e[-w]
}

MNV_LEVELS <- as.character(sapply(strsplit(DINUCLEOTIDES,""), function(d) paste(paste(d, collapse=""), outer(setdiff(NUCLEOTIDES,d[1]), setdiff(NUCLEOTIDES,d[2]), paste, sep=""), sep=">")))


toNCBI <- function(gr, adjust.genome=FALSE) {
    map <- mapSeqlevels(seqlevels(gr), "NCBI")
    grMapped <- renameSeqlevels(gr, map)
    if (adjust.genome==TRUE) {
      genome(grMapped) <- "GRCh37"
    }
    grMapped
}

toUCSC <- function(gr, adjust.genome=FALSE) {
    map <- mapSeqlevels(seqlevels(gr), "UCSC")
    grMapped <- renameSeqlevels(gr, map)
    if (adjust.genome==TRUE) {
      genome(grMapped) <- "hg19"
    }
    grMapped
}

getTrinucleotideContext <- function(vcf, genomePath) {
    t <- scanFa(file=genomePath, resize(granges(vcf), 3, fix="center"))
    stopifnot(length(t) == length(granges(vcf)))
    as.character(t)
}


getTrinucleotideSubs <- function(vcf, tnc) {
    s <- paste0(substr(tnc, 1, 1),"[", ref(vcf), ">", unlist(alt(vcf)),"]", substr(tnc, 3, 3))
    s <- factor(s, levels=SUB)
}


getStrandOrientation <- function(vcf, orientation) {
    h <- findOverlaps(vcf, orientation, ignore.strand=T)
    stopifnot(length(h) == length(granges(vcf))) # sanityCheck
    factor((strand(orientation)[subjectHits(h)]), levels=ORI)
}


getNucleosomeState <- function(vcf, consensusState) {
    h <- findOverlaps(vcf, consensusState, ignore.strand=T)
    stopifnot(length(h) == length(granges(vcf))) # sanityCheck
    factor(mcols(consensusState[subjectHits(h)])$state, levels=c(0,1,2,3))
}


getChromatinState <- function(vcf, consensusState) {
    h <- findOverlaps(vcf, consensusState, ignore.strand=T)
    stopifnot(length(h) == length(granges(vcf))) # sanityCheck
    factor(mcols(consensusState[subjectHits(h)])$state, levels=CHROMATINSTATES)
}


kataegis <- function(vcf, p=1e-3, q=0.1, r=100) {
    d <-  diff(start(vcf))
    w <- d > 1 & diff(as.numeric(seqnames(vcf))) == 0
    #p <- 1e-3 # P N>Kat
    #q <- 0.05 # P Kat>N
    P <- matrix(c(1-p,p,q, 1-q), ncol=2, byrow=TRUE) # Transition matrix
    p0 <- c(1,0)
    s <- c(mean(d[w]), r)
    dw <- d[w]
    l <- length(dw)
    T1 <- T2 <- matrix(0,ncol=l, nrow=2)
    T1[,1] <- log(c(q/(q+p), p/(q+p)))
    lP <- log(P)
    dg <- rbind(dgeom(dw, prob=1/s[1], log=TRUE), dgeom(dw, prob=1/s[2], log=TRUE))
    for(i in 2:l){
        x <- ((T1[,i-1] + lP) + dg[,i])
        T2[1,i] <- (x[1,1] < x[2,1])+1
        T2[2,i] <- (x[1,2] < x[2,2])+1
        T1[1,i] <- x[T2[1,i],1]
        T1[2,i] <- x[T2[2,i],2]
    }
    z <- numeric(l)
    z[l] <- 1
    for(i in l:2){
        z[i-1] <- T2[z[i],i]
    }
    k <- numeric(nrow(vcf))
    k[-1][w][-1] <- z[-l]-1
    k[-nrow(vcf)][w][-1] <- (z[-l]-1) | k[-nrow(vcf)][w][-1]
    #k[-1][w] <- z-1
    #k[-nrow(vcf)][w] <- (k[-nrow(vcf)][w]) | (z-1)
    kf <- factor(k, levels=0:1)
    return(kf)
}

getClustering <- function(vcf){
    tryCatch(
        kataegis(vcf),
        error=function(error_message) {
            return(factor(rep(0, length(vcf)),  levels=0:1))
        }
    )
}


getMNV <- function(vcf){
    if (length(vcf)==0) {
        return(table(factor(NULL, levels=MNV_LEVELS), useNA='a'))
    } else {
        u <- reduce(granges(vcf[which(isMNV(vcf))]))
        u <- u[width(u)>1]
        if(length(u)==0) return(table(factor(NULL, levels=MNV_LEVELS), useNA='a'))
        r <- as.character(ref(vcf)[vcf %over% u])
        h <- subjectHits(findOverlaps(granges(vcf),u))
        a <- as.character(unlist(alt(vcf))[vcf %over% u])
        rr <- sapply(split(r, h), paste, collapse="")
        aa <- sapply(split(a, h), paste, collapse="")
        w <- which(!rr %in% DINUCLEOTIDES)
        rr[w] <- as.character(reverseComplement(DNAStringSet(rr[w])))
        aa[w] <- as.character(reverseComplement(DNAStringSet(aa[w])))
        t <- table(factor(paste(rr,aa,sep=">"), levels=MNV_LEVELS), useNA='a')
        names(t)[length(t)] <- "MNV(other)"
        return(t)
    }
}

isMNV <- function(vcf) {
  d <- diff(start(vcf)) == 1
  w <- c(FALSE, d) | c(d, FALSE)
  return(w)
}


getIndels <- function(vcf){
    r <- ref(vcf)
    r <- substr(r,2,width(r))
    a <- unlist(alt(vcf))
    a <- substr(a,2,width(a))
    c <- r != "" & a != ""
    w <- which(!r %in% c("C","T", DINUCLEOTIDES))
    r[w] <- as.character(reverseComplement(DNAStringSet(r[w])))
    w <- which(width(r)>2)
    l <- cut(width(r)[w],breaks=c(2:9, seq(10,100,10), Inf))
    levels(l) <- c(3:10, levels(l)[-(1:8)])
    r[w] <- as.character(l)
    w <- which(!a %in% c("C","T", DINUCLEOTIDES))
    a[w] <- as.character(reverseComplement(DNAStringSet(a[w])))
    w <- which(width(a)>2)
    l <- cut(width(a)[w],breaks=c(2:9, seq(10,100,10), Inf))
    levels(l) <- c(3:10, levels(l)[-(1:8)])
    a[w] <- as.character(l)
    indel <- as.character(factor((r=="") + 2*(a==""), levels=0:3, labels=c("indel","ins","del","")))
    indel[indel=="ins"] <- paste("ins", a[indel=="ins"], sep="")
    indel[indel=="del"] <- paste("del", r[indel=="del"], sep="")
    i <- c(t(outer(c("del","ins"), c("C","T", DINUCLEOTIDES, levels(l)), paste, sep="")), "indel")
    t <- table(factor(indel, levels=i), useNA='a')
    names(t)[length(t)] <- "indel(other)"
    t
}


## Helper functions data generation functions
getGenome <- function(genomePath) {
  require(Rsamtools)
  genome <- scanFa(genomePath, as="DNAStringSet")
  genome <- genome[grep("^\\d+|X|Y", names(genome))]
  names(genome) <- gsub("\\s.*", "", names(genome))
  return(genome)
}

getGenCode <- function(genCodePath) {
  genCode <- import(genCodePath, format="gtf")
  genCode <- subset(genCode, type=="gene")
  genCode <- toNCBI(genCode)
  seqlevels(genCode, pruning.mode="coarse") <- c(1:22, "X", "Y")
  genCode$symbol <- genCode$gene_name
  genCode
}

getStandardChromosomes <- function(genomePath) {
  require(Rsamtools)
  r <- scanFaIndex(genomePath)[1:24]
  seqlevels(r) <- c(1:22, "X", "Y")
  #genome(r) <- "GRCh37"
  return(r)
}


# functions to generate data
finiteDiff <- function(x) c(0, diff(x, lag=2), 0)

getChromHmm <- function(chrmmHMMData) {
    chrmmHMMList <- GRangesList(chrmmHMMData)
    standardChromosomes <- toUCSC(getStandardChromosomes())
    allChromHMMSegments <- disjoin(unlist(chrmmHMMList))
    allGenomeSegments <- disjoin(unlist(GRangesList(allChromHMMSegments, standardChromosomes)))

    mcols(allGenomeSegments)$TssA <- 0
    mcols(allGenomeSegments)$TssAFlnk <- 0
    mcols(allGenomeSegments)$TxFlnk <- 0
    mcols(allGenomeSegments)$Tx <- 0
    mcols(allGenomeSegments)$TxWk <- 0
    mcols(allGenomeSegments)$EnhG <- 0
    mcols(allGenomeSegments)$Enh <- 0
    mcols(allGenomeSegments)$ZNFRpts <- 0
    mcols(allGenomeSegments)$Het <- 0
    mcols(allGenomeSegments)$TssBiv <- 0
    mcols(allGenomeSegments)$BivFlnk <- 0
    mcols(allGenomeSegments)$EnhBiv <- 0
    mcols(allGenomeSegments)$ReprPC <- 0
    mcols(allGenomeSegments)$ReprPCWk <- 0
    mcols(allGenomeSegments)$Quies <- 0

    for (i in seq(1, length(chrmmHMMData))){
        for (type in sort(unique(mcols(chrmmHMMData[[i]])$name))) {
            col <- as.numeric(strsplit(type, "_")[[1]][1])
            hits <- findOverlaps(allGenomeSegments, chrmmHMMData[[i]][mcols(chrmmHMMData[[i]])$name == type])
            mcols(allGenomeSegments[queryHits(hits)])[, col] <- mcols(allGenomeSegments[queryHits(hits)])[, col] + 1
        }
    }
    return(allGenomeSegments)
}

getRepTimeNew <- function(reduced=T) {
  # make all calculations
  repTime <- granges(repTimeGm12878)
  mcols(repTime) <- DataFrame(mg12878=repTimeGm12878$score, hela=repTimeHela$score, k562=repTimeK569$score, huvec=repTimeHuvec$score, hepg2=repTimeHepg2$score)
  dt <- sapply(mcols(repTime), finiteDiff)  # compute derivative with finite differences
  m <- rowMeans(dt)
  s <- rowMeans(sqrt((dt-m)^2))
  # assign strands
  strand(repTime)[which(m<0)] <- "-"
  strand(repTime)[which(m>0)] <- "+"
  strand(repTime)[which(abs(m)-2*s < 0)] <- "*"
  #strand(repTime)[which(abs(tab$dt)-2*tab$sd<0)] <- "*"
  # set up the whole genome
  #hg19 <- getStandardChromosomes()
  repTime<-sortSeqlevels(repTime)
  repTimeNCBI <- toNCBI(repTime)
  GRCh37 <- getStandardChromosomes()
  noRep <- setdiff(GRCh37, repTimeNCBI, ignore.strand=T)
  allRep <- sort(c(noRep, granges(repTimeNCBI)), ignore.strand=T)

  if (reduced==TRUE) {
    allRep <- reduce(allRep)
    allRep <- sort(allRep, ignore.strand=T)
  }
  return(allRep)
}


getTransStrand <- function(reduced=TRUE) {

  genCodeRed <- reduce(genCode) # reduce length 48075
  # genes on the positive strand
  genCodePos <- subset(genCodeRed, strand(genCodeRed) == "+") # 23992
  # genes on the negative strand
  genCodeNeg <- subset(genCodeRed, strand(genCodeRed) == "-") # 24083
  # now do the disjoin
  genCodeDis <- disjoin(genCodeRed, ignore.strand=T) # 62065

  # now assign disjoined elements
  posHits <- queryHits(findOverlaps(genCodeDis, genCodePos))
  negHits <- queryHits(findOverlaps(genCodeDis, genCodeNeg))
  bothHits <- base::intersect(posHits, negHits)
  onlyPos <- base::setdiff(posHits, bothHits)
  onlyNeg <- base::setdiff(negHits, bothHits)

  # sanity checks
  stopifnot(length(base::intersect(bothHits, onlyPos)) == 0)
  stopifnot(length(base::intersect(bothHits, onlyNeg)) == 0)
  stopifnot(length(base::intersect(onlyPos, onlyNeg)) == 0)
  stopifnot(length(onlyPos) + length(onlyNeg) + length(bothHits) == length(genCodeDis))

  strand(genCodeDis[onlyPos]) <- "+"
  strand(genCodeDis[onlyNeg]) <- "-"
  strand(genCodeDis[bothHits]) <- "*"

  # define standard chromosomes
  GRCh37 <- getStandardChromosomes()

  # get all non-gene associated regions
  noGenes <- setdiff(GRCh37, genCodeDis, ignore.strand=T)
  allRegions <- sort(c(genCodeDis, noGenes), ignore.strand=T)

  #mcols(allRegions)$DIR <- 0
  #allRegions[strand(allRegions)=="+"]$DIR <- 1
  #allRegions[strand(allRegions)=="-"]$DIR <- -1

  if (reduced==TRUE) {
  	allRegions <- reduce(allRegions)
  	allRegions <- sort(allRegions, ignore.strand=T)
  }

  return(allRegions)
}


getConsensusStates <- function(clu, cutoff) {
    EPMatrix <- as.matrix(mcols(EP))
    EPNormed <- EPMatrix / rowSums(EPMatrix)
    consensusArray <-  parRapply(clu, EPNormed, function(x) {
        maxState <- which.max(x)
        if (length(maxState) == 0) {
            return(0)
        }
        maxStateProb <- x[maxState]
        if (maxStateProb > cutoff) {
            return(maxState)
        } else {
            return(0)
        }})

    EPRanges <- granges(EP)
    mcols(EPRanges)$state <- consensusArray
    EPRanges
}


getNucleosomeData <- function(){
    nucleosomes <- read.table(nucleosomePath)
    NucRanges = GRanges(seqnames=nucleosomes$V1,
        ranges=IRanges(start=nucleosomes$V2, end=nucleosomes$V3),
        strand="+")
    NucRanges <- reduce(NucRanges)
    NucRanges <- NucRanges[!width(NucRanges)>2] # get rid of strange points
    NucleosomeCovered <- NucRanges+73
    NucleosomeCovered <- reduce(NucleosomeCovered) # genomic regions within nucleosomes
    NucleosomeCovered <- toNCBI(NucleosomeCovered)
    GRCh37 <- getStandardChromosomes()
    NotNucleosomeCovered <- setdiff(GRCh37, NucleosomeCovered, ignore.strand=T)
    NucleosomePartition <- sort(c(NotNucleosomeCovered, NucleosomeCovered), ignore.strand=T)
    mcols(NucleosomePartition)$state <- 0
    mcols(NucleosomePartition[strand(NucleosomePartition) == "+"])$state <- 1

    return(NucleosomePartition)
}


tncToPyrimidine <- function(vcf){
	a <- unlist(alt(vcf))
	r <- ref(vcf)
	tnc <- DNAStringSet(info(vcf)$TNC)
	rc <- grepl("A|G", r)
	tnc[rc] <- reverseComplement(tnc[rc])
	a[rc] <- reverseComplement(a[rc])
	t <- paste0(substr(tnc,1,1), "[",substr(tnc,2,2), ">",a, "]", substr(tnc,3,3))
	n <- c("A","C","G","T")
	f <- paste0(rep(n, each=4), "[", rep(c("C","T"), each=96/2), ">", c(rep(c("A","G","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
	factor(t, levels=f)
}
