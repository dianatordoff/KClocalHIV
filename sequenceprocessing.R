#load packages
library(seqinr)

#Load data
fastafile <- read.fasta(file = "sequences.fasta")
metadata <- read.csv("metadata.csv")

#Keep PR/RT (protease/reverse transcriptase)
seqPRRT  <- fastafile[c(which(names(fastafile) %in% metadata[metadata$type=="PR/RT",]$seqID))]

#Trim Sequences Lengths
#Drop sequences <917 bp:
lengthPRRT <- as.data.frame(cbind(getLength(seqPRRT), names(seqPRRT)))
names(lengthPRRT) <- c("l", "seqID")
lengthPRRT$l <- as.numeric(lengthPRRT$l)
summary(lengthPRRT$l)
seqPRRT  <- seqPRRT[c(which(names(seqPRRT) %in% lengthPRRT[lengthPRRT$l>=917,]$seqID))] 

#and trim remaining sequences to 917 bp:
trimmedPRRT <- seqPRRT 
for (i in 1:length(seqPRRT)) {
  trimmedPRRT[[i]] <- seqPRRT[[i]][0:917]
}

#Limit to King County Only 
KCtrimmedPRRT  <- trimmedPRRT[c(which(names(trimmedPRRT) %in% metadata[metadata$KC==1,]$seqID))]
write.fasta(KCtrimmedPRRT, names=names(KCtrimmedPRRT), file.out="PRRT.trimmed.KC.fasta")

#created alignment using MAFFT, see https://mafft.cbrc.jp/alignment/software/
#note: alignment subsequently visually inspected for accuracy in Geneious, see https://www.geneious.com/ 

#alternatively, use DECIPHER for multiple sequence alignment
#see https://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.pdf
library(DECIPHER)
dna <- readDNAStringSet("PRRT.trimmed.KC.fasta")
alignment <- AlignSeqs(dna)
writeXStringSet(alignment,"alignment.fasta")

