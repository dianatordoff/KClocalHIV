#load packages
library(seqinr)
library(ape)
library(igraph)
library(purrr)
library(tibble)
library(Hmisc)
library(memisc)
library(stringr)

#alignment <- read.fasta("PRRT.trimmed.KC.aligned.fasta")
dna <- read.dna("PRRT.trimmed.KC.aligned.fasta", format="fasta")
metadata <- read.csv("metadata.csv")


#Use Tamura-Nei (TN93) pairwise distance, R function from "ape" package
class(dna) #must be a DNAbin object

cluster_by_TN93 <- function(dna, min_size = 2, threshold = 0.02) {
  TN93.dist     <- ape::dist.dna(dna, model="TN93", pairwise.deletion = TRUE, as.matrix=TRUE)
  adj_mat       <- TN93.dist <= threshold
  cluster_graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected", diag=FALSE)
  components    <- igraph::components(cluster_graph)
  retain_assignments <- purrr::keep(components$membership, ~ components$csize[.x] >= min_size)
  clusters           <- tibble::tibble(sequence_name = names(retain_assignments), cluster = retain_assignments)
  return(clusters)
}

TN93.02 <- cluster_by_TN93(dna, min_size = 2, threshold = 0.02) 
names(TN93.02) <- c("seqID", "clusterID")
describe(TN93.02$clusterID)

#Note: Clusters were identified based on individual sequences. 
#      Therefore, we must account for clusters composed soley of the same individual, 
#      and define clusters based on 2 or more unique individuals. 

#deduplicate
seqID_to_newnum <- metadata[metadata$seqID %in% names(alignment), c("seqID", "newnum")]
TN93.02$seqID <- with(TN93.02, str_remove_all(seqID , "[']"))
TN93.02 <- merge(TN93.02, seqID_to_newnum, by="seqID", all.x=TRUE)

#create list of newnums with >1 clusterID
TN93.02$count <-1
unique <- with(TN93.02, aggregate(count, by=list(newnum, clusterID), FUN=sum))
names(unique) <- c("newnum", "cluster", "N")
unique$dup <- duplicated(unique$newnum)
unique$dup <- (unique$newnum %in% unique[unique$dup==TRUE, "newnum"])
unique <- subset(unique, dup==TRUE, select=c("newnum", "cluster"))

#create cross walk for combing clusters
unique$n <- ave(unique$cluster, unique$newnum, FUN = seq_along) # create index!
unique <- reshape(unique, direction="wide", timevar = "n", idvar = "newnum") # reshape wide!
l <- list(unique[!is.na(unique$cluster.2),c("cluster.1", "cluster.2")], 
          unique[!is.na(unique$cluster.3),c("cluster.1", "cluster.3")])
xwalk <- rbindlist(l, use.names=FALSE)
names(xwalk) <- c("updated", "clusterID")

#update clusters
TN93.02 <- merge(TN93.02, xwalk, by="clusterID", all.x=TRUE )
TN93.02$final <- with(TN93.02, ifelse(is.na(updated), clusterID, updated))
TN93.02 <- TN93.02[,c("newnum", "final")]
names(TN93.02) <- c("newnum", "clusterID")

#check
describe(TN93.02$newnum) 
TN93.02 <- TN93.02[,1:2]
rm(unique, l, xwalk)

#create variable is.clustered and save
TN93.02$is.clustered.tn02 <- 1
write.csv(TN93.02, "TN93.02.csv")
