setwd("C:/Users/kamoo/OneDrive/Desktop/Karlis/Systems Biology/Network Biology/Project")
library(diffusr)
library(UniprotR)
library(igraph)
library(CINNA)
library(RCy3)
#library(pracma)
library(gplots)
library(arsenal)

#Import Data
#First column = source
#Second column = target
#Third column = score (e.g. 700 = 0.7 confidence)

dat <- read.table("243230.protein.links.v11.0.txt",stringsAsFactors = F)
colnames(dat) <- dat[1,]
dat <- dat[2:520605,]



#Get UniprotIDs from STRING identifiers
#UniprotID <- ConvertID(dat[,1],ID_from = "STRING_ID", ID_to = "ACC+ID")

#write.table(dat[,1],"SourceProteins.txt",col.names=F, row.names=F,quote=F)
#write.table(dat[,2],"TargetProteins.txt",col.names=F, row.names=F,quote=F)



# Connect to cytoscape and import full D. radiodurans network

cytoscapePing ()
#networkSuid = getNetworkSuid()
#setCurrentNetwork(network="Deinococcus_radiodurans_FULL")
ppra <- "243230.DR_A0346"
recA <- "243230.DR_2340"
uvrC <- "243230.DR_1354"
#nodes <- c(ppra,recA,uvrC)

# clearSelection()
# copyVisualStyle("default", "diffusion")
# setVisualStyle("diffusion")
# deleteStyleMapping("diffusion","NODE_FILL_COLOR")
# CyNodes<- selectNodes(nodes, by.col = "shared name")


# recA random walk | PPI with 0.95 confidence

N950 <- createIgraphFromNetwork("950_Deinococcus_radiodurans.sif")
N950_giant <- giant_component_extract(N950)
createNetworkFromIgraph(N950_giant[[1]], title="N950_giant")
Adj_N950 <- as_adjacency_matrix(N950_giant[[1]])
Adj_N950<- as.matrix(Adj_N950)




p0 <- rep(0,length(V(N950_giant[[1]])))
p0[which(vertex_attr(N950_giant[[1]], "name")==recA)] <- 1
graph <- Adj_N950

pt <- random.walk(p0, graph, r = 0.2)
Pt_prob_recA <- as.data.frame(pt[[1]])
row.names(Pt_prob_recA) <- vertex_attr(N950_giant[[1]], "name")
Pt_prob_sorted_recA <- sort(Pt_prob_recA, decreasing = TRUE)


ht <- heat.diffusion(p0, graph)
ht_prob_recA <- as.data.frame(ht)
row.names(ht_prob_recA) <- vertex_attr(N950_giant[[1]], "name")
ht_prob_recA_sorted <- sort(ht_prob_recA, decreasing = TRUE)

# uvrC random walk | PPI with 0.90 confidence

N900 <- createIgraphFromNetwork("900_Deinococcus_radiodurans.sif")
N900_giant <- giant_component_extract(N900)
createNetworkFromIgraph(N900_giant[[1]], title="N900_giant")
Adj_N900 <- as_adjacency_matrix(N900_giant[[1]])
Adj_N900<- as.matrix(Adj_N900)




p1 <- rep(0,length(V(N900_giant[[1]])))
p1[which(vertex_attr(N900_giant[[1]], "name")==uvrC)] <- 1
graph_900 <- Adj_N900

pt_uvrC <- random.walk(p1, graph_900, r = 0.2)
Pt_prob_uvrC <- as.data.frame(pt_uvrC[[1]])
row.names(Pt_prob_uvrC) <- vertex_attr(N900_giant[[1]], "name")
Pt_prob_sorted_uvrC <- sort(Pt_prob_uvrC, decreasing = TRUE)

ht_uvrC <- heat.diffusion(p1, graph_900)
ht_prob_uvrC <- as.data.frame(ht_uvrC)
row.names(ht_prob_uvrC) <- vertex_attr(N900_giant[[1]], "name")
ht_prob_uvrC_sorted <- sort(ht_prob_uvrC, decreasing = TRUE)

# PprA random walk | PPI with 0.70 confidence

N700 <- createIgraphFromNetwork("700_Deinococcus_radiodurans.sif")
N700_giant <- giant_component_extract(N700)
createNetworkFromIgraph(N700_giant[[1]], title="N700_giant")
Adj_N700 <- as_adjacency_matrix(N700_giant[[1]])
Adj_N700<- as.matrix(Adj_N700)




p2 <- rep(0,length(V(N700_giant[[1]])))
p2[which(vertex_attr(N700_giant[[1]], "name")==ppra)] <- 1
graph_700 <- Adj_N700

pt_ppra <- random.walk(p2, graph_700, r = 0.2)
Pt_prob_ppra <- as.data.frame(pt_ppra[[1]])
row.names(Pt_prob_ppra) <- vertex_attr(N700_giant[[1]], "name")
Pt_prob_sorted_ppra <- sort(Pt_prob_ppra, decreasing = TRUE)

ht_ppra <- heat.diffusion(p2, graph_700)
ht_prob_ppra <- as.data.frame(ht_ppra)
row.names(ht_prob_ppra) <- vertex_attr(N700_giant[[1]], "name")
ht_prob_ppra_sorted <- sort(ht_prob_ppra, decreasing = TRUE)


# Match Random Walk output to Heat Diffusion output

a <- match(rownames(Pt_prob_sorted_ppra), rownames(ht_prob_ppra_sorted))
b <- match(rownames(Pt_prob_sorted_recA), rownames(ht_prob_recA_sorted))
c <- match(rownames(Pt_prob_sorted_uvrC), rownames(ht_prob_uvrC_sorted))


# Sort them (by Random Walk) into data frame

Ppra_sorted <- as.data.frame(Pt_prob_sorted_ppra)
Ppra_sorted$heat_diff <- ht_prob_ppra_sorted[a,]
colnames(Ppra_sorted) <- c('RandomWalk', "HeatDiff")

RecA_sorted <- as.data.frame(Pt_prob_sorted_recA)
RecA_sorted$heat_diff <- ht_prob_recA_sorted[b,]
colnames(RecA_sorted) <- c('RandomWalk', "HeatDiff")

UvrC_sorted <- as.data.frame(Pt_prob_sorted_uvrC)
UvrC_sorted$heat_diff <- ht_prob_uvrC_sorted[c,]
colnames(UvrC_sorted) <- c('RandomWalk', "HeatDiff")


# Take top 15 proteins from Random walk; Add the rank of these proteins in heat diffusion

Ppra_top15 <- Ppra_sorted[1:15,]
RecA_top15 <- RecA_sorted[1:15,]
UvrC_top15 <- UvrC_sorted[1:15,]


Ppra_top15$HD_Rank <- a[1:15]        # Rank meaning - row name indicates protein in Random walk 
RecA_top15$HD_Rank <- b[1:15]        # Rank number indicates which position it is in the Heat Diff rank
UvrC_top15$HD_Rank <- c[1:15]



# get uniprot identifiers for top proteins
Uni_ppra <- ConvertID(rownames(Pt_prob_sorted_ppra)[1:15],ID_from = "STRING_ID", ID_to = "ID")  
Uni_recA <- ConvertID(rownames(Pt_prob_sorted_recA)[1:15],ID_from = "STRING_ID", ID_to = "ID")
Uni_uvrC <- ConvertID(rownames(Pt_prob_sorted_uvrC)[1:15],ID_from = "STRING_ID", ID_to = "ID")



Ppra_top15$UniprotID <- gsub("\\_DEIRA*", "", Uni_ppra[,2])
RecA_top15$UniprotID <- gsub("\\_DEIRA*", "", Uni_recA[,2])
UvrC_top15$UniprotID <- gsub("\\_DEIRA*", "", Uni_uvrC[,2])

Ppra_top15$KEGGID <- ConvertID(Uni_ppra[,2], ID_from = "ID", ID_to = "KEGG_ID")[,2]
RecA_top15$KEGGID <- ConvertID(Uni_recA[,2], ID_from = "ID", ID_to = "KEGG_ID")[,2]
UvrC_top15$KEGGID <- ConvertID(Uni_uvrC[,2], ID_from = "ID", ID_to = "KEGG_ID")[,2]


Ppra_top15$STRINGID <- rownames(Ppra_top15)
RecA_top15$STRINGID <- rownames(RecA_top15)
UvrC_top15$STRINGID <- rownames(UvrC_top15)


loadTableData(RecA_top15,data.key.column = "STRINGID", network = "N950_giant")
loadTableData(UvrC_top15,data.key.column = "STRINGID", network = "N900_giant")
loadTableData(Ppra_top15,data.key.column = "STRINGID", network = "N700_giant")


write.table(Ppra_top15$STRINGID, "Ppra_top15_STRING.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(UvrC_top15), "UvrC_top15_STRING.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(RecA_top15$STRINGID, "RecA_top15_STRING_bare.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

ppra_nodes <- selectNodes(rownames(Ppra_top15), network = "N700_giant")
createSubnetwork(ppra_nodes, network = "N700_giant")
############################


looped_pt <- as.data.frame(rep(0, length(p0)))
for (i in 1:9){
  
  looped_pt[,i] <- random.walk(p0, graph, r = 0.1*i)
  
}

row.names(looped_pt) <- vertex_attr(N950_giant[[1]], "name")
colnames(looped_pt) <- seq(0.1,0.9, 0.1)



plot(sort(looped_pt[,2], decreasing = TRUE)[1:30], type = "l", col = 1, 
     main = "RecA random walk with different restart probabilities",
     xlab = "Rank",
     ylab = "Probability",
     ylim = c(0,0.9))

lines(ht_prob_recA_sorted, col = 2)


legend(x = "topright",legend = seq(0.1,0.9, 0.1), col = 1:9, fill = 1:9)
for (i in 2:9){
  
  lines(sort(looped_pt[,i], decreasing = TRUE)[1:30], col = i)
  Sys.sleep(1)
}


##
#PprA
plot(Pt_prob_sorted_ppra[1:15,], type = "p", col = 1, 
     main = "Ppra RWR and HD probabilities",
     xlab = "Rank",
     ylab = "Probability",
     ylim = c(0,0.9))

points(ht_prob_ppra_sorted[1:15,], col = 2)
legend(x = "topright",legend = c("R Walk", "H Diff"), col = 1:2, fill = 1:2)

#RecA
plot(Pt_prob_sorted_recA[1:15,], type = "p", col = 1, 
     main = "RecA RWR and HD probabilities",
     xlab = "Rank",
     ylab = "Probability",
     ylim = c(0,0.9))
text(Pt_prob_sorted_recA[1:4,], labels = RecA_top15$UniprotID[1:4],cex=0.9, font=2, pos=1)
text(x = 5, y = 0.074480969, labels = RecA_top15$UniprotID[5],cex=0.9, font=2, pos=3)
points(ht_prob_recA_sorted[1:15,], col = 2)
legend(x = "topright",legend = c("RWR", "HD"), col = 1:2, fill = 1:2)

#UvrC
plot(Pt_prob_sorted_uvrC[1:15,], type = "p", col = 1, 
     main = "UvrC RWR and HD probabilities",
     xlab = "Rank",
     ylab = "Probability",
     ylim = c(0,0.9))

points(ht_prob_uvrC_sorted[1:15,], col = 2)
legend(x = "topright",legend = c("R Walk", "H Diff"), col = 1:2, fill = 1:2)

################ Plot data

# Decide on color palette
palf <- colorRampPalette(c("dark red","gray80"))

# Plot all probabilities
plot(Pt_prob_sorted_ppra, main = "Distribution of stationary probabilities; PprA", xlab = "Probability")
plot(Pt_prob_sorted_recA, main = "Distribution of stationary probabilities; RecA", xlab = "Probability")
plot(Pt_prob_sorted_uvrC, main = "Distribution of stationary probabilities; UvrC", xlab = "Probability")


# Plot head of each random walk outcome

plot(Pt_prob_sorted_ppra[1:30,], 
     main = "Top 30 Distribution of stationary probabilities; PprA", 
     ylab = "Probability",
     xlab = "Rank", 
     col = palf(30), 
     type = "b")

plot(Pt_prob_sorted_recA[1:30,], 
     main = "Top 30 Distribution of stationary probabilities; RecA", 
     ylab = "Probability",
     xlab = "Rank", 
     col = palf(30), 
     type = "b")
plot(Pt_prob_sorted_uvrC[1:30,], 
     main = "Top 30 Distribution of stationary probabilities; UvrC", 
     ylab = "Probability",
     xlab = "Rank", 
     col = palf(30), 
     type = "b")




write.table(Pt_prob_sorted_ppra[1:33], "")



