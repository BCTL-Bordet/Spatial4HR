# 09_WGCNA.R


#######################################################
# Weighted Gene Correlation Network Analysis (WGCNA) 
#######################################################

library(here)   
library(cowplot)
library(ggdendro)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(knitr)
library(igraph)
library(ggrepel)
library(WGCNA)
library(magrittr)
library(DESeq2)
library(readxl)
library(pgirmess)
library(dplyr)
library(tibble)


#-----------------------------------------------------------------
# I. Network analysis of ductal breast cancer expression data
#-----------------------------------------------------------------

##### 1. Data input and cleaning

options(stringsAsFactors = FALSE)

# ### 1.a Loading expression data
pb_ductals_not_norm <- readRDS('/Users/bengisukarakose/Desktop/final_scripts/00_data/pb_ductals_not_norm.RDS')
pb_ductals_not_norm <- as.matrix(pb_ductals_not_norm)

#filter >>> delete rows which have more than 50% of the different column values as 0
zero_counts <- rowSums(pb_ductals_not_norm == 0)
zero_proportions <- zero_counts / ncol(pb_ductals_not_norm)
pb_ductals_not_norm_filtered <- pb_ductals_not_norm[zero_proportions <= 0.5, ]

#varianceStabilizingTransformation to make the dataset ready for WGCNA
ductal_data <- varianceStabilizingTransformation(pb_ductals_not_norm_filtered)

dim(ductal_data)

datExpr0_ductal = as.data.frame(t(ductal_data))
# datExpr0_ductal is our dataset. columns are genes, rows are patients. 


### 1.b Checking data for excessive missing values and identification of outlier microarray samples

gsg = goodSamplesGenes(datExpr0_ductal, verbose = 3)
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the 
# offending genes and samples from the data:

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are 
# any obvious outliers.

sampleTree = hclust(dist(datExpr0_ductal), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
## no outlier

datExpr <- datExpr0_ductal


### 1.c Loading clinical trait data

# We now read in the trait data and match the samples for which they were measured to the expression samples.

# traitData = read_excel("~/Desktop/bc_cox_model/DUCTAL_CORRECTED.xlsx")
# names(traitData)[2] <- 'name'
# clin_survival <- read.delim("~/Desktop/bc_cox_model/clin_surv_all.txt")
# 
# traitData <- left_join(clin_survival, traitData, by='name')


traitData <- read.csv("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt", sep="")

traitData$KI67_SCORES <- ifelse(traitData$KI67_CATEGORIES == "≤10", 1,
                                ifelse(traitData$KI67_CATEGORIES == "10-20", 2, 
                                       ifelse(traitData$KI67_CATEGORIES == ">20", 3, NA)))

traitData$HER2_SCORES <- ifelse(traitData$HER2_IHC == "1+", 1,
                                ifelse(traitData$HER2_IHC == "0", 0, NA))


# remove columns that hold information we do not need.
allTraits <- traitData[, c(
  "name",
  "AGE",
  "Multifocal",
  "SIZE_MM",
  "GRADE",
  "ER_SCORE",
  "PR_SCORE",
  "Relapse",
  "DISTANT_RELAPSE",
  "OS_status",
  "KI67_SCORES",
  "HER2_SCORES"
)]

allTraits$PR_SCORE <- as.numeric(allTraits$PR_SCORE)
allTraits$OS_status <- as.numeric(allTraits$OS_status)

dim(allTraits)
names(allTraits)

head(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
ductalSamples = rownames(datExpr)
traitRows = match(ductalSamples, allTraits$name)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# We now have the expression data in the variable datExpr, and the corresponding clinical traits in the
# variable datTraits. Before we continue with network construction and module detection, we visualize
# how the clinical traits relate to the sample dendrogram.

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# save(datExpr, datTraits, file = "ductal-01-dataInput.RData")
# 
# 
# write.table(datExpr, file = "datExpr.txt", append = FALSE, quote = F, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE)
# 
# write.table(datTraits, file = "datTraits.txt", append = FALSE, quote = F, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE)


#-----------------------------------------------------------------
# 2.a Automatic network construction and module detection
#-----------------------------------------------------------------

## Extended Data Fig. 8a


# 2 Automatic construction of the gene network and identification of modules
# 2.a Automatic network construction and module detection
# 2.a.1 Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Figure: Analysis of network topology for various soft-thresholding powers. The left panel shows the 
# scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel 
# displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis).

# We choose the power 4, which is the lowest power for which the scale-free topology fit
# index curve flattens out upon reaching a high value (in this case, roughly 0.95).



## 2.a.2 One-step network construction and module detection
# Constructing the gene network and identifying modules 

net = blockwiseModules(datExpr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "breastcaTOM",
                       verbose = 3,
                       maxBlockSize = 17000)
# saveRDS(net, 'net.RDS')
# net$colors contains the module assignment, and net$MEs contains the module eigengenes of the modules.

table(net$colors)

{ #print message 
  print(paste("After filtering a total of", length(which(net$goodGenes==TRUE)),"genes were included in module identification."))
  print(paste("These were clustered into modules by splitting data into", length(net$dendrograms),"blocks."))
  print(paste("A total of", length(table(net$colors))-1,"modules were identified (plus an additional unassigned module (grey))."))
}

bwModuleColors = labels2colors(net$colors)

#create text files with gene and module data
probes <- colnames(datExpr)
bwnet.colors <- c("grey",unique(bwModuleColors)[- which(unique(bwModuleColors)=="grey")])
key    <- 	paste0("M",0:(length(bwnet.colors)-1))
color.key <- cbind(bwnet.colors, key)
colnames(color.key)[1] <- "bwModuleColors"

module.gene.df <- cbind(probes, bwModuleColors,net$colors)
module.gene.df <- merge(module.gene.df,color.key)
module.gene.df <- module.gene.df[,c(2,1,3,4)]
names(module.gene.df) <- c("gene","color","original.number","module")

module.gene.ordered.df <- module.gene.df[order(module.gene.df[,4]),]
# write.table(module.gene.df, file = "module_gene.txt",row.names=F, sep="\t", quote=F)
# write.table(module.gene.ordered.df, file = "module_gene_ordered.txt",row.names=F, sep="\t", quote=F)


# there are 37 modules, labeled 1 through 37 in order of descending size, with sizes 
# ranging from 2416 to 36 genes. The label 0 is reserved for genes outside of all modules. 

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#  Extended Data Fig. 8c - Clustering dendrogram of genes, with dissimilarity based on topological overlap, 
# together with assigned module colors.




# We now save the module assignment and module eigengene information necessary for subsequent analysis.
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "ductal-02-networkConstruction-auto.RData")



#-----------------------------------------------------------------------------
# 3. Relating modules to external information and identifying important genes
#-----------------------------------------------------------------------------

#### 3 Relating modules to external clinical traits 
### 3.a Quantifying module–trait associations

# In this analysis we would like to identify modules that are significantly associated with the 
# measured clinical traits. Since we already have a summary profile (eigengene) for each module, we 
# simply correlate eigengenes with external traits and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

library(pgirmess)
# write.delim(moduleTraitCor, 'moduleTraitCor.txt')
# write.delim(moduleTraitPvalue, 'moduleTraitPvalue.txt')


## Calculate correlations between Module Eigengenes and metadata

# Calculate Module Eigengenes (MEs) and investigate correlations between the different MEs. 
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = bwModuleColors, excludeGrey = TRUE,)
MEs = MEList$eigengenes

# Calculate correlations between the modules
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# (i) Plot hierarchical tree of modules with color names (bwModuleColors) or...
library(ggdendro)
ggdendrogram(METree,rotate = F)

# .. (ii) Plot hierarchical tree of modules with altered module names: M0,M1,M2 etc
MEDiss.module <- MEDiss
colnames.MEDiss.color.key.index <- match(substr(colnames(MEDiss),3,100), color.key[,1])
colnames(MEDiss.module) <- color.key[colnames.MEDiss.color.key.index,2]
rownames(MEDiss.module) <- color.key[colnames.MEDiss.color.key.index,2]
METree.module = hclust(as.dist(MEDiss.module), method = "average");
ggdendrogram(METree.module,rotate = F)

# Reorder MEs so that correlated MEs are side by side
MEs = orderMEs(MEs)

#rename the MEs from their color names to the numbers used in the manuscript
rownames(MEs)=rownames(datTraits)
namekey=color.key

colnames(MEs)=gsub("ME","",colnames(MEs))
colnames(MEs)=namekey[match(colnames(MEs),namekey[,1]),2]

#write out the MEs
# write.table(MEs,file="MEs.txt", quote=F, sep="\t", row.names=TRUE, col.names=TRUE)

#correlate MEs with one another - across all samples
nSamples=nrow(MEs)
betMEcors = cor(MEs)
betMEps = corPvalueStudent(betMEcors, nSamples)
betMEpadj=matrix(p.adjust(betMEps,method="fdr"),ncol= nrow(betMEps),byrow = F)
colnames(betMEpadj)=colnames(betMEps)
rownames(betMEpadj)=rownames(betMEps)
betMEpadj=ifelse(betMEpadj<0.05,1,0)

#simpHeat function
simpHeat=function(cordat,pdat,pal){
  roword=hclust(dist(cordat))
  colord=hclust(dist(t(cordat)))
  cordat=cordat[roword$order,colord$order]
  pdat=pdat[roword$order,colord$order]
  melt=melt(cordat)
  pmelt=melt(pdat)
  pmelt$value[pmelt$value==0]=NA
  colnames(melt)=c("Var1","Var2","Correlation")
  ggplot(melt,aes(x=Var2,y=Var1,fill=Correlation))+geom_tile()+
    geom_rect(data=pmelt, size=1.0, fill=NA, aes(color=as.factor(pmelt$value), xmin=as.numeric(pmelt$Var2)-0.5,xmax=as.numeric(pmelt$Var2)+0.5,
                                                 ymin=as.numeric(pmelt$Var1)-0.5,ymax=as.numeric(pmelt$Var1)+0.5))+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
    scale_fill_gradientn(colors=pal)+scale_color_manual(values=c("black",NA),labels=NULL,guide=F,na.value=NA)
}


library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tibble)

#colors for heatmaps
colorer=colorRampPalette(rev(brewer.pal(11,"RdBu")))

#visualise correlations between MEs in heatmap
simpHeat(betMEcors,betMEpadj,colorer(11)) #correlation between modules

# Extended Data Fig. 8d - Pairwise correlations between WGCNA module eigengenes


## Calculate correlations between Module Eigengenes and metadata
#Find correlation between module eigengenes and metadata traits to identify modules associated with risk."))

#correlate MEs vs numerical versions of traits - samples with rank of risk prognosis
#hier pas al dan niet distinctie tussen poor en good nemen!!!
remove.intcls = "no" #"yes" or "no"
# if(remove.intcls=="yes"){
#   #remove intermediate prognosis
#   datTraits.processed <- datTraits %>% rownames_to_column("sample_id") %>%  filter(Cls.Intermediate != 1) %>% select(-Cls.Intermediate)
# }


if(remove.intcls=="no"){
  #alternatively, do not remove intermediate prognosis!
  datTraits.processed <- datTraits %>% rownames_to_column("sample_id") #only add sample_id
}

# create numeric datTraits
numTraitData=data.frame(lapply(datTraits.processed,as.numeric)) %>% dplyr::select(-sample_id)

# calculate correlations between MEs and trait data
moduleTraitCor = cor(MEs[datTraits.processed$sample_id,], numTraitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits.processed))
padj=matrix(p.adjust(moduleTraitPvalue,method="fdr"),ncol= ncol(moduleTraitCor),byrow = F)
colnames(padj)=colnames(numTraitData)
rownames(padj)=rownames(moduleTraitPvalue)


#(ii) show with all traits
padj.class=ifelse(padj<0.05,1,0)
simpHeat(t(moduleTraitCor),t(padj.class),colorer(11))
#ggsave("ME_trait_correlation_all.tiff")
# Extended Data Fig. 8e - Associations between WGCNA module eigengenes and clinicopathologic variables.


# write.csv(geneModuleMembership, file = "geneModuleMembership.txt")
# write.csv(MMPvalue, file = "MMPvalue.txt")


#GO enrichment via gprofiler2

run.GO = "yes" # "yes" or "no"


if(run.GO == "yes"){
  {
    library(gprofiler2)
    
    module.vector <- key
    module.fdr.df <- c()
    module.bonferroni.df <- c()
    for(i in seq_along(module.vector)){
      print(paste0("calculating module ",module.vector[i]))
      
      # create genes vector for module
      temp.genes <- module.gene.df$gene[module.gene.df$module %in% module.vector[i]]
      
      # Run gprofiler2 GO:BP, and create module list and df
      gost.temp.BP.fdr		 <- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:BP"), correction_method="fdr") #fdr correction
      gost.temp.BP.bonferroni  <- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:BP")) #bonferroni correcion as standard
      gost.temp.CC.fdr 		<- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:CC"), correction_method="fdr") #fdr correction
      gost.temp.CC.bonferroni <- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:CC")) #bonferroni correcion as standard
      gost.temp.MF.fdr 		<- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:MF"), correction_method="fdr") #fdr correction
      gost.temp.MF.bonferroni <- gost(temp.genes, organism = "hsapiens", user_threshold =0.05, sources = c("GO:MF")) #bonferroni correcion as standard
      
      if(length(gost.temp.BP.fdr) > 0)		{temp.BP.fdr 		<-	cbind(rep(module.vector[i],nrow(gost.temp.BP.fdr$result)),  gost.temp.BP.fdr$result[,c(10,9,11,3,1,2,4:8,12:14)])	; colnames(temp.BP.fdr)[1] <- "module"                             }
      if(length(gost.temp.CC.fdr) > 0)		{temp.CC.fdr 		<-	cbind(rep(module.vector[i],nrow(gost.temp.CC.fdr$result)),  gost.temp.CC.fdr$result[,c(10,9,11,3,1,2,4:8,12:14)])	; colnames(temp.CC.fdr)[1] <- "module"                             }
      if(length(gost.temp.MF.fdr) > 0)		{temp.MF.fdr 		<-	cbind(rep(module.vector[i],nrow(gost.temp.MF.fdr$result)),  gost.temp.MF.fdr$result[,c(10,9,11,3,1,2,4:8,12:14)])	; colnames(temp.MF.fdr)[1] <- "module"                             }
      if(length(gost.temp.BP.bonferroni) > 0)	{temp.BP.bonferroni 	<-	cbind(rep(module.vector[i],nrow(gost.temp.BP.bonferroni$result)),  gost.temp.BP.bonferroni$result[,c(10,9,11,3,1,2,4:8,12:14)])    ; colnames(temp.BP.bonferroni)[1] <- "module"   }
      if(length(gost.temp.CC.bonferroni) > 0)	{temp.CC.bonferroni 	<-	cbind(rep(module.vector[i],nrow(gost.temp.CC.bonferroni$result)),  gost.temp.CC.bonferroni$result[,c(10,9,11,3,1,2,4:8,12:14)])    ; colnames(temp.CC.bonferroni)[1] <- "module"   }
      if(length(gost.temp.MF.bonferroni) > 0)	{temp.MF.bonferroni 	<-	cbind(rep(module.vector[i],nrow(gost.temp.MF.bonferroni$result)),  gost.temp.MF.bonferroni$result[,c(10,9,11,3,1,2,4:8,12:14)])    ; colnames(temp.MF.bonferroni)[1] <- "module"   }
      
      if(length(gost.temp.BP.fdr) > 0)		{module.fdr.df <- rbind(module.fdr.df,temp.BP.fdr)}
      if(length(gost.temp.CC.fdr) > 0)		{module.fdr.df <- rbind(module.fdr.df,temp.CC.fdr)}
      if(length(gost.temp.MF.fdr) > 0)		{module.fdr.df <- rbind(module.fdr.df,temp.MF.fdr)}
      if(length(gost.temp.BP.bonferroni) > 0)	{module.bonferroni.df <- rbind(module.bonferroni.df,temp.BP.bonferroni)}
      if(length(gost.temp.CC.bonferroni) > 0)	{module.bonferroni.df <- rbind(module.bonferroni.df,temp.CC.bonferroni)}
      if(length(gost.temp.MF.bonferroni) > 0)	{module.bonferroni.df <- rbind(module.bonferroni.df,temp.MF.bonferroni)}
      
    }
  }
}


# # write GO enrichment output to text file
# write.table(as.matrix(module.fdr.df), file = "GO_enrichment_FDR.txt",row.names=F, sep="\t", quote=F)
# write.table(as.matrix(module.bonferroni.df), file = "GO_enrichment_bonferroni.txt",row.names=F, sep="\t", quote=F)


##### 5 Visualization of networks within R 


#############################################
# Fig. 4a - Create network figure from the modules 
#############################################

#Compute the TOM from your expression data
TOM <- TOMsimilarityFromExpr(
  datExpr,
  power = 4,
  networkType = "unsigned",
  TOMType = "signed",
  corType = "pearson")

# save network
# saveRDS(TOM,"TOM_power_4.Rds")

# load TOM
#TOM <-readRDS("TOM_power_4.Rds")
TOM <- as.matrix(TOM)			

#Extract module labels
bwnet <- net
moduleLabels <- bwnet$colors
moduleColors <- labels2colors(moduleLabels)		

#Extract Edge Data from TOM
threshold <- 0.05 	# Set a threshold for significant TOM interactions

#Convert TOM into an adjacency matrix for igraph
adjacencyMatrix <- TOM > threshold
adjacencyMatrix[lower.tri(adjacencyMatrix, diag = TRUE)] <- 0  # Remove lower triangle and diagonal

#Keep only genes in modules, excluding uncorrelated genes in 'grey' module M0.
adjacencyMatrix <- adjacencyMatrix[bwnet$colors>0,bwnet$colors>0]

#Convert to edge list
keep.genes <- names(bwnet$colors[bwnet$colors > 0]) #remove genes from M0
edgeList <- which(adjacencyMatrix == TRUE, arr.ind = TRUE)
edge.data <- data.frame(
  from = keep.genes[edgeList[, 1]],
  to = keep.genes[edgeList[, 2]]
)
print(paste("There are", nrow(edge.data), "edges"))

#Extract node data
node.data <- data.frame(
  module.colors = moduleColors[bwnet$colors>0],
  gene = keep.genes
)

#create igraph object
# install.packages('influential')
library(influential)
library(igraph)
g <- graph_from_data_frame(d = edge.data, directed=FALSE)

# Extract a graph with only the connected vertices
connected_vertices <- V(g)[degree(g) > 0]
g_connected <- induced_subgraph(g, connected_vertices)

# Adding community information to nodes
node.data.select <- node.data[node.data$gene %in% V(g_connected)$name ,]
node.data.select.reorder <-   node.data.select[ match(names(V(g_connected)), node.data.select$gene),  ]
if(!all(node.data.select.reorder$gene == names(V(g_connected)))){print("ERROR: order of genes in node.data.select do not match those of vertex names in g_connected")}

V(g_connected)$community <- as.factor(node.data.select.reorder$module.colors)
community.color <- data.frame(V(g_connected)$community)
color.key <- as.data.frame(color.key)
community.color$community <- color.key$key[match(community.color[,1], color.key$bwModuleColors)]
V(g_connected)$community <- community.color$community

communities.vector <- unique(V(g_connected)$community)

#create list of subgraphs (one for each module) 
subgraph.list <- list()
for(i in seq_along(communities.vector)){
  subgraph.list[[i]] <- induced_subgraph(g_connected, V(g_connected)[V(g_connected)$community == communities.vector[i]])
  names(subgraph.list)[[i]] <- communities.vector[i]
}   
subgraph.list.backup <- subgraph.list 

# now we will plot the subgraphs
# we will first group them together based on the module correlation plot (see: 'Create module correlation heatmap')

# plot correlation of eigenvalues to see which modules to group together
p=simpHeat(betMEcors,betMEpadj,colorer(11)) #correlation between modules
p

# From our module correlation heatmap we can see 5 groups: We will group the subgraphs accordingly, and merge them into a single figure
#select groups to plot together
names(subgraph.list.backup) #use these elements 

group.1 <- paste0("M", c(4,6,23,27,17,3,11))
group.2 <- paste0("M", c(21,30,20,31,19,13,14,9,28,35,5,7,25))
group.3 <- paste0("M", c(16,24,10,2))
group.4 <- paste0("M", c(8,18,12))
group.5 <- paste0("M", c(1,15,32,37))

#now you will have to put each group somewhere on the plot
#for each group we set a position on the x axis (offset.x) and y axis (offset.y) for the layout
layout.1.offset.x <- 0
layout.1.offset.y <- 0

layout.2.offset.x <- -500
layout.2.offset.y <- 200

layout.3.offset.x <- 350
layout.3.offset.y <- -100

layout.4.offset.x <- -500
layout.4.offset.y <- -50 

layout.5.offset.x <- 350
layout.5.offset.y <- 200 

#set jitter to spread the nodes out
set.jitter.x=40
set.jitter.y=25

#calculate coordinates and store them in node.data.select
#group 1
subgraph.list <- subgraph.list.backup[group.1]
{
  g_joined <- disjoint_union(subgraph.list)
  layouts <- lapply(subgraph.list, layout_with_kk)
  layouts <- lapply(layouts, function(coords) cbind(jitter(coords[, 1], amount = set.jitter.x), jitter(coords[, 2], amount = set.jitter.y)))
  lay <- merge_coords(subgraph.list, layouts)
  lay[,1] <- lay[,1] + layout.1.offset.x
  lay[,2] <- lay[,2] + layout.1.offset.y
  
  #extract coordinates 
  node.data.select <- data.frame(cbind(V(g_joined), V(g_joined)$community, lay[, 1], lay[, 2]))
  node.data.select$gene <- rownames(node.data.select) 
  colnames(node.data.select) <- c("number","module","pos.x","pos.y","gene")
  node.data.select.1 <- node.data.select
}

#group 2
subgraph.list <- subgraph.list.backup[group.2]
{
  g_joined <- disjoint_union(subgraph.list)
  layouts <- lapply(subgraph.list, layout_with_kk)
  layouts <- lapply(layouts, function(coords) cbind(jitter(coords[, 1], amount = set.jitter.x), jitter(coords[, 2], amount = set.jitter.y)))
  lay <- merge_coords(subgraph.list, layouts)
  lay[,1] <- lay[,1] + layout.2.offset.x
  lay[,2] <- lay[,2] + layout.2.offset.y
  
  #extract coordinates 
  node.data.select <- data.frame(cbind(V(g_joined), V(g_joined)$community, lay[, 1], lay[, 2]))
  node.data.select$gene <- rownames(node.data.select) 
  colnames(node.data.select) <- c("number","module","pos.x","pos.y","gene")
  node.data.select.2 <- node.data.select
} 

#group 3
subgraph.list <- subgraph.list.backup[group.3]
{
  g_joined <- disjoint_union(subgraph.list)
  layouts <- lapply(subgraph.list, layout_with_kk)
  layouts <- lapply(layouts, function(coords) cbind(jitter(coords[, 1], amount = set.jitter.x), jitter(coords[, 2], amount = set.jitter.y)))
  lay <- merge_coords(subgraph.list, layouts)
  lay[,1] <- lay[,1] + layout.3.offset.x
  lay[,2] <- lay[,2] + layout.3.offset.y
  
  #extract coordinates 
  node.data.select <- data.frame(cbind(V(g_joined), V(g_joined)$community, lay[, 1], lay[, 2]))
  node.data.select$gene <- rownames(node.data.select) 
  colnames(node.data.select) <- c("number","module","pos.x","pos.y","gene")
  node.data.select.3 <- node.data.select
}

#group 4
subgraph.list <- subgraph.list.backup[group.4]
{
  g_joined <- disjoint_union(subgraph.list)
  layouts <- lapply(subgraph.list, layout_with_kk)
  layouts <- lapply(layouts, function(coords) cbind(jitter(coords[, 1], amount = set.jitter.x), jitter(coords[, 2], amount = set.jitter.y)))
  lay <- merge_coords(subgraph.list, layouts)
  lay[,1] <- lay[,1] + layout.4.offset.x
  lay[,2] <- lay[,2] + layout.4.offset.y
  
  #extract coordinates 
  node.data.select <- data.frame(cbind(V(g_joined), V(g_joined)$community, lay[, 1], lay[, 2]))
  node.data.select$gene <- rownames(node.data.select) 
  colnames(node.data.select) <- c("number","module","pos.x","pos.y","gene")
  node.data.select.4 <- node.data.select
}

#group 5
subgraph.list <- subgraph.list.backup[group.5]
{
  g_joined <- disjoint_union(subgraph.list)
  layouts <- lapply(subgraph.list, layout_with_kk)
  layouts <- lapply(layouts, function(coords) cbind(jitter(coords[, 1], amount = set.jitter.x), jitter(coords[, 2], amount = set.jitter.y)))
  lay <- merge_coords(subgraph.list, layouts)
  lay[,1] <- lay[,1] + layout.5.offset.x
  lay[,2] <- lay[,2] + layout.5.offset.y
  
  #extract coordinates 
  node.data.select <- data.frame(cbind(V(g_joined), V(g_joined)$community, lay[, 1], lay[, 2]))
  node.data.select$gene <- rownames(node.data.select) 
  colnames(node.data.select) <- c("number","module","pos.x","pos.y","gene")
  node.data.select.5 <- node.data.select
} 

#store all layouts in node.data.select
node.data.select <- rbind(node.data.select.1, node.data.select.2, node.data.select.3, node.data.select.4, node.data.select.5)

# Create a data frame of edges with coordinates for ggplot2
edge_data_for_plot <- data.frame(
  x = node.data.select$pos.x[match(edge.data$from, node.data.select$gene)],
  y = node.data.select$pos.y[match(edge.data$from, node.data.select$gene)],
  xend = node.data.select$pos.x[match(edge.data$to, node.data.select$gene)],
  yend = node.data.select$pos.y[match(edge.data$to, node.data.select$gene)]
)

# Create the network plot using ggplot2
node.data.select <- merge(node.data.select, color.key, by.x = "module", by.y = "key", all.x = TRUE)
node.data.select$pos.x <- as.numeric(node.data.select$pos.x)
node.data.select$pos.y <- as.numeric(node.data.select$pos.y)

#plot figure without edges (fast but incomplete)
p=ggplot() + 
  geom_point(data = node.data.select, aes(x = pos.x, y = pos.y, color = module), size = 2) +
  scale_color_manual(values = setNames(unique(node.data.select$bwModuleColors), unique(node.data.select$module))) +
  labs(color = 'Module') # Change legend title
p 

#plot figure without edges with labels (so you see which module is where)
library(ggrepel)
module_labels <- node.data.select %>%
  group_by(module) %>%
  summarise(pos.x = mean(pos.x), pos.y = mean(pos.y)) # Use mean position to label the module centrally

p=ggplot() + 
  geom_point(data = node.data.select, aes(x = pos.x, y = pos.y, color = module), size = 2) +
  scale_color_manual(values = setNames(unique(node.data.select$bwModuleColors), unique(node.data.select$module))) +
  labs(color = 'Module')+ # Change legend title
  geom_label_repel(data = module_labels, aes(x = pos.x, y = pos.y, label = module), 
                   size = 4, fontface = "bold", box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf)  
p #fig. 4a


#plot without label (will take a long time to plot so just save directly)
p1=ggplot()
p1=p1+
  geom_segment(mapping = aes(x = x, y = y, xend = xend, yend = yend),color = "#CCCCCC", size = 0.01, data = edge_data_for_plot)+ # draw a straight line
  geom_point(mapping = aes(x = pos.x, y = pos.y, color = module, alpha=0.7),size = 3, data = node.data.select)+ # add point
  scale_color_manual(values = setNames(unique(node.data.select$bwModuleColors), unique(node.data.select$module))) +
  scale_size(range = c(0, 6) * 2)+ # specifies the minimum and maximum size 
  theme_void()	    

# #save 
# ggsave(p1, filename = "output/plots/wgnca_big_nolabel.jpg"  , width = 8, height=6)
 
#plot with modules labeled
p2=p1+geom_label_repel(data = module_labels, aes(x = pos.x, y = pos.y, label = module),size = 4, fontface = "bold", box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf)
 
# #save 
# ggsave(p2, filename = "output/plots/wgnca_big_withlabel.jpg"  , width = 8, height=6)

p1
p2




### -------------------------------------------
### WGCNA – CARD correlations
### Clustered heatmap + FDR correction
### -------------------------------------------

library(qusage)
library(Seurat)
library(Hmisc)
library(dplyr)
library(pheatmap)

### -----------------------------
### Load data
### -----------------------------

object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data


### -----------------------------
### Compute module scores per spot
### -----------------------------

data_matrix <- object_merged@assays[["SCT"]]@data

pb_all_modules <- read.gmt("~/Desktop/final_scripts/00_data/module_genes.gmt")
names(pb_all_modules) <- trimws(names(pb_all_modules))
pb_all_modules <- lapply(pb_all_modules, trimws)

# Select modules
pb_significant_modules <- pb_all_modules[c("M4", "M14")]

module_sig <- cbind(object_merged$barcode, object_merged$orig.ident)

for (i in 1:length(pb_significant_modules)) {
  
  module_gene_vector <- pb_significant_modules[[i]]
  
  module_subset_matrix <- data_matrix[rownames(data_matrix) %in% module_gene_vector, ]
  
  module_gene_exp <- as.data.frame(colMeans(module_subset_matrix))
  
  module_sig <- cbind(module_sig, module_gene_exp)
  colnames(module_sig)[i+2] <- names(pb_significant_modules[i])
}

colnames(module_sig)[1] <- 'barcode'
colnames(module_sig)[2] <- 'orig.ident'

### -----------------------------
### Merge with metadata
### -----------------------------

full_package_metadata <- cbind(ST_ubermeta, module_sig)

### -----------------------------
### Correlation + p-values
### -----------------------------

celltype_pattern <- paste(
  c(
    "^adipocyte$",
    "^Endothelial",
    "^CAFs",
    "^PVL",
    "^B\\.cells",
    "^T_cells",
    "^Myeloid",
    "^Epithelial$",
    "^Plasmablasts$"
  ),
  collapse = "|"
)

celltype_matrix <- as.matrix(
  full_package_metadata[ , grepl(celltype_pattern, colnames(full_package_metadata)),
    drop = FALSE])

module_matrix <- as.matrix(full_package_metadata[, c("M4", "M14"), drop = FALSE])

cor_results <- rcorr(celltype_matrix, module_matrix, type = "pearson")

# Extract correlation coefficients
cor_matrix <- cor_results$r[1:ncol(celltype_matrix),
                            (ncol(celltype_matrix)+1):
                              (ncol(celltype_matrix)+ncol(module_matrix))]

# Extract p-values
p_matrix <- cor_results$P[1:ncol(celltype_matrix),
                          (ncol(celltype_matrix)+1):
                            (ncol(celltype_matrix)+ncol(module_matrix))]

# FDR correction (BH)
p_adj_matrix <- matrix(p.adjust(p_matrix, method = "BH"),
                       nrow = nrow(p_matrix),
                       ncol = ncol(p_matrix))

rownames(p_adj_matrix) <- rownames(cor_matrix)
colnames(p_adj_matrix) <- colnames(cor_matrix)

### -----------------------------
### Select Top 15 cell types
### (by maximum absolute correlation across modules)
### -----------------------------

cor_df <- as.data.frame(cor_matrix)
cor_df$celltype <- rownames(cor_df)

cor_df$max_abs_cor <- apply(abs(cor_matrix), 1, max)

top15 <- cor_df %>%
  arrange(desc(max_abs_cor)) %>%
  slice(1:15)

top_celltypes <- top15$celltype

cor_top  <- cor_matrix[top_celltypes, , drop = FALSE]
fdr_top  <- p_adj_matrix[top_celltypes, , drop = FALSE]

### -----------------------------
### Symmetric color scaling
### -----------------------------

max_abs <- max(abs(cor_top))
breaks <- seq(-max_abs, max_abs, length.out = 101)

### -----------------------------
### Clustered heatmap (design preserved)
### -----------------------------

pheatmap(cor_top,
         color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = F,
         display_numbers = round(cor_top, 2),
         fontsize_number = 9,
         border_color = NA)

## Fig. 4d - Spot-level Spearman correlations between module signature scores and deconvolved cell-type proportions. 



### -------------------------------------------
### Histology annotations vs WGCNA modules
### Clustered heatmap + FDR correction
### -------------------------------------------

library(Hmisc)
library(dplyr)
library(pheatmap)

### -----------------------------
### Extract matrices
### -----------------------------

annotation_names <- c(
  "Tumor",
  "Necrosis",
  "Fat_tissue",
  "High_TILs_stroma",
  "Cellular_stroma",
  "Acellular_stroma",
  "Vessels",
  "Canal_galactophore",
  "In_situ",
  "Nerve",
  "Lymphocyte"
)

annotation_matrix <- as.matrix(
  full_package_metadata[, annotation_names, drop = FALSE])

module_matrix <- as.matrix(full_package_metadata[, c("M4", "M14"), drop = FALSE])

colnames(annotation_matrix) <- gsub("_", " ", colnames(annotation_matrix))
colnames(annotation_matrix)[colnames(annotation_matrix) == "Canal galactophore"] <- "Normal breast"


### -----------------------------
### Correlation + p-values
### -----------------------------

cor_results <- rcorr(annotation_matrix, module_matrix, type = "pearson")

# Extract correlation coefficients
cor_matrix <- cor_results$r[1:ncol(annotation_matrix),
                            (ncol(annotation_matrix)+1):
                              (ncol(annotation_matrix)+ncol(module_matrix))]

# Extract p-values
p_matrix <- cor_results$P[1:ncol(annotation_matrix),
                          (ncol(annotation_matrix)+1):
                            (ncol(annotation_matrix)+ncol(module_matrix))]

# BH FDR correction
p_adj_matrix <- matrix(p.adjust(p_matrix, method = "BH"),
                       nrow = nrow(p_matrix),
                       ncol = ncol(p_matrix))

rownames(p_adj_matrix) <- rownames(cor_matrix)
colnames(p_adj_matrix) <- colnames(cor_matrix)

### -----------------------------
### Symmetric color scaling
### -----------------------------

max_abs <- max(abs(cor_matrix))
breaks <- seq(-max_abs, max_abs, length.out = 101)

### -----------------------------
### Clustered heatmap
### -----------------------------

pheatmap(cor_matrix,
         color = colorRampPalette(c(c("#4575b4", "white", "#d73027")))(100),
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = F,
         display_numbers = round(cor_matrix, 2),
         fontsize_number = 9,
         border_color = NA)
# save 5.5x6
## Extended Data Fig. 8f





####################################
## Fig. 4e,f - Spatial projection of module signature scores and representative 
## cell-type features in two tumors (ST131 and ST103). 
####################################

#### MODULE SIGNATURES ####

library(Seurat)
library(SeuratObject)
library(semla)
library(hdf5r)
library(stringr)
library(data.table)
library(Polychrome)
library(dplyr)
library(magrittr)
library(tibble)
library(patchwork)
library(parallel)
library(qusage)


root <- "~/Desktop/ductals_st_data/spaceranger/ST"
annot_classes <- c(
  "Tumor", "Necrosis", "Fat_tissue", "High_TILs_stroma",
  "Cellular_stroma", "Acellular_stroma", "Vessels",
  "Artefact", "Canal_galactophore", "Nodule_lymphoid",
  "In_situ", "Nerve", "Lymphocyte", "Hole",
  "Microcalcification", "Out", "Apocrine metaplasia"
)
info_colnames <- c(
  "barcode", "in_tissue", "array_row", "array_col",
  "col_pxl", "row_pxl"
)
coords_colnames <- c(
  "barcode", "in_tissue", "coord1", "coord2",
  "pxl_row_in_fullres", "pxl_col_in_fullres"
)

read_samples <- function(st_id) {
  # Load sample data paths
  samples <- paste0(root, st_id, "/outs/filtered_feature_bc_matrix.h5")
  spotfiles <- paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv")
  imgs <- paste0(root, st_id, "/outs/spatial/tissue_hires_image.png")
  json <- paste0(root, st_id, "/outs/spatial/scalefactors_json.json")
  info_table <- as.data.frame(cbind(samples, spotfiles, imgs, json))
  
  # Verify that the image file exists
  if (!file.exists(imgs)) {
    stop(paste("Image file does not exist at:", imgs))
  }
  
  # Load the Visium data
  st_sample <- ReadVisiumData(
    info_table,
    assay = "Spatial",
    min.cells = 5,
    min.features = 200
  )
  st_sample@meta.data$orig.ident <- st_id
  
  # Load annotations and coordinates
  annotations <- fread(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list_annotation.csv")
  )
  coordinates <- read.csv(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv"),
    header = FALSE
  )
  
  # Process metadata
  st_sample@meta.data$barcode <- rownames(st_sample@meta.data)
  st_spots <- rownames(st_sample@meta.data)
  
  annotations <- set_colnames(annotations, c(info_colnames, annot_classes)) %>%
    filter(barcode %in% st_spots) %>%
    dplyr::select(all_of(c(annot_classes, "barcode"))) %>%
    arrange(match(barcode, st_spots))
  
  coordinates <- set_colnames(coordinates, coords_colnames) %>%
    filter(barcode %in% st_spots) %>%
    arrange(match(barcode, st_spots))
  
  st_sample@meta.data <- Reduce(
    function(x, y) merge(x, y, by = "barcode", all = TRUE),
    list(st_sample@meta.data, annotations, coordinates)
  )
  rownames(st_sample@meta.data) <- st_sample@meta.data$barcode
  
  # Filter spots
  st_sample <- SubsetSTData(st_sample, Hole < 0.3)
  st_sample <- SubsetSTData(st_sample, Artefact < 0.3)
  st_sample <- SubsetSTData(st_sample, Out < 0.3)
  
  # Filter genes
  genes <- rownames(st_sample)
  non_meta_genes <- genes[!(grepl("RPL", genes) | grepl("RPS", genes) | grepl("MT-", genes) | grepl("MTRNR", genes))]
  st_sample <- SubsetSTData(st_sample, features = non_meta_genes)
  
  # Load images with error handling
  
  st_sample <- LoadImages(st_sample, verbose = TRUE, time.resolve = FALSE)
  
  if (is.null(st_sample)) {
    stop("st_sample is NULL after attempting to load images.")
  }
  return(st_sample)
}



############# st 131

st_id <- "131" # Replace with your actual sample ID
st_object <- read_samples(st_id)

st_object <- LoadImages(st_object, verbose = FALSE)
cols <- RColorBrewer::brewer.pal(11, "Spectral") |> rev()

columns_to_add <- c(
  "seurat_clusters",
  "adipocyte",
  "Endothelial.ACKR1",
  "Endothelial.RGS5",
  "Endothelial.CXCL12",
  "CAFs.MSC.iCAF.like.s1",
  "CAFs.MSC.iCAF.like.s2",
  "CAFs.myCAF.like.s4",
  "PVL.Immature.s1",
  "Endothelial.Lymphatic.LYVE1",
  "B.cells.Memory",
  "T_cells_c4_CD8._ZFP36",
  "T_cells_c6_IFIT1",
  "T_cells_c7_CD8._IFNG",
  "T_cells_c8_CD8._LAG3",
  "T_cells_c0_CD4._CCR7",
  "T_cells_c1_CD4._IL7R",
  "T_cells_c2_CD4._T.regs_FOXP3",
  "T_cells_c3_CD4._Tfh_CXCL13",
  "T_cells_c9_NK_cells_AREG",
  "T_cells_c10_NKT_cells_FCGR3A",
  "Myeloid_c10_Macrophage_1_EGR1",
  "Myeloid_c12_Monocyte_1_IL1B",
  "Myeloid_c1_LAM1_FABP5",
  "Myeloid_c8_Monocyte_2_S100A9",
  "Epithelial",
  "CAFs.Transitioning.s3",
  "CAFs.myCAF.like.s5",
  "PVL.Differentiated.s3",
  "PVL_Immature.s2",
  "B.cells.Naive",
  "Plasmablasts",
  "Myeloid_c2_LAM2_APOE",
  "Myeloid_c9_Macrophage_2_CXCL10",
  "Myeloid_c11_cDC2_CD1C",
  "Myeloid_c4_DCs_pDC_IRF7",
  "Myeloid_c3_cDC1_CLEC9A",
  "Myeloid_c0_DC_LAMP3",
  "T_cells_c5_CD8._GZMK",
  "Myeloid_c7_Monocyte_3_FCGR3A",
  "M4",
  "M14"
)

st131_meta <- full_package_metadata[full_package_metadata$orig.ident == '131',]

table(st131_meta$barcode == st_object@meta.data$barcode)


# Remove old versions
st_object@meta.data <- st_object@meta.data[
  , !colnames(st_object@meta.data) %in% columns_to_add,
  drop = FALSE
]

# Add updated versions
st_object@meta.data <- cbind(
  st_object@meta.data,
  st131_meta[, columns_to_add, drop = FALSE]
)


# M4 - IFN module
colnames(st_object@meta.data)

MapFeatures(st_object,
            features = "M4",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Myeloid_c9_Macrophage_2_CXCL10",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Myeloid_c1_LAM1_FABP5",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "M14",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "T_cells_c8_CD8._LAG3",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "M14",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Epithelial",
            # image_use = 'raw',
            colors = cols)



############# st 103

st_id <- "103" 
## replaced ST131 with ST103

st_object <- read_samples(st_id)

st_object <- LoadImages(st_object, verbose = FALSE)
cols <- RColorBrewer::brewer.pal(11, "Spectral") |> rev()

st103_meta <- full_package_metadata[full_package_metadata$orig.ident == '103',]

table(st103_meta$barcode == st_object@meta.data$barcode)


# Remove old versions
st_object@meta.data <- st_object@meta.data[
  , !colnames(st_object@meta.data) %in% columns_to_add,
  drop = FALSE
]

# Add updated versions
st_object@meta.data <- cbind(
  st_object@meta.data,
  st103_meta[, columns_to_add, drop = FALSE]
)


# M4 - IFN module
colnames(st_object@meta.data)

MapFeatures(st_object,
            features = "M4",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Myeloid_c9_Macrophage_2_CXCL10",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Myeloid_c1_LAM1_FABP5",
            # image_use = 'raw',
            colors = cols)


MapFeatures(st_object,
            features = "M14",
            # image_use = 'raw',
            colors = cols)

MapFeatures(st_object,
            features = "Epithelial",
            # image_use = 'raw',
            colors = cols)





# ============================================================
# Fig. 4c - Representative GO enrichment plot
# ============================================================

library(dplyr)
library(ggplot2)
library(stringr)

# ------------------------------------------------------------
# 1. Read enrichment results
# ------------------------------------------------------------

GO_enrichment_bonferroni <- read.delim(
  "~/Desktop/final_scripts/00_data/GO_enrichment_bonferroni.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. Clean enrichment results
# ------------------------------------------------------------

module_enrichment <- GO_enrichment_bonferroni %>%
  mutate(
    module = trimws(module),
    significant = toupper(as.character(significant)) == "TRUE",
    p_value = as.numeric(p_value),
    log_pval = -log10(
      pmax(p_value, .Machine$double.xmin)
    )
  ) %>%
  filter(
    module %in% c("M4", "M14"),
    significant,
    source == "GO:BP"
  )

# ------------------------------------------------------------
# 3. Select representative GO terms
# ------------------------------------------------------------

representative_terms <- c(
  
  # M4: interferon signaling
  "GO:0051607",  # defense response to virus
  "GO:0034340",  # response to type I interferon
  "GO:0140888",  # interferon-mediated signaling pathway
  "GO:0045087",  # innate immune response
  "GO:0048525",  # negative regulation of viral process
  "GO:0045071",  # negative regulation of viral genome replication
  "GO:0019885",  # antigen processing and presentation via MHC class I
  
  # M14: proliferation and cell cycle
  "GO:0007049",  # cell cycle
  "GO:0000278",  # mitotic cell cycle
  "GO:0007059",  # chromosome segregation
  "GO:0051301",  # cell division
  "GO:0000280",  # nuclear division
  "GO:0006260",  # DNA replication
  "GO:0000819"   # sister chromatid segregation
)

modules_plot <- module_enrichment %>%
  filter(term_id %in% representative_terms) %>%
  distinct(module, term_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 4. Set module order
# M14 is shown first on the x-axis, then M4
# ------------------------------------------------------------

modules_plot <- modules_plot %>%
  mutate(
    module = factor(
      module,
      levels = c("M14", "M4")
    )
  )

# ------------------------------------------------------------
# 5. Clean GO labels
# ------------------------------------------------------------

modules_plot <- modules_plot %>%
  mutate(
    term_label = term_name %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

# ------------------------------------------------------------
# 6. Define the exact biological order of terms
# Bottom-to-top order in the plot
# ------------------------------------------------------------

term_order <- c(
  
  # M14 terms: bottom section
  "Cell cycle",
  "Mitotic cell cycle",
  "Chromosome segregation",
  "Cell division",
  "Nuclear division",
  "Dna replication",
  "Sister chromatid segregation",
  
  # M4 terms: top section
  "Defense response to virus",
  "Response to type i interferon",
  "Interferon-mediated signaling pathway",
  "Innate immune response",
  "Negative regulation of viral process",
  "Negative regulation of viral genome replication",
  "Antigen processing and presentation of endogenous peptide antigen via mhc class i"
)

# Keep only terms actually found in the enrichment table
term_order <- term_order[
  term_order %in% modules_plot$term_label
]

modules_plot <- modules_plot %>%
  mutate(
    term_label = factor(
      term_label,
      levels = term_order
    )
  )

# ------------------------------------------------------------
# 7. Optional nicer display labels
# ------------------------------------------------------------

display_labels <- c(
  "Cell cycle" =
    "Cell cycle",
  
  "Mitotic cell cycle" =
    "Mitotic cell cycle",
  
  "Chromosome segregation" =
    "Chromosome segregation",
  
  "Cell division" =
    "Cell division",
  
  "Nuclear division" =
    "Nuclear division",
  
  "Dna replication" =
    "DNA replication",
  
  "Sister chromatid segregation" =
    "Sister chromatid segregation",
  
  "Defense response to virus" =
    "Defense response to virus",
  
  "Response to type i interferon" =
    "Response to type I interferon",
  
  "Interferon-mediated signaling pathway" =
    "Interferon-mediated signaling pathway",
  
  "Innate immune response" =
    "Innate immune response",
  
  "Negative regulation of viral process" =
    "Negative regulation of viral process",
  
  "Negative regulation of viral genome replication" =
    "Negative regulation of viral genome replication",
  
  "Antigen processing and presentation of endogenous peptide antigen via mhc class i" =
    "Antigen processing and presentation via MHC class I"
)

# ------------------------------------------------------------
# 8. Create plot in the style of the previous figure
# ------------------------------------------------------------

go_plot <- ggplot(
  modules_plot,
  aes(
    x = module,
    y = term_label,
    size = intersection_size,
    color = module
  )
) +
  geom_point(
    alpha = 0.9
  ) +
  
  scale_size_continuous(
    range = c(3, 9)
  ) +
  
  scale_color_manual(
    values = c(
      "M14" = "#E41A1C",
      "M4" = "#984EA3"
    )
  ) +
  
  scale_y_discrete(
    labels = display_labels,
    drop = FALSE
  ) +
  
  scale_x_discrete(
    labels = c(
      "M14" = "M14",
      "M4" = "M4"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    size = "Intersecting genes",
    color = "Module"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11,
      color = "grey30"
    ),
    
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      color = "grey30"
    ),
    
    panel.grid.major.x = element_line(
      color = "grey90"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 15,
      l = 15
    )
  )

print(go_plot)



####################################
## Fig. 4b - Association between module activity and RFS
####################################

library(qusage)
library(GSEABase)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)

counts_pseudobulk <- readRDS("~/Desktop/final_scripts/00_data/ductal_pb_rpm_log.RDS")
genesets <- read.gmt(file="/Users/bengisukarakose/Desktop/final_scripts/00_data/module_genes.gmt")

ductal_meta <- read.csv("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt", sep="")

# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_pseudobulk),   # qua ci va l'input
                    kcdf="Gaussian", geneSets = genesets)

# Assuming `normalized_counts_filtered` is your normalized counts matrix and `gene_sets` is your gene set list
gsva_res_filtered <- gsva(params)

module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == ductal_meta$name)

pseudobulk_meta <- cbind(ductal_meta, module_sig)

colnames(pseudobulk_meta) <- trimws(colnames(pseudobulk_meta))
module_matrix <- as.matrix(
  pseudobulk_meta[, paste0("M", 1:38), drop = FALSE]
)


## read forest plot (00_forestplot.R)

# forest plots
allForest(module_matrix, 
          y = Surv(time = pseudobulk_meta$time, event = pseudobulk_meta$status), fdr = T, new_page = T)


allForest(module_matrix, 
          y = Surv(time = pseudobulk_meta$time, event = pseudobulk_meta$status), fdr = T, new_page = T)
## significant ones





#############################################
## Extended Data Fig. 8g - External validation of 
## relapse-associated modules in the METABRIC cohort. 
#############################################

library(readxl)
library(qusage)
library(GSEABase)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)

counts_metabric <- read.delim("~/Desktop/bc_cox_model/2_metabric/counts_metabric.gct")
counts_metabric <- counts_metabric[,-2]
rownames(counts_metabric) <- counts_metabric$NAME
counts_metabric <- counts_metabric[,-1]

data_clinical_patient <- read.delim("~/Desktop/bc_cox_model/2_metabric/data_clinical_patient.txt", comment.char="#")
data_clinical_sample <- read.delim("~/Desktop/bc_cox_model/2_metabric/data_clinical_sample.txt", comment.char="#")

table(data_clinical_patient$PATIENT_ID == data_clinical_sample$PATIENT_ID)

metabric_meta <- cbind(data_clinical_patient, data_clinical_sample)
colnames(counts_metabric) = sub("\\.", "-", colnames(counts_metabric))
metabric_meta = metabric_meta[metabric_meta$PATIENT_ID %in% colnames(counts_metabric),]

genesets <-  read.gmt("~/Desktop/final_scripts/00_data/module_genes.gmt")
names(genesets) <- trimws(names(genesets))
genesets <- lapply(genesets, trimws)

genesets_sig <- genesets[c("M4", "M14")]


# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_metabric),   
                    kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == metabric_meta$PATIENT_ID)

metabric_meta <- cbind(metabric_meta, module_sig)

metabric_meta$OS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$OS_STATUS))
metabric_meta$RFS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$RFS_STATUS))


# forest plots
allForest(metabric_meta[,c('M4', 'M14')], 
          y = Surv(time = metabric_meta$RFS_MONTHS, event = metabric_meta$RFS_STATUS), fdr = T, new_page = T)

allForest(metabric_meta[,c('M4', 'M14')], 
          y = Surv(time = metabric_meta$OS_MONTHS, event = metabric_meta$OS_STATUS), fdr = T, new_page = T)




#############################################
## Extended Data Fig. 8h - External validation of 
## relapse-associated modules in the SCAN-B cohort. 
#############################################

library(readxl)
library(qusage)
library(GSEABase)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)

counts_scanb <- read.delim("~/Desktop/bc_cox_model/3_scanb/counts_scanb.gct")
counts_scanb <- counts_scanb[,-2]
rownames(counts_scanb) <- counts_scanb$NAME
counts_scanb <- counts_scanb[,-1]
clin_rev_scanb <- read.delim("~/Desktop/bc_cox_model/3_scanb/clin_rev_scanb.txt")
clin_all_scanb <- read_excel("~/Desktop/bc_cox_model/3_scanb/Supplementary Data Table 1 - 2023-01-13.xlsx")

# Filter clin_all_scanb to keep only rows with GEX.assay in colnames(counts_scanb)
clin_filtered <- clin_all_scanb %>%
  filter(GEX.assay %in% colnames(counts_scanb))

genesets <-  read.gmt("~/Desktop/final_scripts/00_data/module_genes.gmt")
names(genesets) <- trimws(names(genesets))
genesets <- lapply(genesets, trimws)

genesets_sig <- genesets[c("M4", "M14")]


# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_scanb),   
                    kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == clin_filtered$GEX.assay)
scanb_big_meta <- cbind(clin_filtered, module_sig)

scanb_big_meta$DRFi_days = as.numeric(scanb_big_meta$DRFi_days)
scanb_big_meta$DRFi_event = as.numeric(scanb_big_meta$DRFi_event)

scanb_big_meta$OS_days = as.numeric(scanb_big_meta$OS_days)
scanb_big_meta$OS_event = as.numeric(scanb_big_meta$OS_event)

scanb_big_meta$RFi_days = as.numeric(scanb_big_meta$RFi_days)
scanb_big_meta$RFi_event = as.numeric(scanb_big_meta$RFi_event)

scanb_big_meta$BCFi_days = as.numeric(scanb_big_meta$BCFi_days)
scanb_big_meta$BCFi_event = as.numeric(scanb_big_meta$BCFi_event)


# forest plots

allForest(scanb_big_meta[,c('M4', 'M14')],
          y = Surv(time = scanb_big_meta$DRFi_days, event = scanb_big_meta$DRFi_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c('M4', 'M14')],
          y = Surv(time = scanb_big_meta$OS_days, event = scanb_big_meta$OS_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c('M4', 'M14')],
          y = Surv(time = scanb_big_meta$RFi_days, event = scanb_big_meta$RFi_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c('M4', 'M14')],
          y = Surv(time = scanb_big_meta$BCFi_days, event = scanb_big_meta$BCFi_event), fdr = T, new_page = T)














