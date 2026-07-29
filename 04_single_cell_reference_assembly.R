####################################
## 04_single_cell_reference_assembly.R
####################################

library(Seurat)
library(harmony)
library(ggsci)
library(ggplot2)
library(patchwork)

#retrieve single cell data 

########################
# Part 1 - Wu et al.
########################

# Load the Wu2021 datasets
# Wu SZ, Nature Genetics 2021, A single-cell and spatially resolved atlas of human breast cancers


###run only first time###	  ###First we need to separate the samples for integration, because everything is in one big matrix
###run only first time###	
###run only first time###	library(Seurat)
###run only first time###		setwd("E:/VIROLOGY/scRNAseq")
###run only first time###		
###run only first time###	#data retrieved from https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-download
###run only first time###	SO <- Read10X(data.dir = "scRNAseq_data/Wu2021/BrCa_Atlas_Count_out")	
###run only first time###	
###run only first time###	meta.data <- read.table(file="scRNAseq_data/Wu2021/Whole_miniatlas_meta.csv", sep=",", header=TRUE, row.names=1)
###run only first time###	meta.data <- meta.data[-1,]
###run only first time###	
###run only first time###	PIDS <- names(table(meta.data$Patient))
###run only first time###	
###run only first time###	SO.list <- list()
###run only first time###	for(i in seq_along(PIDS)){
###run only first time###		select.cells <- rownames(meta.data)[meta.data$Patient %in% PIDS[i]]
###run only first time###		SO.list[[i]] <- SO[,colnames(SO) %in% select.cells]
###run only first time###		saveRDS(SO.list[[i]], file=paste0("scRNAseq_data/Wu2021/",PIDS[i],".rds")) #save each expression matrix as rds file
###run only first time###	}


#0. settings
set.assay = "RNA"         # "RNA" ; for now only RNA available
feat.min.cells = 3        # 3 works well;  delete features (gene) expressed in less than this number of cells (for each sample).
set.mito = 10             # set maximal mitochondrial content
min.nCount_RNA = 1500     # 1500 works well; minimum RNA counts per cell; cells with less than this are removed 
max.nCount_RNA = 99999    # 99999 works well; maximum RNA counts per cell; cells with more than this are removed 
min.nFeature_RNA = 200    # 200 generally ok; minimim amount of RNA features (i.e. genes) per cell. Will remove cells with low number of genes
n.variablefeatures = 2000 # 2000 is ok for harmony integration 

visualization.algorithm = "UMAP" # "UMAP" ; for now only "UMAP" available
visualization.dimensions = 30    # 30 works well; uses the first x dimensions for the dimensionality reduction plot
PCA.dimensions = 30              # 30 works well, used the first x dimensions for PCA reduction	

FindCluster.resolution = 2.0     # 2.0 for harmony is a good start for large data sets. 0.4 is a good start for small data sets. Increase for more clusters, decrease for less clusters

seed = 18101983 # set seed for reproducibility
set.seed(seed)

#Create settings report
timestamp <- format(Sys.time(), "%Y-%m-%d_%Hh%M" )
seurat.version <- sessionInfo()$otherPkgs$Seurat[[2]]

report.1 <- paste("Seurat v.",seurat.version,"//", timestamp)
report.1 <- rbind(report.1, paste("n.variablefeatures =",n.variablefeatures))
report.1 <- rbind(report.1, paste("set.seed =",seed,"//","set.assay =",set.assay,"//","set.mito =",set.mito,"//","feat.min.cells = ",feat.min.cells))
report.1 <- rbind(report.1, paste("min.nFeature_RNA =",min.nFeature_RNA,"//","min.nCount_RNA =",min.nCount_RNA,"//","max.nCount_RNA =",max.nCount_RNA))
report.1 <- rbind(report.1, paste("visualization.dimensions =",visualization.dimensions,"//","PCA.dimensions =",PCA.dimensions,"//","FindCluster.resolution =",FindCluster.resolution))
report.1 <- rbind(report.1, paste("visualization.algorithm = ", visualization.algorithm))


#1. load single cell data 
#Create sample name vector and data directory vector
sample.name.vector  <- c("CID3586", "CID3838", "CID3921","CID3941","CID3946","CID3948","CID3963","CID4040","CID4066","CID4067","CID4290A","CID4398","CID44041","CID4461","CID4463","CID4465","CID4471","CID4495","CID44971","CID44991","CID4513","CID4515","CID45171","CID4523","CID4530N","CID4535") #this will be orig.ident in meta.data
data.dir.vector     <- c("CID3586", "CID3838", "CID3921","CID3941","CID3946","CID3948","CID3963","CID4040","CID4066","CID4067","CID4290A","CID4398","CID44041","CID4461","CID4463","CID4465","CID4471","CID4495","CID44971","CID44991","CID4513","CID4515","CID45171","CID4523","CID4530N","CID4535") # directories where the data is to be found
root.dir.path 		<- "scRNAseq_data/Wu2021" #path to main directory that contains data 
SO.list <- list()

#create seurat object list with selected samples, process them and 
for(i in seq_along(sample.name.vector)){
  print(paste("processing sample",i,":",sample.name.vector[i]))
  #SO.list[[i]] <- Read10X(data.dir = eval(paste0(root.dir.path,"/",data.dir.vector[i])))
  #SO.list[[i]] <- Read10X(data.dir = eval(paste0(root.dir.path,"/",data.dir.vector[i])),  gene.column=1)
  SO.list[[i]] <- readRDS(eval(paste0(root.dir.path,"/",data.dir.vector[i],".rds")))
  SO.list[[i]] <- CreateSeuratObject(counts = SO.list[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
  SO.list[[i]][["percent.mt"]]    <- PercentageFeatureSet(SO.list[[i]], pattern = "^MT-")
  #SO.list[[i]] <- RenameCells(SO.list[[i]], add.cell.id = data.dir.vector[i])
  SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > min.nFeature_RNA & percent.mt < set.mito & nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)			
}

#add meta.data (you can add as many meta data as you want.
#SO.list[[1]][["my_category"]] <- "treated"     #sets sample 1
#SO.list[[2]][["my_category"]] <- "untreated"	#sets sample 2 


#integrate using Harmony 
#if one sample, run the non-integrated pipeline	
if(length(SO.list)==1){		
  SO.integrated <- SO.list[[1]]
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunUMAP(SO.integrated, dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}else{
  #if more than one sample integrate using Harmony 
  SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
  SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}

#create quality control plot
plot.qc <- list()
plot.qc$qc.1 <- print(VlnPlot(SO.integrated, "nFeature_RNA", group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.2 <- print(VlnPlot(SO.integrated, "nCount_RNA"  , group.by="orig.ident")   + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.3 <- print(VlnPlot(SO.integrated, "percent.mt"  , group.by="orig.ident")   + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.4 <- print(VlnPlot(SO.integrated, "nFeature_RNA", group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10'))
plot.qc$qc.5 <- print(VlnPlot(SO.integrated, "nCount_RNA"  , group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10'))

#create report plots
report.1 <- rbind(report.1,paste("total cells: ",ncol(SO.integrated)))
report.1 <- rbind(report.1,paste("cells per sample: ",paste(unname(table(SO.integrated$orig.ident))[match( sample.name.vector, names(table(SO.integrated$orig.ident)))], collapse=",")))
report.1 <- cbind(report.1,"\n") #prepare for plotting as ggtitle
report.1 <- paste(t(report.1),collapse =" ") #prepare for plotting as ggtitle

plot.empty <- ggplot() + theme_void()
plot.report.1 <- print(plot.empty + ggtitle(report.1) + theme(plot.title=element_text(family='', face='plain', size=8)))
plot.report.2 <- print(plot.empty)
plot.report.3 <- print(plot.empty)


#plot orig.ident and seurat clusters
plot.orig.ident <- DimPlot(SO.integrated, group.by = "orig.ident") +ggtitle ("showing orig.ident" )                        
plot.ident   <- DimPlot(SO.integrated, group.by = "seurat_clusters", label=TRUE, label.size=7) +ggtitle ("showing seurat_clusters" ) 

source("algorithms\\wrap_plots.R", echo=TRUE, max.deparse.length = 100000)

#add original meta data by Wu et al.
meta.data <- read.table(file="scRNAseq_data/Wu2021/Whole_miniatlas_meta.csv", sep=",", header=TRUE, row.names=1)
meta.data <- meta.data[-1,]

meta.data <- meta.data[rownames(meta.data) %in% colnames(SO.integrated)  ,] #select meta data of cells that are in the data set
if(!is.element(FALSE,c(rownames(meta.data) %in% colnames(SO.integrated)))){print("every cell has meta data: ok you can proceed")}
if( is.element(FALSE,c(rownames(meta.data) %in% colnames(SO.integrated)))){print("not every cell has meta data!!!")}

meta.data <- meta.data[,!colnames(meta.data) %in% c("nCount_RNA","nFeature_RNA","Percent_mito")]
meta.data <- merge(meta.data, SO.integrated@meta.data, by="row.names", all=TRUE)
rownames(meta.data) <- meta.data$Row.names
meta.data <- meta.data[, !colnames(meta.data) %in% "Row.names"] 
meta.data <- meta.data[rownames(SO.integrated@meta.data),]  #order 

SO.integrated@meta.data <- meta.data

#Join Layers 
SO.integrated <- JoinLayers(SO.integrated)		

#save SO.integrated
#saveRDS(SO.integrated, file="saveRDS/Wu_2021_seuratV5.Rds")



########################
# Part 2 - Tang et al.
########################

library(Seurat)
library(harmony)
library(ggsci)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

#0. settings
set.assay = "RNA"         # "RNA" ; for now only RNA available
feat.min.cells = 3        # 3 works well;  delete features (gene) expressed in less than this number of cells (for each sample).
set.mito = 10             # set maximal mitochondrial content
min.nCount_RNA = 1500     # 1500 works well; minimum RNA counts per cell; cells with less than this are removed 
max.nCount_RNA = 99999    # 99999 works well; maximum RNA counts per cell; cells with more than this are removed 
min.nFeature_RNA = 200    # 200 generally ok; minimim amount of RNA features (i.e. genes) per cell. Will remove cells with low number of genes
n.variablefeatures = 2000 # 2000 is ok for harmony integration 

visualization.algorithm = "UMAP" # "UMAP" ; for now only "UMAP" available
visualization.dimensions = 30    # 30 works well; uses the first x dimensions for the dimensionality reduction plot
PCA.dimensions = 30              # 30 works well, used the first x dimensions for PCA reduction	

FindCluster.resolution = 2.0     # 2.0 for harmony is a good start for large data sets. 0.4 is a good start for small data sets. Increase for more clusters, decrease for less clusters

seed = 18101983 # set seed for reproducibility
set.seed(seed)

#Create settings report
timestamp <- format(Sys.time(), "%Y-%m-%d_%Hh%M" )
seurat.version <- sessionInfo()$otherPkgs$Seurat[[2]]

report.1 <- paste("Seurat v.",seurat.version,"//", timestamp)
report.1 <- rbind(report.1, paste("n.variablefeatures =",n.variablefeatures))
report.1 <- rbind(report.1, paste("set.seed =",seed,"//","set.assay =",set.assay,"//","set.mito =",set.mito,"//","feat.min.cells = ",feat.min.cells))
report.1 <- rbind(report.1, paste("min.nFeature_RNA =",min.nFeature_RNA,"//","min.nCount_RNA =",min.nCount_RNA,"//","max.nCount_RNA =",max.nCount_RNA))
report.1 <- rbind(report.1, paste("visualization.dimensions =",visualization.dimensions,"//","PCA.dimensions =",PCA.dimensions,"//","FindCluster.resolution =",FindCluster.resolution))
report.1 <- rbind(report.1, paste("visualization.algorithm = ", visualization.algorithm))


#1. load single cell data 
#Create sample name vector and data directory vector
sample.name.vector  <- c("S53_T","S76_T","S81_T") #this will be orig.ident in meta.data
data.dir.vector     <- c("S53_T","S76_T","S81_T") # directories where the data is to be found
root.dir.path 		<- "scRNAseq_data/Tang_2023_FrontiersImmunol" #path to main directory that contains data 
SO.list <- list()

#create seurat object list with selected samples, process them and 
for(i in seq_along(sample.name.vector)){
  print(paste("processing sample",i,":",sample.name.vector[i]))
  SO.list[[i]] <- Read10X(data.dir = eval(paste0(root.dir.path,"/",data.dir.vector[i])))
  #SO.list[[i]] <- Read10X(data.dir = eval(paste0(root.dir.path,"/",data.dir.vector[i])),  gene.column=1)
  #SO.list[[i]] <- readRDS(eval(paste0(root.dir.path,"/",data.dir.vector[i],".rds")))
  SO.list[[i]] <- CreateSeuratObject(counts = SO.list[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
  SO.list[[i]][["percent.mt"]]    <- PercentageFeatureSet(SO.list[[i]], pattern = "^MT-")
  SO.list[[i]] <- RenameCells(SO.list[[i]], add.cell.id = data.dir.vector[i])
  SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > min.nFeature_RNA & percent.mt < set.mito & nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)			
}

#add meta.data (you can add as many meta data as you want.
#SO.list[[1]][["my_category"]] <- "treated"     #sets sample 1
#SO.list[[2]][["my_category"]] <- "untreated"	#sets sample 2 


#integrate using Harmony 
#if one sample, run the non-integrated pipeline	
if(length(SO.list)==1){		
  SO.integrated <- SO.list[[1]]
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunUMAP(SO.integrated, dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}else{
  #if more than one sample integrate using Harmony 
  SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
  SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}

#create quality control plot
plot.qc <- list()
plot.qc$qc.1 <- print(VlnPlot(SO.integrated, "nFeature_RNA", group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.2 <- print(VlnPlot(SO.integrated, "nCount_RNA"  , group.by="orig.ident")   + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.3 <- print(VlnPlot(SO.integrated, "percent.mt"  , group.by="orig.ident")   + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none'))
plot.qc$qc.4 <- print(VlnPlot(SO.integrated, "nFeature_RNA", group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10'))
plot.qc$qc.5 <- print(VlnPlot(SO.integrated, "nCount_RNA"  , group.by="orig.ident") + stat_summary(fun.y = median, geom='point', size = 10, colour = "red", shape = 95) + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10'))

#create report plots
report.1 <- rbind(report.1,paste("total cells: ",ncol(SO.integrated)))
report.1 <- rbind(report.1,paste("cells per sample: ",paste(unname(table(SO.integrated$orig.ident))[match( sample.name.vector, names(table(SO.integrated$orig.ident)))], collapse=",")))
report.1 <- cbind(report.1,"\n") #prepare for plotting as ggtitle
report.1 <- paste(t(report.1),collapse =" ") #prepare for plotting as ggtitle

plot.empty <- ggplot() + theme_void()
plot.report.1 <- print(plot.empty + ggtitle(report.1) + theme(plot.title=element_text(family='', face='plain', size=8)))
plot.report.2 <- print(plot.empty)
plot.report.3 <- print(plot.empty)


#plot orig.ident and seurat clusters
plot.orig.ident <- DimPlot(SO.integrated, group.by = "orig.ident") +ggtitle ("showing orig.ident" )                        
plot.ident   <- DimPlot(SO.integrated, group.by = "seurat_clusters", label=TRUE, label.size=7) +ggtitle ("showing seurat_clusters" ) 

source("algorithms\\wrap_plots.R", echo=TRUE, max.deparse.length = 100000)


#Join Layers 
SO.integrated <- JoinLayers(SO.integrated)		

#save SO.integrated
#saveRDS(SO.integrated, file="saveRDS/Tang_2023_seuratV510.Rds")
#load SO.integrated
#SO.integrated <- readRDS(file="saveRDS/Tang_2023_seuratV510.Rds")


##############################################################
# FindMarkers 
##############################################################

###4.1 settings
min.pct.expression = 0.25 #standard setting: 0.25
min.logfc = 0.25 #0.25 is standard
p.val.cutoff <- (1/10^3) #(1/10^3) is standard, use (1/10^0) to ignore
DefaultAssay(SO.integrated) <- "RNA" # "RNA" or "SCSIGN" is also possible if single cell signatures are calculated, or "ADT" for citeSeq

###4.2 run algorithm

library(future)
plan("multicore", workers = 4)

cluster.names <- unique(Idents(SO.integrated))[order(unique(Idents(SO.integrated)))]
markers.fm.list <- list()
for (i in 1:length(cluster.names)) {
  print(paste0("calculating markers for cluster ",cluster.names[i],". Total: ",length(cluster.names)," clusters"))
  markers.fm.list[[i]] <- FindMarkers(SO.integrated, ident.1 = cluster.names[i], min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
  markers.fm.list[[i]] <- markers.fm.list[[i]][which(markers.fm.list[[i]]$p_val_adj < (p.val.cutoff) ),] #select genes with p_val_adj > p.val.cutoff setting
}

#markers.fm.list <- markers.fm.list[-which(sapply( markers.fm.list, function(x) length(rownames(x))) == 0 ) ]  #run this line if an element of the marker list is zero
for (i in 1:length(markers.fm.list)){
  print(DotPlot(SO.integrated, cols = "RdYlBu",features = rev( rownames(markers.fm.list[[i]])[1:40]) ) + RotatedAxis() )
  print(paste0(cluster.names[i]))
  #readline(prompt="Press [enter] to continue")
}	

library(RColorBrewer)
FeaturePlot(SO.integrated, features = feat) &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))



#############################################################################################
# Manually annotate clusters
# 
# Here you can computationally isolate clusters, explore their identity and  
# 	manually save the annotations, and then inject them back in the original seurat object
#
#############################################################################################

#set working directory
#setwd("C:/my_dir")

#1. set annotations 
annotations.list <- list()
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"       #celltype 0     
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"     
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =4 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =5 )  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"  #celltype 5         
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =6 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"     
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =7 )  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =8 )  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =10)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"   #celltype 10
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =9 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =11)  ;names(annotations.list)[length(annotations.list)] <-   "ILC3"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =12)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =13)  ;names(annotations.list)[length(annotations.list)] <-   "neutrophil"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =14)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =15)  ;names(annotations.list)[length(annotations.list)] <-   "B cell"  #celltype 15               
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =16)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =17)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =18)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =19)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =20)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"      #celltype 20                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =21)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =22)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =23)  ;names(annotations.list)[length(annotations.list)] <-   "plasma cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =24)  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =25)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"    #celltype 25                    
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =26)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =27)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =28)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =29)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =30)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"  #celltype 30               
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =31)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =32)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =33)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =34)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =35)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"               #celltype 35  
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =36)  ;names(annotations.list)[length(annotations.list)] <-   "mast cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =37)  ;names(annotations.list)[length(annotations.list)] <-   "neutrophil"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =38)  ;names(annotations.list)[length(annotations.list)] <-   "stromal cell"                 
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =39)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =40)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =41)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 

#split cluster 34 and identify adipocytes
SO.34 <- subset(SO.integrated, idents=34)
SO.34 <- FindClusters(SO.34, resolution = 0.2) #hoe hoger res, hoe meer clusters
DimPlot(SO.34, label=TRUE, label.size=8)

#add annotations
annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.34, idents =1)  ;names(annotations.list)[length(annotations.list)] <-   "adipocyte"                 


#2. Inject new annotations in Seurat object
for(i in seq_along(annotations.list)){
  Idents(SO.integrated, cells = annotations.list[[i]])  <- names(annotations.list)[[i]]
}

#3. save in meta.data$annotation
SO.integrated[["annotation"]] <- Idents(SO.integrated)
SO.integrated[["celltype_major"]] <- Idents(SO.integrated)
SO.integrated[["celltype_minor"]] <- Idents(SO.integrated)
SO.integrated[["celltype_subset"]] <- Idents(SO.integrated)
DimPlot(SO.integrated, label=T, label.size=5)

#save 
#saveRDS(SO.integrated, file="saveRDS/Tang_2023_seuratV510_annotated_adipocyte.Rds")      #save
#load
#SO.integrated <- readRDS(file="saveRDS/Tang_2023_seuratV510_annotated_adipocyte.Rds")

#markers of adipocytes
cluster.names <- unique(Idents(SO.integrated))[order(unique(Idents(SO.integrated)))]
i=1
min.pct.expression = 0.25 #standard setting: 0.25
min.logfc = 0.25 #0.25 is standard
p.val.cutoff <- (1/10^3) #(1/10^3) is standard, use (1/10^0) to ignore
DefaultAssay(SO.integrated) <- "RNA" # "RNA" or "SCSIGN" is also possible if single cell signatures are calculated, or "ADT" for citeSeq
print(paste0("calculating markers for cluster ",cluster.names[i]))
markers.adipo <- FindMarkers(SO.integrated, ident.1 = cluster.names[i], min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)

#plot 
print(DotPlot(SO.integrated, cols = "RdYlBu",features = rev( rownames(markers.adipo)[1:25]) ) + RotatedAxis() )

p1 = FeaturePlot(SO.integrated, features = "PLIN1") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p2 = FeaturePlot(SO.integrated, features = "PLIN4") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p3 = FeaturePlot(SO.integrated, features = "ADIPOQ") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p4 = FeaturePlot(SO.integrated, features = "CIDEC") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

p = wrap_plots(list(p1,p2,p3,p4), ncol=2)
p

p.CPA3 = FeaturePlot(SO.integrated, features = "CPA3") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p.KIT = FeaturePlot(SO.integrated, features = "KIT") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

p.TPSAB1 = FeaturePlot(SO.integrated, features = "TPSAB1") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p.TPSB2 = FeaturePlot(SO.integrated, features = "TPSB2") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

#markers of 34
min.pct.expression = 0.25 #standard setting: 0.25
min.logfc = 0.25 #0.25 is standard
p.val.cutoff <- (1/10^3) #(1/10^3) is standard, use (1/10^0) to ignore
DefaultAssay(SO.integrated) <- "RNA" # "RNA" or "SCSIGN" is also possible if single cell signatures are calculated, or "ADT" for citeSeq
markers.34 <- FindMarkers(SO.integrated, ident.1 = '34', min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)

#plot 
print(DotPlot(SO.integrated, cols = "RdYlBu",features = rev( rownames(markers.34)[1:40]) ) + RotatedAxis() )
p.CPA3 = FeaturePlot(SO.integrated, features = "CPA3") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p.KIT = FeaturePlot(SO.integrated, features = "KIT") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

p.TPSAB1 = FeaturePlot(SO.integrated, features = "TPSAB1") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p.TPSB2 = FeaturePlot(SO.integrated, features = "TPSB2") &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))


#save raw read count matrices of adipocytes per patient as rds file
SO.adipo <- subset(SO.integrated, idents="adipocyte")
SO.adipo.S53_T <- subset(SO.adipo, orig.ident == 'S53_T')
SO.adipo.S76_T <- subset(SO.adipo, orig.ident == 'S76_T')
SO.adipo.S81_T <- subset(SO.adipo, orig.ident == 'S81_T')

S53_T.adipo.rawcounts <- GetAssayData(SO.adipo.S53_T, layer='counts' )		
S76_T.adipo.rawcounts <- GetAssayData(SO.adipo.S76_T, layer='counts' )		
S81_T.adipo.rawcounts <- GetAssayData(SO.adipo.S81_T, layer='counts' )		

#save 
#saveRDS(S53_T.adipo.rawcounts, file="scRNAseq_data/Tang_2023_FrontiersImmunol/adipocyte_rawreadcounts_RDS/S53_T_adipo_rawcounts.Rds")      #save
#saveRDS(S76_T.adipo.rawcounts, file="scRNAseq_data/Tang_2023_FrontiersImmunol/adipocyte_rawreadcounts_RDS/S76_T_adipo_rawcounts.Rds")      #save
#saveRDS(S81_T.adipo.rawcounts, file="scRNAseq_data/Tang_2023_FrontiersImmunol/adipocyte_rawreadcounts_RDS/S81_T_adipo_rawcounts.Rds")      #save


#save adipocyte seurat object
#saveRDS(SO.adipo, file="saveRDS/Tang_2023_isolated_adipocytes_seuratV510.Rds")
#load adipocyte seurat object
#SO.adipo <- readRDS(file="saveRDS/Tang_2023_isolated_adipocytes_seuratV510.Rds")



########################
# Part 3 - Integration
########################

############################################################################################################################################
#integrate adipocytes from Tang 2023 Frontiers Immunol with Wu 2021 Nature Genetics
#for pipeline of adipocytes: E:/VIROLOGY/scRNAseq/scRNAseq_data/Tang_2023_FrontiersImmunol/scRNAseq_Tang_2023_FrontiersImmunol_SeuratV5.R
############################################################################################################################################

library(Seurat)
library(harmony)
library(ggsci)
library(ggplot2)
library(patchwork)


#load data sets
SO.tang <- readRDS(file="saveRDS/Tang_2023_isolated_adipocytes_seuratV510.Rds")
SO.wu   <- readRDS(file="saveRDS/Wu_2021_seuratV5.Rds")

#Trim each sample to raw RNA counts only
tang.sample.vector <- as.character(unique(SO.tang@meta.data$orig.ident))
n.tang <- length(tang.sample.vector)
SO.list <- list()

for (i in seq_along(tang.sample.vector)) {
  print(paste("trimming sample",i))
  SO.list[[i]] <- subset(SO.tang, subset = orig.ident == tang.sample.vector[i]) 
  DefaultAssay(SO.list[[i]]) <- "RNA"
  SO.list[[i]] <- DietSeurat(SO.list[[i]] , assays = "RNA")
}

wu.sample.vector <- as.character(unique(SO.wu@meta.data$orig.ident))
for (i in seq_along(wu.sample.vector)) {
  print(paste("trimming sample",i+n.tang))
  SO.list[[i+n.tang]] <- subset(SO.wu, subset = orig.ident == wu.sample.vector[i]) 
  DefaultAssay(SO.list[[i+n.tang]]) <- "RNA"
  SO.list[[i+n.tang]] <- DietSeurat(SO.list[[i+n.tang]] , assays = "RNA")
}  


#remove empty elements from the list, and count cells
SO.list <- SO.list[!sapply(SO.list, is.null)]
cells.per.sample <- paste(unlist(sapply(SO.list,ncol)),collapse=",")

#re-integrate and process
#if one sample, run the non-integrated pipeline	
n.variablefeatures = 2000
visualization.dimensions = 30
PCA.dimensions = 30

if(length(SO.list)==1){		
  SO.integrated <- SO.list[[1]]
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunUMAP(SO.integrated, dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}else{
  #if more than one sample integrate using Harmony 
  SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
  rm(SO.list)
  SO.integrated <- NormalizeData(SO.integrated)
  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
  SO.integrated <- ScaleData(SO.integrated)
  SO.integrated <- RunPCA(SO.integrated)
  SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
  SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:visualization.dimensions)
  SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:PCA.dimensions) 
  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
}


#1. create plots
plot.orig.ident <- DimPlot(SO.integrated, group.by = "orig.ident") + ggtitle(paste("harmony-reintegrated"))
plot.ident   <- DimPlot(SO.integrated,  label=TRUE, label.size=7)  + ggtitle(paste0("reintegrated, del clust ",paste0(delete.clusters, collapse = ",")))   
source("algorithms\\wrap_plots.R", echo=TRUE, max.deparse.length = 100000)

#join layers
SO.integrated <- JoinLayers(SO.integrated)		


#save
#saveRDS(SO.integrated, file="saveRDS/Wu_Tang_isolated_adipocytes_seuratV510.Rds")
#SO.integrated <- readRDS(file="saveRDS/Wu_Tang_isolated_adipocytes_seuratV510.Rds")


