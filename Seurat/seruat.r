library(Seurat)
library(dplyr)
library(argparse)
parser = ArgumentParser()
parser$add_argument("--matrices", help="/path/to/cellranger/output/filtered_gene_bc_matrices")
parser$add_argument("--gene_bar", help="the gene_bar.csv produced by cellranger")
parser$add_argument("--sample", help="the sample of project")

args <- parser$parse_args()
str(args)


matrices = args$matrices
gene_bar= args$gene_bar

if (!is.null(matrices) & is.null(gene_bar)){
	matrices = args$matrices
	pbmc.data <- Read10X(data.dir = matrices) 
}else if(is.null(matrices) & !is.null(gene_bar)){
	pbmc.data <- read.csv(gene_bar,header=T,row.names=1,sep=',')
}

#else if{
#	print "You should provide matrices or gene_bar.csv" 
#}
	
sample=args$sample


#if (!is.na(matrices) && is.na(gene_bar)){
#pbmc.data <- Read10X(data.dir = matrices)
#elif(!is.na(gene_bar) && is.na(matrix)){
#elif(is.na(gene_bar) && is.na(matrices)){
#	print "You should provide gene_bar.csv or martix"
#	}

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 400,project = sample)
umi.data<-cbind(rownames(pbmc@meta.data),pbmc@meta.data)
colnames(umi.data)[1]<-"Barcode"
write.table(umi.data[1:3],file=paste(sample,'_gene_umi.csv',sep=''),col.names=T,row.names=F,sep=',',quote=F)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pdf(paste(sample,'_Gene_UMI_mito.pdf',sep=''))

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
    low.thresholds = c(400, -Inf), high.thresholds = c(2500, 0.05))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
    scale.factor = 10000)

pdf(paste(sample,'_Dispersion',sep=''))

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

dev.off()

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 5)

pca.data<-Seurat::GetCellEmbeddings(object = pbmc, reduction.type = "pca", dims.use = 1:10)
pca.data<-cbind(rownames(pca.data),pca.data)
colnames(pca.data)[1]<-"Barcode"
write.table(pca.data,file=paste(sample,'_pca.csv',sep=''),col.names=T,row.names=F,sep=',',quote=F)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

##cluster resulat
write.table(pbmc@ident,file=paste(sample,'_seruat_cluster.csv',sep=""),row.names=T,col.names=c('Barcode,cluster'),quote=F,sep=',')

##TSNE 
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
tsne.data<-pbmc@dr$tsne@cell.embeddings
tsne.data<-cbind(rownames(tsne.data),tsne.data)
colnames(tsne.data)[1]<-"Barcode"
write.table(tsne.data,file=paste(sample,'_TSNE.csv',sep=""),row.names=F,col.names=T,sep=',',quote=F)
pdf(paste(sample,'_tsne.pdf',sep=''))

TSNEPlot(object = pbmc)
dev.off()

###cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
####write.table(cluster1.markers,file='cluster_1.csv',quote=F,row.names=T,col.names=T,sep='\t')

###differ Analysis
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

write.table(pbmc.markers,file=paste(sample,'_cluster_diff.csv',sep=''),quote=F,row.names=T,col.names=T,sep='\t')

##top genes 
top_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

##violin
for (each in top_markers$gene){
	p<-VlnPlot(pbmc,each,point.size.use=NA,do.return = TRUE)
	p<-p+xlab("Cluster")+ylab("log(UMI)")
	ggsave(paste(each,'_violin.pdf',sep=''), p, width = 6, height =8)
	}	

##gene tsne
for (each in top_markers$gene){
	P<-FeaturePlot(object = pbmc,each,cols.use = c("grey", "blue"), reduction.use = "tsne")
	ggsave(paste(each,'_tsne.pdf',sep=''), P, width = 6, height = 8)
	}

#pdf('gene_tsne.pdf')
#FeaturePlot(object = pbmc,genes,cols.use = c("grey", "blue"), reduction.use = "tsne")
#dev.off()

##Heatmap
##pdf(paste(sample,'Heatmap.pdf',sep=''))
p<-DoHeatmap(object = pbmc, genes.use = top_markers$gene, slim.col.label = TRUE, remove.key = TRUE)
ggsave(paste(sample,'Heatmap.pdf',sep=''),width = 15, height=10)
dev.off()

saveRDS(pbmc, file = paste(sample,"_final.rds",sep=''))
