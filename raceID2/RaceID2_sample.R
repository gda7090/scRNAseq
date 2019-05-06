## install required packages (only at first time)
#install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit","vegan"))

## load class definition and functions
source("RaceID2_StemID_class.R")
library(argparse)

parser = ArgumentParser()
parser$add_argument("--table", help="input matrix, with barcodes as colnames and genes as rownames",required=TRUE)
parser$add_argument("--multisample",help='if libraries with different complexities are combined and jointly analyzed, or if the analyzed sample comprises cells with highly variable total transcript count (several orders of magnitude), it is strongly advised to eliminate cell-to-cell differences in technical noise by using downsampling instead of normalization in the filterdata method by setting the argument downsample to TRUE.', required=TRUE,default="TRUE")
parser$add_argument("--prefix", help="prefix of output files",required=TRUE)
parser$add_argument("--outdir",help='output folder',required=TRUE,default="output")
args <- parser$parse_args()
str(args)
table=args$table
multisample=args$multisample
prefix=args$prefix
outdir=args$outdir
dir.create(outdir)

## input data
x <- read.csv(table,sep="\t",header=TRUE)
rownames(x) <- x$GENEID
save(list =ls(all=TRUE), file=paste(outdir,"/",prefix,".RData",sep=""))
# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]

## RaceID2
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=TRUE, dsn=1, rseed=17000)
# k-medoids clustering
sc <- clustexp(sc,clustnr=30,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids")
# compute t-SNE map
sc <- comptsne(sc,rseed=15555)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)

## diagnostic plots
# gap statistics: only if do.gap == TRUE
##plotgap(sc)
# plot within-cluster dispersion as a function of the cluster number: only if sat == TRUE
pdf(file=paste(outdir,'/',prefix,"plotsaturation_disp.pdf",sep=""))
plotsaturation(sc,disp=TRUE)
dev.off()
# plot change of the within-cluster dispersion as a function of the cluster number: only if sat == TRUE
pdf(file=paste(outdir,'/',prefix,"plotsaturation.pdf",sep=""))
plotsaturation(sc)
dev.off()
# silhouette of k-medoids clusters
pdf(file=paste(outdir,'/',prefix,"plotsilhouette.pdf",sep=""))
plotsilhouette(sc)
dev.off()
# Jaccard's similarity of k-medoids clusters
pdf(file=paste(outdir,'/',prefix,"plotjaccard.pdf",sep=""))
plotjaccard(sc)
dev.off()
# barchart of outlier probabilities
pdf(file=paste(outdir,'/',prefix,"plotoutlierprobs.pdf",sep=""))
plotoutlierprobs(sc)
dev.off()
# regression of background model
pdf(file=paste(outdir,'/',prefix,"plotbackground.pdf",sep=""))
plotbackground(sc)
dev.off()
# dependence of outlier number on probability threshold (probthr)
pdf(file=paste(outdir,'/',prefix,"plotsensitivity.pdf",sep=""))
plotsensitivity(sc)
dev.off()
# heatmap of k-medoids cluster
pdf(file=paste(outdir,'/',prefix,"clustheatmap.pdf",sep=""))
clustheatmap(sc,final=FALSE,hmethod="single")
dev.off()
# heatmap of final cluster
pdf(file=paste(outdir,'/',prefix,"clustheatmap.pdf",sep=""))
clustheatmap(sc,final=TRUE,hmethod="single")
dev.off()
# highlight k-medoids clusters in t-SNE map
pdf(file=paste(outdir,'/',prefix,"plottsne.pdf",sep=""))
plottsne(sc,final=FALSE)
dev.off()
# highlight final clusters in t-SNE map
pdf(file=paste(outdir,'/',prefix,"plottsne_final.pdf",sep=""))
plottsne(sc,final=TRUE)
dev.off()
# highlight cell labels in t-SNE map
pdf(file=paste(outdir,'/',prefix,"plotlabelstsne.pdf",sep=""))
plotlabelstsne(sc,labels=sub("(\\_\\d+)","",names(sc@ndata)))
dev.off()
# highlight groups of cells by symbols in t-SNE map
pdf(file=paste(outdir,'/',prefix,"plotsymbolstsne.pdf",sep=""))
plotsymbolstsne(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))
dev.off()
# highlight transcirpt counts of a set of genes in t-SNE map, e. g. all Apoa genes
pdf(file=paste(outdir,'/',prefix,"plotexptsne.pdf",sep=""))
g <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", "Apoa5__chr9")
plotexptsne(sc,g,n="Apoa genes",logsc=TRUE)
dev.off()
## identification of marker genes
# differentially regulated genes in each cluster compared to the full ensemble
cdiff <- clustdiffgenes(sc,pvalue=.01)

## write results to text files
# final clusters 
x <- data.frame(CELLID=names(sc@cpart),cluster=sc@cpart)
write.table(x[order(x$cluster,decreasing=FALSE),],paste(outdir,'/',prefix,"cell_clust.xls",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
# differentially expressed genes in cluster
for ( n in names(cdiff) ) write.table(data.frame(GENEID=rownames(cdiff[[n]]),cdiff[[n]]),paste(outdir,'/',prefix,paste("cell_clust_diff_genes",sub("\\.","\\_",n),sep="_"),".xls",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

# differentially expressed genes between two sets of clusters, e. g. cluster 1 and clusters 2,3
d <- diffgenes(sc,cl1=1,cl2=c(2,3),mincount=5)
pdf(file=paste(outdir,'/',prefix,"plotdiffgenes.pdf",sep=""))
plotdiffgenes(d,gene=names(d$z)[1])
dev.off()


