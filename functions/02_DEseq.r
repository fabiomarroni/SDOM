# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--abundance"), type="character", default="",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-C", "--configfile"), type="character", default="", 
              help="Config file [default= %default]", metavar="character"),
  make_option(c("-c", "--condition"), type="character", default="Status", 
              help="Name of column of config file representing the sample condition [default= %default]", metavar="character"),
  make_option(c("-S", "--maxspecies"), type="numeric", default=50, 
              help="Maximum number of species to plot [default= %default]", metavar="character"),
  make_option(c("-r", "--minreads"), type="numeric", default=100,
              help="Minimum number of reads for a species to be retained", metavar="character"),
  make_option(c("-R", "--removeme"), type="character", default="", 
              help="Comma delimited list of samples to remove [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="", 
              help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$abundance)) {
  stop("WARNING: No abundance specified with '-I' flag.")
} else {  cat ("abundance is ", opt$abundance, "\n")
  abundance <- opt$abundance  
  }

if (is.null(opt$configfile)) {
  stop("WARNING: No configfile specified with '-C' flag.")
} else {  cat ("configfile is ", opt$configfile, "\n")
  configfile <- opt$configfile  
  }

if (is.null(opt$condition)) {
  stop("WARNING: No condition specified with '-c' flag.")
} else {  cat ("condition is ", opt$condition, "\n")
  condition <- opt$condition  
  }

if (is.null(opt$maxspecies)) {
  stop("WARNING: No maxspecies specified with '-S' flag.")
} else {  cat ("maxspecies is ", opt$maxspecies, "\n")
  maxspecies <- opt$maxspecies  
  }

if (is.null(opt$outdir)) {
  stop("WARNING: No outdir specified with '-S' flag.")
} else {  cat ("outdir is ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

if (is.null(opt$minreads)) {
  stop("WARNING: No minreads parameter specified with '-r' flag.")
} else {  cat ("minreads is ", opt$minreads, "\n")
  minreads <- opt$minreads  
  }

if (is.null(opt$removeme)) {
  stop("WARNING: No outfile specified with '-T' flag.")
} else {  cat ("removeme is ", opt$removeme, "\n")
  removeme <- opt$removeme  
  }


runDEseq<-function() 
{
    library(DESeq2)
    library(data.table)
	library(ape)
	library("RColorBrewer")
	library("pheatmap")
	library("ggplot2")
	library("ggrepel")
    library(ComplexHeatmap)
    toremove<-unlist(strsplit(removeme,","))
	#Read config file
	configdata<-fread(configfile,header=T,data.table=F)
	configdata<-configdata[!configdata[,1]%in%toremove,]
    sampleCondition<-configdata[,colnames(configdata)%in%condition]
	sampleNames<-configdata[,1]
    #Tranlsate from Italian to English for publication
    sampleCondition[sampleCondition%in%"SANO"]<-"Healthy T0"
    sampleCondition[sampleCondition%in%"MALATO"]<-"Moribund T0"
    sampleCondition[sampleCondition%in%"POST"]<-"T1"
	species<-fread(abundance,data.table=F,header=T)
	#Patch to sum genera with same name (but different taxon ID): this is basically used only for unculutured
	if(length(unique(species$name))<nrow(species))
	{
	species<-aggregate(species[,2:ncol(species)],by=list(species$name),FUN="sum")
	names(species)[1]<-"name"
	}


	row.names(species)<-species$name
	species$name<-species$taxonomy_id<-species$Total<-NULL
    names(species)<-gsub("raw_","",names(species))
	species<-species[,!names(species)%in%toremove,]
	species$Total<-rowSums(species)
	species<-species[species$Total>=minreads,]
	species<-species[order(species$Total,decreasing=TRUE),]
	species$Total<-NULL
	#Create a metadata on the fly
	metadata<-data.frame(Sample=sampleNames,condition=sampleCondition)
	metadata$Sample<-factor(metadata$Sample,levels=metadata$Sample)
	row.names(metadata)<-metadata$Sample
	#I change name to easily use a recycled piece of code
	countdata<-species
	if(sum(rownames(metadata)==names(countdata))!=ncol(countdata)) stop("Names in countdata do not match names in metadata")
	ddsHTSeq <- DESeqDataSetFromMatrix(countData  = countdata,
                                        colData = metadata,
                                        design= ~ condition)

 	vsd <- varianceStabilizingTransformation(ddsHTSeq, fitType="local")
	#This is my horrible customization of the DESeq2 PCA function to manually change shapes
	intgroup<-"condition"
	rv <- rowVars(assay(vsd))
    #In order to keep all the results similar to those produce by DESeq2, I fixed the max number of rows to 500, the default of plotPCA
    select <- order(rv, decreasing = TRUE)[seq_len(min(500,
        length(rv)))]
    pca <- prcomp(t(assay(vsd)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(vsd)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(vsd)[, intgroup,
        drop = FALSE])
    myshapes<-colData(vsd)[, intgroup]
	myshapes<-unique(myshapes)
	myshapes<-seq(0,length(myshapes)-1)
	# myshapes<-letters[seq(1,length(myshapes))]
	# myshapes<-unlist(lapply(as.character(myshapes),utf8ToInt))
	manualcol=c("dodgerblue1","firebrick1","dodgerblue4","firebrick4")
    Condition <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(vsd)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = Condition,
        intgroup.df, name = colnames(vsd))
    
#If confidence we use the stat_ellipse method which computes the 95% confidence ellipse.
#If confidence==FALSE we use geom_mark_ellipse to compute the smallest ellipse encompassing all our points (because with 3 points we cannot compute 95% ellipse)

confidence<-TRUE
if(confidence)
{
pdf(paste0(outdir,"PCA_custom_conf.pdf"))
    tt<-ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "Condition")) +
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
        100), "% variance")) + 
         stat_ellipse(aes(x=PC1,y=PC2,color=Condition, fill=Condition),
                geom="polygon", level=0.95, alpha=0.2) +
            scale_color_manual(values=manualcol) +
            scale_fill_manual(values=manualcol) +
		coord_fixed()
	print(tt)
	dev.off()
}else{
pdf(paste0(outdir,"PCA_custom.pdf"))
    tt<-ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "Condition")) +
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
        100), "% variance")) + 
         ggforce::geom_mark_ellipse(aes(x=PC1,y=PC2,color=Condition, fill=Condition),
                n=1000,tol=0.1, alpha=0.2) +
            scale_color_manual(values=manualcol) +
            scale_fill_manual(values=manualcol) +
		coord_fixed()
	print(tt)
	dev.off()
}
	pdf(paste0(outdir,"PCA.pdf"))
    myplot<-plotPCA(vsd,intgroup="condition")
    print(myplot)
	dev.off()
  somma=rowSums(assay(vsd))
  mapvsd=assay(vsd)
  mapvsd=mapvsd[order(somma, decreasing = T),][1:min(nrow(mapvsd),maxspecies),]
  #Some tricks to have colors of classes in heatmap similar to those in plotPCA. We will need to fix them in a better way
  annotation<-data.frame(metadata$condition)
  names(annotation)[1]<-"Condition"
#  annotation[,1]<-factor(annotation[,1],levels=sort(levels(annotation[,1]),decreasing=T))
  row.names(annotation)<-row.names(metadata)
#  pheatmap(mapvsd, cluster_rows=FALSE, border_color=NA, annotation_col= annotation, fontsize=4, cellwidth=6, cellheight=4, filename=paste0(outdir,"heatmap.pdf"))
  pheatmap(mapvsd, cluster_rows=FALSE, border_color=NA, annotation_col= annotation, fontsize=6, cellwidth=8, cellheight=6, filename=paste0(outdir,"heatmap.pdf"))
	ddsHTSeq<-DESeq(ddsHTSeq)

	write.table(counts(ddsHTSeq,normalize=F),paste0(outdir,"raw_counts.txt"),quote=F,sep="\t",col.names=NA)
	write.table(counts(ddsHTSeq,normalize=T),paste0(outdir,"norm_counts.txt"),quote=F,sep="\t",col.names=NA)
	to.contrast<-unique(sampleCondition)
	compare<-combn(to.contrast,2)
	for(aaa in 1:ncol(compare))
    {
 		res <- results(ddsHTSeq,contrast=c("condition",compare[2,aaa],compare[1,aaa]))
		res<-res[order(res$padj),]
        #We also write expression data of the condition
		getnames<-row.names(colData(ddsHTSeq))[colData(ddsHTSeq)$condition%in%c(compare[2,aaa],compare[1,aaa])]
		getexpr<-counts(ddsHTSeq,normalize=F)[,getnames]
		res<-merge(res,getexpr,by="row.names",sort=F)
		resfile<-paste0(outdir,"DA_",compare[2,aaa],"_vs_",compare[1,aaa],".txt")
		write.table(res,resfile,sep="\t",row.names=F,quote=F)
		mysig<-sum(res$padj<=0.05,na.rm=T)
		cat(compare[2,aaa],"vs",compare[1,aaa],",",mysig,"significant species out of",nrow(res),"\n")
	}



}
runDEseq()
