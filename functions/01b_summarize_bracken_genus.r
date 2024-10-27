# Run with --help flag for help.
# Modified 10/07/2022 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",
              help="Input directory containing kraken reports", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-N", "--nspecies"), type="numeric", default=15, 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-R", "--removeme"), type="character", default="", 
              help="Comma delimited list of samples to remove [default= %default]", metavar="character"),
  make_option(c("-G", "--outgraph"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input folder specified with '-I' flag.")
} else {  cat ("Indir is ", opt$indir, "\n")
  indir <- opt$indir  
  }

if (is.null(opt$nspecies)) {
  stop("WARNING: No nspecies specified with '-N' flag.")
} else {  cat ("nspecies is ", opt$nspecies, "\n")
  nspecies <- opt$nspecies  
  }

if (is.null(opt$removeme)) {
  stop("WARNING: No outfile specified with '-T' flag.")
} else {  cat ("removeme is ", opt$removeme, "\n")
  removeme <- opt$removeme  
  }

if (is.null(opt$outgraph)) {
  stop("WARNING: No outgraph specified with '-t' flag.")
} else {  cat ("outgraph is ", opt$outgraph, "\n")
  outgraph <- opt$outgraph  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No output file specified with '-I' flag.")
} else {  cat ("Output file is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_bracken<-function()
{
library("data.table")
library(ggplot2)
setwd(indir)
fullfiles<-dir(pattern="_G.bracken.txt")
krakfiles<-dir(pattern="kraken.report.txt")
outgraph<-gsub("genus_",paste0("genus_",nspecies,"_"),outgraph)
for(aaa in 1:length(fullfiles))
{
	myunk<-fread(krakfiles[aaa],data.table=F)
	mykperc<-myunk$V1[myunk$V6=="root"]/100
	mydata<-fread(fullfiles[aaa],data.table=F)
	mydata$perc<-100*mykperc*mydata$fraction_total_reads
	myraw<-mydata[,c("name","taxonomy_id","new_est_reads")]
	mydata<-mydata[,c("name","taxonomy_id","perc")]
	setnames(mydata,"perc",paste("perc",gsub("_G.bracken.txt","",fullfiles[aaa]),sep="_"))
	setnames(myraw,"new_est_reads",paste("raw",gsub("_G.bracken.txt","",fullfiles[aaa]),sep="_"))
	if(aaa==1) 
	{
	findata<-mydata
	finrawdata<-myraw
	}
	if(aaa>1) 
	{
	findata<-merge(findata,mydata,by=c("name","taxonomy_id"),all=TRUE)
	finrawdata<-merge(finrawdata,myraw,by=c("name","taxonomy_id"),all=TRUE)
	}
}
findata[is.na(findata)]<-0
findata$Total<-rowSums(findata[,3:ncol(findata)])
findata<-findata[order(findata$Total,decreasing=T),]
finrawdata[is.na(finrawdata)]<-0
finrawdata$Total<-rowSums(finrawdata[,3:ncol(finrawdata)])
finrawdata<-finrawdata[order(finrawdata$Total,decreasing=T),]

#findata<-findata[findata$Total>40,]
write.table(findata,outfile,sep="\t",quote=F,row.names=F)
write.table(finrawdata,gsub("bracken.txt","bracken_raw.txt",outfile),sep="\t",quote=F,row.names=F)
findata$Total<-findata$taxonomy_id<-NULL
#Crazy things with strings to have short names
sample.names<-apply(sapply(strsplit(names(findata), "-"), "[", 3:4),2,paste,collapse="_")
sample.names<-apply(sapply(strsplit(sample.names,"_"),"[",1:2),2,paste,collapse="_")
sample.names<-sample.names[-1]
names(findata)[2:ncol(findata)]<-sample.names

maxperc<-100
toremove<-unlist(strsplit(removeme,","))
findata<-findata[,!names(findata)%in%toremove]
#maxperc<-max(colSums(findata[,2:ncol(findata)]))

#Only plot the top "nspecies" species
findata<-findata[1:nspecies,]

meltgenus<-melt(setDT(findata),id.vars="name",variable.name="Sample")
#I change the order just ofr aesthetic reasons
# meltgenus$Sample<-factor(meltgenus$Sample,,levels=sort(levels(meltgenus$Sample)))



#Color palette (remember that it has to have length equal or greater than the number of rows of genus_perc... here we only selected the top 15 Genera)
my_pal<-c("magenta3","tan1","darkslateblue","black","steelblue","darkorchid4",
                        "dodgerblue4","turquoise1","deepskyblue","forestgreen","orangered","firebrick1",
                        "darkorange","gold","springgreen","orange1","orange4","seagreen1",
                        "darkgreen","yellow4","darkred","cyan4","yellowgreen",
                        "yellow","darkmagenta","chartreuse2","cornflowerblue","darksalmon","deeppink","khaki")
pino<-ggplot(meltgenus, aes(x=Sample, y=value, fill=name)) + 
  # ggtitle(my_pal) +  #Add main title to plot
  coord_cartesian(ylim = c(0, maxperc)) +
  scale_fill_manual(values=my_pal) +
  geom_col(show.legend=TRUE) +
  xlab("\nSample") +
  ylab("Percentage\n") +
  theme(legend.key.size = unit(0.7,"line")) + #Change legend size
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -5, b = 0, l = -10))) + #Move axis title
  theme(axis.title.x = element_text(margin = margin(t = -2, r = 0, b = 0, l = 0))) +
  theme(legend.justification=c(0,0),legend.position="top",legend.title=element_blank()) +
  guides(fill=guide_legend(nrow=ceiling(length(unique(meltgenus$name))/4),byrow=TRUE)) +
  scale_y_continuous(expand = c(0.01,0.01)) + scale_x_discrete(expand = c(0.1,0.1) ) + # Expand to require No extra space around plot
  theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust=0.5),text = element_text(size=10)) #Rotate x axis text
ggsave(outgraph,pino)

}
plot_bracken()