library(phyloseq)
library(ggplot2)
library(vegan)
library(ape)

#this loads all variables until the save command
load("data.RData")
source("functions.R")
#dir.create("plots")


### DATA GLOMMING ###
data.phylum = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.family = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)

data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phylum.prop <- transform_sample_counts(data.phylum, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.family.prop <- transform_sample_counts(data.family, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)

### COORDINATE PLOT
coord_plot(data,rank="Phylum",num=5,col=c(rep("blue",17),rep("orange",10)),lty=c(rep(1,17),rep(2,10)),legnames=c("CEL","POT"),legcol=c("blue","orange"),exclude="Cyanobacteria/Chloroplast",exclude_noname=T)

### RAREFACTION ANALYSIS

set.seed(1)
data.rare<-rarefy_even_depth(data, sample.size=0.9*min(sample_sums(data)), replace=F,rngseed=1)
r<-rarecurve(t(otu_table(data)), step=50,label=T,col=color) 
evalslopes(r,0.01)
legend("bottomright",legend=c("Celiac","Possible"),col=c("red","black"),box.lty=0,ncol=4,lty=1)

########## ECOLOGICAL INDICES ####################################

richness<-estimate_richness(data.prop, split = TRUE, measures = c("Observed","Chao1","Shannon"))

par(mfrow=c(3,3))

plot(1:dim(richness)[1],richness$Observed,col=color,xlab="Samples",ylab="Observed OTU",xaxt="n",main="Richness: OTU number")
order<-order(richness$Observed)
plot(1:dim(richness)[1],richness$Observed[order],col=color[order],xlab="Samples",ylab="Observed OTU",xaxt="n",main="Richness: OTU number")
legend("topleft",legend=c("celiac","potential"),col=c("black","red"),box.lty=0,ncol=4,pch=1)
w<-wilcox.test(richness$Observed~sample_data(data)$Type)
plot(richness$Observed~sample_data(data)$Type,xlab="Type",ylab="Observed OTU",main=paste("Richness: OTU number\nW=",round(w$statistic,2),"P=",round(w$p.value,2)),outline=F)
stripchart(richness$Observed~sample_data(data)$Type,vertical=T,add=T,method="jitter",pch=1)

plot(1:dim(richness)[1],richness$Chao1,col=color,xlab="Samples",ylab="Chao1 OTU",xaxt="n",main="Richness: Chao1")
order<-order(richness$Chao1)
plot(1:dim(richness)[1],richness$Chao1[order],col=color[order],xlab="Samples",ylab="Chao1 OTU",xaxt="n",main="Richness: Chao1")
legend("topleft",legend=c("celiac","potential"),col=c("black","red"),box.lty=0,ncol=4,pch=1)
w<-wilcox.test(richness$Chao1~sample_data(data)$Type)
plot(richness$Chao1~sample_data(data)$Type,xlab="Type",ylab="Richness: Chao1",main=paste("Chao1 index\nW=",round(w$statistic,2),"P=",round(w$p.value,2)),outline=F)
stripchart(richness$Chao1~sample_data(data)$Type,vertical=T,add=T,method="jitter",pch=1)

plot(1:dim(richness)[1],richness$Shannon,col=color,xlab="Samples",ylab="Shannon OTU",xaxt="n",main="Richness: Shannon")
order<-order(richness$Shannon)
plot(1:dim(richness)[1],richness$Shannon[order],col=color[order],xlab="Samples",ylab="Shannon OTU",xaxt="n",main="Richness: Shannon")
legend("topleft",legend=c("celiac","potential"),col=c("black","red"),box.lty=0,ncol=4,pch=1)
w<-wilcox.test(richness$Shannon~sample_data(data)$Type)
plot(richness$Shannon~sample_data(data)$Type,xlab="Type",ylab="Richness: Shannon",main=paste("Shannon index\nW=",round(w$statistic,2),"P=",round(w$p.value,2)),outline=F)
stripchart(richness$Shannon~sample_data(data)$Type,vertical=T,add=T,method="jitter",pch=1)


### filters ###


### DEseq2 ###
library("DESeq2")

data.filt = filter_taxa(data, function(x) mean(x) > 1e-5, TRUE)

DE_analysis<-function(data,outfile,pthr,fcthr) {
	DE<-phyloseq_to_deseq2(data, ~ Type)
	DEres<-DESeq(DE, test="Wald", fitType="parametric")
	r<-results(DEres)
	r<-cbind(r,tax_table(data))
	sel<-!is.na(r$padj) & r$padj<=pthr & abs(r$log2FoldChange)>=fcthr
	res<-r[sel,]

	o<-order(sample_data(data)$Type)
	labels<-paste(sep=" - " ,colnames(otu_table(data))[o],sample_data(data)$Type[o])

	pdf(file=paste0(outfile,".pdf"),width=15,height=10)
	par(mar=c(7,4,4,2))
	for (i in 1:dim(res)[1]) {
		x<-match(rownames(res)[i],rownames(d))
		barplot(otu_table(data)[x,o],las=2,names.arg=labels,main=paste(sep="|",collapse="|",c(row.names(r)[x],as.character(tax_table(data)[x,]))))
	}
	dev.off()
	write.table(res,sep="\t",file=paste0(outfile,".txt"))
	out<-list(DEres,r,res)
	out
}

res_OTU<-DE_analysis(data,out="plots/DE_results_OTU",pthr=0.05,fcthr=1)
res_phylum<-DE_analysis(data.phylum,out="plots/DE_results_phylum",pthr=0.05,fcthr=1)
res_class<-DE_analysis(data.class,out="plots/DE_results_class",pthr=0.05,fcthr=1)
res_order<-DE_analysis(data.order,out="plots/DE_results_order",pthr=0.05,fcthr=1)
res_family<-DE_analysis(data.family,out="plots/DE_results_family",pthr=0.05,fcthr=1)
res_genus<-DE_analysis(data.genus,out="plots/DE_results_genus",pthr=0.05,fcthr=1)

DE<-phyloseq_to_deseq2(data.genus, ~ Type)
vst_genus <- assay(varianceStabilizingTransformation(DE, blind=FALSE))
clr_genus<-microbiome::transform(data.genus, transform = "clr", target = "OTU", shift = 0, scale = 1)


pdf(file="plots/heatmaps.pdf",height=10,width=10)
heatmap(vst_genus,ColSideColors=color,cexCol=1.5)
heatmap(otu_table(clr_genus),ColSideColors=color,cexCol=1.5)
dev.off()


