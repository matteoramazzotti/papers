library(phyloseq)
library(ggplot2)
library(vegan)
library(plyr)
library(ape)

### preprocessing performance by micca
micca<-read.delim("micca_prep.log",sep="\t",header=T,row.names=1)
barplot(as.matrix(t(micca)),beside=T,col=c("black","red","blue","green"),las=2)
legend("top",legend=names(micca),fill=c("black","red","blue","green"))

### this loads data as a phyloseq object, plus ancillary functions used in the analysis
load("data.Rdata")
source("functions.R")

### GLOMMING OTUs TO DIFFRENT RANKS
data.phylum = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.family = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)

### EXPORT ABUNDANCE TABLES
write.table(cbind(as(otu_table(data.phylum),"matrix"),as(tax_table(data.phylum),"matrix")),file="plots/counts_phylum.txt",sep="\t")
write.table(cbind(as(otu_table(data.class),"matrix"),as(tax_table(data.class),"matrix")),file="plots/counts_class.txt",sep="\t")
write.table(cbind(as(otu_table(data.order),"matrix"),as(tax_table(data.order),"matrix")),file="plots/counts_order.txt",sep="\t")
write.table(cbind(as(otu_table(data.family),"matrix"),as(tax_table(data.family),"matrix")),file="plots/counts_family.txt",sep="\t")
write.table(cbind(as(otu_table(data.genus),"matrix"),as(tax_table(data.genus),"matrix")),file="plots/counts_genus.txt",sep="\t")
write.table(cbind(as(otu_table(data),"matrix"),as(tax_table(data),"matrix")),file="plots/counts_otu.txt",sep="\t")

### ECOLOGICAL INDICES
ecodata<-estimate_richness(data, split = TRUE, measures = c("Observed","Chao1","Shannon"))
ecodata<-data.frame(ecodata,eveness(data))
factor<-sample_data(data)$Disease

t<-wilcox.test(ecodata[,2]~factor,paired=T)
boxplot(ecodata[,2]~factor,cex.main=1.5,cex.axis=1.5,main=paste("Chao1 index","\n","paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("I","NI"),cex.axis=1.5)
stripchart(ecodata[,2]~factor,vertical=T,add=T,pch=16,cex=0.8)
c<-ifelse(ecodata[11:20,2]-ecodata[1:10,2]>0,"blue","red")
segments(rep(1,10),ecodata[11:20,2],rep(2,10),ecodata[1:10,2],col=c)
legend("topleft",legend=table(c),text.col=c("blue","red"),bty="n")

t<-wilcox.test(ecodata[,4]~factor,paired=T)
boxplot(ecodata[,4]~factor,cex.main=1.5,cex.axis=1.5,main=paste("Shannon index","\n","paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("I","NI"),cex.axis=1.5)
stripchart(ecodata[,4]~factor,vertical=T,add=T,pch=16,cex=0.8)
c<-ifelse(ecodata[11:20,4]-ecodata[1:10,4]>0,"blue","red")
segments(rep(1,10),ecodata[11:20,4],rep(2,10),ecodata[1:10,4],col=c)
legend("topleft",legend=table(c),text.col=c("blue","red"),bty="n")

t<-wilcox.test(ecodata[,5]~factor,paired=T)
boxplot(ecodata[,5]~factor,cex.main=1.5,cex.axis=1.5,main=paste("Eveness index","\n","paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("I","NI"),cex.axis=1.5)
stripchart(ecodata[,5]~factor,vertical=T,add=T,pch=16,cex=0.8)
c<-ifelse(ecodata[11:20,5]-ecodata[1:10,5]>0,"blue","red")
segments(rep(1,10),ecodata[11:20,5],rep(2,10),ecodata[1:10,5],col=c)
legend("topleft",legend=table(c),text.col=c("blue","red"),bty="n")
dev.off()

