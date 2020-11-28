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

### figure 7
dend<-as.dendrogram(hclust(dist(sqrt(t(otu_table(data.prop))))))
col<-ifelse(sample_data(data.prop)$Disease[o]=="NI","blue","red")
o<-match(labels(dend),sample_names(data.prop))
newlab<-ifelse(o>=10,o-10,o)
labels(dend)<-paste0("S",newlab)
labels_colors(dend)<-col
plot(dend,main=paste0("Complete clustering of Euclidean distance on\n sqrt transformed OTU percent abundance"))
text(18,10,"NI",col="blue",cex=1.2)
text(19,10,"I",col="red",cex=1.2)

pcoa<-pcoa(vegdist(t(otu_table(data.prop))))
plot(pcoa$vectors[,1],pcoa$vectors[,2],pch=".",main=paste0("PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),xlab="PC1",ylab= "PC2")
text(pcoa$vectors[,1],pcoa$vectors[,2],paste0("S",rep(1:10,2)),col=ifelse(sample_data(data.prop)$Disease=="NI","blue","red"))

center1 <- apply(pcoa$vectors[1:10,1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[1:10,1:2])
center2 <- apply(pcoa$vectors[11:20,1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[11:20,1:2])
ellipse(center1, cov_mat1, center.pch=0, col="blue", fill=TRUE, fill.alpha=0.1,lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, col="red", fill=TRUE, fill.alpha=0.1, lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)

center1 <- apply(pcoa$vectors[c(1:3,5:10),1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[c(1:3,5:10),1:2])
center2 <- apply(pcoa$vectors[c(11:13,15:20),1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[c(11:13,15:20),1:2])
ellipse(center1, cov_mat1, center.pch=0, lty=4,col="blue", radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, lty=4,col="red", radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
dev.off()

### 
coord_plot(data,rank="Phylum",num=5,col=c(rep("blue",10),rep("orange",10)),legnames=c("NI","I"),legcol=c("blue","orange"),exclude="Cyanobacteria/Chloroplast",exclude_noname=T)
