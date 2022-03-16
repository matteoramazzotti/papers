library(phyloseq)
library(ape)
### this loads data as a phyloseq object, plus ancillary functions used in the analysis
load("data.RData")
source("functions.R")

data.phylum = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.family = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)

### DATA IS NORMALIZED IN COLUMN PERCENT
data.prop <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
data.phylum.prop <- transform_sample_counts(data.phylum, function(otu) otu/sum(otu)*100)
data.class.prop <- transform_sample_counts(data.class, function(otu) otu/sum(otu)*100)
data.order.prop <- transform_sample_counts(data.order, function(otu) otu/sum(otu)*100)
data.family.prop <- transform_sample_counts(data.family, function(otu) otu/sum(otu)*100)
data.genus.prop <- transform_sample_counts(data.genus, function(otu) otu/sum(otu)*100)

####### FIGURE 1 ##########

data=subset_samples(data,sample_data(data)$SEX!="C")
ecodata<-estimate_richness(data, split = TRUE, measures = c("Observed","Chao1","Shannon"))
ecodata<-data.frame(ecodata,eveness(data),sample_data(data)[,3:5])

### HEALTHY VS DISEASE in MALES AND FEMALES

#split by sex and resort for pairing
ecoM=ecodata[ecodata$SEX=="M",]
ecoM=ecoM[order(ecoM$DISEASE),]
ecoF=ecodata[ecodata$SEX=="F",]
ecoF=ecoF[order(ecoF$DISEASE),]

pdf(file="paper/16-3-2022/data_ecoindexes_box_HvsT.pdf",width=18,height=12)
par(mfrow=c(2,3))

#Chao1 in Males H vs T (col 2)
t<-wilcox.test(ecoM[,2]~factor(ecoM$DISEASE),paired=T)
boxplot(ecoM[,2]~factor(ecoM$DISEASE),cex.main=2,cex.axis=2,main=paste("Male - Chao1 index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("MH","MT"),cex.axis=2)
stripchart(ecoM[,2]~factor(ecoM$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoM[1:10,2]-ecoM[11:20,2]>0,"blue","orange")
segments(rep(1,10),ecoM[1:10,2],rep(2,10),ecoM[11:20,2],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

#Shannon in Males H vs T (col 4)
t<-wilcox.test(ecoM[,4]~factor(ecoM$DISEASE),paired=T)
boxplot(ecoM[,4]~factor(ecoM$DISEASE),cex.main=2,cex.axis=2,main=paste("Male - Shannon index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("MH","MT"),cex.axis=2)
stripchart(ecoM[,4]~factor(ecoM$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoM[1:10,4]-ecoM[11:20,4]>0,"blue","orange")
segments(rep(1,10),ecoM[1:10,4],rep(2,10),ecoM[11:20,4],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

#Eveness in Males H vs T (col 5)
t<-wilcox.test(ecoM[,5]~factor(ecoM$DISEASE),paired=T)
boxplot(ecoM[,5]~factor(ecoM$DISEASE),cex.main=2,cex.axis=2,main=paste("Male - eveness index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("MH","MT"),cex.axis=2)
stripchart(ecoM[,5]~factor(ecoM$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoM[1:10,5]-ecoM[11:20,5]>0,"blue","orange")
segments(rep(1,10),ecoM[1:10,5],rep(2,10),ecoM[11:20,5],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

#Chao1 in Females H vs T
t<-wilcox.test(ecoF[,2]~factor(ecoF$DISEASE),paired=T)
boxplot(ecoF[,2]~factor(ecoF$DISEASE),cex.main=2,cex.axis=2,main=paste("Female - Chao1 index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("FH","FT"),cex.axis=2)
stripchart(ecoF[,2]~factor(ecoF$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoF[1:10,2]-ecoF[11:20,2]>0,"blue","orange")
segments(rep(1,10),ecoF[1:10,2],rep(2,10),ecoF[11:20,2],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

#Shannon in Females H vs T
t<-wilcox.test(ecoF[,4]~factor(ecoF$DISEASE),paired=T)
boxplot(ecoF[,4]~factor(ecoF$DISEASE),cex.main=2,cex.axis=2,main=paste("Female - Shannon index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("FH","FT"),cex.axis=2)
stripchart(ecoF[,4]~factor(ecoF$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoF[1:10,4]-ecoF[11:20,4]>0,"blue","orange")
segments(rep(1,10),ecoF[1:10,4],rep(2,10),ecoF[11:20,4],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

#Eveness in Females H vs T
t<-wilcox.test(ecoF[,5]~factor(ecoF$DISEASE),paired=T)
boxplot(ecoF[,5]~factor(ecoF$DISEASE),cex.main=2,cex.axis=2,main=paste("Female - eveness index","\n","Paired Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n",border="white")
axis(1,at=c(1,2),labels=c("FH","FT"),cex.axis=2)
stripchart(ecoF[,5]~factor(ecoF$DISEASE),vertical=T,add=T,pch=16,cex=1.5)
c<-ifelse(ecoF[1:10,5]-ecoF[11:20,5]>0,"blue","orange")
segments(rep(1,10),ecoF[1:10,5],rep(2,10),ecoF[11:20,5],col=c,lwd=3)
legend("topleft",legend=table(c),text.col=c("blue","orange"),bty="n",cex=3)

dev.off()

### MALE VS FEMALES IN HEALTHY AND DISEASE
ecoH=ecodata[ecodata$DISEASE=="H",]
ecoH=ecoH[order(ecoH$SEX),]
ecoT=ecodata[ecodata$DISEASE=="T",]
ecoT=ecoT[order(ecoT$SEX),]

pdf(file="paper/16-3-2022/data_ecoindexes_box_FvsM.pdf",width=18,height=12)
par(mfrow=c(2,3))
#Chao1 in HEALTHY female vs male (col 2)
t<-wilcox.test(ecoH[,2]~factor(ecoH$SEX))
boxplot(ecoH[,2]~factor(ecoH$SEX),cex.main=2,cex.axis=2,main=paste("Healthy Chao1 index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FH","MH"),cex.axis=2)
stripchart(ecoH[,2]~factor(ecoH$SEX),vertical=T,add=T,pch=16,cex=1.5)

#Shannon in HEALTHY female vs male (col 4)
t<-wilcox.test(ecoH[,4]~factor(ecoH$SEX))
boxplot(ecoH[,4]~factor(ecoH$SEX),cex.main=2,cex.axis=2,main=paste("Healthy - Shannon index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FH","MH"),cex.axis=2)
stripchart(ecoH[,4]~factor(ecoH$SEX),vertical=T,add=T,pch=16,cex=1.5)

#Eveness in HEALTHY female vs male (col 5)
t<-wilcox.test(ecoH[,5]~factor(ecoH$SEX))
boxplot(ecoH[,5]~factor(ecoH$SEX),cex.main=2,cex.axis=2,main=paste("Healthy - eveness index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FH","MH"),cex.axis=2)
stripchart(ecoH[,5]~factor(ecoH$SEX),vertical=T,add=T,pch=16,cex=1.5)

#Chao1 in TUMOR female vs male (col 2)
t<-wilcox.test(ecoT[,2]~factor(ecoT$SEX),paired=T)
boxplot(ecoT[,2]~factor(ecoT$SEX),cex.main=2,cex.axis=2,main=paste("Tumor - Chao1 index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FT","MT"),cex.axis=2)
stripchart(ecoT[,2]~factor(ecoT$SEX),vertical=T,add=T,pch=16,cex=1.5)

#Shannon in TUMOR female vs male 
t<-wilcox.test(ecoT[,4]~factor(ecoT$SEX),paired=T)
boxplot(ecoT[,4]~factor(ecoT$SEX),cex.main=2,cex.axis=2,main=paste("Tumor - Shannon index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FT","MT"),cex.axis=2)
stripchart(ecoT[,4]~factor(ecoT$SEX),vertical=T,add=T,pch=16,cex=1.5)

#Eveness in TUMOR female vs male 
t<-wilcox.test(ecoT[,5]~factor(ecoT$SEX),paired=T)
boxplot(ecoT[,5]~factor(ecoT$SEX),cex.main=2,cex.axis=2,main=paste("Tumor - eveness index","\n","Wilcox test: w=",t$statistic,"p=",round(t$p.value,3)),xaxt="n")
axis(1,at=c(1,2),labels=c("FT","MT"),cex.axis=2)
stripchart(ecoT[,5]~factor(ecoT$SEX),vertical=T,add=T,pch=16,cex=1.5)

dev.off()

######## FIGURE 2 #######

# PCOA
library(vegan)

#datasets
dM=subset_samples(data.prop,sample_data(data.prop)$SEX=="M")
dF=subset_samples(data.prop,sample_data(data.prop)$SEX=="F")
dH=subset_samples(data.prop,sample_data(data.prop)$DISEASE=="H")
dT=subset_samples(data.prop,sample_data(data.prop)$DISEASE=="T")

#permanovas
#Male: disease vs healthy
metadata <- as(sample_data(dM), "data.frame")
dist=phyloseq::distance(dM,method="bray")
adonis(dist ~ DISEASE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#DISEASE    1    0.5045 0.50452  2.6697 0.12916  0.005 **

#Female: disease vs healthy
metadata <- as(sample_data(dF), "data.frame")
dist=phyloseq::distance(dF,method="bray")
adonis(dist ~ DISEASE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#DISEASE    1    0.3711 0.37106  1.4606 0.07505  0.144

#healthy: female vs male
metadata <- as(sample_data(dH), "data.frame")
dist=phyloseq::distance(dH,method="bray")
adonis(dist ~ SEX,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SEX        1    0.6434 0.64341  3.0741 0.14587  0.001 ***

#disease: female vs male
metadata <- as(sample_data(dT), "data.frame")
dist=phyloseq::distance(dT,method="bray")
adonis(dist ~ SEX,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#SEX        1    0.6269 0.62694  2.6822 0.12969  0.025 *

#graphs
pdf(file="paper/16-3-2022/PCoA_new.pdf",width=10,height=5)
par(mfrow=c(1,2))

###################################################
library(car)
#TvsH in Males
pcoa<-pcoa(vegdist(t(sqrt(otu_table(dM)))))
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=ifelse(sample_data(dM)$DISEASE=="T","blue","orange"),
	main=paste0("Males - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
	xlab=paste0("PC1 (",round(pcoa$values$Relative_eig[1]*100,2),"%)"),
	ylab= paste0("PC2 (",round(pcoa$values$Relative_eig[2]*100,2),"%)"),
	cex=1.5,
	pch=16
)
legend("topleft",legend=c("T","H"),pch=16,col=ifelse(sample_data(dM)$DISEASE=="T","blue","orange"),cex=1)
center1 <- apply(pcoa$vectors[sample_data(dM)$DISEASE=="T",1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[sample_data(dM)$DISEASE=="T",1:2])
center2 <- apply(pcoa$vectors[sample_data(dM)$DISEASE=="H",1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[sample_data(dM)$DISEASE=="H",1:2])
ellipse(center1, cov_mat1, center.pch=0, col="blue", fill=TRUE, fill.alpha=0.1,lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, col="orange", fill=TRUE, fill.alpha=0.1, lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)

#TvsH in Females
pcoa<-pcoa(vegdist(t(sqrt(otu_table(dF)))))
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=ifelse(sample_data(dF)$DISEASE=="T","blue","orange"),
	main=paste0("Females - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
	xlab=paste0("PC1 (",round(pcoa$values$Relative_eig[1]*100,2),"%)"),
	ylab= paste0("PC2 (",round(pcoa$values$Relative_eig[2]*100,2),"%)"),
	cex=1.5,
	pch=16,
)
legend("bottomright",legend=c("T","H"),pch=16,col=ifelse(sample_data(dF)$DISEASE=="T","blue","orange"),cex=1)
center1 <- apply(pcoa$vectors[sample_data(dF)$DISEASE=="T",1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[sample_data(dF)$DISEASE=="T",1:2])
center2 <- apply(pcoa$vectors[sample_data(dF)$DISEASE=="H",1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[sample_data(dF)$DISEASE=="H",1:2])
ellipse(center1, cov_mat1, center.pch=0, col="blue", fill=TRUE, fill.alpha=0.1,lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, col="orange", fill=TRUE, fill.alpha=0.1, lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)

###################################################
#MvsF in Tumor
pcoa<-pcoa(vegdist(t(sqrt(otu_table(dT)))))
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=ifelse(sample_data(dT)$SEX=="M","blue","orange"),
	main=paste0("Tumor - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
	xlab=paste0("PC1 (",round(pcoa$values$Relative_eig[1]*100,2),"%)"),
	ylab= paste0("PC2 (",round(pcoa$values$Relative_eig[2]*100,2),"%)"),
	cex=1.5,
	pch=16
)
legend("bottomright",legend=c("M","F"),pch=16,col=ifelse(sample_data(dT)$SEX=="M","blue","orange"),cex=1)
center1 <- apply(pcoa$vectors[sample_data(dT)$SEX=="M",1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[sample_data(dT)$SEX=="M",1:2])
center2 <- apply(pcoa$vectors[sample_data(dT)$SEX=="F",1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[sample_data(dT)$SEX=="F",1:2])
ellipse(center1, cov_mat1, center.pch=0, col="blue", fill=TRUE, fill.alpha=0.1,lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, col="orange", fill=TRUE, fill.alpha=0.1, lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)

#MvsF in Healthy
pcoa<-pcoa(vegdist(t(sqrt(otu_table(dH)))))
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=ifelse(sample_data(dH)$SEX=="M","blue","orange"),
	main=paste0("Tumor - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
	xlab=paste0("PC1 (",round(pcoa$values$Relative_eig[1]*100,2),"%)"),
	ylab= paste0("PC2 (",round(pcoa$values$Relative_eig[2]*100,2),"%)"),
	cex=1.5,
	pch=16
)
legend("topright",legend=c("M","F"),pch=16,col=ifelse(sample_data(dH)$SEX=="M","blue","orange"),cex=1)
center1 <- apply(pcoa$vectors[sample_data(dH)$SEX=="M",1:2], 2, mean)
cov_mat1 <- cov(pcoa$vectors[sample_data(dH)$SEX=="M",1:2])
center2 <- apply(pcoa$vectors[sample_data(dH)$SEX=="F",1:2], 2, mean)
cov_mat2 <- cov(pcoa$vectors[sample_data(dH)$SEX=="F",1:2])
ellipse(center1, cov_mat1, center.pch=0, col="blue", fill=TRUE, fill.alpha=0.1,lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)
ellipse(center2, cov_mat2, center.pch=0, col="orange", fill=TRUE, fill.alpha=0.1, lty=0, radius=sqrt(2 * qf(.95, 2, 9999)),add=T)

dev.off()

pdf(file="paper/Top_Phyla.pdf",width=6,height=5)

dM=subset_samples(data.phylum.prop,sample_data(data.prop)$SEX=="M")
dF=subset_samples(data.phylum.prop,sample_data(data.prop)$SEX=="F")

meansM=apply(otu_table(dM),1,mean)
meansM=meansM[order(meansM,decreasing=T)][1:10]
meansM_lab=as.character(tax_table(dM)[match(names(meansM),rownames(tax_table(dM))),2])
meansM_lab[6]="Unknown"

meansF=apply(otu_table(dF),1,mean)
meansF=meansF[order(meansF,decreasing=T)][1:10]
meansF_lab=as.character(tax_table(dF)[match(names(meansF),rownames(tax_table(dF))),2])
meansF_lab[5]="Unknown"
#F and M most abundante are the same in different order
meansF=meansF[match(meansM_lab,meansF_lab)]

p=as.matrix(cbind(meansM,meansF))
par(mar=c(9,4,4,3))
a=barplot(t(p),beside=T,col=c("blue","orange"),ylab="Abundance %",xaxt="no", main="Top 5 most abundant phyla")
text(cex=1, x=colMeans(a), y=-1, meansM_lab, xpd=TRUE, srt=45,adj=1)
legend("topright",legend=c("Male","Female"),text.col=c("blue","orange"),bty="n",cex=1,fill=c("blue","orange"))
abline(h=0)

dev.off()

### DIFFERENTIAL ANALISIS, PAIRED
load("data.RData")
data=subset_samples(data,sample_data(data)$SEX!="C")

dataF=subset_samples(data,sample_data(data)$SEX=="F")
pdf(file="plots/female_HvsT_paired.pdf",width=12,height=12)
ranks<-c("Phylum","Class","Order","Family","Genus","OTU")
resF<-DEanalysis_pair(dataF,"ID+DISEASE",ranks,0.05,1,0.5,g1=seq(1,20,2),g2=seq(2,20,2,),"female_HvsT_paired")
dev.off()

dataM=subset_samples(data,sample_data(data)$SEX=="M")
pdf(file="plots/male_HvsT_paired.pdf",width=12,height=12)
ranks<-c("Phylum","Class","Order","Family","Genus","OTU")
resM<-DEanalysis_pair(dataM,"ID+DISEASE",ranks,0.05,1,0.5,g1=seq(1,20,2),g2=seq(2,20,2,),"male_HvsT_paired")
dev.off()

### DIFFERENTIAL ANALISIS, UNPAIRED
load("data.RData")
data=subset_samples(data,sample_data(data)$SEX!="C")

dataT=subset_samples(data,sample_data(data)$DISEASE=="T")
ranks<-c("Phylum","Class","Order","Family","Genus","OTU")
pdf(file="plots/T_male_vs_female.pdf",width=12,height=12)
res<-DEanalysis_mono(dataT,"SEX",ranks,0.05,1,0.5,"T_male_vs_female")
dev.off()

ranks<-c("Phylum","Class","Order","Family","Genus","OTU")
dataH=subset_samples(data,sample_data(data)$DISEASE=="H")
pdf(file="plots/H_male_vs_female.pdf",width=12,height=12)
res<-DEanalysis_mono(dataH,"SEX",ranks,0.05,1,0.5,"H_male_vs_female")
dev.off()

########## DE ANALYSIS REPRESENTATION ###########
#to instruct this plot all DE rows produced above must be joined in a long (many rows) table
source("DEtax_Circo.R")
pdf(file="DEtax_circo_out.pdf",width=10,height=10)
DEtax_circo("plots/H_male_vs_female_ok.txt","H_male_vs_female")
DEtax_circo("plots/T_male_vs_female.txt","T_male_vs_female")
DEtax_circo("plots/male_HvsT_paired.txt","male_HvsT_paired")
dev.off()

