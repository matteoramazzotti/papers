library(phyloseq)
library(ape)
library(decontam)
### this loads data as a phyloseq object, plus ancillary functions used in the analysis
load("data.RData")
source("functions.R")
sample_data(data)=read.delim("samples_age_ok.txt",row.names=1)

#screening of the top 10 genera in the samples
#screen = data.frame(rep(0,10))
#for (i in 1:41) {
#screen=cbind(screen,as(tax_table(data)[order(otu_table(data)[,i],decreasing=T),6][1:10],"vector"))
#} 
#colnames(screen)=c("remove",sample_names(data))

#dna=read.delim("DNA")
#sample_data(data)$is.control=c(rep(FALSE,40),TRUE)

#deconmtam step
#contamdf.freq <- isContaminant(data, method="frequency", conc="dna_conc",threshold=0.5)
#contamdf.prev <- isContaminant(data, method="prevalence", neg="is.control",threshold=0.5)
contamdf.comb <- isContaminant(data, method="combined", neg="is.control",conc="DNACONC",threshold=0.5)
decontam=rownames(tax_table(data)[contamdf.comb$contaminant,])

#bowtie human removal
human=as.character(read.delim(pipe("grep '>' bowtie/eucar_contaminant_ASV.fasta | perl -ne 's/>//g;print'"),header=F)[,1])

#final purge
toremove=unique(c(human,decontam))
data=prune_taxa(!taxa_names(data) %in% toremove,data)

#data.freq=prune_taxa(!contamdf.freq$contaminant,data)
#data.prev=prune_taxa(!contamdf.prev$contaminant,data)
#data.comb=prune_taxa(!contamdf.comb$contaminant,data)
#the "combined" mode was selected as the most aggressive
#data=data.comb

#minfreq=0.1
#screen=apply(otu_table(data),1,function(x) (sum(ifelse(x>0,1,0))/nsamples(data))>minfreq)
#data_filt=subset_taxa(data,screen)

data.phylum = tax_glom(data, taxrank = "Phylum", NArm = F)
data.class = tax_glom(data, taxrank = "Class", NArm = F)
data.order = tax_glom(data, taxrank = "Order", NArm = F)
data.family = tax_glom(data, taxrank = "Family", NArm = F)
data.genus = tax_glom(data, taxrank = "Genus", NArm = F)

#contamdf.freq.genus <- isContaminant(data.genus, method="frequency", conc="dna_conc",threshold=0.5)
#contamdf.prev.genus <- isContaminant(data.genus, method="prevalence", neg="is.control",threshold=0.5)
#contamdf.comb.genus <- isContaminant(data.genus, method="combined", neg="is.control",conc="dna_conc",threshold=0.5)
#data.genus.freq=prune_taxa(!contamdf.freq.genus$contaminant,data.genus)
#data.genus.prev=prune_taxa(!contamdf.prev.genus$contaminant,data.genus)
#data.genus.comb=prune_taxa(!contamdf.comb.genus$contaminant,data.genus)

#EVALUATION OF THE DECONTAMINATION
#par(mfrow=c(1,2))
#venn(list("freq.OTU"=tax_table(data)[contamdf.freq$contaminant,], "comb.OTU"=tax_table(data)[contamdf.comb$contaminant,] ,"prev.OTU"=tax_table(data)[contamdf.prev.genus$contaminant,]))
#venn(list("freq.genus"=tax_table(data.genus)[contamdf.freq.genus$contaminant,], "comb.genus"=tax_table(data.genus)[contamdf.comb.genus$contaminant,] ,"prev.genus"=tax_table(data.genus)[contamdf.prev.genus$contaminant,]))


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
ecodata<-data.frame(ecodata,eveness(data),sample_data(data)[,3:7])

### HEALTHY VS DISEASE in MALES AND FEMALES

#split by sex and resort for pairing
ecoM=ecodata[ecodata$SEX=="M",]
ecoM=ecoM[order(ecoM$DISEASE),]
ecoF=ecodata[ecodata$SEX=="F",]
ecoF=ecoF[order(ecoF$DISEASE),]

#pdf(file="paper/16-3-2022/data_ecoindexes_box_HvsT.pdf",width=18,height=12)
pdf(file="paper/21-7-2022/data_ecoindexes_box_HvsT.pdf",width=18,height=12)
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

summary(lm(rank(ecoH[,2])~factor(ecoH$SEX)+ecoH$AGE))
#sex is significant, Age is not
summary(lm(rank(ecoT[,2])~factor(ecoT$SEX)+ecoT$AGE))
#sex and age are not significant

pdf(file="paper/21-7-2022/data_ecoindexes_box_FvsM.pdf",width=18,height=12)
#pdf(file="paper/16-3-2022/data_ecoindexes_box_FvsM.pdf",width=18,height=12)
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
#dM=subset_samples(data.prop,sample_data(data.prop)$SEX=="M" | sample_data(data.prop)$SEX=="C")
#dF=subset_samples(data.prop,sample_data(data.prop)$SEX=="F" | sample_data(data.prop)$SEX=="C")
dM=subset_samples(data.prop,sample_data(data.prop)$SEX=="M")
dF=subset_samples(data.prop,sample_data(data.prop)$SEX=="F")
dH=subset_samples(data.prop,sample_data(data.prop)$DISEASE=="H")
dT=subset_samples(data.prop,sample_data(data.prop)$DISEASE=="T")

#permanovas
#Male: disease vs healthy
metadata <- as(sample_data(dM), "data.frame")
dist=phyloseq::distance(dM,method="bray")
adonis(dist ~ DISEASE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
#DISEASE    1    0.4960 0.49600  2.7829 0.1339  0.004 **
#Residuals 18    3.2082 0.17823         0.8661          
#Total     19    3.7042                 1.0000

adonis(dist ~ DISEASE + AGE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#DISEASE    1    0.4960 0.49600  2.9457 0.13390  0.002 **
#AGE        1    0.3456 0.34562  2.0526 0.09331  0.027 * 
#Residuals 17    2.8625 0.16838         0.77279          
#Total     19    3.7042                 1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Female: disease vs healthy
metadata <- as(sample_data(dF), "data.frame")
dist=phyloseq::distance(dF,method="bray")
adonis(dist ~ DISEASE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#DISEASE    1    0.3672 0.36723   1.567 0.08008  0.139
#Residuals 18    4.2183 0.23435         0.91992       
#Total     19    4.5856                 1.00000
adonis(dist ~ DISEASE + AGE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#DISEASE    1    0.4960 0.49600  2.9457 0.13390  0.002 **
#AGE        1    0.3456 0.34562  2.0526 0.09331  0.016 * 
#Residuals 17    2.8625 0.16838         0.77279          
#Total     19    3.7042                 1.00000 


#healthy: female vs male
metadata <- as(sample_data(dH), "data.frame")
dist=phyloseq::distance(dH,method="bray")
adonis(dist ~ SEX+AGE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#SEX        1    0.5695 0.56949 2.94260 0.14352  0.002 **
#AGE        1    0.1085 0.10852 0.56074 0.02735  0.897   
#Residuals 17    3.2901 0.19353         0.82913          
#Total     19    3.9681                 1.00000          


#disease: female vs male
metadata <- as(sample_data(dT), "data.frame")
dist=phyloseq::distance(dT,method="bray")
adonis(dist ~ SEX + AGE,data = metadata)
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#SEX        1    0.6137 0.61370 2.73086 0.13222  0.020 *
#AGE        1    0.2076 0.20756 0.92359 0.04472  0.415  
#Residuals 17    3.8204 0.22473         0.82307         
#Total     19    4.6416                 1.00000 

###################################################

#graphs
#pdf(file="paper/16-3-2022/PCoA_new.pdf",width=10,height=5)
pdf(file="paper/21-7-2022/PCoA_new.pdf",width=10,height=5)
par(mfrow=c(1,2))

library(car)
par(mfrow=c(1,2))
#TvsH in Males
pcoa<-pcoa(vegdist(t(sqrt(otu_table(dM)))))
cols=ifelse(sample_data(dM)$DISEASE=="T","blue","orange")
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=cols,
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
col=ifelse(sample_data(dF)$DISEASE=="T","blue","orange")
cols[21]="black"
plot(pcoa$vectors[,1],pcoa$vectors[,2],
	col=cols,
	main=paste0("Females - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
	xlab=paste0("PC1 (",round(pcoa$values$Relative_eig[1]*100,2),"%)"),
	ylab= paste0("PC2 (",round(pcoa$values$Relative_eig[2]*100,2),"%)"),
	cex=1.5,
	pch=16,
)
legend("topright",legend=c("T","H"),pch=16,col=ifelse(sample_data(dF)$DISEASE=="T","blue","orange"),cex=1)
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
	main=paste0("Healthy - PCoA of Bray-Curtis Distance\n on sqrt OTU abundance %"),
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

pdf(file="paper/21-7-2022/Top_10_Phyla.pdf",width=10,height=5)
dM=subset_samples(data.phylum.prop,sample_data(data.prop)$SEX=="M")
dF=subset_samples(data.phylum.prop,sample_data(data.prop)$SEX=="F")
meansM=apply(otu_table(dM),1,mean)
meansM=meansM[order(meansM,decreasing=T)][1:10]
meansM_lab=as.character(tax_table(dM)[match(names(meansM),rownames(tax_table(dM))),2])
meansM_lab[5]="Unknown"
meansF=apply(otu_table(dF),1,mean)
meansF=meansF[order(meansF,decreasing=T)][1:10]
meansF_lab=as.character(tax_table(dF)[match(names(meansF),rownames(tax_table(dF))),2])
meansF_lab[5]="Unknown"
meansF_lab[10]="Unknown"
par(mfrow=c(1,2),mar=c(9,4,4,3))
a=barplot(meansF,col="blue",ylab="Abundance %",xaxt="no", main="Top 10 most abundant phyla in females")
text(cex=1, x=a, y=-1, meansF_lab, xpd=TRUE, srt=45,adj=1)
a=barplot(meansM,col="orange",ylab="Abundance %",xaxt="no", main="Top 10 most abundant phyla in males")
text(cex=1, x=a, y=-1, meansM_lab, xpd=TRUE, srt=45,adj=1)
dev.off()


pdf(file="paper/21-7-2022/Top_10_genera.pdf",width=10,height=5)
dM=subset_samples(data.genus.prop,sample_data(data.genus.prop)$SEX=="M")
dF=subset_samples(data.genus.prop,sample_data(data.genus.prop)$SEX=="F")
meansF=apply(otu_table(dF),1,mean)
meansF=meansF[order(meansF,decreasing=T)][1:10]
meansF_lab=as.character(tax_table(dF)[match(names(meansF),rownames(tax_table(dF))),6])
meansF_lab[4]="Unknown"
meansF_lab[8]="Unknown"
meansM=apply(otu_table(dM),1,mean)
meansM=meansM[order(meansM,decreasing=T)][1:10]
meansM_lab=as.character(tax_table(dM)[match(names(meansM),rownames(tax_table(dM))),6])
meansM_lab[3]="Unknown"
par(mfrow=c(1,2),mar=c(9,4,4,3))
a=barplot(meansF,col="blue",ylab="Abundance %",xaxt="no", main="Top 10 most abundant genera in females")
text(cex=1, x=a, y=-1, meansF_lab, xpd=TRUE, srt=45,adj=1)
a=barplot(meansM,col="orange",ylab="Abundance %",xaxt="no", main="Top 10 most abundant genera in males")
text(cex=1, x=a, y=-1, meansM_lab, xpd=TRUE, srt=45,adj=1)
dev.off()


### DIFFERENTIAL ANALISIS, PAIRED
ranks<-c("Phylum","Class","Order","Family","Genus")
#data=subset_samples(data,sample_data(data)$SEX!="C")

dataF=subset_samples(data,sample_data(data)$SEX=="F")
pdf(file="paper/rev2/female_HvsT_paired.pdf",width=12,height=12)
resF<-DEanalysis_pair(dataF,"ID+DISEASE",ranks,0.05,1,0.5,g1=seq(1,20,2),g2=seq(2,20,2),1,"female_HvsT_paired","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls filterin=female_HvsT_paired strip=counts_ source=R out=paper/rev2/female_HvsT_paired.xls")

dataM=subset_samples(data,sample_data(data)$SEX=="M")
pdf(file="paper/rev2/male_HvsT_paired.pdf",width=12,height=12)
ranks<-c("Phylum","Class","Order","Family","Genus")
resM<-DEanalysis_pair(dataM,"ID+DISEASE",ranks,0.05,1,0.5,g1=seq(1,20,2),g2=seq(2,20,2),1,"male_HvsT_paired","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls filterin=male_HvsT_paired strip=counts_ source=R out=paper/rev2/male_HvsT_paired.xls")


### DIFFERENTIAL ANALISIS, UNPAIRED
ranks<-c("Phylum","Class","Order","Family","Genus")

dataT=subset_samples(data,sample_data(data)$DISEASE=="T")
pdf(file="paper/rev2/T_male_vs_female.pdf",width=12,height=12)
resT<-DEanalysis_mono(dataT,"SEX,SEX,M,F",ranks,0.05,1,0.5,0,"T_male_vs_female","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls+given filterin=T_male_vs_female strip=counts_ source=R out=paper/rev2/T_male_vs_female.xls")


dataH=subset_samples(data,sample_data(data)$DISEASE=="H")
pdf(file="paper/rev2/H_male_vs_female.pdf",width=12,height=12)
resH<-DEanalysis_mono(dataH,"SEX,SEX,M,F",ranks,0.05,1,0.5,0,"T_male_vs_female","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls+given filterin=T_male_vs_female strip=counts_ source=R out=paper/rev2/H_male_vs_female.xls")



###### PAPER REVISION #########################################
#A NEW SAMPLE DATA HAS BEEN CEATED WITH A TERMINAL AGE COLUMN

F_age=sample_data(data)$AGE[seq(1,40,by=4)]
M_age=sample_data(data)$AGE[seq(3,40,by=4)]
t.test(F_age,M_age)
sd(F_age)
sd(M_age)
#significant difference !!!!

#source("functions.R")
#THIS FILE NOW CONTAINS A FUNCTIN TO INCLUDE COVARIATES INTO THE DE ANALYSIS 

dataT=subset_samples(data,sample_data(data)$DISEASE=="T")
pdf(file="paper/rev2/T_male_vs_female.pdf",width=6,height=6)
resT<-DEanalysis_mono_corr(dataT,"SEX+AGE,SEX,M,F",ranks,0.05,1,0.5,0,"T_male_vs_female_given_age","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls filterin=T_male_vs_female_given strip=counts_ source=R out=paper/rev2/T_male_vs_female_given_age.xls")

dataH=subset_samples(data,sample_data(data)$DISEASE=="H")
pdf(file="paper/rev2/H_male_vs_female.pdf",width=6,height=6)
resH<-DEanalysis_mono_corr(dataH,"SEX+AGE,SEX,M,F",ranks,0.05,1,0.5,0,"H_male_vs_female_given_age","paper/rev2/")
dev.off()
system("data_xls_inject.pl dir=paper/rev2 filterout=pdf+xls filterin=H_male_vs_female_given strip=counts_ source=R out=paper/rev2/H_male_vs_female_given_age.xls")
#system("rm paper/rev1/*de.txt")


########## DE ANALYSIS REPRESENTATION USING CIRCOTAX ###########
#to instruct this plot all DE rows produced above must be joined in a long (many rows) table
conjoin=function(list) {
	out = list[[1]]
	if (length(list) > 1) {
		for (i in 2:length(list)) {
			out=rbind(out,list[[i]])
			cat("rows:",dim(out)[1],"\n")
		}
	}
	out
}
source("CircoTax3.R")
allF=conjoin(resF)
allM=conjoin(resM)
allT=conjoin(resT)
allH=conjoin(resH)

pdf(file="paper/rev2/AllM_DEtax_circo_out.pdf",width=10,height=10)
CircoTax3(allM,title="Males\nHealthy vs Tumor",ramp=c("orange","white","blue"),tax_col=8:13,fc_col=2,sort="no")
dev.off()
pdf(file="paper/rev2/AllF_DEtax_circo_out.pdf",width=10,height=10)
CircoTax3(allF,title="Females\nHealthy vs Tumor",ramp=c("orange","white","blue"),tax_col=8:13,fc_col=2,sort="no")
dev.off()
pdf(file="paper/rev2/AllT_DEtax_circo_out.pdf",width=10,height=10)
CircoTax3(allT,title="Tumor\nMales vs Females",ramp=c("orange","white","blue"),tax_col=8:13,fc_col=2,sort="no")
dev.off()
pdf(file="paper/rev2/AllH_DEtax_circo_out.pdf",width=10,height=10)
CircoTax3(allH,title="Healthy\nMales vs Females",ramp=c("orange","white","blue"),tax_col=8:13,fc_col=2,sort="no")
dev.off()


