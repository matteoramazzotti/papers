library(phyloseq)
library(DESeq2)

good<-function(data) {
	if (class(data)=="phyloseq") {
		data<-otu_table(data)
	}
	s<-apply(data,2,function(x) 1-(sum(ifelse(x==1,1,0))/sum(ifelse(x>0,x,0))))
	t(s*100)
}

eveness<-function(data) {
	if (class(data)=="phyloseq") {
		data<-otu_table(data)
	}
	richness<-estimate_richness(data, split = TRUE, measures = "Shannon")
	colnames(richness)<-"Eveness"
	otus<-apply(otu_table(data),2,function(x) sum(ifelse(x>0,1,0)))
	richness/log(otus)
}

summarize_phyloseq<-function(data) {
	sums<-rowSums(otu_table(data))
	totr=sum(sums)
	toto=dim(otu_table(data))[1]
	tax<-NULL
	otus<-NULL
	reads<-NULL
	ranks<-c("Phylum","Class","Order","Family","Genus")
	if(ncol(tax_table(data)) == 7) {
		ranks<-c("Phylum","Class","Order","Family","Genus","Species")
	}
	for (i in ranks) {
		tax<-c(tax,length(grep("\\w",perl=T,unique(tax_table(data)[,i]))))
		otus<-c(otus,length(grep("\\w",perl=T,tax_table(data)[,i])))
		reads<-c(reads,sum(sums[grep("\\w",perl=T,tax_table(data)[,i])]))
	}
	data.frame("Rank"=ranks,"Count"=tax,"Reads"=reads,"percReads"=reads/totr*100,"OTU"=otus,"percOTU"=otus/toto*100)
}

data.threshold<-function(data,rank="Phylum",thr=0.5,prop=50,method=c("mean","maj")) {
	if (rank != "OTU") {
		data<-tax_glom(data, taxrank = rank, NArm = F)
	}
	data <- transform_sample_counts(data, function(otu) otu/sum(otu)*100)
	#a filter on low occurrence OTUs
	if(method == "mean") {
		sel<-taxa_names(data)[ifelse(rowMeans(otu_table(data))<=thr,T,F)]
	}
	if(method == "maj") {
		v<-apply(otu_table(data),1,function(x) sum(ifelse(x>thr,1,0)))
		sel<-taxa_names(data)[v>=prop]
	}
	#those not passing the filter are merged together
	data.filt = merge_taxa(data, sel,archetype=1)
	#and assigned the "other" taxonomy
	tax_table(data.filt)[,rank][match(NA,tax_table(data.filt)[,rank])]<-"other"
	data.filt
}
			 
coord_plot<-function(original,rank="Phylum",num=5,col="black",lty=1,legnames=null,legcol=null,exclude=NULL,include=NULL,exclude_noname=F) {
	i_data<-tax_glom(original, taxrank = rank, NArm = F)
	data<-transform_sample_counts(i_data, function(otu) otu/sum(otu)*100)
	d<-otu_table(data)[order(rowMeans(otu_table(data)),decreasing=T),]
	t<-tax_table(data)[order(rowMeans(otu_table(data)),decreasing=T),][,rank]
	if(!is.null(exclude)) {
		todel1<-match(exclude,t)
		d<-d[-todel1,]
		t<-t[-todel1,]
	}
	if(exclude_noname==TRUE) {
		todel2<-match("",t)
		d<-d[-todel2,]
		t<-t[-todel2,]
	}
	d<-d[1:num,]
	t<-t[1:num]
	cat(t,"\n")
	ylim<-c(min(d),max(d))
	for (i in 1:nsamples(data)) {
		if (i == 1) {
			par(mar=c(12,4,4,1))
			plot(1:num,d[,i],axes=F,type="l",ylim=ylim,col=col[i],xlab="",ylab="Abundance %",main=paste("Coordinate plot of top",num,rank))
		} else {
			lines(1:num,d[,i],col=col[i],lty=lty[i])
		}
	}
	if(!is.null(legnames)) {
		legend("topright",legnames,lty=unique(lty),bty="n",col=legcol)
	}
	axis(1,1:num,t,las=2)
	axis(2)
	list(d,t)
}

DEanalysis_mono<-function(data,ranks,pval,log2fc,zeros,name) {
	results<-list()
	for (i in 1:length(ranks)) {
		tot<-length(sample_names(data))
		cat(">>> Processing",ranks[i],"with",tot,"rows\n")
		if(ranks[i] != "OTU") {
			d<-tax_glom(data, taxrank = ranks[i], NArm = T)
		}
		else {
			d<-data
		}
		d.prop<-transform_sample_counts(d, function(otu) otu/sum(otu)*100)
		d<-filter_taxa(d,function(x) sum(ifelse(x==0,1,0))/length(x)< zeros,TRUE)
		d.prop<-filter_taxa(d.prop,function(x) sum(ifelse(x==0,1,0))/length(x)< zeros,TRUE)
		otuDE<-phyloseq_to_deseq2(d, ~ Surgery)
		otuDEres<-DESeq(otuDE, test="Wald", fitType="parametric")
		otuDEresults<-results(otuDEres, cooksCutoff = FALSE, alpha=pval)
		#coeff for shrink is set to the the last constrast, the same used by default in the "results" call
		coeff<-rev(colnames(coef(otuDEres)))[1]
		shrink<-lfcShrink(otuDEres,coef=coeff)
		otuDEresults<-data.frame(otuDEresults,"lfcShrink"=shrink[,2])
		#write.table(file="prova.de.txt",otuDEresults,sep="\t")
		otuDEsig<-otuDEresults[which(otuDEresults$padj <= pval & abs(otuDEresults$lfcShrink) >= log2fc), ]
		cat(dim(otuDEsig)[1],"tax significant, ",coeff, "tested\n")
		if (dim(otuDEsig)[1] > 0) {
			otuDEsig.ok = as(otuDEsig, "data.frame")
			tax<-as(tax_table(d)[rownames(otuDEsig), ], "matrix")
			f<-cbind(otuDEsig.ok,tax)
			write.table(file=paste0("plots/",ranks[i],".",name,".de.txt"),sep="\t",f)
			results[[paste0(name,".",ranks[i])]]<-f
			par(cex.lab=2,cex.main=1.5,cex.axis=1.5,mar=c(5.1,5.1,4.1,2.1))
			for (j in 1:dim(otuDEsig.ok)[1]) {
				#toplot<-otu_table(d)[rownames(otuDEsig.ok)[j],]
				toplot.prop<-otu_table(d.prop)[rownames(otuDEsig.ok)[j],]
				if(ranks[i] != 'OTU') {
					#title<-paste0(name,"\n",ranks[i]," - ",tax[j,ranks[i]],"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					#title<-paste0(ranks[i]," - ",tax[j,ranks[i]],"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					title<-paste0(ranks[i]," ",tax[j,ranks[i]],"\nlog2 FC ",round(otuDEsig.ok[j,"lfcShrink"],3))
				} else {
					#title<-paste0(name,"\n","OTU - ",paste(tax[j,],collapse="|"),"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					title<-paste0("OTU - ",paste(tax[j,],collapse="|"),"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
				}
				boxplot(as.numeric(toplot.prop)~sample_data(d.prop)$Surgery,main=title,ylab="Normalized OTU count")
				stripchart(as.numeric(toplot.prop)~sample_data(d.prop)$Surgery,pch=16,vertical=T,add=T,cex=1.5)
			}
		}
	}
	results
}

DEanalysis_pair<-function(data,ranks,pval,log2fc,zeros,name) {
	results<-list()
	for (i in 1:length(ranks)) {
		tot<-length(sample_names(data))
		cat(">>> Processing",ranks[i],"with",tot,"rows\n")
		if(ranks[i] != "OTU") {
			d<-tax_glom(data, taxrank = ranks[i], NArm = T)
		}
		else {
			d<-data
		}
		d.prop<-transform_sample_counts(d, function(otu) otu/sum(otu)*100)
		d<-filter_taxa(d,function(x) sum(ifelse(x==0,1,0))/length(x)< zeros,TRUE)
		d.prop<-filter_taxa(d.prop,function(x) sum(ifelse(x==0,1,0))/length(x)< zeros,TRUE)
		otuDE<-phyloseq_to_deseq2(d, ~ Patient+Disease)
		otuDEres<-DESeq(otuDE, test="Wald", fitType="parametric")
		otuDEresults<-results(otuDEres, cooksCutoff = FALSE)
		#coeff for shrink is set to the the last constrast, the same used by default in the "results" call
		coeff<-rev(colnames(coef(otuDEres)))[1]
		shrink<-lfcShrink(otuDEres,coef=coeff)
		otuDEresults<-data.frame(otuDEresults,"lfcShrink"=shrink[,2])
		#write.table(file="prova.de.txt",otuDEresults,sep="\t")
		otuDEsig<-otuDEresults[which(otuDEresults$padj <= pval & abs(otuDEresults$lfcShrink) >= log2fc), ]
		cat(dim(otuDEsig)[1],"tax significant, ",coeff, "tested\n")
		if (dim(otuDEsig)[1] > 0) {
			otuDEsig.ok = as(otuDEsig, "data.frame")
			tax<-as(tax_table(d)[rownames(otuDEsig), ], "matrix")
			f<-cbind(otuDEsig.ok,tax)
			write.table(file=paste0("plots/",ranks[i],".",name,".de.txt"),sep="\t",f)
			results[[paste0(name,".",ranks[i])]]<-f
			par(cex.lab=2,cex.main=1.5,cex.axis=1.5,mar=c(5.1,5.1,4.1,2.1))
			for (j in 1:dim(otuDEsig.ok)[1]) {
				#toplot<-otu_table(d)[rownames(otuDEsig.ok)[j],]
				toplot.prop<-otu_table(d.prop)[rownames(otuDEsig.ok)[j],]
				if(ranks[i] != 'OTU') {
					#title<-paste0(name,"\n",ranks[i]," - ",tax[j,ranks[i]],"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					#title<-paste0(ranks[i]," - ",tax[j,ranks[i]],"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					title<-paste0(ranks[i]," ",tax[j,ranks[i]],"\nlog2 FC ",round(otuDEsig.ok[j,"lfcShrink"],3))
				} else {
					#title<-paste0(name,"\n","OTU - ",paste(tax[j,],collapse="|"),"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
					title<-paste0("OTU - ",paste(tax[j,],collapse="|"),"\nlog2FC ",round(otuDEsig.ok[j,"log2FoldChange"],3)," & log2shFC ",round(otuDEsig.ok[j,"lfcShrink"],3))
				}
				boxplot(as.numeric(toplot.prop)~sample_data(d.prop)$Disease,border="white",main=title,ylab="Normalized OTU count")
				stripchart(as.numeric(toplot.prop)~sample_data(d.prop)$Disease,pch=16,vertical=T,add=T)
				dlogs<-log2(toplot.prop[1,11:20]/toplot.prop[1,1:10])
				c<-rep("black",tot/2)
				c<-ifelse(dlogs>=log2fc & dlogs>0,"orange",c)
				c<-ifelse(dlogs<=log2fc & dlogs<0,"blue",c)
				c<-ifelse(is.na(c),"black",c)
				dn<-sum(ifelse(c=="orange",1,0))
				up<-sum(ifelse(c=="blue",1,0))
				segments(rep(1,tot/2),as.numeric(toplot.prop)[11:20],rep(2,tot/2),as.numeric(toplot.prop)[1:10],col=c)
				legend("topleft",legend=c(up,dn),text.col=c("blue","orange"),bty="n",cex=2)
			}
		}
	}
	results
}


