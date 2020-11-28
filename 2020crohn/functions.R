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
