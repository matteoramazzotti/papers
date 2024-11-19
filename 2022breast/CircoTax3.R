library(ggplot2)
CircoTax3=function(file,title="CircoTax plot",ramp=c("orange","white","blue"),tax_col=7:11,fc_col=2,sort=c("no","rank","fc","absfc","alpha"),sort_dir="d") {
  data=file
  sort_dir=ifelse(sort_dir == "d",FALSE,TRUE)
  if(length(tax_col) == 5) {
    gplot_labels= unlist(strsplit("PCOFG",""))
  }
  if(length(tax_col) == 6) { #KPCOFG as in RDP
    gplot_labels= unlist(strsplit("KPCOFG",""))
  }
  if(length(tax_col) == 7) { #KPCOFGS as in silva
    gplot_labels= unlist(strsplit("KPCOFGS",""))
  }
  #build the taxa-related variables (tax index and label)
  #ranks are assumed to be in decreasing order from kinkdom(domain) to species
  #y represent the height (from the center of the circle) of the bars => domain=1 (min) species=7 (max)
  y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
  if(sort[1] == "no") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y
    fc=data[,fc_col]
  }
  if(sort[1] == "rank") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    o=order(y,decreasing=sort_dir)
    synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "alphalin") {
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    o=order(synth,decreasing=sort_dir)
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "alpha") {
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    o=order(labels)
    labels=labels[o]
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "fc") {
    o=order(data[,fc_col],decreasing=sort_dir)
    synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=data[o,fc_col]
  }
  if(sort[1] == "absfc") {
    #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    o=order(abs(data[,fc_col]),decreasing=sort_dir)
    synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
    synth=synth[order(abs(data[,fc_col]),decreasing=sort_dir)]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    fc=data[o,fc_col]
  }
  
  #builds the data.frame for ggplot2 
  df=data.frame("id"=1:dim(data)[1],"name"=labels,"rank"=y,"FC"=fc)
  
  #adds the label angle column
  nbar=dim(df)[1]
  angle=90-360*((1:nbar)-0.5)/nbar
  angle<-ifelse(angle < -90, angle+180, angle-10)
# cat(angle,"\n")
  
  #the plot starts here
  ggplot(df, aes(x = id, y = rank)) +
    ggtitle(title) +
    geom_col(
      aes(fill=FC),
      position = "dodge"
    ) + 
    scale_fill_gradient2(
      low="blue",
      mid="white",
      high="orange"
    ) +
    annotate(
      "text",label=gplot_labels, x=rep(0,length(tax_col)), y=1:length(tax_col),size=4.5, fontface=2
    ) +
    geom_hline(
      yintercept=1:6,
      color="grey"
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.margin = unit(rep(1,4), "cm"),      # Adjust the margin to make in sort labels are not truncated!
      plot.title = element_text(hjust = 0.5,face="bold", size=18)
    ) +
    coord_polar(start = 0, clip="off") +
    geom_text(
      aes(
        x= id, 
        y=7, 
        label=name,
      ), 
      color="black", 
      fontface="bold",
      alpha=0.6,
      size=3,
      angle=angle,
    ) +   
	labs(fill="log2FC")
}
