spec_chr_plot <- function(genedata,mapdata,tissuedata,chr=8,tag_p=3,y_lab="p",name_gene="PENK",y_min=NULL,
                          y_max=NULL,y_ticks=NULL,cut_ctr=2,sig_line1=3,sig_line2=-3,color_tissue=NULL) {
  library(readr)
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(data.table)
  library(magrittr)
  library(tidyverse)
  library(pipeR)
  library(gtools)
  
  g<-genedata
  m<-mapdata
  
  tissue<-tissuedata
  
  m%<>%rename(gene=ENSG)
  m$CHR<-as.integer(m$CHR)
  g1<-g%>%filter(0<P & P<0.1)
  test2<-g%>%filter(P>0.1)%>%group_by(CHR)%>%bind_rows(g1)%>%select(-SNP)%>%
    mutate(Tis="Single SNP",
           BP=BP/1000000,
           P=log10(P))%>%filter(CHR!=23&CHR==chr)
  test<-tissue%>%left_join(m,by='gene')%>%mutate(BP=(START_POS+END_POS)/2,
                                                 P=-log10(pvalue))%>%select(CHR,BP,P,tissue,gene_name)%>%
    filter(CHR==chr)
  
  
  postart<-tissue%>%left_join(m,by='gene')%>%mutate(P=-log10(pvalue),
                                                    START_POS=START_POS/1000000,
                                                    END_POS=END_POS/1000000,
                                                    gene_type=paste(gene_name,tissue, sep = "_"))%>%filter(CHR==chr)%>%
    select(CHR,gene_name,tissue,START_POS,P,gene_type)%>%rename(BP=START_POS)
  
  poend<-tissue%>%left_join(m,by='gene')%>%mutate(P=-log10(pvalue),
                                                  START_POS=START_POS/1000000,
                                                  END_POS=END_POS/1000000,
                                                  gene_type=paste(gene_name,tissue, sep = "_"))%>%filter(CHR==chr)%>%
    select(CHR,gene_name,tissue,END_POS,P,gene_type)%>%rename(BP=END_POS)
  
  popopo<-bind_rows(postart,poend)
  xyn<-as.numeric(popopo[popopo$gene_name==name_gene,"BP"][1,1])
  popopo%<>%filter(BP>xyn- cut_ctr &BP<xyn+ cut_ctr)
  popopo<-bind_rows(postart,poend)%>%filter(BP>xyn-cut_ctr&BP<xyn+cut_ctr)
  potest<-test2%>%filter(BP>xyn- cut_ctr  &BP<xyn+cut_ctr)
  # Add highlight and annotation information
  #mutate( is_highlight=ifelse(SNP %in% mysnps, "yes", "no")) %>%
  #mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))
  
  #df2$CHR<-as.integer(df2$CHR)
  for_tag<-popopo%>%filter(P>tag_p)
  
  if(is.null(color_tissue) == TRUE) {
    #color_tissue =  c("firebrick", "blue", "goldenrod", "#D95F02", "#7570B3", "#E7298A", "forestgreen", "darkgray", "black", "dodgerblue")
    color_tissue =  c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered","blue4","darkturquoise","dodgerblue","cyan","lightsteelblue", 
                      "lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "lightpink1","gray32","grey85","goldenrod4","khaki3", 
                      "lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4","plum","chocolate4","paleturquoise4","palevioletred3","seagreen1",
                      "darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0")
    
  }
  else {
    color_tissue = color_tissue 
  }
  
  # get chromosome center positions for x-axis
  #axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  labels_cat <- c(sort(unique(as.character(test$tissue))))
  
  if(is.null(y_min) == TRUE) {y_min <- round(min(potest$P, na.rm=TRUE) - 1)} else {y_min <- y_min}
  if(is.null(y_max) == TRUE) {y_max <- round(max(popopo$P, na.rm=TRUE) + 1)} else {y_max <- y_max}
  if(is.null(y_ticks) == TRUE) {
    if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)} 
  } else {y_ticks <- y_ticks}
  
  break_length <- round(round(y_max - y_min)/y_ticks)
  
  genom<-ggplot() +theme_bw() +theme(axis.text.x = element_text(size=12, color = 'black'),
                                     axis.text.y = element_text(size=12, color = 'black'),
                                     axis.title.x = element_text(size = 12, face = "bold", color ="black"),
                                     axis.title.y = element_text(size = 12, face = "bold", color ="black"),
                                     axis.ticks.x=element_line())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab(paste("Chromosomes",chr,"(MB)", sep=" ")) +ylab(y_lab)+
    geom_label_repel(data = for_tag, aes(x = BP, y = P, label = gene_name))+
    geom_point(data = potest, aes(x = BP, y = P)) +
    geom_line(data = popopo, aes(x = BP, y = P, group=gene_type,size=1.5,color = factor(tissue)))+
    
    scale_colour_manual(name = "Tissue Type", values = c("black", color_tissue), labels = c("Single SNP", labels_cat))+
    guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) +
    
    #geom_label_repel(data = subset(new_data, log_p > gene_tag_p), aes(x = new_pos, y = log_p, label = Gene)) +
    #geom_label_repel(data = for_tag, aes(x = BP/1000000, y = P, label = factor(gene_name))) +
    geom_hline(aes(yintercept = 0), size = 3) +
    geom_hline(yintercept = sig_line1, size = .5)+
    geom_hline(yintercept = sig_line2, size = .5)+
    scale_y_continuous(breaks=seq(y_min, y_max, break_length)) + 
    expand_limits(y=c(y_min, y_max))
  
  return(genom)
}
