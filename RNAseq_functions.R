###########################################RNAseq_functions.R#########################
my_plot_PCA <- function(dds, set_name,intgroup, be,res_dir=DE_result_dir, out="jpg"){
  
  vstcounts_exp <- varianceStabilizingTransformation (dds, blind=TRUE)
  print(DESeq2::plotPCA(vstcounts_exp, intgroup=intgroup) + 
          labs(title=paste(set_name,'PCA'))+
          geom_label_repel(aes(label = dds@colData@rownames), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50'))
  if (out=="jpg"){
    dev.copy(jpeg,filename=paste0(res_dir,"/PCA_",set_name, ".jpg"), 
             width=1000, height=1200)
  } else{
    dev.copy2pdf(device = quartz, file =paste0(res_dir,"/PCA_",set_name, ".pdf"))
  }
  dev.off()
  
  if (!missing(be)){
    print(paste("be:",be))
    # remove batch effect
    assay(vstcounts_exp) <- 
        limma::removeBatchEffect(assay(vstcounts_exp), vstcounts_exp@colData[,be])
    print(DESeq2::plotPCA(vstcounts_exp, intgroup=intgroup) + 
            labs(title=paste(set_name,'PCA- BE removed')) +
            geom_label_repel(aes(label = dds@colData@rownames), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50'))
    if (out=="jpg"){
      dev.copy(jpeg,filename=paste0(res_dir,"/PCA_nbe_",set_name,".jpg"), 
             width=1000, height=1200)
    } else{
      dev.copy2pdf(device = quartz, file =paste0(res_dir,"/PCA_nbe_",set_name, ".pdf"))
    }
    dev.off()
    return(vstcounts_exp)
  }
}

plot_samples_dist_heatmap <- function(dds, vst_nbe, set_name,res_dir=DE_result_dir){
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  vst <- varianceStabilizingTransformation (dds, blind=FALSE)
  sampleDists <- dist(t(assay(vst)))
  
  pheatmap(as.matrix(sampleDists), border_color = NA,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, main="Samples Distance VST Heatmap")
  dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,"_vst_samples_dists", ".jpg"), 
           width=1000, height=1200)
  dev.off()
  
  # remove batch effect
  if (!missing(vst_nbe)){
    print("vst_nbe")
    sampleDists <- dist(t(assay(vst_nbe)))
    
    pheatmap(as.matrix(sampleDists),border_color = NA,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors, main="Samples Distance VST Heatmap - batch effect removed")
    dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,"_vst_samples_dists_nbe", ".jpg"), 
             width=1000, height=1200)
    dev.off()
    
    # vst is now the normalized log2-transformed data
    pheatmap(cor(assay(vst_nbe)), main="Samples Correlation VST Heatmap - batch effect removed")
    dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,"_vst_samples_cor_nbe", ".jpg"), 
             width=1000, height=1200)
    dev.off()
																				   
    
  }
 # rld <- rlog(dds)
#  sampleDists <- dist( t( assay(rld) ) )
#  pheatmap(as.matrix(sampleDists),
#           clustering_distance_rows=sampleDists,
#           clustering_distance_cols=sampleDists,
#           col=colors, main="Samples Distance RLD Heatmap")
#  dev.copy(jpeg,filename=paste0(DE_result_dir,"/heatmap_",set_name,"_rld_samples_dists", ".jpg"), 
#           width=800, height=1000)
#  dev.off()
  
  
}

plot_genes_heatmap <- function(rld,
                               anno_df = NA,
                               set_name,
                               anno_row,
                               cluster, 
                               cols_dend = TRUE, 
                               res_dir=DE_result_dir,
                               y_name="genes", 
                               gap_c = NULL, 
                               gap_r = NULL,
                               set_colnames=TRUE,
                               special_color,
                               an_colors,
                               scl="row",
                               plot_title,
                               out_format="jpeg"){
  if (missing(plot_title)) plot_title<-set_name
  #colors<-colorRampPalette(rev(heat.colors(50)))(255)
  suerat_palette<- Custom_palette(low = "#313695" , high = "#A50026", mid = "#FFFFBF", k = 50)
  suerat_palette2<- Custom_palette(low = "magenta", high = "yellow", mid = "black", k = 50)
  suerat_palette3<- Custom_palette(low = "white", high = "black", mid = NULL, k = 50)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
  my_colors<-colorRampPalette(rev(brewer.pal(n = 8, name="RdYlBu")))(100)
  
  
  if (!missing(special_color)){colors=special_color}
  
  if (out_format=="pdf"){
    pdf(file=paste0(res_dir,"/heatmap_",set_name,"_", y_name,".pdf"),width=7, height=7)  
    }
  #anno_df <- as.data.frame(colnames(rld))
  if (missing(anno_row)){
    if (!cols_dend){
      print("missing anno_row, cols_dend=F")
      
     # print("anno_row is shitttt ")
      
      print(paste("gap_c:",gap_c,"gap_r",gap_r))
      p_title<-paste(plot_title,y_name,"(",nrow(rld),")")
      pheatmap(rld, fontsize=9,annotation_col=anno_df, scale = scl, border_color = NA,  annotation_colors = list(group=c(BL= '#296de3', YL="#FFFF66")),# NL = "#778899"
               fontsize_number = 0.6 * fontsize, col=my_palette, 
               gaps_col =  gap_c,gaps_row =  gap_r,
               cluster_cols = FALSE, show_colnames= set_colnames,
               main = p_title)
      #heatmap.2(as.matrix(rld),scale="row", 
      #          trace="none", Rowv=T, Colv=F ,col= colors,
      #          main=paste(plot_title,y_name,nrow(rld),")-zscore"))
      
    }else{
      print("missing anno_row, cols_dend=T")
      p_title<-paste(plot_title,y_name,"(",nrow(rld),")")
      pheatmap(rld, fontsize=9,annotation_col=anno_df, scale = scl,
               fontsize_number = 0.6 * fontsize, col=my_colors,cluster_cols = FALSE,
               main = p_title)
      }
  }else if (missing(cluster)){
    print("anno_row, no cluster")
    p_title<-paste(plot_title,y_name,"(",nrow(rld),")")
    pheatmap(rld, annotation_col=anno_df, scale = scl, border_color = NA,
             annotation_row=anno_row, show_colnames= set_colnames,
             fontsize=24,fontsize_row = 8, fontsize_number = 34,fontsize_col = 24,
             col=my_colors, cluster_cols = cols_dend, cluster_rows = F,
             main = p_title)
  }else {
    print("anno_row+cluster")
    if (missing(an_colors)){
      anno_colors = rainbow(length(cluster))
      names(anno_colors) <- cluster
      anno_colors <- list(FlowSOM = anno_colors)}
    else  {
      names(an_colors) <- cluster
      anno_colors <- list(FlowSOM = an_colors)
    }
    
    print(anno_colors)
    
    if (cols_dend){
      print("cols_dend=T")
      p_title<-paste(plot_title,set_name,y_name,"(",nrow(rld),")")
      pheatmap(rld, fontsize=9, border_color = NA,
             scale = scl,annotation_col=anno_df,
             annotation_row=anno_row, cluster_rows = FALSE, annotation_colors = anno_colors,
             fontsize_number = 0.6 * fontsize, col=my_colors, show_colnames= set_colnames,
             main = p_title)
    } else{
      #order by column name b/a
      print("cols_dend=F")
      #rld<-rld[ ,order(substr(names(rld),6,6))]
        p_title<-paste(plot_title,set_name,y_name,"(",nrow(rld),")")
        pheatmap(rld,  border_color = NA,
               scale = scl,annotation_col=anno_df,cluster_cols = FALSE,
               annotation_row=anno_row, cluster_rows = FALSE, annotation_colors = anno_colors,
               fontsize=14,fontsize_row = 8, fontsize_number = 34,fontsize_col = 14,
               gaps_col =  gap_c, gaps_row =  gap_r, # border_color = "black",
               #cellheight = 10, cellwidth = 20,
               col=my_colors, show_colnames= set_colnames,
               main = p_title)
    }
  }
  gene_no <-nrow(rld)
  sample_no <-ncol(rld)
  if (out_format=="jpeg"){
    if (gene_no<20){
      wid<-max((sample_no*25)+str_length(rownames(rld))*3,str_length(p_title))*2
      dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,".jpg"), 
               width=wid, height=(gene_no*20)+150)
    }else if (gene_no<50){
      wid<-max((sample_no*25)+str_length(rownames(rld))*5,str_length(p_title))*2
      dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,".jpg"), 
               width=wid, height=(gene_no*10)+170)
    }else if (gene_no<250){
      wid<-max((sample_no*35)+str_length(rownames(rld))*5,str_length(p_title))*2
      dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,".jpg"), 
               width=wid, height=(gene_no*8)+150)
    }else {
      wid<-max((sample_no*50)+str_length(rownames(rld))*5,str_length(p_title))*4
      dev.copy(jpeg,filename=paste0(res_dir,"/heatmap_",set_name,".jpg"), 
               width=wid, height=(gene_no*5)+150)
    }
  }
  if (out_format!="")   dev.off()
									
		   
}


plot_volcano <- function(res,set_name, res_dir=DE_result_dir, labeled=TRUE, element="gene",
                         lfc=1, taxa=taxa_16S_df, label_list,label_list2, output_format="jpg",vtitle,
                         quantile_val=0){
  #res<-res_treat
  #set_name<-proj
  #lfc=1
  #label_list<-spec_genes
  #res_dir<-DE_result_dir
  if(missing(vtitle)) vtitle<-set_name
  res_4plot <- as.data.frame(res) %>% 
    rownames_to_column(var="gene") %>%    #save rowname as col
    mutate(genelabels = "") %>%           #add empty col
    arrange(pvalue) %>%                   #order by pvalue
    mutate(sig = pvalue <= 0.05 & abs(log2FoldChange) >= lfc) # add logical vector where TRUE values

  # add short name to OTUs
  sig_no <- nrow(subset(res_4plot,sig == TRUE))
#browser()  
  if (!labeled) {
    print("!labeled")
    print(ggplot(res_4plot) +
            geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = sig)) +
            ggtitle(paste(set, set_name,"\n significant","(LFC>=",lfc,")",element,"number:",sig_no)) +
            xlab("log2 fold change") + 
            ylab("-log10 p-value") +
            geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue),label = genelabels)) +
            theme(legend.position = "none",
                  plot.title = element_text(size = rel(1.5), hjust = 0.5),
                  axis.title = element_text(size = rel(1.25)))  )
    dev.copy(jpeg,filename=paste0(res_dir,"/volcano_sig_","_",set_name,".jpg"), 
             width=600, height=600)
    dev.off()
  } else {
    print(paste("no of sig:", sig_no))
    if (missing(label_list)) {
      res_4plot$genelabels <- ifelse(res_4plot$padj<=0.05,res_4plot$gene,NA)
      padj_no<-nrow(subset(res_4plot,!is.na(res_4plot$genelabels)))
      print(paste("no of genes padj<=0.05:",padj_no))
    } else {
      if (missing(label_list2)){
        res_4plot<-res_4plot %>%
                    mutate(genelabels = ifelse(gene %in% label_list,gene,NA)) %>%
                    mutate(sig = ifelse(gene %in% label_list,"spec",sig)) %>%
                    arrange(sig,genelabels)
      } else {
        res_4plot<-res_4plot %>%
                    mutate(genelabels = ifelse(gene %in% label_list,gene,
                                               ifelse(gene %in% label_list2,gene,NA))) %>%
                    mutate(sig = ifelse(gene %in% label_list,"pos",
                                        ifelse(gene %in% label_list2,"neg",sig))) %>%
                    arrange(sig,genelabels)
      }

    }   
#browser()
    gg<-ggplot(res_4plot) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = sig))  + 
      scale_color_manual(values=c("#b5b5b5","#CA2259","#258CB5","#25b531"))+
      #scale_color_manual(values=c("#999999", "#d10420", "#56B4E9"))+
      xlab("log2 fold change") + 
      ylab("-log10 p-value") +
      geom_label_repel(aes(x= log2FoldChange, y= -log10(pvalue),label = genelabels),
                           min.segment.length= 1, max.overlaps=100) +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            text=element_text(size = 16,  family="sans") ) 
      
    if (output_format=="pdf") {
      res_4plot<-dplyr::filter(res_4plot,log2FoldChange!=min(log2FoldChange))
      gg<-ggplot() +
        geom_point(data = res_4plot[res_4plot$sig != 'spec',], aes(x = log2FoldChange, y = -log10(pvalue), colour = sig), size = 2)  + 
        geom_point(data = res_4plot[res_4plot$sig == 'spec',], aes(x = log2FoldChange, y = -log10(pvalue), colour = sig), size = 3)  + 
        scale_color_manual(values=c("#b5b5b5","#CA2259","#258CB5"))+
        xlab("log2 fold change") + 
        ylab("-log10 p-value") + xlim(c(-6,6)) +
        geom_text_repel(data = res_4plot[res_4plot$sig == 'spec',], aes(x = log2FoldChange, y = -log10(pvalue),label = gene), 
                        min.segment.length = unit(0, 'lines'),
                        point.padding = 0.5) +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              text=element_text(size = 16,  family="sans"))
      
      
      ggsave(paste0(res_dir,"/volcano_labeled_sig_",set_name,".pdf"), plot = gg,
             width = 7.5, height = 5, device = "pdf", dpi=300, units="in")
      save(res_4plot, file=paste0(home_dir,"/volcano_labeled_4Shlomi", ".Rdata"))
      
    }  else {
#browser()
      if (missing(label_list)) {
        print(gg +
                ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                              "\n Label on PADJ<=0.05 ",":",padj_no)))
        dev.copy(jpeg,filename=paste0(res_dir,"/volcano_labeled_sig_",set_name, ".jpg"), 
                 width=800, height=1000)
        dev.off()
      }else{
        if (missing(label_list2)) {
          print(gg +
                  ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                                "\n Label on Specific genes (",length(label_list),")")))
        }else {
          print(gg +
                  ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                                "\n Label on positive genes (",length(label_list),"), \nnegative (",length(label_list2),")")))
          
        }
        if (output_format=="tiff"){
          ggsave(paste0(res_dir,"/volcano_labeled_spec_",set_name, ".tif"), width = 6, height = 5, device = "tiff", dpi=300, units="in")
        }else{
          dev.copy(jpeg,filename=paste0(res_dir,"/volcano_labeled_spec_",set_name, ".jpg"), 
                   width=800, height=1000)
        }
          dev.off()
      }
      
      
    }   
    
  }
}



#volcano for pathways

plot_volcano_pathway <- function(res,set_name, res_dir=DE_result_dir, labeled=TRUE, element="gene",
                                 lfc=1, taxa=taxa_16S_df, label_list,label_list2, output_format="jpg",vtitle,
                                 quantile_val=0){
  #res<-res_treat
  #set_name<-proj
  #lfc=1
  #label_list<-spec_genes
  #res_dir<-DE_result_dir
  if(missing(vtitle)) vtitle<-set_name
  res_4plot <- as.data.frame(res) %>% 
    rownames_to_column(var="gene") %>%    #save rowname as col
    mutate(genelabels = "") %>%           #add empty col
    arrange(pvalue) %>%                   #order by pvalue
    mutate(sig = pvalue <= 0.05 & abs(log2FoldChange) >= lfc) %>% # add logical vector where TRUE values
    mutate(sig2=sig)
  

  
  
  num_of_labeled_genes<-sum(tolower(res_4plot$gene) %in% tolower(label_list))
  # add short name to OTUs
  sig_no <- nrow(subset(res_4plot,sig == TRUE))
  #browser()  
  if (!labeled) {
    print("!labeled")
    print(ggplot(res_4plot) +
            geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = sig)) +
            ggtitle(paste(set, "\n", set_name,"\n significant","(LFC>=",lfc,")",element,"number:",sig_no)) +
            xlab("log2 fold change") + 
            ylab("-log10 p-value") +
            geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue),label = genelabels)) +
            theme(legend.position = "none",
                  plot.title = element_text(size = rel(1.5), hjust = 0.5),
                  axis.title = element_text(size = rel(1.25)))  )
    dev.copy(jpeg,filename=paste0(res_dir,"/volcano_sig_","_",set_name,".jpg"), 
             width=600, height=600)
    dev.off()
  } else {
    if (missing(label_list)) {
      res_4plot$genelabels <- ifelse(res_4plot$padj<=0.05,res_4plot$gene,NA)
      padj_no<-nrow(subset(res_4plot,!is.na(res_4plot$genelabels)))
      print(paste("no of genes padj<=0.05:",padj_no))
    } else {
      if (missing(label_list2)){
        res_4plot<-res_4plot %>%
          mutate(genelabels = ifelse(tolower(gene) %in% tolower(label_list),gene,NA)) %>%
          mutate(sig = ifelse(tolower(gene) %in% tolower(label_list),"spec",sig)) %>%
          arrange(sig,genelabels)
      } else {
        res_4plot<-res_4plot %>%
          mutate(genelabels = ifelse(tolower(gene) %in% tolower(label_list),gene,
                                     ifelse(tolower(gene) %in% tolower(label_list2),gene,NA))) %>%
          mutate(sig = ifelse(tolower(gene) %in% tolower(label_list),"pos",
                              ifelse(tolower(gene) %in% tolower(label_list2),"neg",sig))) %>%
          arrange(sig,genelabels)
      }
      
    }   
    #browser()
    labs=ifelse(res_4plot$sig == "spec" & res_4plot$sig2, res_4plot$genelabels, NA)
    #browser()  
    gg<-ggplot(res_4plot) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = sig))  + 
      scale_color_manual(values=c("#b5b5b5","#CA2259","#258CB5","#25b531"))+
      #scale_color_manual(values=c("#999999", "#d10420", "#56B4E9"))+
      xlab("log2 fold change") + 
      ylab("-log10 p-value") +
      geom_label_repel(aes(x= log2FoldChange, y= -log10(pvalue),label = labs),
                       min.segment.length= 1, max.overlaps=100) +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            text=element_text(size = 16,  family="sans") )
    
    
    if (output_format=="pdf") {
      res_4plot<-dplyr::filter(res_4plot,log2FoldChange!=min(log2FoldChange))
      gg<-ggplot() +
        geom_point(data = res_4plot[res_4plot$sig != 'spec',], aes(x = log2FoldChange, y = -log10(pvalue), colour = sig), size = 2)  + 
        geom_point(data = res_4plot[res_4plot$sig == 'spec',], aes(x = log2FoldChange, y = -log10(pvalue), colour = sig), size = 3)  + 
        scale_color_manual(values=c("#b5b5b5","#CA2259","#258CB5"))+
        xlab("log2 fold change") + 
        ylab("-log10 p-value") + xlim(c(-6,6)) +
        geom_text_repel(data = res_4plot[res_4plot$sig == 'spec',], aes(x = log2FoldChange, y = -log10(pvalue),label = gene), 
                        min.segment.length = unit(0, 'lines'),
                        point.padding = 0.5) +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              text=element_text(size = 16,  family="sans"))
      
      
      ggsave(paste0(res_dir,"/volcano_pathway_",vtitle,".pdf"), plot = gg,
             width = 7.5, height = 5, device = "pdf", dpi=300, units="in")
      save(res_4plot, file=paste0(home_dir,"/volcano_labeled_4Shlomi", ".Rdata"))
      
    }  else {
      #browser()
      if (missing(label_list)) {
        print(gg +
                ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                              "\n Label on PADJ<=0.05 ",":",padj_no)))
        dev.copy(jpeg,filename=paste0(res_dir,"/volcano_labeled_sig_",set_name, ".jpg"), 
                 width=800, height=1000)
        dev.off()
      }else{
        if (missing(label_list2)) {
          print(gg +
                  ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                                "\n colored genes in pathway (",num_of_labeled_genes,")"), subtitle = set_name))
        }else {
          print(gg +
                  ggtitle(paste(set, vtitle, "\n significant ","(LFC>=",lfc,") number:",sig_no,
                                "\n Label on positive genes (",num_of_labeled_genes,"), \nnegative (",length(label_list2),")"), subtitle = set_name))
          
        }
        if (output_format=="tiff"){
          ggsave(paste0(res_dir,"/volcano_pathway_",vtitle, ".tif"), width = 6, height = 5, device = "tiff", dpi=300, units="in")
        }else{
          dev.copy(jpeg,filename=paste0(res_dir,"/volcano_pathway_",vtitle, ".jpg"), 
                   width=800, height=1000)
        }
        dev.off()
      }
      
      
    }   
    
  }
}







plot_k <- function(data,set,max_k=40){
  #elbow plot to Determining Optimal Clusters
  
  set.seed(123)
  # Set maximum cluster 
    data_scaled <- scale(data, center = TRUE, scale = TRUE)
  wss <- (nrow(data_scaled)-1)*sum(apply(data_scaled,2,var))
  for (i in 2:max_k) wss[i] <- sum(kmeans(data_scaled,
                                          centers=i)$withinss)
  plot(1:max_k, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  dev.copy(jpeg,filename=paste0(DE_result_dir,"/elbow_",set, ".jpg"), 
           width=1000, height=700)
  dev.off()
  
}

plot_pval_hist <- function(res,set_name){
  
  colori <- c(`not passed`="khaki", `passed`="powderblue")
  
  use <- res$baseMean > metadata(res)$filterThreshold
  h_no_use <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h_use <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)

  barplot(height = rbind(h_no_use$counts, h_use$counts), beside = FALSE,
          col = colori, space = 0, ylab="frequency",
          main = paste("p-values results of DEseq  threshold:",
                       round(metadata(res)$filterThreshold,digit=2),
                       "\n passed:",length(res$pvalue[use]),
                       "not passed" ,length(res$pvalue[!use]) )) 
  text(x = c(0, length(h_no_use$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
  legend("top", fill=rev(colori), legend=rev(names(colori)))
  dev.copy(jpeg,filename=paste0(DE_result_dir,"/pval_hist_",set_name,".jpg"), 
           width=600, height=600)
}

plot_gene_cor <- function(df,set_name,rectan=2,mth="pearson",res_dir=DE_result_dir,element="genes"){
  #df<-spec_mat
  #filter out some not interesting genes
  df<-df %>%
        rownames_to_column("gene") %>%
        dplyr::filter(!str_detect(gene,"^Gm")) %>%
        dplyr::filter(!str_detect(gene,"Rik$")) %>%
        column_to_rownames("gene") %>%
        t %>%
        as.data.frame(stringsAsFactors=FALSE)#%>%
  print(paste("nrow df:",nrow(df)))
  cor_mat_sig <- corr.test(df, method=mth)
  if (rectan==0){
    corrplot::corrplot(cor_mat_sig$r, p.mat = cor_mat_sig$p, method = "square",#type = "lower",
                       insig = "label_sig", sig.level = c(.001, .01, .05),
                       title=paste("\nsignificant",element,"correlation", set_name),
                       pch.cex = .9, pch.col = "white", order="hclust")
    
  }else {
    corrplot::corrplot(cor_mat_sig$r, p.mat = cor_mat_sig$p, method = "square",#type = "lower",
                       insig = "label_sig", sig.level = c(.001, .01, .05),
                       title=paste("\nsignificant",element,"correlation", set_name),
                       pch.cex = .9, pch.col = "white", order="hclust", addrect = rectan)
    }
  dev.copy(jpeg,filename=paste0(res_dir,"/correlation_sig_",element,"_",mth,"_",set_name, ".jpg"), 
           width=ncol(df)*20+200, height=ncol(df)*20+250)
  dev.off()
}

  
RNA_sig_pairs <- function(run=NULL,RNA_thresh,RNA_all_results,use_pval=T,element="gene",labeled=T){
    
  wb_RNA <- openxlsx::createWorkbook()
  
  first_res<-TRUE
  all_pairs_genes<-character()
  for (res in RNA_all_results){
    #res <- RNA_all_results[2]
    print(res)
    res_df<- as.data.frame(get(res))
#browser()    
    ## Volcano plot t
    pv <- plot_volcano(res_df,paste0(set,"_",res),lfc = RNA_thresh, labeled=labeled)
    
    res_padj_sig <- subset(res_df, padj <= 0.05)
    print(paste(res, "padj sig:",nrow(res_padj_sig)))
    
    res_sig <- subset(res_df, pvalue <= 0.05)
    print(paste(res, "pvalue sig:",nrow(res_sig)))
    
    gene_res_padj_sig<-subset(res_padj_sig, abs(log2FoldChange) >=RNA_thresh)
    assign(paste0(res,"_padj_sig.",RNA_thresh), gene_res_padj_sig, envir = .GlobalEnv)
    print(paste("res_padj_sig_lfc.",RNA_thresh,nrow(gene_res_padj_sig)))
    
    # significats threshold 
    gene_res_pval_sig <- subset(res_sig, abs(log2FoldChange) >=RNA_thresh)
    assign(paste0(res,"_sig.",RNA_thresh), gene_res_pval_sig, envir = .GlobalEnv)
    print(paste("res_pvalue_sig_lfc.",RNA_thresh,nrow(gene_res_pval_sig)))
#browser()    
    if (use_pval) res_lfc<-gene_res_pval_sig else res_lfc<-gene_res_padj_sig
    #aggregate all sig genes from all the pairs
    if(!str_detect(res,"all") & str_detect(res,"CTRL")) all_pairs_genes<- unique(c(all_pairs_genes,rownames(res_lfc)))
      
    # Up padj significats threshold 
    gene_res_padj_sig_up <- subset(res_padj_sig, log2FoldChange >=RNA_thresh)
    assign(paste0(res,"_padj_up.",RNA_thresh), gene_res_padj_sig_up, envir = .GlobalEnv)
    print(paste("res_padj_up.",RNA_thresh,nrow(gene_res_padj_sig_up)))
    
    # Down padj significats threshold 
    gene_res_padj_sig_down <- subset(res_padj_sig, log2FoldChange <=-RNA_thresh)
    assign(paste0(res,"_padj_down.",RNA_thresh), gene_res_padj_sig_down, envir = .GlobalEnv)
    print(paste("res_padj_down.",RNA_thresh,nrow(gene_res_padj_sig_down)))
    
    # Up pval significats threshold 
    gene_res_pval_sig_up <- subset(res_sig, log2FoldChange >=RNA_thresh)
    assign(paste0(res,"_up.",RNA_thresh), gene_res_pval_sig_up, envir = .GlobalEnv)
    print(paste("res_pvalue_up.",RNA_thresh,nrow(gene_res_pval_sig_up)))
    
    # Down pval significats threshold 
    gene_res_pval_sig_down <- subset(res_sig, log2FoldChange <=-RNA_thresh)
    assign(paste0(res,"_down.",RNA_thresh), gene_res_pval_sig_down, envir = .GlobalEnv)
    print(paste("res_pvalue_down.",RNA_thresh,nrow(gene_res_pval_sig_down)))
    #save sig.1 to file for correlations with 16S
    #write.table(rownames(get(paste0(res,"_sig.",RNA_thresh))), 
    #            file = paste0(DE_result_dir,"/",set,"_",res,"_sig.",RNA_thresh,".txt"),  
    #            sep = "\t", row.names = F, col.names = FALSE, quote = FALSE)
    
#browser()    
    
    #------------------pairs up & down genes  
    ana_type<-ncol(str_split(res,"_",simplify = T))
    if ((ana_type==4 | ana_type==6) & !str_detect(res,"_y_")){
      st<-str_split_fixed(res,"\\_",4)[,2]
      nd<-str_split_fixed(res,"\\_",4)[,4]
    }else{
      st<-str_split_fixed(res,"\\_",5)[,3]
      nd<-str_split_fixed(res,"\\_",5)[,5]
    }
    #rld_pair_df<- get(paste0("rld_",nd))
    
    dds_count_sig<- dds_count_df %>%
        rownames_to_column("gene") %>%
        dplyr::filter(gene %in% rownames(res_lfc)) %>%
        column_to_rownames("gene")
#browser()    
    # keep the normalyzed counts of significant genes
    if(!str_detect(res,"all")){
      vst_pair_df<-get (str_replace(res,element,"vst"))
      anno_df<-get (str_replace(res,element,"anno"))
      
      #vst_pair_df<-vst_pair_df[, str_starts(colnames(vst_pair_df),st)|
      #                           str_starts(colnames(vst_pair_df),nd)]
      #dds_count_sig<-dds_count_sig[, str_starts(colnames(dds_count_sig),st)|
      #                               str_starts(colnames(dds_count_sig),nd)]
      # only the relevant rows
      #anno_df<-cond_df[str_starts(rownames(cond_df),st)|
      #                   str_starts(rownames(cond_df),nd),,drop = FALSE]
    }else{
      vst_pair_df<-vst_all
      anno_df<-cond_df
    }
    assign(paste0("dds_count_sig_",st,"_vs_",nd), dds_count_sig, envir = .GlobalEnv)
    
    # only the relevant cols
    print(paste("dds_count_sig:",nrow(dds_count_sig),ncol(dds_count_sig)))
    
    res_lfc <- res_lfc %>%
      arrange(lfcSE)
    gap_cols <- match(unique(anno_df[,1]), anno_df[,1])-1
    gap_cols <-gap_cols[2:length(gap_cols)]
    save(dds_count_sig, gene_res_pval_sig,res_lfc,
         file=paste0(home_dir,"/",set,"_",run,"_",res,"_sig.",RNA_thresh,".Rdata"))
    
    #res_sig<-get(paste0(res,"_sig_",thresh))
    sig_rld<-vst_pair_df[rownames(res_lfc),] 
    if (use_pval) {
      up <- rownames(gene_res_pval_sig_up)
      down <- rownames(gene_res_pval_sig_down)
    }else{
      up <- rownames(gene_res_padj_sig_up)
      down <- rownames(gene_res_padj_sig_down)
    }
    assign(paste0(res,"_sig_rld.",RNA_thresh), sig_rld, envir = .GlobalEnv)
    print(paste("sig_rld:",nrow(sig_rld)))
    save(sig_rld,anno_df,gap_cols,vst_pair_df,
         file=paste0(home_dir,"/rld_",set,"_",run,"_",res,"_sig.",RNA_thresh,".Rdata"))
    if (nrow(sig_rld)>0){
      # up genes 
      print(paste("up",length(up)))
      #wb_RNA <- openxlsx::createWorkbook()
      if (length(up))      up_rld<- up_or_down_genes(res ,"up", up,RNA_thresh, vst_pair_df, 
                                                     anno_df, wb=wb_RNA, res_df=res_df)
#browser()      
      # down genes 
      print(paste("down",length(down)))
      if (length(down))      down_rld<- up_or_down_genes(res ,"down", down,RNA_thresh, vst_pair_df, 
                                                         anno_df, wb=wb_RNA, res_df=res_df)
      
      all_sig <- rbind(down_rld,up_rld)
      nrow(all_sig)
      nrow(up_rld)
      nrow(down_rld)
      reg_df<- data.frame("up_down" =c(rep("up",nrow(up_rld)),rep("down",nrow(down_rld)))) %>%
        arrange(up_down)
      rownames(reg_df)<-row.names(all_sig)
#browser()      
      plot_genes_heatmap(all_sig,anno_df,paste0(set,"_",res,"_up_down_sig_",
                if(use_pval) "pval" else "padj",".",RNA_thresh),reg_df)
      
      
      if(str_detect(res,"all")){
        print("mean heatmap")
        all_sig_mean_df <- all_sig %>%
          rownames_to_column("gene")  %>%
          #dplyr::filter(!str_detect(gene,"^Gm" ) &  !str_detect(gene,"Rik$" ))%>%
          pivot_longer(cols=-gene) %>%
          separate(name, c("a", "id"),  "_", extra = "drop") %>%
          dplyr::select(-a) %>%
          group_by(id,gene) %>%
          mutate(mean=mean(value,na.rm=TRUE)) %>%
          dplyr::select(-value) %>%
          distinct()  %>%
          pivot_wider(names_from = id, values_from=mean)   %>%
          column_to_rownames("gene")
        
        gap_rows<-head(as.numeric(cumsum(table(reg_df$up_down))), -1)
        plot_genes_heatmap(rld=all_sig_mean_df,
                           anno_df=NA, #gap_c=gap_cols,
                           cols_dend =FALSE,
                           #scl="none",
                           anno_row=reg_df,
                           gap_r = gap_rows,
                           set_colnames=T,
                           set_name=paste0(set,"_",res,"_up_down_sig_",
                                           if(use_pval) "pval" else "padj",".",RNA_thresh,"_mean"))

      }
    }
    #}
    
  }
  assign(paste0("all_pairs_genes.",RNA_thresh,"_",run), all_pairs_genes, envir = .GlobalEnv)
  openxlsx::saveWorkbook(wb_RNA, 
                         file=paste0(home_dir,"/",set,"_",run,"_sig_gene_annotation_",RNA_thresh,".xlsx"), overwrite = TRUE)
}
  
  
  
compare_sig_genes<-function(set_name=NULL,list1,label_list1){
  pdf(file=paste0(DE_result_dir,"/venn_compare_sig_genes_",set_name, ".pdf"),width=7, height=7) 
  wb_comp <- openxlsx::createWorkbook()
  list2<-list1
  for (res1 in list1){
    list2<-list2[list2 != res1]
    for (res2 in list2){
      if (res1!=res2){
        print(paste(res1,res2))
        #res1<-list1[1]
        #res2<-list1[2]
        par(mfrow=c(2,3))
        res1_sig_df<-as.data.frame(get(paste0(res1,"_sig.",thresh_pair)))
        res2_sig_df<-as.data.frame(get(paste0(res2,"_sig.",thresh_pair)))
        venn(c( setNames(list(rownames(res1_sig_df)),
                         paste0(res1,"_sig.",thresh_pair)), 
                setNames(list(rownames(res2_sig_df)),
                         paste0(res2,"_sig.",thresh_pair)) ))
        
        venn(c( setNames(list(get(paste0(res1,"_up.",thresh_pair))),
                         paste0(res1,"_up.",thresh_pair)), 
                setNames(list((get(paste0(res2,"_up.",thresh_pair)))),
                         paste0(res2,"_up.",thresh_pair)) ))
        venn(c( setNames(list(get(paste0(res1,"_up.",thresh_pair))),
                         paste0(res1,"_up.",thresh_pair)), 
                setNames(list((get(paste0(res2,"_down.",thresh_pair)))),
                         paste0(res2,"_down.",thresh_pair)) ))
        venn(c( setNames(list(get(paste0(res1,"_down.",thresh_pair))),
                         paste0(res1,"_down.",thresh_pair)), 
                setNames(list((get(paste0(res2,"_up.",thresh_pair)))),
                         paste0(res2,"_up.",thresh_pair)) ))
        venn(c( setNames(list(get(paste0(res1,"_down.",thresh_pair))),
                         paste0(res1,"_down.",thresh_pair)), 
                setNames(list((get(paste0(res2,"_down.",thresh_pair)))),
                         paste0(res2,"_down.",thresh_pair)) ))
#browser()        
        res1_df<-as.data.frame(get(res1)) %>%
                  rownames_to_column("gene")
        res2_df<-as.data.frame(get(res2)) %>%
          rownames_to_column("gene")
        res1_df_lfc<- subset(res1_sig_df,select="log2FoldChange")
        res2_df_lfc<- subset(res2_sig_df,select="log2FoldChange")
        res_df_both <- merge(res1_df_lfc, res2_df_lfc, by="row.names") %>%
          column_to_rownames("Row.names") %>%
          setnames(c(paste0("lfc_",res1),paste0("lfc_",res2)))
        if(missing(label_list1))  label_list1<-intersect(rownames(res1_df_lfc),rownames(res2_df_lfc))
        label_list2<-c("null")
        
        lfc_merged<- res1_df %>%
          full_join(res2_df,by="gene") %>%
          dplyr::rename(log2FoldChange1=log2FoldChange.x, 
                        log2FoldChange2=log2FoldChange.y,
                        pval1=pvalue.x,
                        pval2=pvalue.y)  %>%
          mutate_at(vars(!contains("gene")),funs(replace(., is.na(.), 0)))  %>%
          # add sig, "darkorange2" both, green-H, magenta-Drug, grey-none
          mutate(significant=ifelse(abs(log2FoldChange1) >=thresh_pair & abs(log2FoldChange2) >= thresh_pair &
                                      pval2<0.05 & pval1<0.05, "both", 
                                    ifelse(abs(log2FoldChange2) >=thresh_pair & pval2<0.05,eval(res2), 
                                           ifelse(abs(log2FoldChange1) >=thresh_pair & pval1<0.05,eval(res1), "none"))))  %>%
          # add label to both sig
          mutate(genelabels = ifelse(gene %in% label_list1, gene,
                                     ifelse(gene %in% label_list2,
                                            gene,NA))) %>%
          column_to_rownames("gene")
#browser()        
        par(mfrow=c(1,1))
        print(ggscatter(lfc_merged,x="log2FoldChange1", y="log2FoldChange2", 
                        #ylim = c(-5,5),
                        main = paste("Scatter plot LFC of",res1,"and",res2),
                        shape  = 20, 
                        size  = ifelse(lfc_merged$significant=="none",0.2,2.5), 
                        frame = FALSE,
                        #label = "genelabels", repel = TRUE, font.label=c(7, "bold", "black"),  label.rectangle = T,
                        color = "significant", palette = c("darkorange2", "green", "magenta", "grey"))+
                xlab(str_replace(res1,"gene_y_","LFC "))+
                ylab(str_replace(res2,"gene_y_","LFC "))+
                geom_label_repel(aes(label=genelabels),color="black", max.overlaps=Inf, size =2, 
                                 box.padding = unit(0.2, "lines"))+
                geom_hline(yintercept=0) +
                geom_vline(xintercept=0) )#+ 
        
        #plot(res_df_both[,1], res_df_both[,2], 
        #     main = paste("Scatter plot LFC of ",res1,"and LFC of",res2),
        #     xlab = paste("LFC",res1), ylab = paste("LFC",res2),
        #     pch = 20, cex = .1, frame = FALSE)
        #lines(x = c(min(res_df_both[,1]),max(res_df_both[,1])), 
        #      y = c(min(res_df_both[,2]),max(res_df_both[,2])))
        
        
        
        sheet_name<-str_replace_all(paste0(res1,"+",res2),"gene_","")
        addWorksheet(wb_comp,sheet_name )
        comm_list<-intersect(rownames(get(paste0(res1,"_sig.",thresh_pair))),
                             rownames(get(paste0(res2,"_sig.",thresh_pair))))
        anno_gene_list <- anno_genes(comm_list) %>%
          dplyr::rename("gene"="external_gene_name")%>%
          left_join(rownames_to_column(res_df_both,"gene"),by="gene")
        writeData(wb_comp,sheet_name , anno_gene_list, colNames = TRUE, startRow = 1) 
        setColWidths(wb_comp, sheet_name, cols = 1:10, widths = "auto")
        
        sheet_name<-str_replace_all(paste0(res1,"-",res2),"gene_","")
        addWorksheet(wb_comp,sheet_name )
        only_res1<-setdiff(rownames(get(paste0(res1,"_sig.",thresh_pair))),
                             rownames(get(paste0(res2,"_sig.",thresh_pair))))
        anno_gene_list <- anno_genes(only_res1) %>%
          dplyr::rename("gene"="external_gene_name")%>%
          left_join(rownames_to_column(res1_df_lfc,"gene"),by="gene")
        writeData(wb_comp,sheet_name , anno_gene_list, colNames = TRUE, startRow = 1) 
        setColWidths(wb_comp, sheet_name, cols = 1:10, widths = "auto")
#browser()        
        sheet_name<-str_replace_all(paste0(res2,"-",res1),"gene_","")
        addWorksheet(wb_comp,sheet_name )
        only_res2<-setdiff(rownames(get(paste0(res2,"_sig.",thresh_pair))),
                             rownames(get(paste0(res1,"_sig.",thresh_pair))))
        anno_gene_list <- anno_genes(only_res2) %>%
          dplyr::rename("gene"="external_gene_name")%>%
          left_join(rownames_to_column(res2_df_lfc,"gene"),by="gene")
        writeData(wb_comp,sheet_name , anno_gene_list, colNames = TRUE, startRow = 1) 
        setColWidths(wb_comp, sheet_name, cols = 1:10, widths = "auto")
        
      }
    }
  }
  dev.off()
  openxlsx::saveWorkbook(wb_comp, 
                         file=paste0(home_dir,"/",set_name,"_common_sig.",thresh_pair,"_genes.xlsx"), overwrite = TRUE)
}
  
  
  prepare_GSEA_rnk <- function(res_df,rnk_file,human,mouse){
    #res_df<-res_treat_df
    res_df$fcSign <- sign(res_df$log2FoldChange)
    res_df$logP=-log10(res_df$pvalue)
    res_df$metric=res_df$logP/res_df$fcSign
    res_df$mouse_gene<-rownames(res_df)
    print(paste("nrow(res_df):",nrow(res_df)))
    print(res_df)
    
     #convert genes from mouse to human for GSEA
    Mouse2HumanTable <- Mouse2Human(MouseGenes = res_df$mouse_gene,human=human,mouse=mouse)
    print(paste("nrow(Mouse2HumanTable):",nrow(Mouse2HumanTable)))
    res_df$human_gene<- 
      Mouse2HumanTable[match(res_df$mouse_gene, Mouse2HumanTable$MGI), "HGNC"] 
    print(paste("nrow(res_df):",nrow(res_df)))
    
    
    #human_GSEA<- subset(res_df,!(is.na(human_gene)| is.na(metric)), select= c("human_gene", "metric"))
    human_GSEA<- subset(res_df,!is.na(human_gene)& !is.na(metric)& str_length(human_gene)>0, 
                        select= c("human_gene", "metric"))
    #human_GSEA_LFC<- subset(res_df,!is.na(human_gene)& !is.na(log2FoldChange)& str_length(human_gene)>0, 
    #                    select= c("human_gene", "log2FoldChange"))
    print(paste("nrow(human_GSEA):",nrow(human_GSEA), 
                "unique nrow(human_GSEA):", length(unique(human_GSEA$human_gene))))
    #save rnk file for GSEA 
    
    nameF = paste0(home_dir,"/",rnk_file,"_GSEA.rnk")
    write.table(human_GSEA, file = nameF, sep = "\t", 
                row.names = F, col.names = FALSE, quote = FALSE)
    #write.table(human_GSEA_LFC, file = paste0(home_dir,rnk_file,"_LFC.rnk"), sep = "\t", 
    #            row.names = F, col.names = FALSE, quote = FALSE)
    return(human_GSEA)
}
  
Mouse2Human <- function(MouseGenes,human,mouse){
  if (missing(human)){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }
  print(paste("in Mouse2Human MouseGenes:",length(MouseGenes)))
  #browser()
  genesMousetoHuman = getLDS(attributes = c("external_gene_name","mgi_symbol"), 
                             filters = "external_gene_name", 
                             values = MouseGenes , 
                             mart = mouse, 
                             attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
                             martL = human, 
                             uniqueRows = TRUE)
  
  colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Human.Gene_ID", "HGNC")
  print(paste("in Mouse2Human nrow(genesMousetoHuman):",nrow(genesMousetoHuman)))
  
  return(genesMousetoHuman) 
  }
  
  
GSEA_results_file <- function(file_type="tsv",set_name,sub_dir_prefix,fdr=0.05){
    file_names <-list.files(home_dir, 
                            #pattern= glob2rx("*gsea_report_for_na*.xls$"),
                            pattern= glob2rx(paste0("*gsea_report_for_na*.",file_type,"$")),
                            #pattern= glob2rx(paste0("*",set,"*gsea_report_for_na*.tsv$")),
                            recursive = TRUE)
   if (!missing(sub_dir_prefix))       file_names <-file_names[startsWith(file_names, sub_dir_prefix)]
    duplicated(file_names)
    file_names[duplicated(file_names)]
    file_names <-as.data.frame(file_names) #%>%
      #dplyr::filter(str_detect(.,set))
    wb <- createWorkbook()
    

    for (fn in file_names$file_names) {
      #fn<- "CS_vs_CTRL_C2.GseaPreranked.1592844211821/gsea_report_for_na_neg_1592844211821.xls"
      print(fn)
      underscore_loc<-as.numeric(as.data.frame(str_locate_all(fn, "\\."))[1,1])
    
      de_set<-paste(substring(fn,1,underscore_loc-1), ifelse(str_detect(fn,"neg"),"neg","pos"))
      #print(de_set)
      
      GSEA_xls<-as.data.frame(read.delim(paste0(home_dir,"/",fn)))
      
      sheet_name <-de_set

      if (str_length(sheet_name)>31) sheet_name <-str_sub(de_set,9,str_length(de_set)-2)
      addWorksheet(wb, sheet_name)
      writeData(wb,sheet_name , subset(GSEA_xls,FDR.q.val<=fdr), rowNames = TRUE) 
      setColWidths(wb, sheet_name, cols = 1:15, widths = "auto")
      
    }
    
    #browser()
    saveWorkbook(wb, file=paste0(home_dir,"/GSEA_all_results_FDR.",fdr,"_",set,
                                 ifelse(missing(set_name),"",paste0("_",set_name)),".xlsx"), 
                                 overwrite = TRUE)
    
}  
  
  GSEA_terms_4words <- function(file_type="tsv",set_name,sub_dir_prefix){
 #browser()
    GSEA_word_list <- openxlsx::read.xlsx(paste0(home_dir,"/GSEA_word_list.xlsx"))
    file_names <-list.files(home_dir, 
                            pattern= glob2rx(paste0("*gsea_report_for_na*.",file_type,"$")),
                            recursive = TRUE)
    if (!missing(sub_dir_prefix))       file_names <-file_names[startsWith(file_names, sub_dir_prefix)]
    duplicated(file_names)
    file_names[duplicated(file_names)]
    
    # only WT files
    #file_names<- subset(file_names,str_detect(file_names,"WT"))
    term_data<-data.frame("NAME"=character(),"NES"=numeric(),"FDR q-val"=numeric())
    for (fn in file_names) {
      #fn<-file_names[1] 
      print(fn)
      fn_list<-unlist(str_split(fn,"_"))
      #underscore_df<-as.data.frame(str_locate_all(fn, "_"))
      #underscore_loc<-as.numeric(underscore_df[nrow(underscore_df)-5,1])
      #de_set<-substring(fn,9,underscore_loc-1)
      
      de_set<-paste0(fn_list[1],"_",fn_list[2],"_",fn_list[3],"_",fn_list[4])
      print(de_set)
      GSEA_xls<-as.data.frame(read.delim(paste0(home_dir,"/",fn)))
      
      
      for (word in unique(GSEA_word_list$word)){
        #   term<-"GO_CARTILAGE_DEVELOPMENT"
        word_found<-GSEA_xls %>%
          dplyr::filter(str_detect(NAME,word) & FDR.q.val<=0.05) %>% 
          dplyr::select("NAME","NES","FDR.q.val") %>%
          mutate(NES=as.numeric(NES))
        term_data<- rbind(term_data, 
                          cbind(word_found,
                                data.frame(word = rep(word,nrow(word_found))),
                                data.frame(de_set = rep(de_set,nrow(word_found)))))
      }
    }
    nrow(term_data)

    term_data_uni <- aggregate(word~ NAME+NES+FDR.q.val+de_set, term_data, 
                               FUN=paste, collapse= " ")
    nrow(term_data_uni)
 

    # 1 row per term
    term_data_cast <- data.table::dcast(setDT(term_data_uni),NAME+word~de_set,fun.aggregate= mean,
                                        value.var=c("NES","FDR.q.val"))
    nrow(term_data_cast)
    term_data_cast<- term_data_cast %>% replace(is.na(.), NA)
    print(term_data_cast)
    
    save(term_data_uni,term_data_cast, file=paste0(home_dir,"/term_data_4GSEA_",set,
                                                   ifelse(missing(set_name),"",paste0("_",set_name)), ".Rdata"))
    
    #term_data_cast$NAME[duplicated(term_data_cast$NAME)]
    
    sheet_name <-"Terms"
    wb <- createWorkbook()
    
    addWorksheet(wb, sheet_name)
    posStyle <- createStyle(fontColour = "#006100", fgFill = "#C6EFCE")
    to_color<- which( !is.na(term_data_cast$NES_VHOM_vs_CTRL) &
                        !is.na(term_data_cast$NES_VHOS_vs_CTRL) &
                        !is.na(term_data_cast$NES_CS_vs_CTRL))
    writeData(wb,sheet_name , term_data_cast, rowNames = F) 
    addStyle(wb,sheet_name, style =posStyle, rows = to_color+1, 
             cols = 1:ncol(term_data_cast)+1, gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 3:9, widths = "auto")
    saveWorkbook(wb, file=paste0(home_dir,"/GSEA_suggested_terms_",set,
                                 ifelse(missing(set_name),"",paste0("_",set_name)),".xlsx"), 
                 overwrite = TRUE)
    
  }
  
  GSEA_add_categories <- function(de_set_order="all", category_order, per_cat=10, 
                                  w=5200, h=6000,subs,set_name=NULL,comment="",output="jpg"){
    #per_cat=30
    #subs<-c("MOVEMENT","MORPHOGENESIS")
    #subs<-c("GASTRO","MUSCLE")
    #subs=c("IMMUNE","IMMUNE-CYTOKINE")
    #browser()
    load(paste0(home_dir,"/term_data_4GSEA_",set,
                ifelse(missing(set_name),"",paste0("_",set_name)), ".Rdata"),verbose = T)
    print(nrow(term_data_uni))
    #GSEA_category <- read_excel(paste0(home_dir,"GSEA_category2.xlsx"))
    #just for demonstration
    GSEA_category <- read_excel(paste0(home_dir,"/GSEA_terms_category_",set,
                                       ifelse(missing(set_name),"",paste0("_",set_name)),".xlsx"))
    GSEA_category <- GSEA_category %>%
      dplyr::select(NAME,category) %>%
      dplyr::filter(category!="NA")
    
    #GSEA_category$category<- paste0(GSEA_category$category,
    #                               ifelse(!is.na(GSEA_category$subcategory),
    #                                      paste(" -\n",GSEA_category$subcategory),""))
    #category_order<-c(#"Development and organogenesis", 
    #                  "Muscle development and functions",
    #                  "Neuronal functions",
    #                  "Extracellular matrix", 
    #                  "Protein glycosylation",
    #                  "Energy harvest and Metabolism",
    #                  "Protein degradation")
    
    if (!missing(subs)) {
      GSEA_category <- GSEA_category %>%
        dplyr::filter(str_detect(category,paste(subs,collapse = '|'))) %>%
        arrange(category)
    }
    if (missing(category_order)) if (missing(subs)) {
      category_order<-sort(unique(GSEA_category$category))
    } else{
      category_order<-subs
      GSEA_category <- GSEA_category %>%
        mutate(category=factor(category, levels = category_order)) 
    }
    
    print(table(GSEA_category$category)  )
    print(setdiff(GSEA_category$NAME,term_data_uni$NAME))
    term_data_uni<- merge(term_data_uni, GSEA_category, by="NAME") 
      
    
    term_data_cast_cat<- merge(term_data_cast, GSEA_category, by="NAME")
    print(table(term_data_uni$category)  )
    print(table(term_data_cast_cat$category)  )
    
    term_data_uni<-term_data_uni %>% 
      mutate(category=factor(category, levels = category_order)) %>%
      mutate(NAME=as.factor(NAME)) %>%
      mutate(NAME=str_pad(NAME,max(length(NAME)))) %>%
      arrange(category)
  
    if (de_set_order!="all"){
      
      term_data_uni_4plot <- term_data_uni %>%
        dplyr::filter(de_set %in% de_set_order)
      
      
    } else {
      de_set_order<-unique(term_data_uni$de_set)
      assign("term_data_uni_4plot",term_data_uni,envir = .GlobalEnv)
    }
    print(table(term_data_uni_4plot$category)  )
    
    
      
    #add artificial rows to empty de_set
    
    for (cat in unique(term_data_uni_4plot$category)){
      
      
      
      term_cat<-subset(term_data_uni_4plot, category==cat)
      empty_de_set<-setdiff(de_set_order,unique(term_cat$de_set))
      #if (length(empty_de_set)>1){
       # browser()
      #  print(empty_de_set)
      #}
        
      #browser()       
      if (length(empty_de_set)>0)
        if (empty_de_set!="all")
          for (ds in empty_de_set){
            term_data_uni_4plot<-rbind(term_data_uni_4plot,
                                       cbind(NAME=term_cat$NAME[1],NES=0,FDR.q.val=0,de_set=ds,word="NA",category=cat))%>%
              mutate(NES=as.numeric(NES)) %>%
              mutate(FDR.q.val=as.numeric(FDR.q.val)) %>%
              arrange(category)
          }
    }
    
    #browser()
    
    if (!missing(subs)) {
      set_name<-paste0(set_name,"_",paste(subs,collapse = '_'))
    }
    #browser()  
    if (output=="jpg"){
      #browser()  
      pg <- plot_GSEA(term_data_uni_4plot, name_font=50, per_cat=per_cat)
      dev.copy(jpeg,filename=paste0(DE_result_dir,"/GSEA_barplot2_",set,"_",set_name,"_",comment, ".jpg"), 
               width=w, height=h)
      dev.off()
      #ggsave(paste0(DE_result_dir,"/GSEA_barplot2_",set,"_",set_name,"_",comment,".pdf"), 
      #       width = w, height = h, device = "pdf", dpi=300, units="mm")
    }else{
      
      pg <- plot_GSEA(term_data_uni_4plot, name_font=3.5, per_cat=per_cat)
      
      
      dev.copy2pdf(device = quartz, file =paste0(DE_result_dir,"/GSEA_barplot2_",set,"_",set_name,"_",comment,".pdf"),
                   width = 6, height = 6)
      dev.off()
    }
    save(term_data_uni,term_data_cast_cat, file=paste0(home_dir,"/term_data_categ_4GSEA_",set_name,"_",set,".Rdata"))
    #openxlsx::write.xlsx (term_data_cast, paste0(home_dir,"/GSEA_terms_cat_",set,"_all",".xlsx"), row.names = F, col.names = T) 
    
    
    #-------------------------heatmap
  
    
    term_data_cast_4plot<- term_data_cast_cat %>%
      dplyr::filter(category %in% GSEA_category$category)  %>%
      mutate(category=factor(category, levels = category_order)) %>%
      arrange(category) %>%
      dplyr::select(NAME,starts_with("NES"))  %>%
      setNames(str_replace_all(names(.),"NES_", ""))  #%>%
    
    
    if (de_set_order!="all") term_data_cast_4plot<-term_data_cast_4plot%>%
      dplyr::select(NAME,de_set_order)  %>%
      relocate((de_set_order))  #%>%    #rev
    term_data_cast_4plot<-term_data_cast_4plot%>%
      column_to_rownames("NAME")  %>%
      mutate_all(~replace(., is.na(.), 0))
    print(nrow(term_data_cast_4plot))
    
    
    if(ncol(term_data_cast_4plot)>1)  {
      term_data_cast_4plot<-as.data.frame(term_data_cast_4plot[rowSums(term_data_cast_4plot != 0) >0, ] )
      
      print(nrow(term_data_cast_4plot))
      
      GSEA_category<-GSEA_category %>%
        #rownames_to_column("NAME")  %>%
        dplyr::filter(NAME %in% rownames(term_data_cast_4plot))  #%>%
      #column_to_rownames("NAME")  #%>%
      print(nrow(GSEA_category))
      
      #anno_df<- data.frame(NES=colnames(term_data_cast_4plot),anno="")%>%
      #            column_to_rownames("NES")  
      #browser()
      anno_r<-column_to_rownames(GSEA_category,"NAME")
      gap_rows<-head(as.numeric(cumsum(table(anno_r$category))), -1)
      
      plot_genes_heatmap(rld=term_data_cast_4plot,
                         anno_df=NA, #gap_c=gap_cols,
                         cols_dend =FALSE,
                         scl="none",
                         anno_row=anno_r,
                         gap_r = gap_rows,
                         cluster=unique(GSEA_category$category),
                         plot_title= paste("GSEA",proj,set,set_name,comment),
                         set_colnames=T,
                         y_name="Pathways",
                         set_name=paste0("GSEA_",proj,set,"_",set_name,"_",comment))
    }
  }
  
  plot_GSEA <- function(path_data,per_cat=10, name_font=60){

    
    #Loading objects:
     # term_data_uni
    #term_data_cast
    #[1] 1055
    #path_data<-term_data_uni
    #path_data<-term_data_uni_subset
    
    
    gg<-1
    my_plots = list()
    h_all<-numeric()
    cat_uni<-levels(path_data$category)
    
    
    print(paste("no of categiries:",length(levels(path_data$category))))
    print(table(path_data$category)  )
    
    for (cat in cat_uni){
      #cat<-cat_uni[1]
      print(cat)
      #browser()      
      path_data_names<-subset(path_data, category==cat) %>%
        group_by(de_set)  %>%
        arrange(desc(abs(NES))) %>%
        head(per_cat)
      path_data_p<-path_data %>%
        dplyr::filter(category==cat & NAME %in% path_data_names$NAME)
      print(paste("no of NAMEs:",nrow(path_data_p)))
      #browser()
      my_plots[[gg]]<- ggplot(path_data_p,aes(fill=FDR.q.val, y=NES, x=NAME))+
        geom_bar(position="dodge", stat="identity")+
        theme_grey(base_size = name_font) +
        #coord_cartesian(ylim = Ylims) +
        scale_y_continuous(limits=c(ifelse(min(path_data$NES)<0,min(path_data$NES),0), 
                                    ifelse(max(path_data$NES)>0,max(path_data$NES),0))) +
        scale_fill_continuous(limits=c(min(path_data$FDR.q.val), max(path_data$FDR.q.val)))+
        labs(fill = "FDR qvalue")+
        theme(axis.title.y= element_blank(), 
              axis.text.y= element_text(size = name_font,  family="sans"),
              legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank())+#, 
        #axis.line = element_line(colour = "black"))+
        facet_wrap(~de_set, strip.position = "top",drop=FALSE, scales = "fixed",nrow = 1)+
        scale_color_manual(values=c("#5B1031","#455F4A","#7A5E2E","#144756"))+
        coord_flip()
      #last one shall have y axis
      if(cat!= last(cat_uni)){
        my_plots[[gg]]<-my_plots[[gg]]+ theme(axis.ticks.x= element_blank(),
                                              axis.title.x= element_blank(), 
                                              axis.text.x= element_blank())
      } else {
        my_plots[[gg]]<-my_plots[[gg]]+ theme(axis.title.x= element_text(size = name_font),
                                              axis.text.x= element_text(size = name_font))
      }
      # first shall include titles and space for them
      if(cat!= cat_uni[1] ){ #| length(unique(path_data_p$de_set))==1
        my_plots[[gg]]<-my_plots[[gg]]+ theme(strip.background = element_blank() ,
                                              strip.text.x = element_text(size=0))
      } else {
        my_plots[[gg]]<-my_plots[[gg]]+ theme(strip.text.x = element_text(size=name_font))
      }
      
      #add cat as text. gpar Alpha channel for transparency
      ##### CHANGE CATEGORIES STYLE
      my_plots[[gg+1]] <- grobTree(rectGrob(gp=gpar(fill="grey90", alpha=0.5)),
                                   textGrob(cat,gp=gpar(fontsize=name_font+50)))
      
      gg<-gg+2
      h_all<- c(h_all,nrow(path_data_p)/nrow(path_data)*200)
    }
    #calculate plot hight
    print(length(h_all))
    print(h_all)
    
    print(c((length(unique(path_data$de_set))*150)+max(str_length(path_data$NAME),na.rm =T),max(str_length(cat_uni))*10))
    print(wrap_plots(my_plots,ncol=2,heights = round(h_all*35), 
                     widths=c((length(unique(path_data$de_set))*150)+max(str_length(path_data$NAME),na.rm =T),max(str_length(cat_uni))*10)) +
            plot_layout(guides = 'collect')& theme(legend.position = 'bottom'))
    #save(path_data, file=paste0(home_dir,"/GSEA_barplot_4Shlomi", ".Rdata"))
    return(my_plots)
    
  }

  
  up_or_down_genes<- function(res_name ,up_or_down, gene_list,thresh, rld_mtx,anno_df,wb,res_df){
    #gene_list<-up
    #res_name<-res
    #up_or_down<-"up"
 # browser()
    print(length(gene_list))
    gene_list_rld <- rld_mtx%>%
      as.data.frame(stringAsFactor=F)  %>%
      rownames_to_column("gene")  %>%
      dplyr::filter(gene %in% gene_list)  %>%
      column_to_rownames("gene")
    gene_list_rld <- gene_list_rld - rowMeans(gene_list_rld)
#browser()    
    plot_genes_heatmap(gene_list_rld,anno_df,paste0(set,"_",res_name,"_sig_",up_or_down,".",thresh))
    assign(paste0(res_name,"_",up_or_down,"_",thresh),gene_list, envir = .GlobalEnv)
    
    sheet_name <- paste0(res_name,".",up_or_down,".",thresh )
    print(sheet_name)
    if (str_length(sheet_name)>31) sheet_name<-substring(sheet_name,5)
    addWorksheet(wb, sheet_name)
    anno_gene_list <- anno_genes(gene_list) %>%
      dplyr::rename("gene"="external_gene_name")%>%
      left_join(rownames_to_column(res_df,"gene"),by="gene")
    writeData(wb,sheet_name , anno_gene_list, colNames = TRUE, startRow = 1) 
    setColWidths(wb, sheet_name, cols = 1:10, widths = "auto")
    return(gene_list_rld)
  }
  
  anno_genes<- function(gene_list){
   # if (!exists(paste0(home_dir,"../../mouse.Rdata"))) {
    #if (!exists("mouse_ensembl")) mouse_ensembl = useDataset("mmusculus_gene_ensembl",mart=useMart("ensembl"))
    if (!exists("mouse_ensembl"))  {load(paste0(home_dir,"../../mouse.Rdata"),verbose = T)}
    #listAttributes(mart=mouse_ensembl)
#browser()    
  result = tryCatch({
    anno_data<-
        getBM(attributes=c('external_gene_name','description','gene_biotype'), useCache = FALSE,
              filters = 'mgi_symbol', 
              values = gene_list, 
              mart = mouse_ensembl)
  },
  error = function (condition) {
    anno_data<-data.frame(external_gene_name=gene_list)
  })
  
  #return(anno_data)
  return(result)
  }
  
  
  
    
pathway_gene_heatmap <- function(rld_mtx,anno_df,res_name,inc_sig=F,res_sig,res_df,file_type="tsv",sub_dir_prefix,categ,min_genes=3){
    pair_name<-str_replace(res_name,"res_","")
    pair_name<-str_replace(pair_name,"gene_y_","")
    print(pair_name)
    file_names_all <-list.files(home_dir, 
                            pattern= glob2rx(paste0("*.",file_type)), #xls
                            recursive = TRUE)
    if (!missing(sub_dir_prefix))       file_names <-file_names[startsWith(file_names, sub_dir_prefix)]
    if (file_type=="tsv") GSEA_colname<-"SYMBOL" else GSEA_colname<-"PROBE"
#browser()
    #check if pathway is significant
    if (str_detect(pair_name,"all")){
          sig_path<- term_data_cast$NAME
          file_names<-file_names_all[str_detect(file_names_all,paste(sig_path,collapse = "|"))]
    }else {
          sig_path<-str_trim(subset(term_data_uni,str_detect(de_set,pair_name),"NAME")$NAME)
          file_names<-file_names_all[str_detect(file_names_all,paste(sig_path,collapse = "|")) & 
                                       str_detect(file_names_all,pair_name)]
      }
    
    if (!missing(categ)){
      GSEA_category_cat <- read_excel(paste0(home_dir,"GSEA_terms_category_",set,"_",run,".xlsx")) %>%
                      dplyr::select(NAME,category) %>%
                      dplyr::filter(category==categ)
      file_names<-file_names[str_detect(file_names,paste(GSEA_category_cat$category,collapse = "|"))]
    }
    print(paste("no of files:",length(file_names)))
    print(paste("no of sig path:",length(sig_path)))
#browser()    
    categ_genes<-data.frame(pw=character(),gene=character())
    for (fn in file_names){
      # fn<-file_names[1]
      print(fn)
      pathway<-str_replace(basename(fn),".tsv","")
      pw_xls<-as.data.frame(read.delim(paste0(home_dir,fn)))
      print(paste("excel lines:", nrow(pw_xls)))
#browser()
      pw_genes<- subset(pw_xls, CORE.ENRICHMENT=="Yes", eval(GSEA_colname)) #PROBE
      print(paste("excel enriched genes:",nrow(pw_genes)))

      #collect all genes in category
      categ_genes<- unique(rbind(categ_genes,cbind(pathway,pw_genes)))
      
      pw_genes_rld <- subset( rld_mtx,str_to_upper(rownames(rld_mtx)) %in% pw_genes[,eval(GSEA_colname)]) #PROBE
      print(paste("enriched genes found in set:",nrow(pw_genes_rld)))
      pw_genes_sig <- subset( res_sig,str_to_upper(rownames(res_sig)) %in% pw_genes[,eval(GSEA_colname)]) #PROBE
      nrow(pw_genes_sig)
      if (nrow(pw_genes_rld)>min_genes){
          pw_genes_rld<-
            pw_genes_rld[ ,order(match(colnames(pw_genes_rld),rownames(meta))) ]
          if (inc_sig) {
#browser()    
            if(!missing(res_df)){
              pv <- plot_volcano(res_df,
                                 set_name=paste(substr(pathway,1,50),"genes on", res_name,set,proj),
                                 vtitle=paste(substr(pathway,1,50),"genes","\non", res_name,set,proj),
                                 lfc= thresh_pair,label_list = rownames(pw_genes_sig))
            }
            plot_genes_heatmap(pw_genes_rld,anno_df, gap_c=gap_cols,
                               paste(proj,set,res_name," SIG- Pathway",substr(pathway,1,70)),
                               plot_title = paste(proj,set,res_name,"SIG- \nPathway",pathway),
                               cols_dend =FALSE,
                               set_colnames=T)
            if (!str_detect(res_name,"all")) {
#browser()
              tit<-paste(proj,set,res_name,"LFC SIG-Pathway \n",substr(pathway,1,70))
              print(ggplot(rownames_to_column(pw_genes_sig), aes(fill=pvalue, 
                                                                       y=log2FoldChange, 
                                                                       x=reorder(rowname,log2FoldChange)))+
                      geom_bar(position="dodge", stat="identity")+
                      #theme_grey(base_size = 20) +
                      theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"),
                            text=element_text(size = 18,  family="sans") )+ 
                      ggtitle(tit)+
                      theme(axis.title.y= element_blank())+#, 
                      coord_flip())
#browser()
              w<-max((max(pw_genes_sig$log2FoldChange)-min(res_sig$log2FoldChange))*3,str_length(tit)*3) 
              h<-if (nrow(pw_genes_sig)<10) nrow(pw_genes_sig)*20 else nrow(pw_genes_sig)*10
              ggsave(paste0(DE_result_dir,"/barplot_SIG_Pathway_",substr(pathway,1,70),"_",proj,set,res_name,".jpg"), 
                     width = w, height = h,  dpi=300, units="mm",limitsize = FALSE)
            }
          }else{
            plot_genes_heatmap(pw_genes_rld,anno_df, gap_c=gap_cols,
                               paste(proj,set,res_name,"- Pathway",substr(pathway,1,70)),
                               plot_title = paste(proj,set,res_name,"- \nPathway",pathway),
                               cols_dend =FALSE,
                               set_colnames=T)
          }
      }else {print("too little genes")}
    }
    save(categ_genes,file=paste0(home_dir,"/Top_sig_pathways_genes_4corr_",set,".Rdata"))
  }
  
  
  my_dimReduction<- function(data, markers = NULL,method = "tsne", tsneSeed = 42){
    data <- as.matrix(data)
    rnames <- row.names(data)
    out_dim <- 2
    if (is.null(markers)) {markers <- colnames(data)}
    marker_filtered_data <- data[, markers]
    print(ncol(marker_filtered_data))
    cat("  Running t-SNE...with seed", tsneSeed)
    set.seed(tsneSeed) # Set a seed if you want reproducible results
    tsne_out <- Rtsne(marker_filtered_data, initial_dims = ncol(marker_filtered_data), 
                      dims = 2, 
                      check_duplicates = FALSE, 
                      pca = TRUE)
    mapped <- tsne_out$Y
    if(!is.null(mapped)){
      if(ncol(mapped) < out_dim){
        out_dim <- ncol(mapped)
        message("Run ",method," for dimensional reduction, out dimension coerced to ",out_dim)
      }
      mapped <- mapped[ ,c(1:out_dim)]
      colnames(mapped) <- paste(method, c(1:out_dim), sep = "_")
      rownames(mapped) <- rnames
    }
    cat("  DONE\n")
    return(mapped)
    
  }
  
  
  cluster_function <- function(df, df_tsne, k, set,clusterColorName="Set1"){
    #df <- fcs_markers_immune4cluster_equal
    #df_tsne <- fcs_tsne_immune
    #k<-selected_k
    #set<-paste0("immune_",threshold,"_",cells_number)
    
    set.seed(123)
    cluster_FlowSOM <- my_cluster(xdata = df, ydata= df_tsne, method = "FlowSOM", FlowSOM_k = k)
    
    data_FlowSOM <- as.data.frame(cbind(as.matrix(df), df_tsne, FlowSOM=cluster_FlowSOM))
    print(paste("clusterColorName:",clusterColorName))
    #pdf(file=paste0(result_dir,"/FlowSOM_cluster_marker_",set,"_", k,"_",day,"_",clusterColorName, ".pdf"),width=7, height=7)  
    #  print(my_clusterPlot(data=data_FlowSOM, xlab="tsne_1", ylab="tsne_2", cluster="FlowSOM",
    #                       sampleLabel = FALSE, addLabel = TRUE, 
    #                       title = paste("FlowSOM marker",set, "k=",k)))
    print(my_clusterPlot(data=data_FlowSOM, xlab="tsne_1", ylab="tsne_2", cluster="FlowSOM",
                           sampleLabel = FALSE, #addLabel = FALSE, 
                           title = paste("FlowSOM marker",set, "k=",k),
                           clusterColorName=clusterColorName))
  
    dev.copy(jpeg,filename=paste0(DE_result_dir,"/FlowSOM_cluster_marker_",set,"_", k,"_", ".jpg"), 
             width=1000, height=700)
    dev.off()
    return(data_FlowSOM)
  }
  
  my_cluster<- function(xdata = NULL, ydata= NULL, method= c("Rphenograph", "FlowSOM"),
                        FlowSOM_k=30, Rphenograph_k=30, flowSeed = 123){
 
    xdata <- as.matrix(xdata)
    print(method)
    method = match.arg(method)
    if(method == "NULL"){
      return(NULL)
    }
    switch(method, 
           Rphenograph = {
             cat("  Running PhenoGraph...")
             clusters <- as.numeric(membership(Rphenograph(xdata, k=Rphenograph_k)))
           },
           FlowSOM = {
             cat("  Running FlowSOM...")
             set.seed(flowSeed)
             clusters <- FlowSOM_integrate2cytofkit(xdata, FlowSOM_k, flowSeed = flowSeed)
           })
    
    if( length(clusters) != ifelse(is.null(ydata), nrow(xdata), nrow(ydata)) ){
      message("Cluster is not complete, cluster failed, try other cluster method(s)!")
      return(NULL)
    }else{
      if(!is.null(xdata) && !is.null(row.names(xdata))){
        names(clusters) <- row.names(xdata)
      }else if(!is.null(ydata) && !is.null(row.names(ydata))){
        names(clusters) <- row.names(ydata)
      }
      cat(" DONE!\n")
      return(clusters)
    }
    print(length(cluster))
    return(cluster)
    #}
  }
  
  FlowSOM_integrate2cytofkit <- function(xdata, k, flowSeed = NULL, ...){
    cat("    Building SOM...\n")
    xdata <- as.matrix(xdata)
    
    ord <- tryCatch({
      map <- SOM(xdata, silent = TRUE,  ...) #
      cat("    Meta clustering to", k, "clusters...\n")
      metaClusters <- suppressMessages(metaClustering_consensus(map$codes, k = k, seed = flowSeed))
      cluster <- metaClusters[map$mapping[,1]]
    }, error=function(cond) {
      message("Run FlowSOM failed")
      message("Here's the error message:")
      message(cond)
      return(NULL)
    }) 
    
    if(is.null(ord)){
      cluster <- NULL
    }else{
      if(length(ord) != nrow(xdata)){
        message("Run FlowSOM failed!")
        return(NULL)
      }
      cluster <- ord
    }
    
    return(cluster)
  }
  
  my_clusterPlot <- function(data, xlab, ylab, cluster, sample, title = "cluster", 
                             type = 1, point_size = NULL, addLabel=TRUE, 
                             labelSize=8, sampleLabel=TRUE, 
                             labelRepel = FALSE, fixCoord=TRUE, clusterColorName="Set1") {
    if(!is.data.frame(data))
      data <- as.data.frame(data)
    
    if(missing(sample)){
      sample <- "sample"
      data$sample <- "data"
    }
    
    paraCheck <- c(xlab, ylab, cluster, sample) %in% colnames(data)
    if(any(!paraCheck)){
      stop("Undefined parameters found: ",
           c(xlab, ylab, cluster, sample)[!paraCheck])
    }
    
    data[[cluster]] <- as.factor(data[[cluster]])
    data[[sample]] <- as.factor(data[[sample]])
    cluster_num <- length(unique(data[[cluster]]))
    sample_num <- length(unique(data[[sample]]))
    col_legend_row <- ceiling(cluster_num/8)
    size_legend_row <- ceiling(sample_num/4)
    grid_row_num <- round(sqrt(sample_num))
    if (sample_num >= 8) {
      shape_value <- LETTERS[1:sample_num]
    } else {
      shape_value <- c(1:sample_num) + 15
    }
    if (is.null(point_size)) {
      point_size <- ifelse(nrow(data) > 10000, 0.85, 3)
    }
    if (clusterColorName !="hila"){
      getPalette <- colorRampPalette(brewer.pal(9, clusterColorName))
      clusterColor <- getPalette(max(cluster_num))
      
    } else {
      clusterColor <- c("#4a8ff7","#f56733","#fca030", "#f7e045", "#bed160", "#808080", "#8bb559","#72a37f", "#549e88", "#60b1d1", "#325363", "#035496",  "#325cbf", "#6e5399", "#775c82", "#bf91b7", "#e8b7d1", "#bababa", "#e590f0", "#8c6f4f", "#999267", "#64566e")       }

    if(type == 1){
      if(sampleLabel){
        cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster, shape = sample)) + 
          geom_point(size = point_size) + scale_shape_manual(values = shape_value) + 
          scale_colour_manual(values = clusterColor) + 
          xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Scatter Plot", sep = " ")) + 
          theme_bw() + theme(legend.position = "bottom") + 
          theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)), 
                 shape = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
      }else{
        cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster)) + 
          geom_point(size = point_size) + scale_shape_manual(values = shape_value) + 
          scale_colour_manual(values = clusterColor) + 
          xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Scatter Plot", sep = " ")) +
          theme_classic()+ 
          #theme_bw() + 
          theme(legend.position = "bottom") + 
          theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #,
          #      axis.line = element_line(colour = "black")) +
          guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4))) #
      }
      
      if(addLabel){
        edata <- data[ ,c(xlab, ylab, cluster)]
        colnames(edata) <- c('x', "y", "z")
        center <- aggregate(cbind(x,y) ~ z, data = edata, median)
        
        if(labelRepel && !sampleLabel){
          cp <- cp + geom_text_repel(data=center, aes_string(x = "x", y = "y", label = "z"),
                                     size = labelSize, fontface = 'bold', color = "black",
                                     box.padding = unit(0.5, 'lines'),
                                     point.padding = unit(1.6, 'lines'),
                                     segment.color = '#555555',
                                     segment.size = 0.5,
                                     arrow = arrow(length = unit(0.02, 'npc')))
        }else{
          cp <- cp + annotate("text", label = center[,1], x=center[,2], y = center[,3],
                              size = labelSize, colour = "black")
        }
      }
      
    }else if (type == 2){
      cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster)) +
        facet_wrap(~sample, nrow = grid_row_num, scales = "fixed") + 
        geom_point(size = point_size - 0.05 * sample_num) + 
        scale_colour_manual(values = clusterColor) + 
        xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Grid Plot", sep = " ")) + 
        theme_bw() + theme(legend.position = "bottom") + 
        theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
        guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)), 
               shape = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
      
      if(addLabel){
        edata <- data[ ,c(xlab, ylab, cluster)]
        colnames(edata) <- c('x', "y", "z")
        center <- aggregate(cbind(x,y) ~ z, data = edata, median) 
        
        if(labelRepel && !sampleLabel){
          cp <- cp + geom_text_repel(data=center, aes_string(x = "x", y = "y", label = "z"),
                                     size = labelSize, fontface = 'bold', color = "black",
                                     box.padding = 0.7,
                                     point.padding = unit(1.6, 'lines'),
                                     segment.color = '#555555',
                                     segment.size = 0.5,
                                     arrow = arrow(length = unit(0.02, 'npc')))
        }else{
          cp <- cp + geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), 
                               size = labelSize, colour = "black") 
        }
      }
      
    }else{ 
      stop("Undefined type, only 1 or 2.") 
    }
    
    if(fixCoord){
      cp <- cp + coord_fixed()
    }
    
																	   
																					  
										   
																	
																										   
																				   
   
   
										   
										   
    return(cp)
    
  }
  
 
  Custom_palette <- function(
    low = "white",
    high = "red",
    mid = NULL,
    k = 50
  ) {
    low <- col2rgb(col = low) / 255
    high <- col2rgb(col = high) / 255
    if (is.null(x = mid)) {
      r <- seq(from = low[1], to = high[1], len = k)
      g <- seq(from = low[2], to = high[2], len = k)
      b <- seq(from = low[3], to = high[3], len = k)
    } else {
      k2 <- round(x = k / 2)
      mid <- col2rgb(col = mid) / 255
      r <- c(
        seq(from = low[1], to = mid[1], len = k2),
        seq(from = mid[1], to = high[1], len = k2)
      )
      g <- c(
        seq(from = low[2], to = mid[2], len = k2),
        seq(from = mid[2], to = high[2],len = k2)
      )
      b <- c(
        seq(from = low[3], to = mid[3], len = k2),
        seq(from = mid[3], to = high[3], len = k2)
      )
    }
    return(rgb(red = r, green = g, blue = b))
  }
  
  appendRData <- function(robj, filename) {
    
    tmpEnv <- new.env()
    
    if (!file.exists(filename)) warning("'file' does not exist. A new file will be created.")
    savedObjects <- load(file = filename, envir = tmpEnv, verbose = TRUE)
    
    # quick check for name collisions
    stopifnot(!(deparse(substitute(robj)) %in% savedObjects))
    
    save(list = c(savedObjects, deparse(substitute(robj))),
         file = filename,
         envir = tmpEnv)
  }