# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import cowplot
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import Seurat
## usethis namespace: end
NULL


###################################
## colour by cell type
###################################


#' Visualisation to colour the UMAP points by cell type label
#'
#' Colour UMAP by Cell Type.
#' @param seurat_object object with UMAP coordinates
#' @param cell_type_col column name of cell type
#' @param cell_cell_connections data.frame structured like the example “cell_cell_connections”.
#' @param ligand_cell_type do you want to limit the curves drawn to only starting at one cell type of interest? For example only cDC2 cells. Defaults to NA
#' @param self_interactions do you want to see loops for when cells are interacting with themselves? Defaults to F
#' @param label_cell_type do you want text labels for the cell types? Defaults to T
#' @keywords UMAP
#' @export
#' @examples
#' plot_gene_UMAP_color_bycelltype()
#'
plot_gene_UMAP_color_bycelltype<-function(seurat_object, cell_type_col, cell_cell_connections, ligand_cell_type = NA, self_interactions = F, label_cell_type = T){
  # extract UMAP coordinates
  umap_mat<-as.data.frame(Embeddings(object = seurat_object, reduction = "umap"))#
  umap_mat$cell<-rownames(umap_mat)
  meta<-seurat_object@meta.data
  meta$cell<-rownames(meta)
  plt<-merge(meta, umap_mat, by="cell")

  if("umap_1"%in%colnames(plt)){colnames(plt)[which(colnames(plt)=="umap_1")]<-"UMAP_1"}
  if("umap_2"%in%colnames(plt)){colnames(plt)[which(colnames(plt)=="umap_2")]<-"UMAP_2"}


  ## only ligands in one cell type
  if(is.na(ligand_cell_type)){}else{cell_cell_connections<-cell_cell_connections[which(cell_cell_connections$Cell1==ligand_cell_type),]}

  ## Cell type centroids
  plt_median<-plt %>% group_by(!!sym(cell_type_col)) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
  plt_median<-as.data.frame(plt_median)

  adjust_lab_by<-diff(range(plt_median$mean_umap2))*0.05

  cell_cell_connections_median<-merge(cell_cell_connections,plt_median, by.x="Cell1", by.y=cell_type_col)
  colnames(cell_cell_connections_median)[which(colnames(cell_cell_connections_median)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")

  cell_cell_connections_median<-merge(cell_cell_connections_median,plt_median, by.x="Cell2", by.y=cell_type_col)
  colnames(cell_cell_connections_median)[which(colnames(cell_cell_connections_median)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


  ## remove self interactions or if keeping draw circles at the points
  if(self_interactions==F){
    cell_cell_connections_median<-do.call(rbind,lapply(1:nrow(cell_cell_connections_median), function(x) if(cell_cell_connections_median$Cell1[x]==cell_cell_connections_median$Cell2[x]){}else{cell_cell_connections_median[x,]}))
  }else{
    cell_cell_connections_median_self<-cell_cell_connections_median[which(cell_cell_connections_median$Cell1 == cell_cell_connections_median$Cell2),]
    cell_cell_connections_median<-do.call(rbind,lapply(1:nrow(cell_cell_connections_median), function(x) if(cell_cell_connections_median$Cell1[x]==cell_cell_connections_median$Cell2[x]){}else{cell_cell_connections_median[x,]}))
    ## self loops
    radius <- 0.5

    circle_data<-do.call(rbind, lapply(1:nrow(cell_cell_connections_median_self), function(x){
      center_x <- cell_cell_connections_median_self$Cell1x[x]-(radius)
      center_y <- cell_cell_connections_median_self$Cell1y[x]
      data.frame(
        angle = seq(0, 2 * pi, length.out = 100),
        x = center_x + radius * cos(seq(0, 2 * pi, length.out = 100)),
        y = center_y + radius * sin(seq(0, 2 * pi, length.out = 100)),
        celltype=cell_cell_connections_median_self$Cell2[x])}))
  }


  ## axes
  len_x_bar<-((range(plt$UMAP_1))[2]-(range(plt$UMAP_1))[1])/10
  len_y_bar<-((range(plt$UMAP_2))[2]-(range(plt$UMAP_2))[1])/10
  arr <- list(x = min(plt$UMAP_1), y = min(plt$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

  umap_network<-ggplot() +
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=8,hjust = 0.05),
                       axis.title.y = element_text(size=8,hjust = 0.05,angle = 90),
                       legend.position = "none")+
    geom_point(aes(UMAP_1,UMAP_2), data=plt, size = 0.6, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2, color=!!sym(cell_type_col)), data=plt,size=0.5)+
    annotate("rect", xmin = min(plt$UMAP_1) * 1.1, xmax = max(plt$UMAP_1) * 1.1, ymin = min(plt$UMAP_2) * 1.1, ymax = max(plt$UMAP_2) * 1.1, fill = "white", alpha = 0.75)+
    geom_point(aes(mean_umap1,mean_umap2, fill=!!sym(cell_type_col)), data=plt_median,size=2, shape=21)+
    geom_curve(
      data = cell_cell_connections_median,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y),
      arrow = arrow(length = unit(0.02, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3, #-0.3,
      lineend = "round") +
    annotate("segment",
             x = arr$x, xend = arr$x + c(arr$x_len, 0),
             y = arr$y, yend = arr$y + c(0, arr$y_len), linewidth=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt')))+
    xlab("")+ylab("")+
    annotate("text",  x = arr$x+(arr$x_len/2), y = arr$y-adjust_lab_by, label="UMAP 1", size=3)+
    annotate("text",  x = arr$x-adjust_lab_by, y = arr$y+(arr$y_len/2), label="UMAP 2", size=3, angle = 90)
  umap_network

  if(self_interactions==T){
    umap_network<-umap_network+
      geom_path(aes(x = x, y = y, group=celltype), circle_data)+
      geom_point(aes(x = Cell1x-0.005, y = Cell1y+adjust_lab_by/4), cell_cell_connections_median_self, shape=25, size=1.5, fill="black")}
  if(label_cell_type==T){umap_network<-umap_network+geom_text(aes(mean_umap1,mean_umap2+adjust_lab_by, label=!!sym(cell_type_col)), data=plt_median,size=3)}
  umap_network}







###############
## Color by receptor ligand expression
###############
#' Visualisation to colour the UMAP points by expression of two genes
#'
#' Colour UMAP by expression of two genes. The gene don't have to be receptors and ligands just a typical use case.
#' @param seurat_object object with UMAP coordinates
#' @param cell_type_col column name of cell type
#' @param cell_cell_connections data.frame structured like the example “cell_cell_connections”.
#' @param ligand_cell_type do you want to limit the curves drawn to only starting at one cell type of interest? For example only cDC2 cells. Defaults to NA
#' @param self_interactions do you want to see loops for when cells are interacting with themselves? Defaults to F
#' @param label_cell_type do you want text labels for the cell types? Defaults to T
#' @param receptor gene name of the receptor (or any gene).
#' @param ligand gene name of the ligand (or any gene)
#' @param percentile the percentile cut off to call a gene as “highly expressed” and color a cell as having that gene “on”. Defaults to 0.8
#' @keywords UMAP
#' @export
#' @examples
#' plot_gene_UMAP_exp_colored()
plot_gene_UMAP_exp_colored<-function(seurat_object,cell_type_col, cell_cell_connections,ligand_cell_type=NA, self_interactions=F, label_cell_type=T, receptor, ligand, percentile=0.8){

  # extract UMAP coordinates
  umap_mat<-as.data.frame(Embeddings(object = seurat_object, reduction = "umap"))#
  umap_mat$cell<-rownames(umap_mat)
  meta<-seurat_object@meta.data
  meta$cell<-rownames(meta)
  plt<-merge(meta, umap_mat, by="cell")

  if("umap_1"%in%colnames(plt)){colnames(plt)[which(colnames(plt)=="umap_1")]<-"UMAP_1"}
  if("umap_2"%in%colnames(plt)){colnames(plt)[which(colnames(plt)=="umap_2")]<-"UMAP_2"}

  ## only ligands in one cell type
  if(is.na(ligand_cell_type)){}else{cell_cell_connections<-cell_cell_connections[which(cell_cell_connections$Cell1==ligand_cell_type),]}

  ## Cell type centroids
  plt_median<-plt %>% group_by(!!sym(cell_type_col)) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
  plt_median<-as.data.frame(plt_median)

  adjust_lab_by<-diff(range(plt_median$mean_umap2))*0.05

  cell_cell_connections_median<-merge(cell_cell_connections,plt_median, by.x="Cell1", by.y=cell_type_col)
  colnames(cell_cell_connections_median)[which(colnames(cell_cell_connections_median)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")

  cell_cell_connections_median<-merge(cell_cell_connections_median,plt_median, by.x="Cell2", by.y=cell_type_col)
  colnames(cell_cell_connections_median)[which(colnames(cell_cell_connections_median)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


  ## remove self interactions or if keeping draw circles at the points
  if(self_interactions==F){
    cell_cell_connections_median<-do.call(rbind,lapply(1:nrow(cell_cell_connections_median), function(x) if(cell_cell_connections_median$Cell1[x]==cell_cell_connections_median$Cell2[x]){}else{cell_cell_connections_median[x,]}))
  }else{
    cell_cell_connections_median_self<-cell_cell_connections_median[which(cell_cell_connections_median$Cell1 == cell_cell_connections_median$Cell2),]
    cell_cell_connections_median<-do.call(rbind,lapply(1:nrow(cell_cell_connections_median), function(x) if(cell_cell_connections_median$Cell1[x]==cell_cell_connections_median$Cell2[x]){}else{cell_cell_connections_median[x,]}))
    ## self loops
    radius <- 0.5

    circle_data<-do.call(rbind, lapply(1:nrow(cell_cell_connections_median_self), function(x){
      center_x <- cell_cell_connections_median_self$Cell1x[x]-(radius)
      center_y <- cell_cell_connections_median_self$Cell1y[x]
      data.frame(
        angle = seq(0, 2 * pi, length.out = 100),
        x = center_x + radius * cos(seq(0, 2 * pi, length.out = 100)),
        y = center_y + radius * sin(seq(0, 2 * pi, length.out = 100)),
        celltype=cell_cell_connections_median_self$Cell2[x])}))
  }

  ## axes
  len_x_bar<-((range(plt$UMAP_1))[2]-(range(plt$UMAP_1))[1])/10
  len_y_bar<-((range(plt$UMAP_2))[2]-(range(plt$UMAP_2))[1])/10
  arr <- list(x = min(plt$UMAP_1), y = min(plt$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

  ## gene expression
  gene_exp<-FetchData(seurat_object, vars=c(receptor, ligand))
  gene_exp$cell<-rownames(gene_exp)
  plt<-merge(plt, gene_exp, by='cell')


  ### non zero expression over the 90th percentile
  exp_limit1<-quantile(plt[, which(colnames(plt)==receptor)], percentile)
  exp_limit2<-quantile(plt[, which(colnames(plt)==ligand)], percentile)

  both_gene_over<-plt[which(plt[, which(colnames(plt)==receptor)]>=exp_limit1 &
                              plt[, which(colnames(plt)==receptor)]>0 &
                              plt[, which(colnames(plt)==ligand)]>=exp_limit2 &
                              plt[, which(colnames(plt)==ligand)]>0),]
  both_gene_over$color<-"Both Highly Expressed"

  gene1_over<-plt[which(plt[, which(colnames(plt)==receptor)]>=exp_limit1 &
                          plt[, which(colnames(plt)==receptor)]>0),]
  gene1_over<-gene1_over[which(!(gene1_over$cell%in%both_gene_over$cell)),]
  gene1_over$color<-paste(receptor,"Highly Expressed")

  gene2_over<-plt[which(plt[, which(colnames(plt)==ligand)]>=exp_limit2 &
                          plt[, which(colnames(plt)==ligand)]>0),]
  gene2_over<-gene2_over[which(!(gene2_over$cell%in%both_gene_over$cell)),]
  gene2_over$color<-paste(ligand,"Highly Expressed")

  plt_point_color<-rbind(both_gene_over, gene1_over, gene2_over)
  plt_point_color<-plt_point_color[order(plt_point_color$cell),]


  ## Network UMAP overlay
  umap_network<-ggplot() +
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=8,hjust = 0.05),
                       axis.title.y = element_text(size=8,hjust = 0.05,angle = 90),
                       legend.position = "none")+
    geom_point(aes(UMAP_1,UMAP_2),plt, size = 1, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2),plt, color="grey",size=0.75)+
    geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=0.5)+
    scale_color_manual(values=c("#88a000","#b80783","#03008e"))+
    annotate("rect", xmin = min(plt$UMAP_1) * 1.1, xmax = max(plt$UMAP_1) * 1.1, ymin = min(plt$UMAP_2) * 1.1, ymax = max(plt$UMAP_2) * 1.1, fill = "white", alpha = 0.75)+
    geom_point(aes(mean_umap1,mean_umap2, fill=!!sym(cell_type_col)), data=plt_median,size=2, shape=21)+
    geom_curve(
      data = cell_cell_connections_median,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y),
      arrow = arrow(length = unit(0.02, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3, #-0.3,
      lineend = "round") +
    annotate("segment",
             x = arr$x, xend = arr$x + c(arr$x_len, 0),
             y = arr$y, yend = arr$y + c(0, arr$y_len), linewidth=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt')))+
    xlab("")+ylab("")+
    annotate("text",  x = arr$x+(arr$x_len/2), y = arr$y-adjust_lab_by, label="UMAP 1", size=3)+
    annotate("text",  x = arr$x-adjust_lab_by, y = arr$y+(arr$y_len/2), label="UMAP 2", size=3, angle = 90)
  umap_network

  if(self_interactions==T){
    umap_network<-umap_network+
      geom_path(aes(x = x, y = y, group=celltype), circle_data)+
      geom_point(aes(x = Cell1x-0.005, y = Cell1y+adjust_lab_by/4), cell_cell_connections_median_self, shape=25, size=1.5, fill="black")}
  if(label_cell_type==T){umap_network<-umap_network+geom_text(aes(mean_umap1,mean_umap2+adjust_lab_by, label=!!sym(cell_type_col)), data=plt_median,size=3)}
  umap_network

  ## Cell type node legend
  nice_legend<-get_leg(ggplot()+
                         geom_point(aes(mean_umap1,mean_umap2, fill=!!sym(cell_type_col)),
                                    data=plt_median,size=3, shape=21, color="white")+
                         theme_void()+
                         theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                               plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))+
                         guides(fill=guide_legend(ncol=2)))

  ## gene expression legend
  exp_legend<-get_leg(ggplot()+
                        geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=2, alpha=0.5)+
                        scale_color_manual(values=c("#88a000","#b80783","#03008e"),
                                           name=paste("Gene Expression\n(Higher than ", percentile*100, "th\npercentile for gene",sep=""))+
                        theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                              plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))

  interaction_legend<- ggplot()+
    annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
    annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
    annotate("text", label=ligand, x = 1, y = 1, color="black")+
    annotate("text", label=receptor, x = 2, y = 1,color="black")+
    annotate("segment",
             x = 1.3, xend = 1.7 ,
             y = 1, yend = 1, linewidth=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)

  plot_grid(umap_network,
            plot_grid(interaction_legend,nice_legend, exp_legend,ncol=1,rel_heights = c(0.2,1,0.5)),
            rel_widths = c(2,1))
}



######################
## Format cellphonedb output
######################
#' Restructure output of cellphonedb
#'
#' This function formats cpdb data into a simplified cell-cell connection dataframe
#' @param cpdb_out Output from cellphonedb. The object named something like "statistical_analysis_significant_means.txt"
#' @param receptor gene name of the receptor
#' @param ligand gene name of the ligand
#' @keywords UMAP
#' @export
#' @examples
#' cell_cell_format_cpdb()
#'
cell_cell_format_cpdb<-function(cpdb_out, receptor, ligand){
  cell_cell_connections_plt<-melt(cpdb_out, id=colnames(cpdb_out)[c(1:12)])
  cell_cell_connections_plt$variable<-as.character(cell_cell_connections_plt$variable)

  cell_cell_connections_plt$Cell1<-sapply(1:nrow(cell_cell_connections_plt), function(x) strsplit(cell_cell_connections_plt$variable[x], "[.]")[[1]][1])
  cell_cell_connections_plt$Cell2<-sapply(1:nrow(cell_cell_connections_plt), function(x) strsplit(cell_cell_connections_plt$variable[x], "[.]")[[1]][2])
  cell_cell_connections_plt<-cell_cell_connections_plt[which(!(is.na(cell_cell_connections_plt$value))),]

  cell_cell_connections_plt[which(cell_cell_connections_plt$gene_a==ligand & cell_cell_connections_plt$gene_b==receptor), c("Cell1","Cell2", "gene_a","gene_b","value")]}


######################
## Format CellChat output
######################
#' Restructure output of CellChat
#'
#' This function formats CellChat data into a simplified cell-cell connection dataframe
#' @param cell_chat_out Output from CellChat. The object is the out put of the "subsetCommunication" function
#' @param receptor gene name of the receptor
#' @param ligand gene name of the ligand
#' @keywords UMAP
#' @export
#' @examples
#' cell_cell_format_cellchat()
#'
cell_cell_format_cellchat<-function(cell_chat_out, receptor, ligand){
  colnames(cell_chat_out)[which(colnames(cell_chat_out)=="source")]<-"Cell1"
  colnames(cell_chat_out)[which(colnames(cell_chat_out)=="target")]<-"Cell2"
  colnames(cell_chat_out)[which(colnames(cell_chat_out)=="ligand")]<-"gene_a"
  colnames(cell_chat_out)[which(colnames(cell_chat_out)=="receptor")]<-"gene_b"
  colnames(cell_chat_out)[which(colnames(cell_chat_out)=="prob")]<-"value"
  cell_chat_out[which(cell_chat_out$gene_a==ligand & cell_chat_out$gene_b==receptor), c("Cell1","Cell2", "gene_a","gene_b","value")]}


######################
## grab legened from plot
######################
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
