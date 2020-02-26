striat_query <- function(atlas, MNI_coords, afni.dir, mni_x, mni_y,mni_z) {
  net_arr      <- NULL
  for ( iter in 1:nrow(MNI_coords) ){
    print(iter)
    #coords       <- c(data_struct$striat_samples[iter,]$mni_x, data_struct$striat_samples[iter,]$mni_y, data_struct$striat_samples[iter,]$mni_z)
    use_coords <- MNI_coords[iter,c(mni_x,mni_y,mni_z)]

    coords     <- c(as.numeric(use_coords[mni_x]), as.numeric(use_coords[mni_y]), as.numeric(use_coords[mni_z]))
    
    network      <- query_atlas(atlas, coords, 0, afni_dir)
    write(network, '')
    # Check neighboring voxels
    if ( network == 0 ){
      neigh_arr <- query_atlas(atlas, coords, 1, afni_dir)
      num_zeros <- length(neigh_arr[neigh_arr %in% 0])
      if ( num_zeros == 27 ){
        neigh_arr <- query_atlas(atlas, coords, 2, afni_dir)
        num_zeros <- length(neigh_arr[neigh_arr %in% 0])
        if ( num_zeros == 125 ){
          network <- 0 
        } else {
          use       <- neigh_arr[neigh_arr>0]
          use_table <- summary(as.factor(use))
          network   <- as.integer(names(rev(sort(use_table))))
        }
      } else {
        use       <- neigh_arr[neigh_arr>0]
        use_table <- summary(as.factor(use))
        network   <- as.integer(names(rev(sort(use_table))))
      }
    }
    print(paste(iter,':', network))
    write(as.numeric(use_coords), '')
    net_arr      <- c(net_arr, network[1])
  }
  return(net_arr)
}


query_atlas <- function(atlas, coords, rad, afni.dir) {
  cmd       <- paste(afni.dir, '/3dmaskdump -nbox ', coords[1]+rad, ':', coords[1]-rad, ' ', coords[2]+rad, ':', coords[2]-rad, ' ', coords[3]+rad, ':', coords[3]-rad, ' ', atlas, sep='')
  output    <- system(cmd, intern=TRUE)
  neigh_arr <- NULL
  for ( idx in 1:length(output) ) {
    cur_out   <- output[[idx]]
    split_out <- strsplit(cur_out, split=' ')
    network   <- as.numeric(split_out[[1]][4])
    neigh_arr <- c(neigh_arr, network)
  }
  return(neigh_arr)
}



plot_expression <- function(reg_micro_scale, out_file, reg, reg_df){

  x_min  = floor(min(reg_micro_scale$SST))
  x_max  = ceiling(max(reg_micro_scale$SST))
  y_min  = floor(min(reg_micro_scale$PVALB))
  y_max  = ceiling(max(reg_micro_scale$PVALB))
  xy_min = min(x_min, y_min)
  xy_max = max(x_max, y_max)
  
  # Plot correlation
  CairoPDF(out_file, width=5, height=5)
  
  p = ggplot(reg_micro_scale, aes(SST, PVALB)) + geom_point(size=3, alpha=.5, stroke=0, colour='black') +
    ggtitle(as.character(reg_df[reg])) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black",size=16),
          axis.text.y = element_text(colour="black",size=16),
          axis.title.x = element_text(colour="black",size=16),
          axis.title.y = element_text(colour="black",size=16)) +
    ylab("Normalized PVALB Expression (z)") +
    xlab("Normalized SST Expression (z)") +
    geom_smooth(method='lm', se=FALSE, linetype = "dashed") + 
    expand_limits(x = c(xy_min, xy_max), y = c(xy_min, xy_max)) + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(xy_min, xy_max, by = 1)) + 
    scale_y_continuous(expand = c(0, 0), breaks = seq(xy_min, xy_max, by = 1)) + 
    coord_fixed()
  print(p)
  dev.off()
}




plot_distribution <- function(cor.table, distr.gene, distr.color, compare.gene, compare.color){
  
  #cairo_pdf(out_file, width=4, height=2)

  avg.cors <- rowMeans(cor.table)
  avg.cors <- avg.cors[names(avg.cors) != distr.gene]
  plot.df  <- as.data.frame(avg.cors)
  
  p <- ggplot(plot.df, aes(avg.cors)) +
    geom_density(color="grey20", fill=distr.color, alpha=1) +
    ggtitle(paste0('Avg ', compare.gene,' correlation to ', distr.gene)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=16),
          axis.title.y = element_text(colour="black", size=16)) +
    ylab("Density") +
    xlab(paste0("Avg Cross-Region Correlation to ", distr.gene)) +
    expand_limits(y=c(0,4)) + 
    scale_x_continuous(expand = c(0, 0), limits=c(-.4, .6), breaks = seq(-.4,.6,.2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 4, by = 1))
  for ( i in 1:length(compare.gene) ){
    print(compare.gene[i])
    compare.corr <- plot.df[rownames(plot.df) == compare.gene[i],1]
    p <- p + geom_vline(xintercept = compare.corr, color=compare.color[i], size=1, linetype='longdash') + 
      annotate("text", x=compare.corr, y=1, label= compare.gene[i], color = compare.color[i])
  }
  return(p)
  #dev.off()
}

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}



plot_foci <- function(mni.coords, vertices, foci.names, ref.surface, rgba.color, wb.path, projection.structure, out.name){
  #coords.file <- '~/mni_coords.txt'
  #write.table(file=coords.file, x=mni.coords, row.names = FALSE, col.names = FALSE)
  gifti.file  <- readgii(ref.surface)
  #vertex.file <- '~/vertex_file.txt'
  #system(paste(wb.path, '-surface-closest-vertex', ref.surface, coords.file, vertex.file, sep=' '))
  #system(paste(wb.path, '-gifti-convert ASCII ', ref.surface, gsub('surf.gii', 'txt', ref.surface), sep=' '))
  #input <- read.table(file=vertex.file)
  #vertices <- input$V1
  
  top     <- newXMLNode("FociFile", attrs=c(Version="1.0"))
  metadat <- newXMLNode("MetaData", parent = top)
  md      <- newXMLNode("MD", parent = metadat)
  name    <- newXMLNode("Name", parent = md)
  value   <- newXMLNode("Name", parent = md)
  newXMLCDataNode('Caret-Version', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('5.64', parent = value, doc = NULL, at = NA, sep = "\n")
  
  md2      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md2)
  value      <- newXMLNode("Value", parent = md2)
  newXMLCDataNode('Date', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('2017-05-14T16:03:20', parent = value, doc = NULL, at = NA, sep = "\n")
  
  md3      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md3)
  value      <- newXMLNode("Value", parent = md3)
  newXMLCDataNode('UserName', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('vanessen', parent = value, doc = NULL, at = NA, sep = "\n")
  
  md4      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md4)
  value      <- newXMLNode("Value", parent = md4)
  newXMLCDataNode('comment', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('none', parent = value, doc = NULL, at = NA, sep = "\n")
  
  labels <- newXMLNode("LabelTable", parent = top)
  #base.node <- newXMLNode("Label", parent=labels)
  node.list <- list()
  for (idx in 1:dim(rgba.color)[1]){
    
    # Hacky fix for a single vertex that gives us trouble
    vertex   <- vertices[idx]
    v.idxs   <- which(gifti.file$data$triangle == vertex, arr.ind=TRUE)
    tri.vertices <- gifti.file$data$triangle[v.idxs[,1],]
    tmp.vert     <- tri.vertices[1:2,]
    vert.six <- c(tmp.vert[1,], tmp.vert[2,])
    TriAnat <- NULL
    for (i in 1:6){
      tmp <- gifti.file$data$pointset[vert.six[i],]
      TriAnat <- c(TriAnat, tmp)
    }
    if (length(TriAnat) != 18){
      write(idx,'')
      next
    }
    
    cur.rgba  <- rgba.color[idx,]
    cur.label <- paste('Label_', foci.names[idx], sep='')
    tmp.node <- newXMLNode("Label", attrs=c(Key=idx, Red=cur.rgba$r, Green=cur.rgba$g, Blue=cur.rgba$b, Alpha=cur.rgba$a))
    newXMLCDataNode(cur.label, parent=tmp.node, doc = NULL, at = NA, sep = "\n")
    node.list[[idx]] <- tmp.node
  }
  addChildren(labels, node.list)
  
  foci.list <- list()
  for (idx in 1:dim(rgba.color)[1]){
    cur.label <- paste('Label_', foci.names[idx], sep='')
    cur.mni   <- mni.coords[idx,]
    foci.node <- newXMLNode('Focus', attrs=c(Index=idx))
    name.node <- newXMLNode('Name')
    newXMLCDataNode(cur.label, parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    xyz.node <- newXMLNode('SearchXYZ',parent=foci.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=xyz.node)
    name.node <- newXMLNode('Geography')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Area')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('RegionOfInterest')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Size')
    newXMLCDataNode('0', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Statistic')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Comment')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('ClassName')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsIDNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsRepeatNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsParentFocusBaseID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsVersionNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsMSLID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('AttributeID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('StudyMetaDataLinkSet')
    name.node <- newXMLNode('AttributeID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SurfaceProjectedItem')
    struct.node <- newXMLNode('Structure', parent=name.node)
    newXMLTextNode(projection.structure, parent=struct.node)
    struct.node <- newXMLNode('StereotaxicXYZ', parent=name.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node)
    struct.node <- newXMLNode('VolumeXYZ', parent=name.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node)
    addChildren(foci.node, name.node)
    
    vertex   <- vertices[idx]
    v.idxs   <- which(gifti.file$data$triangle == vertex, arr.ind=TRUE)
    tri.vertices <- gifti.file$data$triangle[v.idxs[,1],]
    tmp.vert     <- tri.vertices[1:2,]
    vert.six <- c(tmp.vert[1,], tmp.vert[2,])
    TriAnat <- NULL
    for (i in 1:6){
      tmp <- gifti.file$data$pointset[vert.six[i],]
      TriAnat <- c(TriAnat, tmp)
    }
    if (length(TriAnat) != 18){
      write(idx,'')
      next
    }
    ve.node <- newXMLNode('VanEssenProjection')
    
    struct.node <- newXMLNode('DR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)
    struct.node <- newXMLNode('TriAnatomical', parent=ve.node)
    newXMLTextNode(paste(TriAnat, collapse=' '), parent=struct.node)
    struct.node <- newXMLNode('ThetaR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)
    struct.node <- newXMLNode('PhiR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)    
    struct.node <- newXMLNode('TriVertices', parent=ve.node)
    newXMLTextNode(paste(vert.six, collapse=' '), parent=struct.node)   
    two.vert <- names(rev(sort(table(vert.six))))[1:2]
    struct.node <- newXMLNode('Vertex', parent=ve.node)
    newXMLTextNode(paste(two.vert, collapse=' '), parent=struct.node) 
    struct.node <- newXMLNode('VertexAnatomical', parent=ve.node)
    tmp1 <- gifti.file$data$pointset[as.numeric(two.vert[1]),]
    tmp2 <- gifti.file$data$pointset[as.numeric(two.vert[2]),]
    newXMLTextNode(paste(c(tmp2, tmp1), collapse=' '), parent=struct.node) 
    struct.node <- newXMLNode('PosAnatomical', parent=ve.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node) 
    struct.node <- newXMLNode('FracRI', parent=ve.node)
    newXMLTextNode(5, parent=struct.node) 
    struct.node <- newXMLNode('FracRJ', parent=ve.node)
    newXMLTextNode(5, parent=struct.node) 
    
    addChildren(name.node, ve.node)
    addChildren(foci.node, name.node)
    
    foci.list[[idx]] <- foci.node
  }
  addChildren(top, foci.list)
  saveXML(top, file=out.name)
}







