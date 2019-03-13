plot_riail_geno <- function(riail_gt) {
  
  df$index <- dplyr::group_indices(riail_gt)
  
  strain_index <- df$sample
  names(strain_index) <- df$index + 0.5
  
  r_gt <- ggplot(df,  aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
    geom_rect(aes(alpha = low_sites)) +
    scale_alpha_discrete(range = c(1.0, 0.65)) +
    scale_fill_manual(values = strain_colors) +
    facet_grid(.~chrom, scales="free", space="free") +
    scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
    scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
    theme(strip.background = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None")
  return(r_gt)
}



maxlodplot_edit <- function(map){
  
  cidefiner <- function(cis, map) {
    ci_lod <- vapply(1:nrow(map), function(marker) {
      pos <- map$pos[marker]
      chr <- map$chr[marker]
      lod <- map$lod[marker]
      s <- sum(chr == cis$chr & (pos >= cis$ci_l_pos & pos <= cis$ci_r_pos))
      inci <- as.logical(s)
      cilodscore <- ifelse(inci, lod, 0)
      return(cilodscore)
    }, numeric(1))
    map$ci_lod <- ci_lod
    return(map)
  }
  
  map1 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(lod == max(lod))
  
  cis <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  
  if(nrow(cis) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype, y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  
  map1 <- cidefiner(cis, map1)
  
  plot <- ggplot2::ggplot(map1) +
    ggplot2::aes(x = pos/1e6, y = lod)
  
  if(nrow(cis) != 0) {
    plot <- plot + 
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "blue", alpha = 0.5) +
      ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=(1+maxlod)),
                          fill ="red", shape=25, size=2, show_guide = FALSE) +
      ggplot2::geom_text(data = cis,
                         ggplot2::aes(x=pos/1e6,
                                      y=(ifelse(lod < 10, 2+lod, 2+lod)),
                                      label = paste0(100*round(var_exp, digits = 4),"%")),
                         colour = "black", size=3)
  }
  
  plot <- plot + ggplot2::geom_line(size = 0.75, alpha = 1, color = axis_color) +
    ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_discrete(name="Mapping\nIteration")
  return(plot)
}

pxgplot_edit <- function(cross, map, parent="N2xCB4856") {
  peaks <- map %>% 
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  if(nrow(peaks) == 0) {
    stop("No QTL identified")
  }
  
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  
  colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  
  pheno <- cross$pheno %>%
    dplyr::select_(map$trait[1])
  geno <- data.frame(extract_genotype(cross)) %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                             function(marker) {
                                               paste(
                                                 unlist(
                                                   peaks[
                                                     peaks$marker == 
                                                       gsub("\\.",
                                                            "-",
                                                            marker),
                                                     c("chr", "pos")]),
                                                 collapse = ":")
                                             })
  colnames(geno)[ncol(geno)] <- "pheno"
  
  split <- tidyr::gather(geno, marker, genotype, -pheno)
  
  split$genotype <- sapply(split$genotype, function(x){
    if(is.na(x)) {
      return(NA)
    }
    if(parent=="N2xCB4856") {
      if(x == -1) {
        "N2"
      } else {
        "CB4856"
      }
    } else if(parent=="N2xLSJ2") {
      if(x == -1) {
        "N2"
      } else {
        "LSJ2"
      }
    } else if(parent=="AF16xHK104") {
      if(x==-1) {
        "AF16"
      } else {
        "HK104"
      }
    }
  })
  
  split$genotype <- factor(split$genotype, levels = c("N2","CB4856","LSJ2","AF16","HK104"))
  
  ggplot2::ggplot(split) +
    ggbeeswarm::geom_beeswarm(ggplot2::aes(x = genotype, y = pheno), alpha = .4,priority = "density",cex = 1.2) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA,alpha=0.7) +
    ggplot2::scale_fill_manual(values = strain_colors) +
    ggplot2::facet_wrap(~ marker, ncol = 5) +
    ggplot2::ggtitle(peaks$trait[1]) +
    ggplot2::labs(x = "Genotype", y = "Phenotype")
}


cegwas2_manplot <- function(plot_df, 
                            bf_line_color = "gray",
                            eigen_line_color = "gray",
                            eigen_cutoff = independent_test_cutoff) {
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(i) {
    
    bf_cut <- -log10(0.05/length(unique(plot_df$marker)))
    
    plot_df_pr <- plot_df %>%
      dplyr::filter(trait == i,
                    CHROM != "MtDNA") %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      dplyr::mutate(EIGEN_CUTOFF = eigen_cutoff,
                    BF = bf_cut) %>%
      dplyr::mutate(EIGEN_SIG = ifelse(log10p > bf_cut, "1", 
                                       ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )
    
    plot_df_pr  %>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                             "1" = 1,
                                             "2" = 1)) +
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "hotpink3"), 
                         color = "hotpink",linetype = 2, 
                         alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG!="1") %>% na.omit())+
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "blue"), 
                         color = "blue",fill = "cyan",linetype = 2, 
                         alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG=="1") %>% na.omit())+
      ggplot2::geom_hline(ggplot2::aes(yintercept = bf_cut),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = eigen_line_color, 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)), size = 1 ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))),
                    title = plot_traits)
  })
  plots
}
