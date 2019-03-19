try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
albendazole_data <- "201809_rappaport/Data/ben1_data/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")

trait_of_interest <- "q90.TOF"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"

######################################################################################################################## Dose Response

dr_data <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS2_DR_Processed.tsv"))

dr_data$strain <- gsub("Ju", "JU", dr_data$strain)

# length
pr_df <- dr_data%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == trait_of_interest, grepl("alb|DMSO",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition)) %>%
  dplyr::arrange(value) %>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

pr_df%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", 
                            "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = final_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Albendazole (µM)",
       y = "Relative Animal Length") + 
  base_theme+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))


ggsave(filename = "Plots/Albendazole_DoseResponse.png", height = 4.5, width = 9, dpi = 400)
ggsave(filename = "Plots/Albendazole_DoseResponse.pdf", height = 4.5, width = 9, dpi = 400)


######################################################################################################################## GWAS

pr_maps <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait != "norm.n")

manplot_edit(pr_maps)[[1]]+
  base_theme + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill = background_color, color = axis_color))+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(paste0(plot.dir,"snv_mapping.png"), 
       dpi = 600,
       height = 4, 
       width = 12)

q90burden_pr <- data.table::fread(paste0( "~/Dropbox/AndersenLab/LabFolders/Stefan/labmeeting/201809_rappaport/Data/updated_burden/Albendazole_q90.TOF.Skat.assoc"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))%>%
  dplyr::filter(CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(Pvalue)

q90burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(as.numeric(Pvalue)), alpha = significant, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","#D7263D"))+
  scale_alpha_manual(values = c(0.5,1))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  ggplot2::theme(legend.position = "none")+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = axis_color,size = 1.2),
                 legend.position = "none")+
  base_theme + 
  theme(panel.border = element_rect(fill = NA),
        legend.position = "none")+
  labs(x = "Genomic Position (Mb)", 
       y = expression(-log[10](italic(p))))

ggsave(paste0(plot.dir,"abz-q90_snv-indel_skat.png"), 
       dpi = 600,
       height = 4, 
       width = 12)

load( "~/Dropbox/AndersenLab/LabFolders/Stefan/labmeeting/201809_rappaport/Data/updated_burden/outlierRemoved_All_drugs_Significant_genes_SKAT.Rda")

summarized_bz <- skat.assoc.pr %>%
  dplyr::filter(grepl("zole|sulfone", drug),
                NumVar > 5,
                !grepl("cv|iqr|var|f.L|PC", trait))%>%
  dplyr::select(CHROM,POS,Pvalue,drug,trait) %>%
  dplyr::mutate(log10p = -log10(as.numeric(Pvalue))) %>%
  dplyr::ungroup()%>%
  dplyr::arrange(Pvalue)%>%
  dplyr::mutate(tidydrug = ifelse(grepl("sulfone", drug), "Albendazole\nSulfone",
                                  ifelse(grepl("sulfoxide", drug), "Albendazole\nSulfoxide", drug)))

summarized_bz%>%
  ggplot()+
  aes(x = POS/1e6, y = trait, fill = log10p)+
  geom_point(shape = 21, size = 3)+
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = expression(-log[10](italic(p))), breaks = c(7,10,13,16))+
  facet_grid(tidydrug~CHROM, scales = "free", space = "free")+
  ggplot2::theme(legend.position = "none")+
  ggplot2::theme(panel.background = ggplot2::element_rect(color = axis_color,size = 1.2))+
  base_theme + 
  theme(panel.border = element_rect(fill = NA),
        axis.text.y = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(size = 16, angle = 360,hjust = 0))+
  labs(x = "Genomic Position (Mb)", 
       y = "Trait")+
  xlim(c(0,20))

ggsave(paste0(plot.dir,"bz_summarized_skat.png"), 
       dpi = 600,
       height = 10, 
       width = 12)

load( "~/Dropbox/AndersenLab/LabFolders/Stefan/labmeeting/201809_rappaport/Data/updated_burden/ben1Regressed_All_drugs_Significant_genes_SKAT.Rda")

summarized_bz <- skat.assoc.pr %>%
  dplyr::filter(grepl("zole|sulfone", drug),
                NumVar > 5,
                !grepl("cv|iqr|var|f.L|PC", trait))%>%
  dplyr::select(CHROM,POS,Pvalue,drug,trait) %>%
  dplyr::mutate(log10p = -log10(as.numeric(Pvalue))) %>%
  dplyr::ungroup()%>%
  dplyr::arrange(Pvalue)%>%
  dplyr::mutate(tidydrug = ifelse(grepl("sulfone", drug), "Albendazole\nSulfone",
                                  ifelse(grepl("sulfoxide", drug), "Albendazole\nSulfoxide", drug)))

summarized_bz%>%
  ggplot()+
  aes(x = POS/1e6, y = trait, fill = log10p)+
  geom_point(shape = 21, size = 3)+
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = expression(-log[10](italic(p))), breaks = c(7,10,13,16))+
  facet_grid(tidydrug~CHROM, scales = "free", space = "free")+
  ggplot2::theme(legend.position = "none")+
  ggplot2::theme(panel.background = ggplot2::element_rect(color = axis_color,size = 1.2))+
  base_theme + 
  theme(panel.border = element_rect(fill = NA),
        axis.text.y = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(size = 16, angle = 360,hjust = 0))+
  labs(x = "Genomic Position (Mb)", 
       y = "Trait")+
  xlim(c(0,20))

ggsave(paste0(plot.dir,"ben1regressed_bz_summarized_skat.png"), 
       dpi = 600,
       height = 10, 
       width = 12)

modified_ld_plot(plot_df = pr_maps, trait = "q90.TOF")+
  base_theme+
  theme(legend.position = "none")

ggsave(paste0(plot.dir,"marker_peak_ld.png"), 
       dpi = 600,
       height = 6, 
       width = 8)
