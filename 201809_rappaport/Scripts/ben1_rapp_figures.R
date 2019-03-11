# Finalized figures for BEN-1 paper

# calculate ka/ks with http://services.cbu.uib.no/tools/kaks
# use "seqdump.txt" in data/custom_Popgen folder


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# set to location of files
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
main.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/labmeeting/201809_rappaport/"
setwd(main.dir)
data.dir <- paste0("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Raw_data/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")
final.dir <- paste0(main.dir,"Data/ben1_data/Final_Tables/")


source(paste0(script.dir,"ben1_processing_functions.R"))
source(paste0(script.dir,"N2_CB.R"))

trait_of_interest <- "q90.TOF"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"

library(ggbeeswarm)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# Dose response
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

dr_data <- data.table::fread(paste0(final.dir,"TS2_DR_Processed.tsv"))

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
  base_theme



ggsave(paste0(plot.dir,"doseresponse.png"), 
       dpi = 300,
       height = 4.5, 
       width = 9)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# PxG splits for G-P association
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #


load("/Users/stefan/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/2018_analysis/Data/2018_all_traits_combined.Rda")

bz_pheno <- df_merge2[,c("strain","Abamectin_q75.TOF_ctrl-regressed" )] 

colnames(bz_pheno) <- c("strain", "value")

bz_scaled <- bz_pheno%>%
  na.omit() %>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ggplot(bz_scaled)+
  aes(x = final_pheno,stat(count))+
  geom_density(fill = "#C6C6C6") +
  base_theme +
  labs(y = "Trait")+
  scale_x_continuous(breaks = c(0.1,0.9), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"distribution.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

ggplot(bz_scaled)+
  aes(y = final_pheno, x = "REF")+
  geom_beeswarm(fill = "#C6C6C6", shape = 21, size =2) +
  base_theme +
  labs(y = "Trait")+
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"beeswarm_distribution.png"), 
       dpi = 600,
       height = 4, 
       width = 6)


bz_scaled%>%
  dplyr::rowwise()%>%
  dplyr::mutate(nonsig = ifelse(rnorm(n=1,mean = 0, sd = 1)  >0, "A","G")) %>%
  ggplot()+
  aes(x = nonsig, y = final_pheno, fill = nonsig)+
  geom_beeswarm(size = 2, shape = 21) +
  scale_fill_manual(values = c( "#94F1A2","#479B76")) +
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  base_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"non_sig_split.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

bz_scaled%>%
  dplyr::rowwise()%>%
  dplyr::mutate(nonsig = ifelse(rnorm(n=1,mean = 0, sd = 1)  >0, "AA","REF")) %>%
  ggplot()+
  aes(x = nonsig, y = final_pheno, fill = nonsig)+
  geom_beeswarm(size = 2, shape = 21) +
  scale_fill_manual(values = c("#71B2E0", "#C6C6C6")) +
  base_theme +
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"non_sig_split_insertion.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

t1 <-  quantile(bz_scaled$final_pheno, probs = 0.75, na.rm = T)
bz_scaled%>%
  dplyr::rowwise()%>%
  dplyr::mutate(sig = ifelse(final_pheno > t1 & rnorm(n=1,mean = 0, sd = 1)  > -1, "G","T")) %>%
  na.omit()%>%
  ggplot()+
  aes(x = sig, y = final_pheno, fill = sig)+
  geom_beeswarm(size = 2, shape = 21) +
  scale_fill_manual(values = c( "#94F1A2","#479B76")) +
  base_theme +
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"sig_split.png"), 
       dpi = 600,
       height = 4, 
       width = 6)


bz_scaled%>%
  dplyr::rowwise()%>%
  dplyr::mutate(nonsig = ifelse(rnorm(n=1,mean = 0, sd = 1)  >0, "dACG","REF")) %>%
  ggplot()+
  aes(x = nonsig, y = final_pheno, fill = nonsig)+
  geom_beeswarm(size = 2, shape = 21) +
  scale_fill_manual(values = c("#E4B044", "#C6C6C6")) +
  base_theme +
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0(plot.dir,"non_sig_split_deletion.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

t1 <-  quantile(bz_scaled$final_pheno, probs = 0.7, na.rm = T)

bz_pheno_plt <- bz_scaled%>%
  dplyr::rowwise()%>%
  dplyr::mutate(sig = ifelse(final_pheno > t1 & rnorm(n=1,mean = 0, sd = 1)  > 1.3, "ALT1",
                             ifelse(final_pheno > t1 & rnorm(n=1,mean = 0, sd = 1)  > 1.3, "ALT2",
                                    ifelse(final_pheno > t1 & rnorm(n=1,mean = 0, sd = 1)  > 1.3, "ALT3", 
                                           ifelse(final_pheno > t1 & rnorm(n=1,mean = 0, sd = 1)  > 1.3, "ALT4", "REF"))))) %>%
  dplyr::mutate(burden_all = ifelse(sig == "REF", "REF","ALT")) %>%
  dplyr::mutate(burden1 = ifelse(sig %in% c("REF","ALT2","ALT3","ALT4"), "REF","ALT")) %>%
  dplyr::mutate(burden2 = ifelse(sig %in% c("REF","ALT1","ALT3","ALT4"), "REF","ALT")) %>%
  dplyr::mutate(burden3 = ifelse(sig %in% c("REF","ALT1","ALT2","ALT4"), "REF","ALT")) %>%
  dplyr::mutate(burden4 = ifelse(sig %in% c("REF","ALT1","ALT2","ALT3"), "REF","ALT")) %>%
  na.omit()

ggplot(bz_pheno_plt)+
  aes(x = burden_all, y = final_pheno, fill = sig)+
  geom_beeswarm(size = 2, shape = 21,priority='density', dodge.width=0.15) +
  scale_fill_manual(values = c( "#DCA137","#479B76","#BA5122", "#6EB3E4","#C6C6C6")) +
  base_theme +
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggsave(paste0(plot.dir,"burden_sig_split_all.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

ggplot(bz_pheno_plt)+
  aes(x = burden1, y = trait1, fill = sig)+
  geom_beeswarm(size = 3, shape = 21,priority='density') +
  scale_fill_manual(values = c( "#DCA137","#479B76","#BA5122", "#6EB3E4","#C6C6C6")) +
  base_theme +
  labs(y = "Trait")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave(paste0(plot.dir,"burden_sig_split1.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

ggplot(bz_pheno_plt)+
  aes(x = burden2, y = trait1, fill = sig)+
  geom_beeswarm(size = 3, shape = 21,priority='density') +
  scale_fill_manual(values = c( "#DCA137","#479B76","#BA5122", "#6EB3E4","#C6C6C6")) +
  base_theme +
  labs(y = "Trait")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave(paste0(plot.dir,"burden_sig_split2.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

ggplot(bz_pheno_plt)+
  aes(x = burden3, y = trait1, fill = sig)+
  geom_beeswarm(size = 3, shape = 21,priority='density') +
  scale_fill_manual(values = c( "#DCA137","#479B76","#BA5122", "#6EB3E4","#C6C6C6")) +
  base_theme +
  labs(y = "Trait")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave(paste0(plot.dir,"burden_sig_split3.png"), 
       dpi = 600,
       height = 4, 
       width = 6)

ggplot(bz_pheno_plt)+
  aes(x = burden4, y = trait1, fill = sig)+
  geom_beeswarm(size = 3, shape = 21,priority='density') +
  scale_fill_manual(values = c( "#DCA137","#479B76","#BA5122", "#6EB3E4","#C6C6C6")) +
  base_theme +
  labs(y = "Trait")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave(paste0(plot.dir,"burden_sig_split4.png"), 
       dpi = 600,
       height = 4, 
       width = 6)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# 330 distribution
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")

strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude)

ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color=axis_color, fill=background_color, size=0.5)+
  geom_point(data = strains_330, 
             aes(x=as.numeric(long), y=as.numeric(lat)), 
             color = "black",
             fill = "red",
             shape = 21, 
             size = 3) +
  theme_map()+
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        legend.text = element_text(size = 22, family = number_font),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm")) 

ggsave(paste0(plot.dir,"330_world_distribution.pdf"),
       height = 10,
       width = 20)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  FIGURE 1 -PLOT GWAS and BURDEN
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_maps <- data.table::fread(paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
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

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  summarize marker mappings
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
load("/Users/Stefan/Dropbox/AndersenLab/LabFolders/Stefan/labmeeting/201809_rappaport/Data/updated_marker/20180815_gwas_ben-1_regressed-data_BF5-mappings.Rda")

significant <- sig_maps_df %>%
  na.omit()
# save(significant, file = "~/Dropbox/AndersenLab/LabFolders/Steffen/Projects/GWA_anthelmintics (May 2017)/data/signficant.Rda")

uniquemap <- significant %>%
  distinct(trait, log10p, .keep_all = TRUE)

# uniquemap <- uniquemap[-1,]

goodtraits <- uniquemap %>%
  ungroup() %>%
  dplyr::mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[,1]) %>%
  dplyr::mutate(trait = stringr::str_split_fixed(trait, "\\.", 2)[,2]) %>%
  dplyr::filter(!grepl("red|green|yellow|f.|iqr", trait)) %>%
  ungroup() %>%
  dplyr::group_by(CHROM) %>%
  dplyr::arrange(CHROM, startPOS)

goodtraits$POS <- as.numeric(goodtraits$POS)
goodtraits$startPOS <- as.numeric(goodtraits$startPOS)
goodtraits$endPOS <- as.numeric(goodtraits$endPOS)
goodtraits$peakPOS <- as.numeric(goodtraits$peakPOS)
goodtraits$log10p <- as.numeric(goodtraits$log10p)

# save(goodtraits, file = "~/Dropbox/AndersenLab/LabFolders/Steffen/Projects/GWA_anthelmintics (May 2017)/data/goodtraits.Rda")

goodtraits <- goodtraits %>%
  dplyr::mutate(class = NA)

for(i in 1:nrow(goodtraits)) {
  if(goodtraits$condition[i] %in% c("albendazole", "albendazolesulfone", "albendazolesulfoxide", "mebendazole", "thiabendazole", "triclabendazole")) {
    goodtraits$class[i] <- "BZ"
  } else if(goodtraits$condition[i] %in% c("abamectin", "ivermectin")) {
    goodtraits$class[i] <- "ML"
  } else if(goodtraits$condition[i] %in% c("levamisole","pyrantel", "morantel")) {
    goodtraits$class[i] <- "NARA"
  } else if(goodtraits$condition[i] %in% c("emodepside", "diethylcarbamazine", "praziquantel")) {
    goodtraits$class[i] <- "Others"
  }
}

# #Factor drug class
goodtraits$class <- factor(goodtraits$class)
goodtraits <- goodtraits %>%
  arrange(class, condition)
goodtraits$condition <- factor(goodtraits$condition, levels = c("albendazole", "albendazolesulfone", "albendazolesulfoxide", "mebendazole", "thiabendazole", "triclabendazole"),
                               labels = c("albendazole", "albendazolesulfone", "albendazolesulfoxide", "mebendazole", "thiabendazole", "triclabendazole"))

#Set chromosome boundaries
newrows <- goodtraits[1,]
newrows$POS <- as.numeric(newrows$POS)
newrows$startPOS <- as.numeric(newrows$startPOS)
newrows$endPOS <- as.numeric(newrows$endPOS)
newrows$peakPOS <- as.numeric(newrows$peakPOS)
newrows$log10p <- as.numeric(newrows$log10p)

newrows[1,] = c(NA,"I",1,"pc1",1,NA,NA,NA,NA,NA,NA,0,NA,14972282,NA, NA, "albendazole", NA)
newrows[2,] = c(NA,"II",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,15173999,NA, NA, "albendazole", NA)
newrows[3,] = c(NA,"III",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,13829314,NA, NA, "albendazole", NA)
newrows[4,] = c(NA,"IV",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,17450860,NA, NA, "albendazole", NA)
newrows[5,] = c(NA,"V",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,20914693,NA, NA, "albendazole", NA)
newrows[6,] = c(NA,"X",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,17748731,NA, NA, "albendazole", NA)
newrows$POS <- as.numeric(newrows$POS)
newrows$startPOS <- as.numeric(newrows$startPOS)
newrows$endPOS <- as.numeric(newrows$endPOS)
newrows$peakPOS <- as.numeric(newrows$peakPOS)
newrows$log10p <- as.numeric(newrows$log10p)
newrows <- dplyr::ungroup(newrows)

goodtraits%>%
  dplyr::ungroup()%>%
  ggplot()+
  aes(x=POS/1E6, y=trait)+
  theme_bw() +
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = expression(-log[10](italic(p))))+
  scale_color_gradient(high = "#D7263D", low = "#0072B2",
                       name = expression(-log[10](italic(p))))+
  geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = log10p), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = trait, xend = endPOS/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=log10p),colour = "black",size = 4, alpha = 1, shape = 25)+
  xlab("Genomic Position (Mb)") + ylab("") +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
  base_theme + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill = background_color, color = axis_color),
        strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  facet_grid(condition ~ CHROM, scales = "free", space = "free")+
  ggplot2::labs(x = "Genomic Position (Mb)")

ggsave(paste0(plot.dir,"bz_reg_summarized_marker.png"), 
       dpi = 600,
       height = 10, 
       width = 8) 

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  FIGURE 2A - BAR PLOT OF BEN-1 VARIATION
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")

gwa_mappings <- data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait!="norm.n")%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")

gwa_mappings$snpMarker <- gsub("_",":", gwa_mappings$snpMarker)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")
gwa_mappings <- gwa_mappings%>%
  dplyr::mutate(length = end-start)%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), "A", 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))




gwa_mappings_bar <- gwa_mappings %>%
  dplyr::distinct(strain, value, marker, ben1_prediction)%>%
  dplyr::arrange(marker)%>%
  dplyr::distinct(strain, value, .keep_all=T)%>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ben1_all <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon\nInsertion"), ben1_prediction, "A"))

ben1_all$build_GT[is.na(ben1_all$build_GT)] <- "A"

number_font <- "Itim"
axes_text_size <- 20
axes_title_font <- "Montserrat ExtraBold"
axes_title_size <- 18
title_size <- 20


allgray <- c("Deletion" = "#999999","Insertion" = "#999999","Inversion" = "#999999","Stop Gained" = "#999999",
                    "Transposon\nInsertion" = "#999999","Missense" = "#999999","Splice Donor" = "#999999","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = allgray) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_nocolor.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = colors) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))
  
ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_all_colors.png"), 
       dpi = 300,
       height = 6, 
       width = 10)
 

build_missense <- c("Deletion" = "#999999","Insertion" = "#999999","Inversion" = "#999999","Stop Gained" = "#999999",
            "Transposon\nInsertion" = "#999999","Missense" = "#009E73","Splice Donor" = "#999999","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_missense) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_missense.png"), 
       dpi = 300,
       height = 6, 
       width = 10)


build_m_sp <- c("Deletion" = "#999999","Insertion" = "#999999","Inversion" = "#999999","Stop Gained" = "#999999",
                     "Transposon\nInsertion" = "#999999","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_m_sp) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_miss_splice.png"), 
       dpi = 300,
       height = 6, 
       width = 10)


build_miss_sp_st <- c("Deletion" = "#999999","Insertion" = "#999999","Inversion" = "#999999","Stop Gained" = "#F0E442",
                    "Transposon\nInsertion" = "#999999","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_miss_sp_st) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_miss_splice_stop.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

build_miss_sp_st_del <- c("Deletion" = "#E69F00","Insertion" = "#999999","Inversion" = "#999999","Stop Gained" = "#F0E442",
                    "Transposon\nInsertion" = "#999999","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_miss_sp_st_del) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_miss_splice_stop_del.png"), 
       dpi = 300,
       height = 6, 
       width = 10)


build_miss_sp_st_del_ins <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#999999","Stop Gained" = "#F0E442",
                    "Transposon\nInsertion" = "#999999","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_miss_sp_st_del_ins) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_miss_splice_stop_del_ins.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

build_miss_sp_st_del_ins_te <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#999999","Stop Gained" = "#F0E442",
                              "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = build_miss_sp_st_del_ins_te) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_miss_splice_stop_del_ins_te.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_variation_bar <- plot_bar(ben1_all, pt_colors = colors) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ben1_variation_bar

ggsave(plot = ben1_variation_bar,
       paste0(plot.dir,"ben1_resistance_barplot_all_colors.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE 2B - BEN-1 HAPLOTYPE GENE MODEL PLOT
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# THINGS TO MANUALLY ADJUST IN AN SVG SOFTWARE
# COLORS OF VARIANTS
# GROUP THESE STRAINS TOGETHER - ED3005 JU1808 MY18 MY795 - BC MY795 HAS A VARIANT UNRELIABLE CALLED

# # # LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")


for_plot <- data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait!="norm.n")%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  # dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  na.omit()%>%
  dplyr::filter(peak_id==1)%>%
  dplyr::select(strain, phenotype = value) %>%
  dplyr::left_join(.,ben1_variants,by="strain")


# construct ben-1 model
ben1_gene_info <- data.frame(feature = c("5'UTR", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon", "3'UTR"),
                             start = c(3541630,3541595,3541432,3540411,3540065,3540004,3539841,3539784,3539310,3538501,3538325),
                             stop = c(3541596,3541433,3540412,3540064,3540005,3539842,3539785,3539311,3538502,3538326,3538293),
                             color = c("blue","orange","gray","orange","gray","orange","gray","orange","gray","orange","blue"),
                             min_y = min(for_plot$phenotype)-10,
                             max_y = max(for_plot$phenotype)+10,
                             gene = "ben-1",
                             strand = "Watson")%>%
  dplyr::mutate(size = start-stop)


plot_df<- for_plot%>%
  dplyr::filter(GT != "HET", !is.na(GT), GT == "ALT")%>%
  dplyr::filter(!grepl("intron_variant_MODIFIER|modif|synon|del_3542405_3542407",marker,ignore.case = T))%>%
  dplyr::arrange(start)%>%
  dplyr::select(strain, phenotype, marker)%>%
  dplyr::mutate(GT = 1)%>%
  tidyr::spread(marker, GT)

plot_df[is.na(plot_df)] <- 0

ben1Cluster <- kmeans(plot_df[, 4:ncol(plot_df)], 26, nstart = 20)
max(ben1Cluster$cluster)

ben1haps <- data.frame(cluster= ben1Cluster$cluster, plot_df)%>%
  dplyr::group_by(cluster)%>%
  dplyr::mutate(mean_pheno = mean(phenotype))%>%
  dplyr::ungroup()%>%
  dplyr::arrange((mean_pheno))%>%
  dplyr::mutate(f_clust = factor(cluster, levels = unique(cluster), ordered = T))%>%
  tidyr::gather(marker, GT, -strain, -phenotype, -cluster, -f_clust,-mean_pheno)

marker_starts <- dplyr::select(for_plot,marker, start,end)%>%
  dplyr::distinct()

marker_starts$marker <- gsub(">","\\.",marker_starts$marker)
marker_starts$marker <- gsub("-","\\.",marker_starts$marker)
marker_starts$marker <- gsub("\\+","\\.",marker_starts$marker)
marker_starts$marker <- gsub("&","\\.",marker_starts$marker)

ben1haps1 <- dplyr::left_join(ben1haps, marker_starts, by = "marker")%>%
  dplyr::mutate(plot_shape = ifelse(grepl("ins",marker,ignore.case = T), "insertion",
                                    ifelse(grepl("del",marker), "deletion",
                                           ifelse(grepl("missense",marker), "missense",
                                                  ifelse(grepl("stop",marker), "stop_gained",
                                                         ifelse(grepl("spli",marker),"spli",
                                                                ifelse(grepl("MODIF",marker),"mod",
                                                                       ifelse(grepl("inv",marker),"inv", 
                                                                              ifelse(grepl("trans",marker),"trans", "other")))))))))%>%
  dplyr::filter(GT!=0)%>%
  dplyr::group_by(cluster)%>%
  dplyr::mutate(hap_strains= paste(unique(strain), collapse = " "))%>%
  dplyr::distinct(cluster,marker,.keep_all=T)%>%
  dplyr::ungroup()%>%
  dplyr::arrange(mean_pheno)

pt_size <- 5
ben1_high_var<-ggplot(ben1haps1)+
  aes(x=start, 
      y = factor(f_clust,labels=unique(hap_strains),levels=unique(cluster), ordered=T), 
      fill = factor(GT), 
      shape = plot_shape,
      size = plot_shape)+
  scale_y_discrete(labels = scales::wrap_format(25),position = "right")+
  geom_segment(aes(x = as.numeric(start), 
                   xend = as.numeric(end), 
                   yend = factor(f_clust,labels=unique(hap_strains),levels=unique(cluster), ordered=T)), 
               size =5, 
               color = "red")+
  geom_point( color ="black")+
  scale_shape_manual(values = c("insertion" = 73, 
                                "deletion" = 2, 
                                "missense" = 77, 
                                "stop_gained" = 13, 
                                "mod"=96,
                                "spli" = 83,
                                "trans" = 84,
                                "inv" = 37,
                                other= 16))+
  scale_size_manual(values = c("insertion" = pt_size,
                               "deletion" = pt_size,
                               "missense" = pt_size,
                               "stop_gained" = pt_size,
                               "mod"=pt_size,
                               "spli" = pt_size,
                               "trans"=pt_size,
                               "inv" = pt_size,
                               "other"= pt_size))+
  theme_bw()+
  labs(y="",x = "Genomic Position (Mb)")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black",vjust=1),
        plot.title = element_text(size=10, face="bold",vjust=1),
        legend.title=element_blank(),
        legend.position = "none")+
  geom_vline(aes(xintercept=c(3541595)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3541432)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540411)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540065)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540004)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539841)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539784)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539310)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3538501)), linetype="dotdash", alpha = .25, color = "blue") +
  geom_vline(aes(xintercept=c(3538325)), linetype="dotdash", alpha = .25, color = "blue")

missense_df <- ben1_snps%>%
  dplyr::select(CHROM,POS, aa_change,effect)%>%
  dplyr::filter(effect == "missense_variant")%>%
  dplyr::distinct(aa_change, .keep_all =T)%>%
  dplyr::mutate(plot_aa = gsub(pattern = "p\\.","",aa_change))

missense_df$letter <- c("D404N", "M257I", "F200Y", "A185P","S145F","Q131L","E69G","K58E")

missense_df_low <- dplyr::filter(missense_df, letter %in% c("M257I","A185P","Q131L","K58E"))
missense_df_high <- dplyr::filter(missense_df, !(letter %in% c("M257I","A185P","Q131L","K58E")))

ben1_model <- gene_model(ben1_gene_info)+
  ylim(-4,4)+
  geom_segment(aes(x = as.numeric(POS), y = -1, xend = as.numeric(POS), yend = -2),data = missense_df_low, size =2)+
  geom_segment(aes(x = as.numeric(POS), y = 1, xend = as.numeric(POS), yend = 2),data = missense_df_high, size =2)+
  geom_label(aes(x = POS, y = -3, label = letter), data = missense_df_low,size = 2)+
  geom_label(aes(x = POS, y = 3, label = letter), data = missense_df_high,size = 2)


cowplot::ggdraw() +
  cowplot::draw_plot(ben1_high_var, 0, 0, 1, .8)+
  cowplot::draw_label(label = expression(bolditalic('ben-1')), 0.5, y = .92)+
  cowplot::draw_plot(ben1_model, .13, .81, .52, .1)+
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.8), size = 14)

ggsave(paste0(plot.dir,"F2B_Ben1_haps_by_phenotype.png"),
       dpi = 300,
       height = 12, 
       width = 18)

ggsave(paste0(plot.dir,"F2B_Ben1_haps_by_phenotype.pdf"),
       dpi = 300,
       height = 12, 
       width = 18)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S3 - Generated CHIMERA
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# TEXT INFORMATION -  ben-1 gene tajima d with manually curated variants
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # COMBINE SNPS AND INDELS - GET THIS INFORMATION FROM PREVIOUS SECTION

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker)) #  REMOVE MODIFIER (INTRON) and SYNONYMOUS VARIANTS

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # For coding variation
non_syn <- gene_level_TajimasD(ben1_variants)
non_syn
# Wattersons_Theta Average_PWD TajimasD Singletons FuLiF
# 1             5.25        0.39    -2.59         26 -8.09

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # For coding + MODIFIER + SYNONYMOUS variation
syn <- gene_level_TajimasD(ben1_variants)
syn
# Wattersons_Theta Average_PWD TajimasD Singletons FuLiF
# 1            12.64        4.15    -2.02         42 -6.07


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE 3 - ben-1 neutrality statistics with phylo and map
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
### TREE
tree <-  ape::read.tree( file = paste0(final.dir, "genome.tree"))

tree <- ggtree(tree, 
               branch.length = "rate", 
               layout = "circular") +xlim(-.075,.475)

tree

ben1_variation_matrix <- gwa_mappings%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,marker)%>%
  dplyr::distinct()

missing_strains <- data.frame(strain = row.names(cegwas::kinship)[!row.names(cegwas::kinship)%in%ben1_variation_matrix$strain],
                              marker = NA, stringsAsFactors = F)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

ben1_strain_markers <- ben1_variation_matrix %>%
  dplyr::bind_rows(.,missing_strains)%>%
  dplyr::distinct()%>%
  dplyr::filter(!grepl("missense_variant_MODERATE_c.599T>A_353", marker))%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(marker), NA, 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))

ben1_tree <- tree %<+% ben1_strain_markers+ 
  geom_tippoint(aes( color=ben1_prediction), size = 1)+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1")))) +
  theme(plot.background = element_rect(fill = background_color, color = NA))

ben1_tree

ggsave(paste0(plot.dir,"ben1_tree.pdf"),
       height = 15,
       width  = 15,
       dpi = 300)


# # # MAP 
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::select(strain, marker, ben1_prediction, long=longitude, lat=latitude, landscape, substrate)

isolation_info$long <- as.numeric(isolation_info$long)
isolation_info$lat <- as.numeric(isolation_info$lat)

map <- ggplot()+ geom_map(data=world, map=world,
                          aes(x=long, y=lat, map_id=region),
                          color=axis_color, fill=background_color, size= 0.5, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 4, shape = 21)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()+
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        legend.text = element_text(size = 22, family = number_font),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = background_color)) 

map 

ggsave(paste0(plot.dir,"ben-1_world_dist.png"),
       height = 10,
       width  = 20,
       dpi = 300)

### NEUTRALITY STATS

# . . . . 
# . . . . DEFINE SLIDE WINDOWS
# . . . . 

load(paste0(final.dir,"TS10_ben1_popgen.Rda"))

windowStarts <- data.frame(snp_index = 1:length(colnames(s_gen_stats@BIG.BIAL[[1]])),
                           position = colnames(s_gen_stats@BIG.BIAL[[1]]))

slide_index <- cbind(data.frame(lapply(s_gen_stats@SLIDE.POS,  function(x) as.numeric(floor(mean(x))))))%>%
  tidyr::gather(temp, snp_index)%>%
  dplyr::select(-temp)%>%
  dplyr::left_join(., windowStarts, by = "snp_index")

# . . . . 
# . . . . ANALYZYE NEUTRALITY/DIVERSITY STATS
# . . . .


neutrality_ls <- list()
for(popnumber in 1){
  popname <- LETTERS[[popnumber]]
  neutrality_ls[[popnumber]] <- data.frame(PopGenome::get.neutrality(s_gen_stats, theta = T, stats = T)[[popnumber]],
                                           PopGenome::get.diversity(s_gen_stats)[[popnumber]],
                                           PI = s_gen_stats@Pi[,popnumber])%>%
    dplyr::mutate(Population = popname,
                  WindowPosition = slide_index$position)%>%
    dplyr::select(-Pi)
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)
neutrality_df <- neutrality_df[,colSums(is.na(neutrality_df))<nrow(neutrality_df)]
neutrality_df <- tidyr::gather(neutrality_df, statistic, value, -Population, -WindowPosition)


plt_df <- neutrality_df%>%
  dplyr::filter(statistic %in%c("Tajima.D","Fay.Wu.H","Zeng.E","Fu.Li.F", "PI","nuc.diversity.within"))%>%
  dplyr::mutate(f_stats = ifelse(statistic == "Tajima.D", "Tajima's D",
                                 ifelse(statistic == "Fay.Wu.H", "Fay and Wu's H", 
                                        ifelse(statistic == "Fu.Li.F","Fu and Li's F",
                                               ifelse(statistic == "PI", "pi","Zeng's E")))))



plt_df%>%
  dplyr::filter(grepl("Taji",f_stats))%>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = as.numeric(value), fill = f_stats)+
  scale_color_manual(values = "blue")+
  # geom_smooth(size =1.5, span = .065, method = "loess", se = F)+
  geom_point(shape = 21, size =3)+
  scale_fill_manual(values=c(axis_color,"hotpink3",axis_color,highlight_color, "#F9A227","#0072B2"), name = "")+
  theme_bw(18)+
  labs(x = "Genomic Position (Mb)")+
  geom_vline(aes(xintercept=c(3.541628)), linetype="dotdash", alpha = 0.5, color = "black")+
  geom_vline(aes(xintercept=c(3.537688)), linetype="dotdash", alpha = 0.5, color = "black") +
  theme(strip.background = element_blank(),legend.position = "top")+
  guides(colour=guide_legend(override.aes = list(size = 4)))+
  xlim(c(3.53,3.55))+
  base_theme+
  theme(axis.title.y = ggplot2::element_blank())

ggsave(paste0(plot.dir,"ben-1_popgenstats_tajima.png"),
       height = 6,
       width  = 12,
       dpi = 300)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE S4 - Sampling location and statistics
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

phenos <- gwa_mappings %>%
  dplyr::select(strain,value)

ben_by_subs <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, landscape, substrate)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::left_join(.,phenos,by="strain") %>%
  dplyr::distinct(strain, value, .keep_all = T) %>%
  dplyr::mutate(LoF = ifelse(is.na(marker) | strain == "CB4856", "REF","ALT")) %>%
  dplyr::filter(landscape != "None")


ben_by_subs$landscape <- gsub(" zoo","",ben_by_subs$landscape)
ben_by_subs$landscape <- gsub(" ","\n",ben_by_subs$landscape)

sample_loc_plot <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(landscape))+
  geom_bar(position="fill",color = "black")+
  theme_classic(20)+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  scale_fill_brewer(palette = "Set1",name = "Sampling\nLocation")+
  plot.theme

substrate_type <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(substrate))+
  geom_bar(position="fill",color = "black")+
  scale_fill_manual(values=ancestry.colours,name = "Sampling\nSubstrate")+
  theme_classic(20)+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  plot.theme+
  theme(axis.title.y = element_blank())

cowplot::plot_grid(sample_loc_plot,
                   substrate_type, labels = c("A","B"),
                   ncol = 2,
                   align = 'h',
                   rel_widths = c(1,1),
                   label_size = 18)

ggsave(paste0(plot.dir,"FS4_sampling_location_substrate.png"),
       dpi = 300,
       height = 6,
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE 4 - HTA and competition assay of F200Y and DEL alleles
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_df <- data.table::fread(file = paste0(final.dir,"TS13_HTA_ben-1_regressed.tsv"))

main_figure <- pr_df%>%
  dplyr::filter(strain %in% c("N2", "ECA883", "ECA919"))%>%
  dplyr::arrange(phenotype)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(phenotype == min(phenotype), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(phenotype) - phenotype)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

main_figure%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  dplyr::mutate(strain1 = factor(group, levels = c("N2","del","F200Y"),
                                 labels =c("N2","Del","F200Y")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = final_pheno, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 1)+
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  theme_bw()+
  labs( y = "Animal length")+
  base_theme+
  theme(axis.title.x = element_blank(),
        legend.position="none")

ggsave(paste0(plot.dir,"ben1_allele_hta.png"), 
       dpi = 300,
       height = 6, 
       width = 8)

############ Figure 4 B: competition assay plot ########

competition_assay <- data.table::fread(file = paste0(final.dir,"TS14_competition_assay.csv"))%>%
  dplyr::filter(TargetType == "Ch1Unknown")%>%
  dplyr::select(Condition, Replicate, Generation, FractionalAbundance)%>%
  dplyr::mutate(id = paste(Condition,Generation, sep = "-")) %>%
  na.omit()%>% 
  dplyr::group_by(Condition, Generation)%>% 
  dplyr::mutate(Mean =  mean(100-FractionalAbundance)) %>% 
  dplyr::mutate(SND =  sd(100-FractionalAbundance))%>%
  dplyr::ungroup()%>%
  dplyr::select(id, FractionalAbundance)

outliers  <-  competition_assay%>%
  dplyr::group_by(id)%>%
  dplyr::transmute_if(is.numeric, isnt_out_funs)

competition_assay_outliers <- dplyr::bind_cols(competition_assay,outliers)%>%
  tidyr::separate(id, c("Condition", "Generation"), sep = "-")%>%
  # dplyr::mutate(outlier = ifelse(!z | !mad | !tukey, "OUTLIER", "OK"))%>%
  dplyr::ungroup() %>%
  dplyr::filter(!(!z | !mad | !tukey))%>%
  dplyr::group_by(Condition, Generation)%>%
  dplyr::mutate(Mean =  mean(100-FractionalAbundance)) %>% 
  dplyr::mutate(SND =  sd(100-FractionalAbundance))%>%
  tidyr::separate(Condition, into = c("Strain", "Condition"))%>% 
  dplyr::filter(Strain == "N2")

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  # filter(Condition == "DMSO")%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=.1, size = 1)+
  geom_line(aes(linetype=Condition, group=id), size = 2)+
  # geom_point(size = 2, shape =21, color = "black",aes(fill=Strain))+
  scale_color_manual(values = c("#999999","#999999"))+
  scale_fill_manual(values = c("#999999","#999999"))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_y_continuous(breaks = c(0,50,100), labels = c(0, 0.5, 1), limits = c(0,100))

ggsave(paste0(plot.dir,"competetion_N2.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

competition_assay_outliers <- dplyr::bind_cols(competition_assay,outliers)%>%
  tidyr::separate(id, c("Condition", "Generation"), sep = "-")%>%
  # dplyr::mutate(outlier = ifelse(!z | !mad | !tukey, "OUTLIER", "OK"))%>%
  dplyr::ungroup() %>%
  dplyr::filter(!(!z | !mad | !tukey))%>%
  dplyr::group_by(Condition, Generation)%>%
  dplyr::mutate(Mean =  mean(100-FractionalAbundance)) %>% 
  dplyr::mutate(SND =  sd(100-FractionalAbundance))%>%
  tidyr::separate(Condition, into = c("Strain", "Condition"))%>% 
  dplyr::filter(!Strain == "N2")

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  filter(Condition == "DMSO", Strain == "Del")%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=.1, size = 1)+
  geom_line(aes(linetype=Condition, group=id), size = 2)+
  # geom_point(size = 2, shape =21, color = "black",aes(fill=Strain))+
  scale_color_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 24))+
  scale_y_continuous(breaks = c(0,50,100), labels = c(0, 0.5, 1), limits = c(0,100))

ggsave(paste0(plot.dir,"competetion_DMSO_deletion.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  filter(Condition == "DMSO")%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=.1, size = 1)+
  geom_line(aes(linetype=Condition, group=id), size = 2)+
  # geom_point(size = 2, shape =21, color = "black",aes(fill=Strain))+
  scale_color_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 24))+
  scale_y_continuous(breaks = c(0,50,100), labels = c(0, 0.5, 1), limits = c(0,100))

ggsave(paste0(plot.dir,"competetion_DMSO_deletion_swap.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND, alpha=Condition), width=.1,size =1)+
  geom_line(aes(alpha=Condition, group=id), size =2)+
  # geom_point(aes(alpha=Condition))+
  scale_color_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_alpha_manual(values = c(1,0.2))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank())+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 24))+
  scale_y_continuous(breaks = c(0,50,100), labels = c(0, 0.5, 1), limits = c(0,100))

ggsave(paste0(plot.dir,"competetion_alpha.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

g3 <- competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition),
                b_f = 100 - FractionalAbundance)%>%
  dplyr::filter(Generation == 3, Condition == "ABZ")

ggstatsplot::ggbetweenstats(
  data = g3, 
  x = Strain, 
  y = b_f,
  messages = FALSE
)+         
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_color_manual(values = c("#E69F00", "#009E73","#999999"))

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  dplyr::filter(Generation == 3, Condition == "ABZ")%>%
  ggplot(aes(x=Strain,fill=Strain,y = 100 - FractionalAbundance))+
  # geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND, alpha=Condition), width=.1,size =1)+
  geom_boxplot()+
  geom_beeswarm(shape =21, size = 4)+
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_alpha_manual(values = c(1,0.2))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank())+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 24))+
  scale_y_continuous(breaks = c(0,50,75), labels = c(0, 0.5, .75), limits = c(50,75))

ggsave(paste0(plot.dir,"competetion_g3.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE S5 - HTA and competition assay of F200Y and DEL alleles - all strains
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_df <- data.table::fread(file = paste0(final.dir,"TS13_HTA_ben-1_regressed.tsv"))

pr_df_scaled <- pr_df %>%
  dplyr::arrange(phenotype)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(phenotype == min(phenotype), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(phenotype) - phenotype)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

pr_df_scaled%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  ggplot(.) +
  aes(x = factor(strain, levels = c("N2","ECA882","ECA883","ECA884","ECA917","ECA918","ECA919","ECA920","ECA921" )), 
      y = final_pheno, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("cadetblue2","hotpink2","grey75","grey75"))+
  theme_bw()+
  labs( y = "Animal length")+
  plot.theme +
  theme(axis.title.x = element_blank(),
        legend.position="none")

ggsave(paste0(plot.dir,"FS5_CRISPR_strains_HTA_complete.png"), 
       dpi = 300,
       height = 6, 
       width = 12)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE 5 - CHRX manhattan plot and fine mapping
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_resid_maps <- data.table::fread( paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"))
manhattan_plots <- manplot_edit(pr_resid_maps)

manhattan_plots[[1]] +
  base_theme + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill = background_color, color = axis_color))+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(paste0(plot.dir,"chrx_marker_qtl.png"), 
       dpi = 300,
       height = 4, 
       width = 12)


xqtl_fm <- process_correlations(variant_correlation(pr_resid_maps,condition_trait = F))
save(xqtl_fm,file = "updated_finemap_xqtl.rda")

# ,
# effect %in% c("missense_variant",
#               "splice_donor_variant&intron_variant",
#               "stop_gained",
#               "splice_acceptor_variant&intron_variant",
#               "splice_region_variant&non_coding_exon_variant",
#               "missense_variant&splice_region_variant",
#               "splice_region_variant")

fine_map <- xqtl_fm %>%
  dplyr::filter(strain == "EG4725")%>%
  # dplyr::mutate(tidy_effect = factor(effect, levels = c("missense_variant",
  #                                                       "splice_donor_variant&intron_variant",
  #                                                       "stop_gained",
  #                                                       "splice_acceptor_variant&intron_variant",
  #                                                       "splice_region_variant&non_coding_exon_variant",
  #                                                       "missense_variant&splice_region_variant",
  #                                                       "splice_region_variant"),
  #                                    labels = c("Missense",
  #                                               "Splice Donor\nIntron",
  #                                               "Stop Gained",
  #                                               "Splice Acceptor\nIntron",
  #                                               "Splice Region\nExon",
  #                                               "Splice Region\nMissense",
  #                                               "Splice Region")))%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(corrected_spearman_cor_p))+
  geom_point(size =2)+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(bold(-log[10](bolditalic(p)))))+
  scale_color_brewer(palette="Set1", name = "Effect:")+
  scale_fill_brewer(palette="Set1", name = "Effect:")+
  scale_shape_manual(values = c(16:25), name = "Effect:")+
  # facet_grid(.~qtl_name, scales= "free_x")+
  theme_bw(15)+
  plot.theme+
  theme(axis.title.y = element_blank())

plot_grid(ben1_resid_manplot, fine_map, labels = c("A","B"), 
          label_size = 18, ncol = 2, align = 'v', rel_widths = c(2,1))

ggsave(paste0(plot.dir,"F5_Ben-1_variation_regressed_manplot_finemap_",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 6, 
       width = 18)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S6 - CHRX barplot
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
pr_resid_maps <- data.table::fread( paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::left_join(.,ben1_variants, by = "strain")

pr_resid_maps$snpMarker <- gsub("_",":", pr_resid_maps$snpMarker)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

pr_resid_maps <- pr_resid_maps%>%
  dplyr::mutate(length = end-start)%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), "A", 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))
pr_resid_maps_bar <- pr_resid_maps %>%
  dplyr::distinct(strain, value, marker, ben1_prediction)%>%
  dplyr::arrange(marker)%>%
  dplyr::distinct(strain, value, .keep_all=T)%>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ben1_all <- pr_resid_maps_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon\nInsertion"), ben1_prediction, "None"))

ben1_all$build_GT[is.na(ben1_all$build_GT)] <- "None"

ggplot(ben1_all)+
  aes(x=strain2, y=final_pheno, fill = build_GT)+
  geom_bar(stat="identity", color = "black", size = .1) +
  scale_fill_manual(values= colors,
                    name = expression(paste("Variation at ", italic("ben-1"))))+
  labs(x = "Strain", y = paste0("Albendazole Resistance"))+
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ggsave(paste0(plot.dir,"ben1_regressed_barplot_all_colors.png"), 
       dpi = 300,
       height = 6, 
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S6 - CHRX burden
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

benreg_burden <- read.table(paste0(final.dir, "TS21_GWA_ben1_regressed.VariableThresholdPrice_processed.tsv"),header = T)

benreg_burden%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(paste0(plot.dir,"FS6_GWA_Ben1regressed_BURDEN_VTprice_manplot.png"), 
       height = 4, 
       width = 12)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S7 - single marker PxG, colored by ben-1 variants
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")

gwa_mappings <-  data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")

gwa_mappings$snpMarker <- gsub("_",":", gwa_mappings$snpMarker)
colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")
gwa_mappings <- gwa_mappings%>%
  dplyr::mutate(length = end-start)%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), "A", 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))

gwa_mappings <- gwa_mappings %>%
  tidyr::separate(snpMarker, into = c("CHROM","POS"),sep = ":",remove=F,convert=T)%>%
  dplyr::arrange(CHROM,POS)%>%
  dplyr::mutate(marker2 = factor(snpMarker, levels = unique(snpMarker), labels = unique(snpMarker)))%>%
  dplyr::distinct(marker2, value, .keep_all=T) %>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ggplot(gwa_mappings)+
  aes(x = snpGT, y = final_pheno)+
  geom_boxplot(outlier.colour = NA, fill = background_color)+
  facet_grid(.~marker2)+
  geom_beeswarm(shape = 21, color= "black", alpha = 0.8, fill = "#999999", data = dplyr::filter(gwa_mappings, is.na(marker)),cex=1.2)+
  geom_beeswarm(aes(fill=ben1_prediction), data = na.omit(gwa_mappings), shape = 21, color= "black", size = 3,cex=3)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2))+
  base_theme+
  theme(legend.position = "none")+
  labs(x = "SNV Genotype at QTL", y = "Relative BZ resistance")

ggsave(paste0(plot.dir,"marker_qtl_bxplot_colored_by_ben1.png"), 
       height = 4, 
       width = 12)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S8 - CHRX PxG - LD
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
pr_resid_maps <- data.table::fread( paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"))
regress_ben1 <- data.table::fread(paste0(final.dir,"TS16_ben1_regressed_with_ben1_covariate.tsv"))

x_pg_df <- pr_resid_maps %>%
  dplyr::filter(trait == "regressed_trait")%>%
  na.omit()%>%
  dplyr::select(qtl_marker=marker, strain, reg_val = value, qtl_allele = allele )%>%
  dplyr::left_join(regress_ben1,.,by ="strain")%>%
  dplyr::mutate(qtl_marker1 = factor(qtl_marker, levels = unique(qtl_marker)))

x_pg_df%>%
  ggplot()+
  aes(x = factor(qtl_allele, levels = c(-1,1),
                 labels = c("REF", "ALT")), y = reg_val)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.25, aes(fill = factor(ben1_variant_trait)), shape = 21, size =2)+
  scale_fill_manual(values = c("cadetblue3","hotpink3"), labels = c("REF", "ALT"), name=expression(italic(ben-1))) +
  theme_bw(15)+
  facet_grid(.~qtl_marker1)+
  labs(x="QTL allele", y = "Regressed Animal Length")+
  theme(axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"))


ggsave(paste0(plot.dir,"FS8_GWAS_Xqtl_split_colored_by_ben1_variation.png"), 
       dpi = 300,
       height = 8, 
       width = 8)



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S11 - phenotyping WN2002 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# read wild strain HTA data 


assayregressed_WI <- data.table::fread(paste0(final.dir,"TS23_WI_rephenotype_HTA_processed.tsv"))

assayregressed_WI%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  ggplot(.) +
  aes(x = factor(strain, levels = c("N2", "JU2141", "WN2002","JU2581")), 
      y = phenotype, 
      fill=strain) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  theme_bw()+
  labs( y = "Animal length")+
  scale_fill_manual(values=c("hotpink3","cadetblue3","orange","gray75"),name="Strain")+
  plot.theme +
  theme(axis.title.x = element_blank(),
        legend.position="none")

ggsave(paste0(plot.dir,"FS9_WN2002_HTA.png"), 
       dpi = 300,
       height = 6, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S10 - Generated in Jalview
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S11 - Tajima's D for tubulins
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
ys <- c(-2.6, 0.5)
# III:3537688..3541628
st = 3537688
en = 3541628
ben1_d <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5, site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")
ben1td <- ben1_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 3541628/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# III:10740868..10742461
st = 10740868 
en = 10742461
tbb1_d <- tajimas_d_temp(vcf_path =paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")
tbb1td <- tbb1_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 10742461/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# III:4015769..4017643
st = 4015769
en = 4017643

tbb2_d <- tajimas_d_temp(vcf_path = paste0(final.dir),vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb2td <- tbb2_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 4017643/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# X:9434452..9436893
st = 9434452
en = 9436893
tbb4 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb4td <- tbb4[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 9436893/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

#V:12261804..12263368
st = 12261804
en = 12263368

tbb6 <- tajimas_d_temp(vcf_path = paste0(final.dir),vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "V", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb6td <- tbb6[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 12263368/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# X:7774859..7776689
st = 7774859
en = 7776689
mec7 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
mec7td <- mec7[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size =16))+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 7776689/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# combine all plots
plots <- cowplot::plot_grid(ben1td,tbb1td,tbb2td,tbb4td,tbb6td,mec7td, ncol = 1, labels = "AUTO")
plots

# now add the title
title <- ggdraw() + draw_label("Tajima's D for beta-tubulins", fontface='bold', size = 18)
plot_w_title <- plot_grid(title, plots, ncol=1, rel_heights=c(0.05, 1))

ytitle <- ggdraw() + draw_label("Tajima's D", fontface='bold', angle = 90, size = 16)

plot_grid(ytitle, plot_w_title, ncol=2, rel_widths=c(0.025, 1))

ggsave(paste0(plot.dir,"FS11_TajimaDfigure.png"), 
       dpi = 300,
       height = 10, 
       width = 10)

ggsave(paste0(plot.dir,"TajimaDfigure.pdf"), 
       dpi = 300,
       height = 10, 
       width = 10)

