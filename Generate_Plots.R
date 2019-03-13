try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
arsenic_data <- "Arsenic/Data/"
etoposide_data <- "Etoposide/Data/"
ben1_data <- "201809_rappaport/Data/ben1_data/Final_Tables/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")




########################################################################################################################
# INTRO 
########################################################################################################################

######################################################################################################################## RIAIL GENOTYPES

df <- readr::read_tsv("Data/gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "N2", "CB4856")) %>%
  dplyr::mutate(low_sites = ifelse(sites < 100, TRUE, FALSE)) %>%
  dplyr::filter(grepl("QX", sample)) %>%
  dplyr::mutate(num_riail = as.numeric(gsub("QX", "",sample))) %>%
  dplyr::filter(num_riail > 249)

plot_riail_geno(df) +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = background_color),
        panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        axis.text = element_text(family = number_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")))

ggsave(filename = "Plots/RIAIL_Genotypes.png", height = 8, width = 12, dpi = 400)

########################################################################################################################
# ARSENIC
########################################################################################################################

######################################################################################################################## Linkage LOD PLOT

arsenic_linkage <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 11.tsv"))

arsenic_linkage %>%
  dplyr::filter(trait == ".PC1") %>%
  maxlodplot_edit() +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(filename = "Plots/Arsenic_PC1_Linkage.png", height = 4, width = 12, dpi = 400)

######################################################################################################################## Linkage LOD PLOT

arsenic_linkage_pheno <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 8.tsv"))

arsenic_linkage <- arsenic_linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# Plot Linkage Mapping - PxG
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)

pxgplot_edit(arsenic_cross, dplyr::filter(arsenic_linkage, trait == ".PC1"))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = "PC1", x = "RIAIL Genotype at Peak Marker") +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color)) +
  labs(title = "")
 
ggsave(filename = "Plots/Arsenic_PC1_Linkage_PxG.png", height = 6, width = 12, dpi = 400)

######################################################################################################################## NILs

arsenic_nil_pheno <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 13.tsv"))

arsenic_nil_pheno %>%
  dplyr::filter(Trait == "PC1",
                Strain%in%c("N2","CB4856","ECA414","ECA434"))%>%
  dplyr::mutate(strain1 = factor(Strain, levels = c("N2","CB4856","ECA414","ECA434","ECA581",
                                                    "ECA582","ECA589","ECA590","ECA591"),
                                 labels =c("N2",
                                           "CB4856", 
                                           "CB4856>N2\n II:5.75 - 8.02Mb", 
                                           "CB4856>N2\n II:7.83 - 9.66Mb",
                                           "N2(C78S)\nA", "N2(C78S)\nB",
                                           "CB4856(S78C)\nA","CB4856(S78C)\nB",
                                           "CB4856(S78C)\nC")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = -Value, 
      fill=Strain) +
  geom_beeswarm(alpha = .4,priority = "density",cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "#F9A227","CB4856" = "#2790F9",
                               "ECA414" = "#2790F9","ECA434" = "#2790F9"))+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank(),
        plot.title = element_blank()) +
  labs( y = paste0("PC1"))

ggsave(filename = "Plots/Arsenic_PC1_NIL.png", height = 6, width = 12, dpi = 400)

######################################################################################################################## GWA MANHATTAN PLOT

arsenic_gwa <- data.table::fread(glue::glue("{arsenic_data}Figure 2-source data 4.tsv"))
independent_tests <- 500
gwas_fine_mapping <- readr::read_tsv(glue::glue("{arsenic_data}Figure 2-source data 8.zip")) 
# geno_matrix <- readr::read_tsv(glue::glue("{arsenic_data}Figure 2-source data 8.zip"))%>%
#   na.omit()

cegwas2_manplot(arsenic_gwa, eigen_cutoff = -log10(0.05/independent_tests))[[1]] + 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(title = "")

ggsave(filename = "Plots/Arsenic_PC1_GWA.png", height = 4, width = 12, dpi = 400)

######################################################################################################################## PxG

peak_pos <- na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::pull(facet_marker) %>%
  unique()

na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::group_by(allele, facet_marker)%>%
  dplyr::mutate(mean_pheno = mean(as.numeric(value), na.rm = T))%>%
  dplyr::mutate(n2_cb = case_when(
    strain == "N2" ~ "N2",
    strain == "CB4856" ~ "CB4856", 
    TRUE ~ "Other"
  )) %>%
  ggplot()+
  aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
  geom_boxplot(aes(y=as.numeric(value)), alpha = 0.5, fill = "gray70") +
  geom_beeswarm(cex=2,priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21)+
  scale_fill_manual(values=strain_colors)+
  scale_size_manual(values=c(4,4,2))+
  labs(y = "PC1",
       x = glue::glue("Genotype at {peak_pos}")) +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Arsenic_PC1_GWA_PxG.png", height = 4, width = 6, dpi = 400)



########################################################################################################################
# BEN-1
########################################################################################################################

# # # LOAD AND PROCESS INDELS
ben1_variants <- data.table::fread(glue::glue("{ben1_data}TS9_ben-1_variants.tsv"))%>%
  na.omit()

gwa_mappings <-  data.table::fread(file = glue::glue("{ben1_data}TS5_GWA_processed_marker_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")