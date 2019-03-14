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

arsenic_linkage_pheno <- arsenic_linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# Plot Linkage Mapping - PxG
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage_pheno, set = 2)

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
  geom_beeswarm(alpha = point_alpha,
                priority = "density",
                cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227","CB4856" = "#2790F9",
                               "ECA414" = "#2790F9","ECA434" = "#2790F9"))+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank(),
        plot.title = element_blank()) +
  labs( y = paste0("PC 1"))

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

######################################################################################################################## GWA PxG

peak_pos <- na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::pull(facet_marker) %>%
  unique()

pxg_df <- na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::group_by(allele, facet_marker)%>%
  dplyr::mutate(mean_pheno = mean(as.numeric(value), na.rm = T))%>%
  dplyr::mutate(n2_cb = case_when(
    strain == "N2" ~ "N2",
    strain == "CB4856" ~ "CB4856", 
    TRUE ~ "Other"
  )) 

pxg_df %>%
  ggplot()+
  aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
  geom_beeswarm(cex=1.2,
                priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21, 
                alpha = point_alpha,
                data = dplyr::filter(pxg_df, !strain %in% c("CB4856", "N2")))+
  geom_boxplot(aes(y=as.numeric(value)), alpha = boxplot_alpha, fill = "gray70", outlier.colour = NA) +
  geom_beeswarm(cex=1.2,
                priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21, 
                data = dplyr::filter(pxg_df, strain %in% c("CB4856", "N2")))+
  scale_fill_manual(values=strain_colors)+
  scale_size_manual(values=c(point_highlight_size,point_highlight_size,point_size))+
  labs(y = "PC 1",
       x = glue::glue("Genotype at {peak_pos}")) +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Arsenic_PC1_GWA_PxG.png", height = 6, width = 8, dpi = 400)

######################################################################################################################## GWA PEAK LD
gm <- readr::read_tsv(glue::glue("{arsenic_data}Figure 2-source data 5.tsv"))
arsenic_gwa <- data.table::fread(glue::glue("{arsenic_data}Figure 2-source data 4.tsv"))

LD_output <- Plot_Peak_LD(arsenic_gwa, gm)

LD_output[[1]] + 
  base_theme + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = "Plots/Arsenic_PC1_GWA_PeakLD.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic_PC1_GWA_PeakLD.pdf", height = 8, width = 12, dpi = 400)

######################################################################################################################## SWAP


arsenic_swap <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 13.tsv"))

arsenic_swap%>%
  dplyr::filter(Trait == "PC1",
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(Strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = -Value, 
      fill= Strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA581" = "gray50","ECA590" = "gray50"))+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("PC 1"))

ggsave(filename = "Plots/Arsenic_PC1_SWAP.png", height = 6, width = 10, dpi = 400)

######################################################################################################################## Metabolites ratio in arsenic


metabolites <- data.table::fread(glue::glue("{arsenic_data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{arsenic_data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value)

arsenic100_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Arsenic")%>%
  dplyr::left_join(.,controls_new, by = c("Strain", "FA", "Replicate")) %>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::filter(FA %in% c("17_ratio", "15_ratio"))%>%
  dplyr::select(Strain, FA, control_value, delta_control)

arsenic100_new$FA <- gsub("17_ratio", "C17iso/C17n", arsenic100_new$FA)
arsenic100_new$FA <- gsub("15_ratio", "C15iso/C15n", arsenic100_new$FA)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("rC15", "rC17"), concentration != "200") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$compound <- gsub("rC15", "C15iso/C15n", rC15$compound)
rC15$compound <- gsub("rC17", "C17iso/C17n", rC15$compound)
rC15$concentration <- gsub("Mock", "Water", rC15$concentration)
rC15$concentration <- gsub("100", "Arsenic", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)

arsenic100 <- dplyr::filter(rC15, concentration == "Arsenic")%>%
  dplyr::left_join(.,controls, by = c("strain", "compound","replicate"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(Strain = strain, FA = compound, control_value, delta_control) 

complete_arsenic <- dplyr::bind_rows(arsenic100,arsenic100_new) %>%
  dplyr::group_by(FA) %>%
  dplyr::mutate(scale_delta = scale(delta_control)) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_scale_value = mean(scale_delta),
                sd_scale_value = sd(scale_delta),
                mean_value = mean(delta_control),
                sd_value = sd(delta_control),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), 
                    ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs(y = "Arsenic - Control")

ggsave(filename = "Plots/Arsenic_PC1_ISOratio_Arsenic.png", height = 10, width = 10, dpi = 400)

######################################################################################################################## Metabolites straight in arsenic
controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value)

arsenic100_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Arsenic")%>%
  dplyr::left_join(.,controls_new, by = c("Strain", "FA", "Replicate")) %>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::filter(FA %in% c("C17n", "C15n"))%>%
  dplyr::select(Strain, FA, control_value, delta_control)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("C17n", "C15n"), concentration != "200") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$concentration <- gsub("Mock", "Water", rC15$concentration)
rC15$concentration <- gsub("100", "Arsenic", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)

arsenic100 <- dplyr::filter(rC15, concentration == "Arsenic")%>%
  dplyr::left_join(.,controls, by = c("strain", "compound","replicate"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(Strain = strain, FA = compound, control_value, delta_control) 

complete_arsenic <- dplyr::bind_rows(arsenic100,arsenic100_new) %>%
  dplyr::group_by(FA) %>%
  dplyr::mutate(scale_delta = scale(delta_control)) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_scale_value = mean(scale_delta),
                sd_scale_value = sd(scale_delta),
                mean_value = mean(delta_control),
                sd_value = sd(delta_control),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Arsenic - Control")

ggsave(filename = "Plots/Arsenic_PC1_StraightChain_Arsenic.png", height = 10, width = 10, dpi = 400)

######################################################################################################################## Metabolites ratios in control

controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value) %>%
  dplyr::filter(FA %in% c("15_ratio", "17_ratio"))%>%
  dplyr::select(Strain, FA, control_value)

controls_new$FA <- gsub("17_ratio", "C17iso/C17n", controls_new$FA)
controls_new$FA <- gsub("15_ratio", "C15iso/C15n", controls_new$FA)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("rC17", "rC15"), concentration == "Mock") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$compound <- gsub("rC15", "C15iso/C15n", rC15$compound)
rC15$compound <- gsub("rC17", "C17iso/C17n", rC15$compound)
rC15$concentration <- gsub("Mock", "Water", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)%>%
  dplyr::ungroup()%>%
  dplyr::select(Strain = strain, FA = compound, control_value)


complete_arsenic <- dplyr::bind_rows(controls,controls_new) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_value = mean(control_value),
                sd_value = sd(control_value),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Metabolite Levels")

ggsave(filename = "Plots/Arsenic_PC1_ISO_Control.png", height = 10, width = 10, dpi = 400)

######################################################################################################################## RESCUE
arsenic_rescue <- data.table::fread(glue::glue("{arsenic_data}Figure 4-source data 5.tsv"))

rescue_pheno_pr <- arsenic_rescue %>%
  dplyr::group_by(Strain, Condition, Trait)%>%
  dplyr::mutate(mph = median(Value),
                sph = sd(Value))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(Value >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(Value <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")%>%
  dplyr::ungroup()%>%
  dplyr::select(-cut_l, -cut_h,-flag_l,-flag_h,-sph,-mph)

mc <- "64"
cc <- c("Arsenic",
        "ArsenicC15ISO64")

boxplot_plt(df = rescue_pheno_pr,
            trt = "PC1",
            cond = cc,
            fancy_name = paste0("PC 1"),
            strains = c("N2",
                        "CB4856",
                        "ECA581",
                        "ECA590"),
            fancy_strains = c("Bristol",
                              "Hawaii",
                              "Bristol\n(C78S)",
                              "Hawaii\n(S78C)"),
            ordered_conditions = c("EtOH", 
                                   "C15ISO12",
                                   "C15ISO24", 
                                   "C15ISO48",
                                   "C15ISO64", 
                                   "C15ISO100", 
                                   'Arsenic', 
                                   "ArsenicC15ISO12",
                                   "ArsenicC15ISO24",
                                   "ArsenicC15ISO48",
                                   "ArsenicC15ISO64",
                                   "ArsenicC15ISO100"),
            fancy_ordered_conditions = c("EtOH", 
                                         "C15ISO\n12",
                                         "C15ISO\n24", 
                                         "C15ISO\n48",
                                         "C15ISO\n64", 
                                         "C15ISO\n100", 
                                         'Arsenic', 
                                         "Arsenic\nC15ISO\n12",
                                         "Arsenic\nC15ISO\n24",
                                         "Arsenic\nC15ISO\n48",
                                         "Arsenic\nC15ISO\n64",
                                         "Arsenic\nC15ISO\n100"),
            r_conc = mc)+ 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) 

ggsave(filename = "Plots/Arsenic_PC1_Rescue.png", height = 8, width = 14, dpi = 400)

######################################################################################################################## PopGene

load("Data/Arsenic_II_7598325-8210489_Diversity_Statistics.Rda")

qtl_start <- 7598325
qtl_end <- 8210489
dbt_start <- 7942463
dbt_end <- 7945206

td_df <- neutrality_df %>%
  dplyr::filter(statistic == "Tajima.D") 

td_df %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value)+
  annotate("rect",xmin=qtl_start/1e6, 
           xmax=qtl_end/1e6, 
           ymin=-Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="cyan") +
  geom_point(size = point_size, alpha = point_alpha)+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Tajima's D")

load("Data/Arsenic_test_interval_II_7598325-8210489_Diversity_Statistics.Rda")

td_df <- neutrality_df %>%
  dplyr::filter(statistic %in% c("Fay.Wu.H","Zeng.E","Tajima.D")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value))

td_df %>%
  dplyr::filter(WindowPosition > qtl_start) %>%
  dplyr::filter(WindowPosition < qtl_end) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  annotate("rect",
           xmin=dbt_start/1e6, 
           xmax=dbt_end/1e6, 
           ymin=-Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="cyan") +
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = col_blind_colors)+
  facet_grid(statistic~., scales = "free")+
  base_theme +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)")

td_df %>%
  dplyr::filter(WindowPosition > qtl_start) %>%
  dplyr::filter(WindowPosition < qtl_end) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  annotate("rect",
           xmin=dbt_start/1e6, 
           xmax=dbt_end/1e6, 
           ymin=-Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="cyan") +
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = col_blind_colors)+
  facet_grid(statistic~., scales = "free")+
  base_theme +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)")+
  xlim(c((dbt_start-50000)/1e6, (dbt_end+50000)/1e6))




