library(tidyverse)
library(extrafont)
library(ggthemes)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

source("N2_CB.R")

theme_set(theme_bw())

sth <- readr::read_csv("../Data/Soil_transmitted_Helminths.csv")
colnames(sth) <- gsub(" ","_",colnames(sth))
lf <- readr::read_csv("../Data/Lymphatic_filariasis.csv")
colnames(lf) <- gsub(" ","_",colnames(lf))
country_pops <- readr::read_csv("../Data/country_populations")

# /*
#   =======================================
#   ~ > *                             * < ~
#   ~ ~ > *                         * < ~ ~
#   ~ ~ ~ > * PLOT DISTRIBUTIONS  * < ~ ~ ~ 
#   ~ ~ > *                         * < ~ ~ 
#   ~ > *                             * < ~
#   =======================================
# */
world_map <- map_data("world")
world_map <- world_map[world_map$region != "Antarctica",] # intercourse antarctica

base_pt <- ggplot()+
  geom_map(data = world_map, aes(long,lat, map_id=region), 
           map = world_map,
           fill = background_color,
           color = axis_color)+
  expand_limits(x = world_map$long, y = world_map$lat) 

yr <- 2016

sth_fill <- sth %>%
  dplyr::filter(Year == yr) %>%
  dplyr::select(fill_c = Population_requiring_PC_for_STH_SAC,
                Country) %>%
  dplyr::mutate(fill_c1 = as.numeric(gsub(",","",fill_c)))%>%
  dplyr::left_join(.,country_pops, by = "Country") %>%
  dplyr::filter(Year == yr) %>%
  dplyr::mutate(rate = fill_c1/Value*100)

base_pt + geom_map(data=sth_fill, 
                   map=world_map,
                   aes(fill=rate, map_id=Country), size =2) +
  theme_map()+
  scale_fill_gradient(low = "#0072B2",
                      high = "#F0E442",
                      na.value = background_color,
                      name = "Percent Infected",
                      breaks = c(0, 10, 20, 30),
                      guide = guide_colorbar(title.vjust = 0, 
                                             label.vjust = 1))+
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        legend.background =  element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        legend.text = element_text(size = 22, family = number_font),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm")) 

ggsave("../Plots/soil_helminth_2016_distribution.png", height = 10, width = 20)


overtime <- sth %>%
  dplyr::select(fill_c = Population_requiring_PC_for_STH_SAC,
                Country) %>%
  dplyr::mutate(fill_c1 = as.numeric(gsub(",","",fill_c)))%>%
  dplyr::left_join(.,country_pops, by = "Country") %>%
  dplyr::mutate(rate = fill_c1/Value*100)

overtime%>%
  dplyr::filter(rate<100)%>%
  na.omit()%>%
  ggplot()+
  aes(x = Year, y = rate)+
  geom_smooth(color = highlight_color)+
  ylim(0,100)+
  base_theme+
  labs(y = "Rate")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave("../Plots/soil_helminth_overtime.png",
       dpi = 300,
       height = 6, width = 8)


drugs_used <- sth %>%
  dplyr::select(Country, Year, 
                n_req=Population_requiring_PC_for_STH_SAC,
                n_targetd=Reported_number_of_SAC_treated,
                drugs = Drug_used_SAC)%>%
  dplyr::filter(n_req != "-",
                n_targetd != "-")%>%
  tidyr::separate(drugs, c("drugs", "round"),sep = " - ")%>%
  dplyr::mutate(alb = ifelse(grepl("Alb|Mbd", drugs), "BZ", "other"))%>%
  # dplyr::filter(alb == "BZ") %>%
  dplyr::group_by(Year, alb) %>%
  dplyr::summarise(n_targeted = sum(as.numeric(gsub(",","",n_targetd))),
                   n_req = sum(as.numeric(gsub(",","",n_req)), na.rm = T))

drugs_used%>%
  ggplot()+
  geom_smooth(se=F, color = highlight_color,
              aes(x = Year, y = n_targeted),
              data = dplyr::filter(drugs_used, alb == "BZ"))+
  geom_smooth(se=F, color = axis_color,
              aes(x = Year, y = n_req),
              data = dplyr::filter(drugs_used, alb == "BZ"))+
  # geom_point(color = axis_color,
  #             aes(x = Year, y = n_targeted),
  #             data = dplyr::filter(drugs_used, alb == "other", Year == 2016))+
  base_theme+
  labs(y = "People treated (Millions)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = axis_color))+
  scale_y_continuous(breaks=c(0e8,6e8,12e8),
                     labels=c("0", "600","1200"),
                     limits = c(0,1.4e9))+
  scale_x_continuous(breaks = c(2002, 2009, 2016),
                     limits = c(2002,2017))

ggsave("../Plots/soil_helminth_ABZ_use_overtime.png",
       dpi = 300,
       height = 6, width = 8)


# discrepancy between population data and reporting by WHO
# lf_fill <- lf %>%
#   dplyr::filter(Year == 2016) %>%
#   dplyr::select(fill_c = Population_requiring_PC_for_LF,
#                 Country) %>%
#   dplyr::mutate(fill_c1 = as.numeric(gsub(",","",fill_c)))%>%
#   dplyr::left_join(.,country_pops, by = "Country") %>%
#   dplyr::filter(Year == 2016) %>%
#   dplyr::mutate(rate = fill_c1/Value*100)
# 
# base_pt + geom_map(data=lf_fill, 
#                    map=world_map,
#                    aes(fill=rate, map_id=Country), size =1) +
#   theme_map()+
#   scale_fill_gradient(low = axis_color,
#                       high = highlight_color,
#                       na.value = background_color,
#                       name = "Percent Infected",
#                       breaks = c(0, 10, 20, 30),
#                       guide = guide_colorbar(title.vjust = 0, 
#                                              label.vjust = 1))+
#   theme(panel.background = element_rect(fill = background_color, colour = NA),
#         text = element_text(family = axes_title_font, size = axes_text_size),
#         legend.text = element_text(size = 18, family = number_font),
#         legend.title = element_text(size = 22),
#         legend.key.size = unit(1, "cm")) 
# 
# ggsave("../Plots/soil_helminth_2016_distribution.pdf", height = 10, width = 20)

# geographic distribution
ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="black", fill="gray85", size=0.05)+
  geom_point(data = strains_tim_hawaii, aes(x=long, y=lat),color = "black",fill="red",shape=21, size =3) +
  theme_map()+
  scale_fill_manual(values=strain_islands)+
  coord_cartesian(xlim = c(-161,-155), ylim = c(19,22.3))+
  theme(legend.position = "none")
