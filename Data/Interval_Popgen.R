#!/usr/bin/env Rscript
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

require(PopGenome)
require(WhopGenome)
require(data.table)
require(tidyverse)

# 1 - chromosome
# 2 - start of window
# 3 - end of window
# 4 - VCF file
# 5 - GFF file
# 6 - Save Label
# & - Strain file
# test args
# args <- c("II","7598325","8210489","Ce330_annotated.vcf.gz","WS245_exons.gff", "Arsenic", "249_samples.txt")
# Running script
# Rscript --vanilla Interval_Popgen.R II 10050022 12062611 Ce330_annotated.vcf.gz WS245_exons.gff Etoposide 249_samples.txt
args <- commandArgs(TRUE)

system(glue::glue("echo Initializing PopGenome Parameters"))

chr1 <- c(1,15072434)
chr2 <- c(1,15279421)
chr3 <- c(1,13783801)
chr4 <- c(1,17493829)
chr5 <- c(1,20924180)
chr6 <- c(1,17718942)

chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6)
chroms <- c("I","II","III","IV","V","X")

ANALYSIS_CHROM <- as.character(args[1])
CHROM_START <- chr.lengths[which(chroms == as.character(args[1]))][[1]][1]
CHROM_END <- chr.lengths[which(chroms == as.character(args[1]))][[1]][2]

SLIDE_DISTANCE <- 100
WINDOW_SIZE <-  10000

REGION_START <- as.numeric(args[2])
REGION_END <- as.numeric(args[3])

OUTGROUP <- "XZ1516"

system(glue::glue("echo Done Initializing PopGenome Parameters - WindowSize = {WINDOW_SIZE}, StepSize = {SLIDE_DISTANCE}, Whole Population, Chromosome = {ANALYSIS_CHROM}"))

system(glue::glue("bcftools view -S {args[7]} -r {ANALYSIS_CHROM}:{args[2]}-{args[3]} -v snps {args[4]} | sed 's/0\\/0/0|0/g' | sed 's/1\\/1/1|1/g' | sed 's/0\\/1/1|1/g' | sed 's/1\\/0/1|1/g' | sed 's/.\\/./.|./g' | bcftools view -Oz -o tempRegion.vcf.gz"))
system(glue::glue("tabix tempRegion.vcf.gz"))
system(glue::glue("bcftools query -l tempRegion.vcf.gz > samples.txt"))

samples <- data.table::fread("samples.txt",header = F) %>%
  dplyr::pull(V1)

POPGENOME_VCF <- "tempRegion.vcf.gz"
POPGENOME_GFF <- args[5]

system(glue::glue("echo PopGenome - Reading VCF file {ANALYSIS_CHROM}")) 

vcf_handle <- WhopGenome::vcf_open(POPGENOME_VCF)

GENOME_OBJECT <- Whop_readVCF(
  vcf_handle, 
  numcols = 100, 
  tid = ANALYSIS_CHROM, 
  from = REGION_START,
  to = REGION_END)

system(glue::glue("echo PopGenome - Setting Outgroup and Defining Window Size"))

GENOME_OBJECT <- PopGenome::set.populations(GENOME_OBJECT, list(samples), diploid = FALSE)

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- sliding.window.transform(GENOME_OBJECT, WINDOW_SIZE, SLIDE_DISTANCE, type=2)

system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

save(GENOME_OBJECT, file = glue::glue("{args[6]}_{args[1]}_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Window Bins"))

windowStarts <- data.frame(snp_index = 1:length(colnames(GENOME_OBJECT@BIG.BIAL[[1]])),
                           position = colnames(GENOME_OBJECT@BIG.BIAL[[1]]))

slide_index <- cbind(data.frame(lapply(GENOME_OBJECT@SLIDE.POS,  function(x) as.numeric(floor(mean(x)))))) %>%
  tidyr::gather(temp, snp_index) %>%
  dplyr::select(-temp) %>%
  dplyr::left_join(., windowStarts, by = "snp_index")

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Neutrality and Diversity Stats"))

neutrality_df <- data.frame(Tajima.D = c(GENOME_OBJECT@Tajima.D),
                            n.segregating.sites = c(GENOME_OBJECT@n.segregating.sites),
                            Fu.Li.F = c(GENOME_OBJECT@Fu.Li.F),
                            Fu.Li.D = c(GENOME_OBJECT@Fu.Li.D),
                            Fay.Wu.H = c(GENOME_OBJECT@Fay.Wu.H),
                            Zeng.E = c(GENOME_OBJECT@Zeng.E),
                            theta_Tajima = c(GENOME_OBJECT@theta_Tajima),
                            theta_Watterson = c(GENOME_OBJECT@theta_Watterson),
                            theta_Achaz.Watterson = c(GENOME_OBJECT@theta_Achaz.Watterson),
                            theta_Achaz.Tajima = c(GENOME_OBJECT@theta_Achaz.Tajima),
                            theta_Fay.Wu = c(GENOME_OBJECT@theta_Fay.Wu),
                            theta_Zeng = c(GENOME_OBJECT@theta_Zeng),
                            nuc.diversity.within = c(GENOME_OBJECT@nuc.diversity.within),
                            PI = c(GENOME_OBJECT@Pi),
                            theta_Fu.Li = c(GENOME_OBJECT@theta_Fu.Li),
                            hap.diversity.within = c(GENOME_OBJECT@hap.diversity.within)) %>%
  dplyr::mutate(Population = "WHOLE_POPULATION",
                WindowPosition = slide_index$position) %>%
  tidyr::gather(statistic, value, -Population, -WindowPosition)

save(neutrality_df, file = glue::glue("{args[6]}_{args[1]}_{args[2]}-{args[3]}_Diversity_Statistics.Rda"))

system("rm tempRegion.vcf.gz")
system("rm tempRegion.vcf.gz.tbi")
system("rm samples.txt")
