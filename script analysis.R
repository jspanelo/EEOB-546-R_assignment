install.packages("plotly")
install.packages("htmlwidgets")
library(tidyverse)
#Built in functions
source("functions_R_assignment.R")
## Data Inspection

# Broad inspection of the files

# loop for going through the snp and genotype files for getting basic info. Those are the only txt files in the wd
for( i in list.files(pattern = '*.txt')){ # match all text files to iterate through them
  #The output of this loop is some information as a string, easy to read.
      print(paste0('The file ', i, # print() statement for getting an output after the loop; paste0() to built the string   
                   # read.table() to read the appropriate file in each step of the loop
                   # nrow() and ncol() for counting rows and columns respectively in each step of the loop
                   '  is ', nrow(read.table(i, sep = '\t')), ' rows long and ', 
                   ncol(read.table(i, sep = '\t')), ' columns wide.'))
}

### Loading the files as genotypes ans snp_pos ### 

snp_pos <- read_delim("snp_position.txt", delim = '\t')
head(snp_pos)

# summary for snp_pos
snp_pos %>% filter(Position < 10e1000) %>% 
  group_by(Chromosome = as.numeric(Chromosome)) %>%  
  summarise(SNPs = n(),
            First_Pos_Mb = (min(as.double(Position)/1000000)), 
            Last_Pos_Mb = (max(as.double(Position)/1000000)),
            Coverage_Mb = Last_Pos_Mb - First_Pos_Mb)

# Some extra information about different factors in snp_pos

str(snp_pos)

Info <- NULL

for (i in 1:ncol(snp_pos)) {
  rbind(Info, c(names(snp_pos)[i], 
        length(levels(as.factor(snp_pos[[i]]))))) -> Info
  Info %>% as.table()
  }

as_tibble(Info) %>% rename(., Feature = "V1", `Counts or Levels` = "V2")

# 6 out of 15 variables are numerical (double), 9 out of 15 are characters

# There are 983 SNPs_IDs, 941 Positions (suggesting duplicated SNP_IDs) and 12 Chromosomes listed, which are supposed to be 10, and numerical 
            
# Checking levels and counts in Chromosomes by creating a table for chromosomes

as.data.frame(table(snp_pos$Chromosome, deparse.level = 2), 
              responseName = 'Counts')

# After checking "Chromosome", multiple and unknown chromosomes needs to be filtered out for the analysis
# There are 10 chromosomes and two sets of unassigned SNPs (multiple and unknown, 33 SNPs)
# The largest number of SNPs are on chr1 (155), the smallest on chr10 (53)

# Load fang_et_al_genotypes.txt and get some information

genotypes <- read_delim("fang_et_al_genotypes.txt", delim = '\t')

Info_G <- NULL

for (i in 1:ncol(genotypes)) {
  rbind(Info_G, c(names(genotypes)[i],
    length(levels(as.factor(genotypes[[i]]))))) -> Info_G
  Info_G %>% as.table()
  }
# Summarizyng the data based on the number of counts of each feature
as_tibble(Info_G) %>% rename(., Feature = "V1", Counts = "V2")

# Selecting only groups belonging to teosinte and maize, for extracting some data to be used afterwards
# Counts of sample size and number of SNPs for both teosinte and maize

as_tibble(cbind(c('Teosinte','Maize'), rbind(dim(genotypes[genotypes$Group %in% c('ZMPBA','ZMPIL','ZMPJA'),]), 
                                             dim(genotypes[genotypes$Group %in% c('ZMMIL','ZMMLR','ZMMMR'),])))) %>% 
  rename(Species = 'V1', SNPs = "V3", Samples = "V2")

# The counts are 986 SNPs for both sets, with 975 samples from Teosinte and 1573 from Maize

## Data Processing
#Create a directory for future outputs
dir.create('./output')
dir.create('./plots')

#Filtering teosinte's group

genotypes %>% 
  filter(Group %in% c('ZMPBA','ZMPIL','ZMPJA')) %>% 
  select(-c(2:3)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() -> teosinte
names(teosinte) <- teosinte[1,]
teosinte <- teosinte[-1,]
teosinte <- rename(teosinte, SNP_ID = "Sample_ID")

#Filtering maize's group

genotypes %>% 
  filter(Group %in% c('ZMMIL','ZMMLR','ZMMMR')) %>% 
  select(-c(2:3)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() -> maize
names(maize) <- maize[1,]
maize <- maize[-1,]
maize <- rename(maize, SNP_ID = "Sample_ID")


# Merging genotypic and SNP data from teosinte

left_join(teosinte, select(snp_pos, c('SNP_ID', 'Chromosome', 'Position')), by = "SNP_ID") %>%
  select(c('SNP_ID', 'Chromosome', 'Position'), everything()) %>% 
  filter(Chromosome %in% c(1:10) & Position < 10e1000) %>% droplevels() -> teosinte

# Merging genotypic and SNP data from maize

left_join(maize, select(snp_pos, c('SNP_ID', 'Chromosome', 'Position')), by = "SNP_ID") %>%
    select(c('SNP_ID', 'Chromosome', 'Position'), everything()) %>% 
    filter(Chromosome %in% c(1:10) & Position < 10e1000) %>% droplevels() -> maize

## Files generation. Created using built-in functions under functions_R_assignment.R
# Files generation

for (i in 1:10) {
  write_chrom(maize, i, 'maize')
  write_chrom(teosinte, i, 'teo')
  write_chrom_reverse(maize, i, 'maize')
  write_chrom_reverse(teosinte, i, 'teo')
}

#### Plot SNPs ####

# SNPs Counts #

# SNPs Counts from snp_pos file. The the filter() is for getting rid of unknown and multiple snps 

snp_pos %>% 
  filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Chromosome))) +
  geom_bar(fill = 'orange', color = 'darkred') + 
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3) +
  scale_x_continuous(breaks = 1:10) +
  theme_replace() +
  ggtitle("SNPs count by Chromosome") +
  ylab('Number of SNPs') +
  xlab('Chromosome') 
ggsave('./plots/SNPs_count.pdf')

# SNPs distribution. Position is divided by 1000000 just for getting numbers in MegaBases

snp_pos %>% filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Position)/1000000)) +
  geom_histogram(aes(y = ..density..), color = 'orange', fill = "orange", alpha = 0.4, bins = 20) + 
  geom_density(aes(as.double(Position)/1000000), color = "darkred") + 
  facet_wrap(~ as.double(Chromosome), scales = "free_x") +
  theme_replace() +
  ggtitle("SNPs distribution by Chromosome") +
  xlab('Genome position (Mb)') +
  ylab('SNP density')
ggsave(paste0("./plots/SNP_distribution_by_chrom.pdf"))


# Homozygotes, Heterozygotes and Missing data by sample
# I. Wrangling data by sample
geno_long <- 
genotypes %>% select(-JG_OTU, -Group) %>%   
  pivot_longer(!Sample_ID) %>% 
  mutate(Locus = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote')))  

### II. Plotting data by sample. Plots generated using plotly (interactive) in order to make easier to visualize individual values

color_plots <- c("#009E73",  "#D55E00", "#999999") # Colorblind-friendly colors for plots

plotly::ggplotly(
  geno_long %>% group_by(Sample_ID) %>%  count(Locus) %>% 
    ggplot(aes(fill = Locus, y = n, x = Sample_ID)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_plots) +
    ggtitle("Proportion of Homozygotes, Heterozygotes \nand Missing Data by Genotype ") +
    ylab('Proportion') +
    theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
) -> P
  
setwd('./plots')
htmlwidgets::saveWidget(P, "Proportions_by_genotype.html")
setwd('../')

# Homozygotes, Heterozygotes and Missing data by Group
# I. Wrangling data by group

groups_long <- 
  genotypes %>% select(-JG_OTU, -Sample_ID) %>%   
  pivot_longer(!Group) %>% 
  mutate(Locus = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote'))) 

# II. Plotting data by group
  
groups_long %>% group_by(Group) %>%  count(Locus) %>% 
    ggplot(aes(fill = Locus, y = n, x = Group)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_plots) +
    ggtitle("Proportion of Homozygotes, Heterozygotes and Missing Data \nby Group ") +
    ylab('Proportion') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('./plots/Proportions_by_group.pdf')


### Own plot: Proportion of nucleotides for homozygotic snps for each group ###

color_plots <- c("#009E73", "#999999","#D55E00", "#545454")

groups_long %>% filter(Locus == "Homozygote") %>% group_by(Group) %>%  count(value) %>% 
  ggplot(aes(fill = value, y = n, x = Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_plots) +
  ggtitle("Proportion of nucleotides for Homozygotic sites in each group") +
  ylab('Proportion') +
  theme(axis.text.x = element_text(angle = 90))
ggsave('./plots/own_plot.pdf')
