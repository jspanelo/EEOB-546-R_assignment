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

genotypes <- read_delim("fang_et_al_genotypes.txt", delim = '\t')
snp_pos <- read_delim("snp_position.txt", delim = '\t')

# summary for snp_pos
snp_pos %>% filter(Position < 10e1000) %>% 
  group_by(Chromosome = as.numeric(Chromosome)) %>%  
  summarise(SNPs = n(),
            First_Pos_Mb = (min(as.double(Position)/1000000)), 
            Last_Pos_Mb = (max(as.double(Position)/1000000)),
            Coverage_Mb = Last_Pos_Mb - First_Pos_Mb)

# Some extra information about different factors in snp_pos

for (i in 1:ncol(snp_pos)) {
  print(names(snp_pos)[i])
  print(length(levels(as.factor(snp_pos[[i]]))))
}

levels(as.factor(snp_pos$Chromosome))
# Digging into Chromosome, multiple and unknown chromosomes needs to be filtered out

for (i in 1:3) {
  print(names(genotypes)[i])
  print(length(levels(as.factor(genotypes[[i]]))))
}

# Selecting only groups belonging to teosinte and maize

as_tibble(cbind(c('Teosinte','Maize'), rbind(dim(genotypes[genotypes$Group %in% c('ZMPBA','ZMPIL','ZMPJA'),]), 
                                             dim(genotypes[genotypes$Group %in% c('ZMMIL','ZMMLR','ZMMMR'),])))) %>% 
  rename(Species = 'V1', SNPs = "V3", Samples = "V2")


## Data Processing
dir.create('./output')


#Filtering teosinte

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

#Filtering maize

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


# Merge genotypic and SNP data

left_join(teosinte, select(snp_pos, c('SNP_ID', 'Chromosome', 'Position')), by = "SNP_ID") %>%
  select(c('SNP_ID', 'Chromosome', 'Position'), everything()) %>% 
  filter(Chromosome %in% c(1:10) & Position < 10e1000) %>% droplevels() -> teosinte

left_join(maize, select(snp_pos, c('SNP_ID', 'Chromosome', 'Position')), by = "SNP_ID") %>%
    select(c('SNP_ID', 'Chromosome', 'Position'), everything()) %>% 
    filter(Chromosome %in% c(1:10) & Position < 10e1000) %>% droplevels() -> maize

# Files generation

for (i in 1:10) {
  write_chrom(maize, i, 'maize')
  write_chrom(teosinte, i, 'teo')
  write_chrom_reverse(maize, i, 'maize')
  write_chrom_reverse(teosinte, i, 'teo')
}

#### Plot SNPs ####

# SNPs Counts #

snp_pos %>% 
  filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Chromosome))) +
  geom_bar(fill = 'orange', color = 'darkred') + 
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1) +
  scale_x_continuous(breaks = 1:10) +
  theme_replace() +
  ggtitle("SNPs count by Chromosome") +
  ylab('Number of SNPs') +
  xlab('Chromosome') 
ggsave('./output/SNPs_count.jpg')

# SNPs distribution #

snp_pos %>% filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Position)/1000000)) +
  geom_histogram(aes(y = ..density..), color = 'orange', fill = "orange", alpha = 0.4, bins = 20) + 
  geom_density(aes(as.double(Position)/1000000), color = "darkred") + 
  facet_wrap(~ as.double(Chromosome), scales = "free") +
  theme_replace() +
  ggtitle("SNPs distribution by Chromosome") +
  xlab('Genome position (Mb)') +
  ylab('SNP density')
ggsave(paste0("./output/SNP_distribution_by_chrom.jpg"))

### Wrangling data by sample ###

geno_long <- 
genotypes %>% select(-JG_OTU, -Group) %>%   
  pivot_longer(!Sample_ID) %>% 
  mutate(Locus = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote')))  

### Plotting by Sample ###

color_plots <- c("#009E73",  "#D55E00", "#999999")

geno_long %>% group_by(Sample_ID) %>%  count(Locus) %>% 
  ggplot(aes(fill = Locus, y = n, x = Sample_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_plots) +
  ggtitle("Proportion of Homozygotes, Heterozygotes and Missing Data by Genotype ") +
  ylab('Proportion') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('./output/Proportions_by_genotype.jpg')

### Wrangling Data by Group ###

groups_long <- 
  genotypes %>% select(-JG_OTU, -Sample_ID) %>%   
  pivot_longer(!Group) %>% 
  mutate(Locus = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote')))  

#### Plot by group ####
  
groups_long %>% group_by(Group) %>%  count(Locus) %>% 
    ggplot(aes(fill = Locus, y = n, x = Group)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_plots) +
    ggtitle("Proportion of Homozygotes, Heterozygotes and Missing Data by Group ") +
    ylab('Proportion') 
ggsave('./output/Proportions_by_group.jpg')


### Plotting by Sample ###

color_plots <- c("#009E73", "#999999","#D55E00", "#545454")

groups_long %>% filter(Locus == "Homozygote") %>% group_by(Group) %>%  count(value) %>% 
  ggplot(aes(fill = value, y = n, x = Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_plots) +
  ggtitle("Proportion of nucleotides for Homozygotic sites in each group") +
  ylab('Proportion') +
  theme_bw()
ggsave('./output/own_plot.jpg')
