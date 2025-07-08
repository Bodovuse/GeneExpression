# to manipulate gene expression data. 
# breast cancer data set 30 pairs of tumor and normal samples
# setwd("~/Projects/R/RFiles/GeneExpression")

library(dplyr)
library(tidyverse)
library(GEOquery)
install.packages("BiocManager")
BiocManager::install("GEOquery")

#read data
dat <- read.csv(file = "GSE183947_fpkm.csv")

#get metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

view(metadata)

#select columns
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ","", tissue)) %>%#remove the tissue and metastasis prefix 
  mutate(metastasis = gsub("metastasis: ","", metastasis)) 

head(dat)

#reshaping data
dat.long <- dat %>%
  rename(gene =X) %>%
  gather(key = 'sample', value = 'FPKM', -gene)#convert wide to long format, exclude the gene column
  
#join dat.long with metadata.modified
dat.long <- left_join(dat.long, metadata.modified, by = c("sample"="description")) 

#data exploration
#looking at just BRCA1 BRCA2 gene mutation
dat.long %>%
  filter(gene == 'BRCA1'| gene == 'BRCA2' ) %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM) %>%
  head()

#visualising data

#bar plot
#compare the expression of samples for BRCA1 gene mutation
bar <- dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = sample, y = FPKM, fill = tissue)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle=90,size=7.5))

#desnity
density <- dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.5)
  
#box plot
box <- dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = tissue)) +
  geom_boxplot() 

#violin
violin <- dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = tissue)) +
  geom_violin()

#scatterplot
scatter <- dat.long %>%
  filter(gene == 'BRCA1'| gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM)%>%#want brca1 and 2 as own columns
  ggplot(., aes(x = BRCA1, y = BRCA2, colour = tissue)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

#heatmap
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

heat <- dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = sample, y = gene, fill = FPKM)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle=90,size=7.5)) +
  scale_fill_gradient(low = 'white', high = 'purple')

#save all graphs made
ggsave(heat, filename = 'heatmap_1.png', width = 10, height = 8)
ggsave(box, filename = 'boxplot_1.png', width = 10, height = 8)
ggsave(bar, filename = 'bargraph_1.png', width = 10, height = 8)
ggsave(density, filename = 'density_1.png', width = 10, height = 8)
ggsave(scatter, filename = 'scatter_1.png', width = 10, height = 8)
ggsave(violin, filename = 'violin_1.png', width = 10, height = 8)
