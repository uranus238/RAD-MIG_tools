# vcf2treemix

# dependency: vcfR, dplyr
# Usage: Rscript vcf2treemix.R <vcf_file> <popmap> <out_file>
# e.g. Rscript ./script/vcf2treemix.R ./output/tools/stacks_treemix/neutral.p.snps.vcf ./data/popmap/popmap_treemix.txt ./output/tools/stacks_treemix/test.treemix

suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(dplyr))

# get arguments from shell
args <- commandArgs(TRUE)

vcf_file <- args[1]
popmap <- args[2]
out_file <- args[3]

popmap <- read.table(popmap, header = FALSE, stringsAsFactors = FALSE)
pop <- unique(popmap$V2)
pop_ind <- list()

for (i in 1:length(pop)) {
pop_ind[pop[i]] <- popmap %>%
  filter(popmap$V2 == pop[i]) %>%
  select(V1)
}

vcf <- read.vcfR(vcf_file)
gt <- vcf@gt

ch1 <- gt[, -1] %>% substring(1, 1)
ch2 <- gt[, -1] %>% substring(3, 3)

ch1_0 <- ifelse(is.na(ch1), FALSE, ch1 == "0")
ch1_1 <- ifelse(is.na(ch1), FALSE, ch1 == "1")

ch2_0 <- ifelse(is.na(ch2), FALSE, ch2 == "0")
ch2_1 <- ifelse(is.na(ch2), FALSE, ch2 == "1")

gt_0 <- ch1_0 + ch2_0
gt_1 <- ch1_1 + ch2_1

count0 <- list()
count1 <- list()

for (i in 1:length(pop)) {
count0[pop[i]] <- gt_0 %>%
  as_tibble() %>%
  select(pop_ind[pop[i]] %>% unlist()) %>% 
  mutate(count = rowSums(.)) %>%
  select(count)

count1[pop[i]] <- gt_1 %>%
  as_tibble() %>%
  select(pop_ind[pop[i]] %>% unlist()) %>% 
  mutate(count = rowSums(.)) %>%
  select(count)
}

count_pop <- list()
treemix <- c()
for (i in 1:length(pop)) {
count_pop[[i]] <- paste(count0[pop[i]] %>% unlist() %>% as.matrix(), count1[pop[i]] %>% unlist() %>% as.matrix(), sep = ",")
treemix <- cbind(treemix, count_pop[[i]])
}
treemix <- rbind(pop, treemix)

treemix %>%
  write.table(file = out_file,
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
