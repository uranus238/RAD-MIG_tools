########################
library(vcfR)
vcf <- read.vcfR("example.vcf")

chrom <- create.chromR(vcf)
chromoqc(chrom)

vcf@gt
colnames(vcf@gt)

vcf_e <- vcf@gt %>%
  as.tibble() %>%
  select(-FORMAT)
vcf_e

########################
# https://github.com/shenglin-liu/vcf2sfs を参考にした
# 最も単純な1DSFSをvcfとpopmapから作成
rm(list = ls())
library(tidyverse)

# 関数の定義
calc_1dsfs <- function(pop_for_calc = "xx",
                       vcf_file = "example.vcf",
                       popmap_file = "example_popmap.txt") {
  vcf.gt <- read.table(vcf_file, sep = "\t", stringsAsFactors = FALSE) %>%
    select(-c(1:9)) %>%
    as.matrix()
  popmap <- read.table(popmap_file, sep = "\t", stringsAsFactors = FALSE) %>%
    select(2)
  index <- popmap %>%
    mutate(id = 1:nrow(popmap)) %>%
    filter(V2 %in% c(pop_for_calc)) %>%
    select(id) %>%
    c() %>%
    unlist(use.names = FALSE)
  nrow.vcf <- nrow(vcf.gt)
  ncol.vcf <- ncol(vcf.gt)
  chrom1 <- substring(vcf.gt, 1, 1)
  chrom2 <- substring(vcf.gt, 3, 3)
  chrom <- matrix(as.integer(chrom1) + as.integer(chrom2),
                  nrow = nrow.vcf,
                  ncol = ncol.vcf)
  chrom <- chrom[, index]
  count_per_snp <- rowSums(chrom, na.rm = TRUE)
  sfs <- tibble()
  for (Si in 0:(ncol(chrom) * 2)) {
    sfs[1, Si + 1] <- Si
    sfs[2, Si + 1] <- sum(count_per_snp == Si)  
  }
  return(sfs)
}

# 実行
calc_1dsfs(pop_for_calc = "Kib",
           vcf_file = "example.vcf",
           popmap_file = "example_popmap.txt")

##########################################
# 1dsfsを図示
# 1dsfsをプロットする関数の定義
# calc_1dsfs関数を使用するのでよみこんでおく必要がある
# 理論曲線はWF平衡集団で得られる期待値
# thetaは観測値から絶対値法を用いて推定（とくに根拠はない）

plot_1dsfs <- function(pop_for_calc = "Kib",
                       vcf_file = "example.vcf",
                       popmap_file = "example_popmap.txt",
                       plot_expected_value = TRUE) {
  sfs <- calc_1dsfs(pop_for_calc = pop_for_calc,
                    vcf_file = vcf_file,
                    popmap_file = popmap_file)
  sfs_t <- t(sfs)
  colnames(sfs_t) <- c("S", "Counts")
  
  i <- sfs_t[2:nrow(sfs_t), 1]
  Si_obs <- sfs_t[2:nrow(sfs_t), 2]
  
  # 関数を定義（理論値-観察値の絶対値を返す）
  # E(Si) = theta / i (for i in 1:n-1)
  # Fu (1995), Griffiths and Tavare (1998)
  f <- function(theta) {
    Si_exp <- theta / i
    return(sum(abs(Si_obs - Si_exp)))
  }
  
  # 定義した関数の返り値を最小にするパラメータthetaの探索(上とあわせていわゆる絶対値法の実装)
  # optim(par = 1, # initial values for the parameters to be optimized over
  #      fn = f) # a function to be minimized
  theta <- optim(1, f)$par[1] # the best set of parameters found
  
  smax <- nrow(sfs_t) - 1
  Si_exp <- tibble(x = c(1:(nrow(sfs_t) - 2)),
                 y = theta / x)
  
  # プロット
  
  p <- sfs_t %>%
    as.tibble() %>%
    ggplot() +
    geom_bar(aes(x = S, y = Counts),
             stat = "identity")
  # p
  
  p_no_invariants <- sfs_t %>%
    as.tibble() %>%
    filter(S != 0) %>%
    ggplot() +
    geom_bar(aes(x = S, y = Counts),
             stat = "identity") +
  #   if (plot_expected_value) {
  # p_no_invariants +
  #       geom_line(data = Si_exp,
  #                 aes(x = x, y = y)) +
  #       geom_point(data = Si_exp,
  #              aes(x = x, y = y))
  #       } else {}
  # p_no_invariants +
    geom_line(data = Si_exp,
              aes(x = x, y = y)) +
              geom_point(data = Si_exp,
                         aes(x = x, y = y)) +
    scale_x_continuous(breaks = seq(1, smax - 1, 1),
                       labels = seq(1, smax - 1, 1)) +
    labs(x = "Number of derived alleles in Population",
         y = "Number of sites") +
    theme_bw()
  
  p_no_invariants
}

plot_1dsfs(pop_for_calc = "Kib",
           vcf_file = "example.vcf",
           popmap_file = "example_popmap.txt",
           plot_expected_value = TRUE)

#######################################################

### fitting

Si_exp <- tibble(x = c(1:(nrow(sfs_t) - 2)),
                 y = theta / x)
Si_exp
theta <- 58


i <- sfs_t[2:33, 1]
Si_obs <- sfs_t[2:33, 2]

# 関数を定義（理論値-観察値の絶対値を返す）
# E(Si) = theta / i (for i in 1:n-1)
# Fu (1995), Griffiths and Tavare (1998)
f <- function(theta) {
  Si_exp <- theta / i
  return(sum(abs(Si_obs - Si_exp)))
}

# 定義した関数の返り値を最小にするパラメータthetaの探索(上とあわせていわゆる絶対値法の実装)
optim(par = 1, # initial values for the parameters to be optimized over
      fn = f) # a function to be minimized
theta <- optim(1, f)$par[1] # the best set of parameters found
theta

Si_e <- tibble(x = c(1:31),
               y = theta / x)
Si_e


# プロット
sfs
sfs_t <- t(sfs)
sfs_t
colnames(sfs_t) <- c("S", "Counts")
sfs_t 

p <- sfs_t %>%
  as.tibble() %>%
  ggplot() +
  geom_bar(aes(x = S, y = Counts),
           stat = "identity")
p

p_no_invariants <- sfs_t %>%
  as.tibble() %>%
  filter(S != 0) %>%
  ggplot() +
  geom_bar(aes(x = S, y = Counts),
           stat = "identity") +
  geom_line(data = si_e,
            aes(x = x, y = y)) +
  geom_point(data = si_e,
             aes(x = x, y = y)) +
  scale_x_continuous(breaks = seq(1, ncol(chrom) * 2, 1),
                     labels = seq(1, ncol(chrom) * 2, 1)) +
  labs(x = "Number of derived alleles in Population", y = "Number of sites") +
  theme_bw()
p_no_invariants
###############






popname_1dsfs <- "Kib" # population name for calculate 1D SFS

vcf.gt <- read.table("example.vcf", sep = "\t", stringsAsFactors = FALSE) %>%
  select(-c(1:9)) %>%
  as.matrix()

popmap <- read.table("example_popmap.txt", sep = "\t", stringsAsFactors = FALSE) %>%
  select(2)

nrow.vcf <- nrow(vcf.gt)
ncol.vcf <- ncol(vcf.gt)

chrom1 <- substring(vcf.gt, 1, 1)
chrom2 <- substring(vcf.gt, 3, 3)
chrom <- matrix(as.integer(chrom1) + as.integer(chrom2),
                nrow = nrow.vcf,
                ncol = ncol.vcf)

# gt <- list(popmap = popmap, genotype = chrom)

# 指定した集団に対応する個体のindex（行名）を抽出
index <- popmap %>%
  mutate(id = 1:nrow(popmap)) %>%
  filter(V2 %in% c(popname_1dsfs)) %>%
  select(id) %>%
  c() %>%
  unlist(use.names = FALSE)
index

# popmap
# popmap <- popmap[index, ]
# popmap
# 
# popmap <- popmap %>%
#   filter(V2 %in% c("Kib"))
# popmap

chrom <- chrom[, index]
chrom

 x <- rowSums(chrom, na.rm = TRUE)
 x

sfs <- tibble()
for (i in 0:(ncol(chrom) * 2)) {
sfs[1, i + 1] <- i
sfs[2, i + 1] <- sum(x == i)  
}
sfs

sfs_t <- t(sfs)
sfs_t
colnames(sfs_t) <- c("S", "Counts")
sfs_t 

p <- sfs_t %>%
  as.tibble() %>%
  ggplot() +
  geom_bar(aes(x = S, y = Counts),
           stat = "identity")
p

p_no_invariants <- sfs_t %>%
  as.tibble() %>%
  filter(S != 0) %>%
  ggplot() +
  geom_bar(aes(x = S, y = Counts),
           stat = "identity") +
  geom_line(data = si_e,
            aes(x = x, y = y)) +
  geom_point(data = si_e,
            aes(x = x, y = y)) +
  scale_x_continuous(breaks = seq(1, ncol(chrom) * 2, 1),
                     labels = seq(1, ncol(chrom) * 2, 1)) +
  labs(x = "Number of derived alleles in Population", y = "Number of sites") +
  theme_bw()
p_no_invariants

############################

for (i in 0:(ncol(chrom) * 2)) {
  print(i)
}


########################
# vcf2sfs_customed.R
library(tidyverse)

# function vcf2gt
# read vcf file and popmap file.
vcf.gt <- read.table("example.vcf", sep = "\t", stringsAsFactors = FALSE) %>%
  select(-c(1:9)) %>%
  as.matrix()

popmap <- read.table("example_popmap.txt", sep = "\t", stringsAsFactors = FALSE) %>%
  select(2)

nrow.vcf <- nrow(vcf.gt)
ncol.vcf <- ncol(vcf.gt)

chrom1 <- substring(vcf.gt, 1, 1)
chrom2 <- substring(vcf.gt, 3, 3)
chrom <- matrix(as.integer(chrom1) + as.integer(chrom2),
                nrow = nrow.vcf,
                ncol = ncol.vcf)

gt <- list(popmap = popmap, genotype = chrom)

# function choose.pops
# choose populations (subsetting the gt object by populations).
# pops: a character or integer vector; IDs of the chosen populations.

popmap <- gt$popmap # なくてもよい（直接読みこんだものと変わらない）
chrom <- gt$genotype # なくてもよい


##################################### 検討
popmap
chrom
chrom_t <- t(chrom) # 行列入れ替え
chrom_t

dim(popmap)
pop_chrom <- cbind(popmap, chrom_t)
pop_chrom
########################################


# 指定した集団に対応する個体のindex（行名）を抽出
index <- popmap %>%
  mutate(id = 1:nrow(popmap)) %>%
  filter(V2 %in% c("Kib")) %>%
  select(id) %>%
  c() %>%
  unlist(use.names = FALSE)
index

popmap
popmap <- popmap[index, ]
popmap

popmap <- popmap %>%
  filter(V2 %in% c("Kib"))
popmap
chrom <- chrom[, index]
chrom

x <- rowSums(chrom, na.rm = TRUE)
table(x)

#############################



# SFS based on raw count.
cnt <- matrix(0, nrow.vcf, n.pop)
cnt
ext <- list()
for(i in 1:n.pop) {
  index <- which(popmap == pops[i])
  cnt[, i] <- rowSums(chrom[,index], na.rm = TRUE)
  ext <- c(ext, list(0:ns.chr[i]))
}
ext<-as.matrix(expand.grid(ext))
cnt<-data.frame(rbind(cnt,ext))
sfs.raw<-table(cnt)-1
names(dimnames(sfs.raw))<-pops

sfs.raw















# function samSize
# Calculate sample sizes of the populations
popmap <- gt$popmap
chrom <- gt$genotype

pops <- unique(popmap)

sampleSizes <- popmap %>%
  group_by(V2) %>%
  summarize(n = n())
sampleSizes

# function minSamSize
# Calclate minimum sample sizes of the populations accounting for the missing values. ???
sampleSizes

# function viewMissing
## view the distribution of the missing values.
# help decide wheter some individulas or SNPs should be filtered out.

missing <- is.na(chrom)
missing
layout(matrix(c(1,2,1,3), 2,2))
image(x = 1:nrow(missing),
      y = 1:ncol(missing),
      z = missing,
      xlab = paste("SNPs: ", dim(missing)[1], sep=""),
      ylab = paste("Individuals: ", dim(missing)[2], sep=""),
      col = c("olivedrab3","firebrick3"), main = "Missing values (red)")
plot(sort(colSums(missing)),
     main="By individual",
     xlab="Individuals",
     ylab="Number of missing values")
plot(table(rowSums(missing)),
     main="By SNP",
     xlab="Number of missing values",
     ylab="Number of SNPs")

# filtering 集団ごとに何個体以上
chrom
missingVal <- is.na(chrom)
# mdelete.snp <- which(rowSums)

# function gt2sfs.raw
# generate a SFS from the gt object
# it will output a sfs based on row count without accounting for the missing values.
# pops: a character or integer vector; IDs of populations to be included in the SFS.

# n.pop <- length(pops)
n.pop <- nrow(pops)
n.pop
pops
# Number of chromosomes.
ns.chr <- sapply(pops, function(x) {sum(popmap==x)}) * 2

ns.chr <- sampleSize * 2

sampleSizes<-sapply(pops,function(x){sum(popmap==x)})
names(sampleSizes)<-pops

popmap
length(popmap)
chrom
x <- rowSums(chrom, na.rm = TRUE)
x

cnt <- c()
for(i in 0:length(popmap) * 2) {
  cnt[i] <- sum(x == i) 
}
cnt[0]

sum(x == 0)
table(x)

for(i in 1:length(popmap) * 2) {
  cnt[i] <-
}


# SFS based on raw count.
cnt <- matrix(0, nrow.vcf, n.pop)
cnt
ext <- list()
for(i in 1:n.pop) {
  index <- which(popmap == pops[i])
  cnt[, i] <- rowSums(chrom[,index], na.rm = TRUE)
  ext <- c(ext, list(0:ns.chr[i]))
}
ext<-as.matrix(expand.grid(ext))
cnt<-data.frame(rbind(cnt,ext))
sfs.raw<-table(cnt)-1
names(dimnames(sfs.raw))<-pops

sfs.raw


#### test original script
mygt <- vcf2gt("example.vcf", "example_popmap.txt")
mygt

mysfs1 <- gt2sfs.raw(mygt, "Kib")
mysfs1

mysfs2 <- gt2sfs.raw(mygt, c("Kap1", "Nor"))
mysfs2



