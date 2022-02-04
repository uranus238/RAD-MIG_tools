# calculate sequence similariy

library(Rphylip)
library(adegenet)
library(tidyverse)

dat <- ape::read.dna(file = "./stacks_mu_ind/populations.all.phylip")
dat
class(dat)
str(dat)
names(dat)

# calculate number of substitutions
dis_n <- ape::dist.dna(dat, model = "N", pairwise.deletion = TRUE)
dis_n

# scaled by sequence length
dis_raw <- ape::dist.dna(dat, model = "raw", pairwise.deletion = TRUE)
dis_raw

con <- dis_raw %>%
  as.matrix() %>%
  as_tibble()
con

con_n <- dis_n %>%
  as.matrix() %>%
  as_tibble()
con_n

samples <- names(con) %>% str_replace(pattern = "\t", "")


# plot heatmap
p <- con %>%
  mutate(x = samples) %>%
  pivot_longer(cols = names(con),
               names_to = "y",
               values_to = "dissimilarity") %>%
  mutate(y = y %>% str_replace(pattern = "\t", "")) %>% # tidy data
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = dissimilarity)) +
  scale_fill_gradient(low = "cornsilk",
                      high = "brown3") +
  theme_bw()
p

p <- con_n %>%
  mutate(x = samples) %>%
  pivot_longer(cols = names(con_n),
               names_to = "y",
               values_to = "number of substitutions") %>%
  mutate(y = y %>% str_replace(pattern = "\t", "")) %>% # tidy data
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = `number of substitutions`)) +
  scale_fill_gradient(low = "cornsilk",
                      high = "brown3") +
  theme_bw()
p
