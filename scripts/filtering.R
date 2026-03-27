#!/usr/bin/env Rscript

library(tidyverse)
library(readr)
library(dplyr)

# Per site

# quality 

var_qual <- read_delim("Stats/quality_site.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "pink", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 10000) 

summary(var_qual$qual)

# depth

var_depth <- read_delim("Stats/depth_site.ldepth.mean", delim = "\t",
            col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "plum", colour = "black", alpha = 0.3)
b + theme_light() + xlim(0, 15) 

summary(var_depth$mean_depth) # min=5 max=27

# missing data 

var_miss <- read_delim("Stats/missing_site.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "yellow", colour = "black", alpha = 0.3)
c + theme_light() + xlim(0, 0.1)

summary(var_miss$fmiss) # min 5% or 10%

# Maf
var_freq <- read_delim("Stats/max_al.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

maf_al <- ggplot(var_freq, aes(maf)) + geom_density(fill = "darkseagreen4", colour = "black", alpha = 0.3)
maf_al + theme_light()

summary(var_freq$maf) # >=0.03


# Per individual 

# depth
ind_depth <- read_delim("Stats/depth_ind.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

d <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "darkgreen", colour = "black", alpha = 0.3)
d + theme_light()

summary(ind_depth$depth)

# missing data 

ind_miss  <- read_delim("Stats/missing_ind.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

e <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "magenta4", colour = "black", alpha = 0.3)
e + theme_light()

summary(ind_miss$fmiss)

# Quality per individual

# Read without the column names and turn missing values into NA
gq <- read_delim("quality_ind.lqual", delim = "\t", col_names = FALSE, na = ".")
# Make the data into a long form 
gq_long <- gq %>%
  select(-X1, -X2) %>%  # Remove chrom and pos
  # Take all the individuals and stack them into the column "individual" 
  pivot_longer(cols = everything(), names_to = "individual", values_to = "GQ") %>%
  drop_na(GQ)

f <- ggplot(gq_long, aes(x = GQ)) + geom_density(fill = "deeppink4", colour = "black", alpha = 0.3)
f + theme_light() +  xlim(0, 60)

summary(gq_long$GQ) # >=10

