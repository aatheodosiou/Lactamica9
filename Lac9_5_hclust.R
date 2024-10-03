# Load packages, functions and data ---------------------------------------

#First load relevant packages (installed using install.packages)

library(tidyverse)
library(phyloseq)
library(here)
library(decontam)
library(knitr)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(ggrepel)
library(fpc)
library(ggdendro)
library(dendextend)
library(rstatix)
library(tableone)
library(MVN)

#Load additional functions (created by Dr Wouter de Steenhuijsen Piters - https://gitlab.com/wsteenhu/muis_trx)

source(here("src/load_dada.R"))
source(here("src/utils.R"))

# set paths
knitr::opts_knit$set(root.dir=".", aliases=c(h = "fig.height", w = "fig.width", ow = "out.width"))

#Define phylo & BC

phylo <- readRDS(here("results/ps_final.Rds"))
phylo_RA <- readRDS(here("results/ps_final_RA.Rds"))
phylo_meta <- data.frame(sample_data(phylo))

ASV_RA <- as.data.frame(t(as(otu_table(phylo_RA), "matrix")))
braycurtis <- vegdist(ASV_RA, method="bray")
braycurtis.matrix <- as.matrix(braycurtis)

# Functions for visualising indices

#to visualise all four indices
clust_all_function <- function(dist, k, method = "complete") {
  stats_k <- c(k=k, 
               cluster.stats(dist, cutree(hclust(dist, method=method), k))[
                 c("ch", "avg.silwidth","within.cluster.ss","dunn")])
  return(stats_k)
}

#to only visualise CH and average silhouette width -> present this in paper/thesis figures
clust_chsw_function <- function(dist, k, method = "complete") {
  stats_k <- c(k=k, 
               cluster.stats(dist, cutree(hclust(dist, method=method), k))[
                 c("ch", "avg.silwidth")])
  return(stats_k)
}

#Test for clusters

#CH index is ratio of between-cluster variance to within-cluster variance. 
#SW measures cohesion (similarity to own cluster) vs separation (compared to other clusters)
#i.e. peaks indicate best compromise in terms of defining clusters
#Complete linkage - distance b/w 2 clusters = max distance b/w any 2 points in cluster; 
#VS average linkage - average distance between all pairs of points
#Complete more sensitive to noise/outliers & chaining effect (clusters connected linearly)

n <- 30 #how many clusters to test

#complete linkage
clust_complete_all <- lapply(2:n, clust_chsw_function, dist=braycurtis, method="complete") %>% 
  bind_rows %>%
  setNames(c("k", "Calinski-Harabasz", "Silhouette width"))

plot_clust_complete_all <- clust_complete_all %>%
  gather(index, value, -k) %>%
  ggplot(aes(x=k, y=value)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(index), scales="free_y")

#average linkage
clust_average_all <- lapply(2:n, clust_chsw_function, dist=braycurtis, method="average") %>% 
  bind_rows %>%
  setNames(c("k", "Calinski-Harabasz", "Silhouette width"))

plot_clust_average_all <- clust_average_all %>%
  gather(index, value, -k) %>%
  ggplot(aes(x=k, y=value)) +
  geom_point() +
  geom_line() +
  ylab("Value") + xlab("Number of clusters") +
  facet_wrap(vars(index), scales="free_y") +
  theme_bw(base_size=13) #19 clusters?

#Heatmap-cluster plot (min_size = min number samples in each cluster)

set.seed(100)
hm_plot_INMN <- 
  cluster_plot_hm_sampletype2_RA(subset_samples(phylo_RA, sample_type %in% c("MN","IN")),k=8,n_otus=20,min_size=3, method = "average")
hm_plot_ISMS <- 
  cluster_plot_hm_sampletype2_RA(subset_samples(phylo_RA, sample_type %in% c("MS","IS")),k=11,n_otus=25,min_size=3, method = "average")
hm_plot_MOMSBMIS <- 
  cluster_plot_hm_sampletype2_RA(subset_samples(phylo_RA, sample_type %in% c("MO","MS","BM","IS")),k=13,n_otus=25,min_size=3, method = "average")
hm_plot_allniches <- 
  cluster_plot_hm_sampletype2_RA(phylo_RA,k=19,n_otus=30,min_size=3, method = "average")

