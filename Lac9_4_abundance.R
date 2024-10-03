# Load packages & functions -----------------------------------------------

library(tidyverse)
library(phyloseq)
library(here)
library(decontam)
library(knitr)
library(ggpubr)
library(cowplot)
library(data.table)
library(ape)
library(vegan)
library(RColorBrewer)
library(ggrepel)
library(Maaslin2)

#Load additional functions (created by Dr Wouter de Steenhuijsen Piters - https://gitlab.com/wsteenhu/muis_trx)

source(here("src/load_dada.R"))
source(here("src/utils.R"))

# set paths
knitr::opts_knit$set(root.dir=".", aliases=c(h = "fig.height", w = "fig.width", ow = "out.width"))

#Load phylo
phylo <- readRDS(here("results/ps_final.Rds"))
phylo_RA <- readRDS(here("results/ps_final_RA.Rds"))
phylo_meta <- data.frame(sample_data(phylo))

#Phylum-level phylo

phylo_phylum <- tax_glom(phylo_RA, taxrank = "phylum") %>% prune_taxa(taxa_sums(.)>0, .) 
taxa_names(phylo_phylum) <- tax_table(phylo_phylum)[,"phylum"]

phylo_phylum %>% merge_samples(group = "type", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame() #Phyla across all samples & visits

phylum_visit <- phylo_phylum %>% merge_samples(group = "visit", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame()
rowSums(phylum_visit) #% reads per visit belonging to each phylum (i.e. rowSums=100)

phylum_type <- phylo_phylum %>% merge_samples(group = "sample_type", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame()
rowSums(phylum_type) #% reads per visit belonging to each phylum (i.e. rowSums=100)

#Now look at genera

phylo_genus <- tax_glom(phylo_RA, taxrank = "genus") %>% prune_taxa(taxa_sums(.)>0, .) 
taxa_names(phylo_genus) <- tax_table(phylo_genus)[,"genus"]
ntaxa(phylo_genus) #423 genera

genus_all <- phylo_genus %>% merge_samples(group = "type", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame() %>% t() #Genera across all samples & visits

genus_visit <- phylo_genus %>% merge_samples(group = "visit", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame()
rowSums(genus_visit) #% reads per visit belonging to each phylum (i.e. rowSums=100)

genus_type <- phylo_genus %>% merge_samples(group = "sample_type", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame()
rowSums(genus_type) #% reads per visit belonging to each phylum (i.e. rowSums=100)

genus_type_transposed <- t(genus_type)

#And now by ASVs

ASV_all <- phylo %>% merge_samples(group = "type", fun=mean) %>% 
  transform_sample_counts(., function(x) 100*x/sum(x)) %>% 
  otu_table %>% as.data.frame() %>% t() #ASVs across all samples & visits
nrow(ASV_all)

# Stacked bar charts (without using create_bar function) --------------------

#First make long dataframe: extract OTU table, sort by decreasing RA, then pivot 
#to give one row per ASV per sample

stack.phylo <- t(otu_table(phylo_RA)[order(rowSums(otu_table(phylo_RA)), decreasing=T),][1:35, ]) %>% 
  data.frame(., Residuals = 1-rowSums(.), check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "ASV", values_to="RA") %>% 
  mutate(ASV = fct_rev(fct_inorder(ASV)))

#Then make df from phylo_RA, adding a column called "sample_id"
phylo_RA_df <- data.frame(sample_data(phylo_RA))
phylo_RA_df$sample_id <- rownames(phylo_RA_df)

stack.phylo <- left_join(phylo_RA_df, 
                         stack.phylo, by = c("sample_id"))

stack_colours <- c("grey60","forestgreen","pink","deepskyblue","green3","darkred","forestgreen","deepskyblue",
                   "lightgoldenrod1","orange","navyblue","orange","lightgoldenrod1","yellow","lightgoldenrod1","red","gold",
                   "lightgoldenrod1","forestgreen","lightgoldenrod1","navyblue","red","gold","gold","indianred1","navyblue",
                   "darkred","lightblue1","deepskyblue","turquoise2","gold","navyblue","orange","dodgerblue3","orange","navyblue")

bar_phylo <- stack.phylo %>%
  ggplot(aes(x = participant, y = RA, fill = ASV)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = stack_colours) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition (ASV)") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

bar_phylo_NP <- stack.phylo %>%
  filter(sample_type %in% c("MN", "IN")) %>% 
  ggplot(aes(x = participant, y = RA, fill = ASV)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = stack_colours) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

png(file=here("plots", "stacked_taxonomy_NP.png"), width = 18, height = 10, units = "in",res=300)
bar_phylo_NP
dev.off()

bar_phylo_ISBM <- stack.phylo %>%
  filter(sample_type %in% c("IS","BM")) %>% 
  ggplot(aes(x = participant, y = RA, fill = ASV)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = stack_colours) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

bar_phylo_MSCS <- stack.phylo %>%
  filter(sample_type %in% c("MS","CS")) %>% 
  ggplot(aes(x = participant, y = RA, fill = ASV)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = stack_colours) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

bar_phylo_mother_inoc_MO <- stack.phylo %>%
  filter(sample_type %in% c("MO")) %>% 
  ggplot(aes(x = participant, y = RA, fill = ASV)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = stack_colours) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("MO taxonomic composition (inoculated TRUE/FALSE)") +
  facet_grid(visit ~ inoculated, scales = "free_x", space = "free_x")

#Now do stacked bar charts by genus

phylo_genus_df <- data.frame(sample_data(phylo_genus))
phylo_genus_df$sample_id <- rownames(phylo_genus_df)

stack.phylo.genus <- t(otu_table(phylo_genus)[order(rowSums(otu_table(phylo_genus)), decreasing=T),][1:15, ]) %>% 
  data.frame(., Residuals = 1-rowSums(.), check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "genus", values_to="RA") %>% 
  mutate(genus = fct_rev(fct_inorder(genus)))
stack.phylo.genus <- left_join(phylo_genus_df, 
                               stack.phylo.genus, by = "sample_id") 
bar_genus <- stack.phylo.genus %>%
  ggplot(aes(x = participant, y = RA, fill = genus)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = 
                      c("grey60","pink","yellow","green3","lightblue1","indianred1","turquoise2"
                        ,"red","darkred","deepskyblue","forestgreen","lightgoldenrod1","dodgerblue3",
                        "gold","orange","navyblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo.genus) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition (genus)") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

#And now phylum!

phylo_phylum_df <- data.frame(sample_data(phylo_phylum))
phylo_phylum_df$sample_id <- rownames(phylo_phylum_df)

stack.phylo.phylum <- t(otu_table(phylo_phylum)[order(rowSums(otu_table(phylo_phylum)), decreasing=T),][1:5, ]) %>% 
  data.frame(., Residuals = 1-rowSums(.), check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "phylum", values_to="RA") %>% 
  mutate(phylum = fct_rev(fct_inorder(phylum)))
stack.phylo.phylum <- left_join(phylo_phylum_df, 
                                stack.phylo.phylum, by = "sample_id") 

bar_phylum <- stack.phylo.phylum %>%
  ggplot(aes(x = participant, y = RA, fill = phylum)) + 
  geom_bar(stat = "identity", position = "stack", colour = "white") + 
  scale_fill_manual(values = c("grey60","pink","forestgreen","red","gold","navyblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  scale_x_discrete(name="Participant", label=stack.phylo.phylum) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Relative abundance") + ggtitle("Taxonomic composition (phylum)") +
  facet_grid(visit ~ sample_type, scales = "free_x", space = "free_x")

# Maaslin2 - one per niche, multivariable ---------------------------------

set.seed(123)
maaslin_MO_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="MO" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="MO" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "MO_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

set.seed(123)
maaslin_MS_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="MS" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="MS" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "MS_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

set.seed(123)
maaslin_MN_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="MN" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="MN" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "MN_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

set.seed(123)
maaslin_BM_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="BM" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="BM" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "BM_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

set.seed(123)
maaslin_IS_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="IS" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="IS" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "IS_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

set.seed(123)
maaslin_IN_multivar <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type=="IN" & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type=="IN" & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "IN_multivar"),
  random_effects = c("participant"),
  fixed_effects = c("visit","inoculated","sib_5yr","m_nlac_col","abx_mum_todate"))

phylo_RA %>%
  prune_taxa("Corynebacterium_33", .) %>%
  subset_samples(sample_type == "MN") %>%
  ps_to_df() %>%
  pivot_longer(names_to = "ASV", values_to = "RA", "Corynebacterium_33") %>%
  select(sample_id2, inoculated, ASV, RA) %>%
  ggplot(aes(x = inoculated, y = RA, color = inoculated)) +
  geom_jitter(width = 0.5, size = 2.5, alpha = 0.6) +
  geom_boxplot(alpha = 0.5, fill = "lightblue", outlier.colour = NA) +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.title.y = ggtext::element_markdown()) +
  scale_color_manual(values = c("TRUE" = "deepskyblue", "FALSE" = "lightpink")) +
  labs(y = "Relative abundance", x = "Inoculated") +
  theme(legend.position = "bottom")

# Maaslin2 - comparing niches ---------------------------------------------

set.seed(123)
maaslin_INvIS <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type %in% c("IN","IS") & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type %in% c("IN","IS") & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "INvIS"),
  random_effects = c("id_v"),
  reference=c("sample_type","IS"),
  fixed_effects = c("sample_type"))

maaslin_INvIS_bar <- as.data.frame(maaslin_INvIS$results) %>%
  arrange(coef) %>%
  mutate(
    group = if_else(coef < 0, "IS", "IN"),
    feature_new = format_OTU(feature),
    feature_new = factor(feature_new, levels = feature_new)
  ) %>%
  ggplot(aes(x = coef, y = feature_new, fill = group)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(y = feature_new, xmin = coef - stderr, xmax = coef + stderr),
                width = 0.05, colour = "black", alpha = 0.9, size = 0.3) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 12)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-6.5, 6.5)) +
  geom_vline(xintercept = 0, col = "black") +
  xlab("Standardised β with SE") + ylab("") +
  ggtitle("Maaslin - IN vs IS")

set.seed(123)
maaslin_INvMN <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type %in% c("IN","MN") & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type %in% c("IN","MN") & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "INvMN"),
  random_effects = c("id_v"),
  reference=c("sample_type","MN"),
  fixed_effects = c("sample_type"))

maaslin_INvMN_bar <- as.data.frame(maaslin_INvMN$results) %>%
  arrange(coef) %>%
  mutate(
    group = if_else(coef < 0, "MN", "IN"),
    feature_new = format_OTU(feature),
    feature_new = factor(feature_new, levels = feature_new)
  ) %>%
  ggplot(aes(x = coef, y = feature_new, fill = group)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(y = feature_new, xmin = coef - stderr, xmax = coef + stderr),
                width = 0.05, colour = "black", alpha = 0.9, size = 0.3) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 12)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-6.5, 6.5)) +
  geom_vline(xintercept = 0, col = "black") +
  xlab("Standardised β with SE") + ylab("") +
  ggtitle("Maaslin - IN vs MN")

set.seed(123)
maaslin_ISvMS <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type %in% c("MS","IS") & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type %in% c("MS","IS") & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "ISvMS"),
  random_effects = c("id_v"),
  reference=c("sample_type","MS"),
  fixed_effects = c("sample_type"))

maaslin_ISvMS_bar <- as.data.frame(maaslin_ISvMS$results) %>%
  arrange(coef) %>%
  mutate(
    group = if_else(coef < 0, "MS", "IS"),
    feature_new = format_OTU(feature),
    feature_new = factor(feature_new, levels = feature_new)
  ) %>%
  ggplot(aes(x = coef, y = feature_new, fill = group)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(y = feature_new, xmin = coef - stderr, xmax = coef + stderr),
                width = 0.05, colour = "black", alpha = 0.9, size = 0.3) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 12)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-6.5, 6.5)) +
  geom_vline(xintercept = 0, col = "black") +
  xlab("Standardised β with SE") + ylab("") +
  ggtitle("Maaslin - IS vs MS")

set.seed(123)
maaslin_ISvBM <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type %in% c("BM","IS") & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type %in% c("BM","IS") & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "ISvBM"),
  random_effects = c("id_v"),
  reference=c("sample_type","BM"),
  fixed_effects = c("sample_type"))

maaslin_ISvBM_bar <- as.data.frame(maaslin_ISvBM$results) %>%
  arrange(coef) %>%
  mutate(
    group = if_else(coef < 0, "BM", "IS"),
    feature_new = format_OTU(feature),
    feature_new = factor(feature_new, levels = feature_new)
  ) %>%
  ggplot(aes(x = coef, y = feature_new, fill = group)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(y = feature_new, xmin = coef - stderr, xmax = coef + stderr),
                width = 0.05, colour = "black", alpha = 0.9, size = 0.3) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 12)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-6.5, 6.5)) +
  geom_vline(xintercept = 0, col = "black") +
  xlab("Standardised β with SE") + ylab("") +
  ggtitle("Maaslin - IS vs BM")

set.seed(123)
maaslin_MOvMS <- Maaslin2(
  input_data = t(data.frame(otu_table(subset_samples(phylo_RA, sample_type %in% c("MO","MS") & visit!="V1")))),
  input_metadata = data.frame(sample_data(subset_samples(phylo_RA, sample_type %in% c("MO","MS") & visit!="V1"))),
  normalization = "NONE",  # already RA
  min_abundance = 0.05,
  min_prevalence = 0.05,
  output = here("maaslin_results", "MOvMS"),
  random_effects = c("id_v"),
  reference=c("sample_type","MO"),
  fixed_effects = c("sample_type"))

maaslin_MOvMS_bar <- as.data.frame(maaslin_MOvMS$results) %>%
  arrange(coef) %>%
  mutate(
    group = if_else(coef < 0, "MO", "MS"),
    feature_new = format_OTU(feature),
    feature_new = factor(feature_new, levels = feature_new)
  ) %>%
  ggplot(aes(x = coef, y = feature_new, fill = group)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(y = feature_new, xmin = coef - stderr, xmax = coef + stderr),
                width = 0.05, colour = "black", alpha = 0.9, size = 0.3) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = ggtext::element_markdown(size = 12)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-6.5, 6.5)) +
  geom_vline(xintercept = 0, col = "black") +
  xlab("Standardised β with SE") + ylab("") +
  ggtitle("Maaslin - MO vs MS")

maaslin_MOvMS$results %>%
  filter(qval < 0.1)

