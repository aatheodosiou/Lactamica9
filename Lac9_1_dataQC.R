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

#Load additional functions (created by Dr Wouter de Steenhuijsen Piters - https://gitlab.com/wsteenhu/muis_trx)

source(here("src/load_dada.R"))
source(here("src/utils.R"))

# set paths
knitr::opts_knit$set(root.dir=".", aliases=c(h = "fig.height", w = "fig.width", ow = "out.width"))

#Load DADA2 output ASV table (seqtab.rds)

seqtab <- readRDS(here("output/seqtab.rds"))

#Load taxonomy table

read_tax <- function(path) {
  read.table(path, stringsAsFactors = F) %>% 
    as.matrix()
}

taxa <- read_tax(here("taxonomy/assigntaxonomy.tsv"))

#Load participant metadata and inspect matrix

meta <- read_tsv(here("metadata/meta.tsv")) %>% 
  column_to_rownames("sample_id")

# Make and inspect phyloseq object ----------------------------------------

#Create a raw read phyloseq object

ps <- list()
ps$raw <- create_phylo(seqtab, taxa, meta)

#Create phylo function removes chlor, mitoch, archaea & euk; and also samples with NO ASVs

#Add raw read count to data

sample_data(ps$raw)$raw_reads <- sample_sums(ps$raw)
readcount_vector <- c(sample_data(ps$raw)$raw_reads) %>% sort() #Vector of read counts

#Convert raw reads to relative abundance (RA)

ps$RA <- ps$raw %>% to_RA 

# Plot mocks and blanks ---------------------------------------------------

box_reads_per_sample <- data.frame(sample_data(ps$raw)) %>%
  subset(.,ext_run != "NA") %>%
  ggplot(aes(x = type, y = raw_reads)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha = 0.3, width = 0.2) +
  labs(x="Sample type", y="Raw read count", title="Reads per sample") +
  theme(legend.position = "none")

#Plot mocks - RA for top 20 taxa for isolation mocks

bar_iso_mock_RA <- subset_samples(ps$RA, mock_type == "iso_mock") %>% 
  create_bar(., n=12, ncol_legend=1) + 
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) + 
  labs(title = "Taxonomic composition of isolation mocks")

bar_iso_mocks_raw_count <- subset_samples(ps$RA, mock_type == "iso_mock") %>%
  meta_to_df() %>% rownames_to_column(var="sample_id") %>% 
  ggplot(aes(x = sample_id, y = raw_reads)) +
  geom_col() +
  geom_text(aes(label = scales::comma(raw_reads, accuracy = 1)),  angle = 90, 
            color = "white", hjust = 2, size = 2.5, nudge_y = -0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Raw read count")

cow_iso_mocks <- cowplot::plot_grid(bar_iso_mock_RA,
                                    bar_iso_mocks_raw_count,
                                    align = "v", axis = "lr", nrow = 2, 
                                    rel_heights = c(0.8, 0.2))

#Plot iso blanks

bar_iso_blanks_RA <- subset_samples(ps$RA, blank_type == "iso_blank") %>%
  create_bar(., n = 20, ncol_legend = 1) +
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) +
  labs(title = "Taxonomic composition of isolation blanks")

bar_iso_blanks_raw_count <- subset_samples(ps$RA, blank_type == "iso_blank") %>%
  meta_to_df() %>% rownames_to_column(var="sample_id") %>% 
  ggplot(aes(x = sample_id, y = raw_reads)) +
  geom_col() +
  geom_text(aes(label = scales::comma(raw_reads, accuracy = 1)),  angle = 90, 
            color = "white", hjust = 2, size = 2.5, nudge_y = -0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "", y = "Raw read count")

cow_iso_blanks <- cowplot::plot_grid(bar_iso_blanks_RA, 
                                     bar_iso_blanks_raw_count,
                                     align = "v", axis = "lr", nrow = 2, 
                                     rel_heights = c(0.8, 0.2))

plotgrid_mock.blank.taxa <- plot_grid(cow_iso_blanks, cow_iso_mocks, 
                                      ncol = 2, nrow = 1)

# Define colours ----------------------------------------------------------

sample_types <- c("MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")
sample_type_colors <- c("firebrick2", "darkorange2", "gold2", "springgreen3","dodgerblue2", "purple2","violet","gray40")
names(sample_type_colors) <- sample_types
sample_type_colors <- sample_type_colors[sample_types]

visits <- c("V1","V3","V4","V5","V6","blank")
visit_colors <- c("deeppink3","yellow2","palegreen3","skyblue2","mediumpurple","gray40")
names(visit_colors) <- visits
visit_colors <- visit_colors[visits]

participant_types <- c("mother","infant","sibling")
participant_type_colors <- c("lightpink1","darkseagreen2","goldenrod1")
names(participant_type_colors) <- participant_types
participant_type_colors <- participant_type_colors[participant_types]

sampletype_visits <- c("V3.MO","V4.MO","V5.MO","V6.MO", "V3.MS","V4.MS","V5.MS","V6.MS",
                       "V3.MN","V4.MN","V5.MN","V6.MN","V3.BM","V4.BM","V5.BM","V6.BM",
                       "V3.IS","V4.IS","V5.IS","V6.IS","V3.IN","V4.IN","V5.IN","V6.IN", "V6.CS")
sampletype_visit_colors <- c("firebrick1","firebrick2","firebrick3","firebrick4", 
                             "orange1","orange2","darkorange2","darkorange3","yellow1","yellow2","yellow3","yellow4",
                             "palegreen","green","green3","green4","deepskyblue","deepskyblue3","dodgerblue3","blue3",
                             "orchid1","mediumorchid2","magenta3","purple3","lightpink1")
names(sampletype_visit_colors) <- sampletype_visits
sampletype_visit_colors <- sampletype_visit_colors[sampletype_visits]

sampletype_visit_colors_INV3 <- c("firebrick2","firebrick2","firebrick2","firebrick2", 
                                  "darkorange2","darkorange2","darkorange2","darkorange2",
                                  "gold2","gold2","gold2","gold2",
                                  "springgreen3","springgreen3","springgreen3","springgreen3",
                                  "dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2",
                                  "magenta","purple2","purple2","purple2","lightpink1")
names(sampletype_visit_colors_INV3) <- sampletype_visits
sampletype_visit_colors_INV3 <- sampletype_visit_colors_INV3[sampletype_visits]

sampletype_TF <- c("TRUE.MO","FALSE.MO","TRUE.MS","FALSE.MS","TRUE.MN","FALSE.MN",
                   "TRUE.BM","FALSE.BM","TRUE.IS","FALSE.IS","TRUE.IN","FALSE.IN")
sampletype_TF_colors <- c("firebrick1","firebrick3","orange1","darkorange3",
                          "yellow1","yellow3","palegreen","green4","deepskyblue","blue3",
                          "orchid1","purple3")
names(sampletype_TF_colors) <- sampletype_TF
sampletype_TF_colors <- sampletype_TF_colors[sampletype_TF]

# Inspect by biomass & libary size, & prepare for decontam ----------------------------------------------------

#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#Filter out mocks & PCR blanks to create input for decontam

ps_raw_decontam_full <- ps$raw %>%
  subset_samples(type == "sample" | blank_type == "iso_blank") %>%
  prune_taxa(!taxa_sums(.) == 0, .) #Gets rid of any taxa with 0 readcount

ps_raw_samples <- ps$raw %>%
  subset_samples(type == "sample") %>%
  prune_taxa(!taxa_sums(.) == 0, .) #Another with just samples, for plots

#Note about boxplots: box gives IQR (bounded by 25th & 75th C); midpoint is median;
#Whisker extends from min to max; outlier defined as >1.5IQR either side of median

#Biomass vs ext run, visit, & sample type

box_conc_extrun <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  ggplot(aes(x = ext_run, y = qpcr_pgul)) +
  geom_jitter(width=0.2, fill="white", size=2.5, alpha=0.8) +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank()) +
  labs(y = "DNA concentration (pg/ul)", x = "DNA extraction run", title = "Extraction run") +
  geom_hline(yintercept = 0.3, linetype = "dashed") # Cut-off 0.3 chosen to highlight

box_conc_visit_all <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  mutate(visit = fct_explicit_na(visit, "blank") %>% 
           fct_relevel(., "V1", "V3", "V4", "V5", "V6", "blank")) %>%
  ggplot(aes(x = visit, y = qpcr_pgul, color=visit)) +
  geom_jitter(width=0.2, fill="white", size=2.5, alpha=0.6) +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank()) +
  labs(y = "DNA concentration (pg/ul)", x="Visit", title = "Visit") +
  geom_hline(yintercept = 0.3, linetype="dashed")

box_conc_extrun_visit <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  mutate(visit = fct_explicit_na(visit, "blank") %>% 
           fct_relevel(., "V1","V3","V4","V5","V6", "blank")) %>% 
  ggplot(aes(x = ext_run, y = qpcr_pgul, color = visit)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "DNA concentration (pg/ul)", x = "DNA extraction run", title = "Extraction run and visit") +
  geom_hline(yintercept = 0.3, linetype = "dashed") # Cut-off 0.3 chosen to highlight

box_conc_sampletype <- ps_raw_samples %>%
  meta_to_df() %>%
  mutate(sample_type = factor(sample_type, levels = c("MO", "MS", "MN", "BM", "IS", "IN", "CS"))) %>%
  ggplot(aes(x = sample_type, y = qpcr_pgul, color=sample_type)) +
  geom_jitter(aes(fill = sample_type), width = 0.2, size = 2.5, alpha = 0.8) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1), 
                limits = c(0.05, 100000)) +  
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "DNA concentration (pg/ul)", x="Sample type")

box_conc_extrun_sampletype <- ps_raw_samples %>%
  meta_to_df() %>%
  mutate(sample_type = factor(sample_type, levels = c("MO", "MS", "MN", "BM", "IS", "IN", "CS"))) %>%
  ggplot(aes(x = ext_run, y = qpcr_pgul, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, width = 0.2) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "DNA concentration (pg/ul)", x = "DNA extraction run")

box_conc_visit_infant <- ps_raw_samples %>%
  subset_samples(., sample_type %in% c("IN", "IS")) %>% 
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "IS", "IN")) %>% 
  ggplot(aes(x = visit, y = qpcr_pgul, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, position = position_jitterdodge(dodge.width = 0.8)) +
  stat_summary(fun = median, geom = "line", aes(group = sample_type),
               size = 1, alpha = 0.6, position = position_dodge(0.8)) + # Calculate median within ggplot
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1), 
                limits = c(0.05, 10000)) +  
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
  labs(y = "DNA concentration (pg/ul)", x = "Visit")

box_conc_visit_mother <- ps_raw_samples %>%
  subset_samples(., sample_type %in% c("BM", "MO", "MS", "MN")) %>% 
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM")) %>% 
  ggplot(aes(x = visit, y = qpcr_pgul, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, position = position_jitterdodge(dodge.width = 0.8)) +
  stat_summary(fun = median, geom = "line", aes(group = sample_type),
               size = 1, alpha = 0.6, position = position_dodge(0.8)) + # Calculate median within ggplot
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_log10(labels = function(x) scales::comma(x, accuracy = 0.1), 
                limits = c(0.05, 10000)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
  labs(y = "DNA concentration (pg/ul)", x = "Visit") 

plotgrid_DNAconc <- plot_grid(box_conc_extrun_sampletype,box_conc_sampletype,
                              box_conc_visit_infant, box_conc_visit_mother,
                              ncol = 2, nrow = 2, labels="auto")

#Plot based on library size (reads per sample) in decontam input (excl all mocks & PCR blanks)

plot_reads_per_sample <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  arrange(raw_reads) %>%
  mutate(index = row_number()) %>%
  ggplot(aes(x = index, y = raw_reads, color = type)) +
  geom_point(alpha = 0.5, size = 2, shape = 16, stroke=1) +  
  scale_color_manual(values = c("blank" = "grey40", "sample" = "turquoise"), 
                     guide = guide_legend(title = "Type")) +
  geom_hline(yintercept = 5000, linetype="dashed") +
  labs(y = "Raw reads per sample", x = "Samples", 
       title = "Reads per sample (ascending)")

box_reads_Illumina <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>% 
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x = illumina_run, y = raw_reads, color=sample_type)) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  geom_boxplot(fill = "grey", color = "grey40", alpha = 0.2, outlier.colour = NA) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x="Illumina run", title="Illumina run")

box_reads_extrun <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  ggplot(aes(x = ext_run, y = raw_reads)) +
  geom_jitter(width=0.2, fill="white", size=2.5, alpha=0.6) +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank()) +
  labs(y = "Raw reads per sample", x = "DNA extraction run", title = "Extraction run")

box_reads_extrun_sampletype <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>% 
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>% 
  ggplot(aes(x = ext_run, y = raw_reads, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x = "DNA extraction run", title = "Extraction run")

box_reads_visit_all <- ps_raw_decontam_full %>%
  meta_to_df() %>%
  mutate(visit = fct_explicit_na(visit, "blank") %>% 
           fct_relevel(., "V1", "V3", "V4", "V5", "V6", "blank")) %>%
  ggplot(aes(x = visit, y = raw_reads)) +
  geom_jitter(width=0.2, fill="white", size=2.5, alpha=0.6) +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank()) +
  labs(y = "Raw reads per sample", x="Visit", title = "Visit")

plot_reads.vs.conc <- subset_samples(ps_raw_decontam_full, type=="sample") %>% 
  meta_to_df() %>% 
  mutate(sample_type = factor(sample_type, levels = c("MO", "MS", "MN", "BM", "IS", "IN", "CS"))) %>% 
  ggplot(aes(x = qpcr_pgul, y = raw_reads, color = sample_type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x = "DNA concentration (pg/ul)")

box_reads_sampletype <- subset_samples(ps_raw_decontam_full, type=="sample") %>%
  meta_to_df() %>%
  mutate(sample_type = factor(sample_type, levels = c("MO", "MS", "MN", "BM", "IS", "IN", "CS"))) %>%
  ggplot(aes(x = sample_type, y = raw_reads, color=sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, width = 0.2) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_continuous(limits=c(0,65000)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title=element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x="Sample type")

box_reads_visit_infant <- ps_raw_decontam_full %>%
  subset_samples(., sample_type %in% c("IN", "IS")) %>% 
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "IS", "IN")) %>% 
  ggplot(aes(x = visit, y = raw_reads, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, position = position_jitterdodge(dodge.width = 0.8)) +
  stat_summary(fun = median, geom = "line", aes(group = sample_type),
               size = 1, alpha = 0.6, position = position_dodge(0.5)) + # Calculate median within ggplot
  scale_color_manual(values = sample_type_colors) + # Color for data points
  scale_y_continuous(limits=c(0,65000)) +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x = "Visit")

box_reads_visit_mother <- ps_raw_decontam_full %>%
  subset_samples(., sample_type %in% c("BM", "MO", "MS", "MN")) %>% 
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM")) %>% 
  ggplot(aes(x = visit, y = raw_reads, color = sample_type)) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  geom_jitter(alpha = 0.8, position = position_jitterdodge(dodge.width = 0.8)) +
  stat_summary(fun = median, geom = "line", aes(group = sample_type),
               size = 1, alpha = 0.6, position = position_dodge(0.5)) + # Calculate median within ggplot
  scale_color_manual(values = sample_type_colors) + # Color for data points
  theme(axis.title.x = element_text(), axis.text.x = element_text(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = "none") +
  labs(y = "Raw reads per sample", x = "Visit")

plotgrid_rawreads <- plot_grid(plot_reads.vs.conc, box_reads_sampletype,
                               box_reads_visit_infant, box_reads_visit_mother,
                               ncol=2, nrow=2, labels="auto")

# Run decontam - combined model -------------------------------------------

sl(contam, 
   isContaminant(ps_raw_decontam_full, conc = "qpcr_pgul", neg = "is_neg", 
                 method = "combined", threshold = 0.10), 
   dir_path = here(), overwrite = TRUE)

#Create sub-table of contaminant = TRUE

contam_TRUE <- contam %>% filter(contaminant == TRUE)

#Histogram of decontam score against ASV count, to see if threshold appropriate.
#Should be bimodal, with clear separation between contam and non-contam

prevalencePalette <- c("1-2" = "#edf8e9", "3-5" = "#bae4b3", "6-10" = "#74c476", 
                       "11-20" = "#41ab5d", "21-50" = "#238b45", "51-100" = "#006d2c", 
                       ">100" = "#00441b")
contam$prev_range <- cut(contam$prev, 
                         breaks = c(1, 2, 5, 10, 20, 50, 100, 600), 
                         labels = c("1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))
table(contam$prev_range)

hist_comb <- ggplot(contam, aes(x=p, fill=contaminant)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='Decontam score', y='ASV count') +
  scale_fill_manual(values = prevalencePalette)

ggplot(contam, aes(x=p, fill=contaminant)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='decontam score', y='number of ASVs')

hist_prev <- ggplot(contam, aes(x=p.prev, fill=prev_range)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='decontam score', y='number of ASVs') 

hist_freq <- ggplot(contam, aes(x=p.freq, fill=prev_range)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='decontam score', y='number of ASVs')

#Play around with different thresholds:
sl(contam0.01, 
   isContaminant(ps_raw_decontam_full, conc = "qpcr_pgul", neg = "is_neg", 
                 method = "combined", threshold = 0.01), 
   dir_path = here(), overwrite = TRUE)

contam0.01 %>% filter(contaminant == TRUE) %>% count()
ggplot(contam0.01, aes(x=p, fill=contaminant)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='decontam score', y='number of ASVs')

sl(contam0.5, 
   isContaminant(ps_raw_decontam_full, conc = "qpcr_pgul", neg = "is_neg", 
                 method = "combined", threshold = 0.5), 
   dir_path = here(), overwrite = TRUE)

contam0.5 %>% filter(contaminant == TRUE) %>% count()
ggplot(contam0.5, aes(x=p, fill=contaminant)) +
  geom_histogram(binwidth = 0.005, na.rm = TRUE) + 
  labs(x='decontam score', y='number of ASVs')

#Inspect contam

table(contam$contaminant)

which(contam$contaminant) #Lists contaminant ASVs by abundnace ranking

taxa_contam.byprev_TRUE <- arrange(contam_TRUE, desc(prev)) %>% 
  rownames() #Arrange contamdf.comb output by desc prev & save taxa as vector

taxa_contam.byprev_FALSE <- contam %>% 
  filter(contaminant == FALSE) %>% arrange(., desc(prev)) %>% rownames()

taxa_contam.byfreq_TRUE <- arrange(contam_TRUE, desc(freq)) %>% 
  rownames() #Arrange contamdf.comb output by desc freq & save taxa as vector

taxa_contam.byfreq_FALSE <- contam %>% 
  filter(contaminant == FALSE) %>% arrange(., desc(freq)) %>% rownames()

#Plot contaminants: freq vs DNA conc (combined model)

plot_contam_freq.vs.conc.all <- plot_frequency(ps_raw_decontam_full, 
                                               taxa_contam.byfreq_TRUE, conc="qpcr_pgul") + 
  xlab("DNA concentration (pg/ul)") #Use this to plot ALL contam on one plot

#Use this to plot groups of contam (by desc freq) on one plot

plot_freq_function <- function(taxa) {
  plot_frequency(ps_raw_decontam_full, taxa, conc="qpcr_pgul") + 
    xlab("DNA concentration (pg/ul)")} 

#Use this to plot groups of contam (by desc freq) on one plot (1:20, 21:40, etc)

plot_freq_function(taxa_contam.byfreq_TRUE[1:20]) #Etc...

#Now plot contam prev in samples vs blanks (comb model)

ps_decontam_prev <- transform_sample_counts(ps_raw_decontam_full, 
                                            function(abund) 1*(abund>0))

contam.prev_sample.vs.blank <- 
  data.frame(prev.sample = taxa_sums(prune_samples(sample_data(ps_decontam_prev)$type == "sample", ps_decontam_prev)), 
             prev.blank = taxa_sums(prune_samples(sample_data(ps_decontam_prev)$type == "blank", ps_decontam_prev)), 
             contaminant=contam$contaminant)

contam.prev_sample.vs.blank$prop.sample <- 
  contam.prev_sample.vs.blank$prev.sample / nrow(subset(sample_data(ps_decontam_prev), type=="sample"))
contam.prev_sample.vs.blank$prop.blank <- 
  contam.prev_sample.vs.blank$prev.blank / nrow(subset(sample_data(ps_decontam_prev), type=="blank"))

plot_contam.prev_sample.vs.blank <- ggplot(data=contam.prev_sample.vs.blank, 
                                           aes(x=prop.blank, y=prop.sample, colour=contaminant)) +
  geom_point() + xlab("Prevalence - Blanks") + ylab("Prevalence - Samples") +
  ggtitle("Prevalence of contaminants in samples vs blanks")

# Review results of comb model ---------------------------------

#Compare RA in blanks vs sample for vector of ?contam taxa

plot_RA_sample.vs.blank_function <- function(taxa_vector, phyloseq_object) {
  phyloseq_object %>%
    to_RA() %>%
    prune_taxa(taxa_vector, .) %>%
    ps_to_df() %>%
    pivot_longer(names_to = "ASV", values_to = "RA", head(taxa_vector, n=1):tail(taxa_vector, n=1)) %>%
    ggplot(aes(x = type, y = RA)) +
    geom_jitter(width = 0.2, size=2.5, alpha=0.6) +
    geom_boxplot(alpha=0.5, outlier.colour = NA) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.y = ggtext::element_markdown()) +
    facet_wrap(~ASV) + labs(y = "Relative abundance", x="")
}

plot_RA_sample.vs.blank_function(taxa_contam.byfreq_TRUE[1:20], ps_raw_decontam_full)

#Can also use this function to look at RA in blank vs sample for just one taxon:
plot_RA_sample.vs.blank_function("Streptococcus_1", ps_raw_decontam_full)

#Now find list of ASVs where RA in samples is higher than RA in blanks:

list_RA_sample.vs.blank_function <- function(taxa_vector, phyloseq_object) {
  df <- phyloseq_object %>%
    to_RA() %>%
    prune_taxa(taxa_vector, .) %>%
    ps_to_df() %>%
    pivot_longer(names_to = "ASV", values_to = "RA", head(taxa_vector, n=1):tail(taxa_vector, n=1))
  
  mean_RA_samples <- df %>%
    filter(type == "sample") %>%
    group_by(ASV) %>%
    summarise(mean_RA = mean(RA))
  
  mean_RA_blanks <- df %>%
    filter(type == "blank") %>%
    group_by(ASV) %>%
    summarise(mean_RA = mean(RA))
  
  RA_samples.higherthan.blanks_ASVs <- mean_RA_samples %>%
    filter(mean_RA > mean_RA_blanks$mean_RA) %>%
    pull(ASV)
  
  return(RA_samples.higherthan.blanks_ASVs)
}

RA_samples.higherthan.blanks_ASVs <- list_RA_sample.vs.blank_function(taxa_contam.byfreq_TRUE, ps_raw_decontam_full)
RA_samples.higherthan.blanks_ASVs #Contains 21 ASVs

#Now can examine these in more detail:

plot_freq_function(RA_samples.higherthan.blanks_ASVs) 

#Examine RA vs DNA conc plots again. Any ?false pos (NOT contam) - inspect individually.

# Manually check any individual ASV ?contam ----------------------------------

# Insert ASV name into plotcontam_grid_function bottom line!

plot_contam_grid_function <- function(taxa) {
  
  contam_prev_data <- data.frame(
    prev.sample = 
      taxa_sums(prune_samples(sample_data(ps_decontam_prev)$type == "sample", 
                              ps_decontam_prev))[taxa]/nrow(subset(sample_data(ps_decontam_prev), type=="sample")),
    prev.blank = 
      taxa_sums(prune_samples(sample_data(ps_decontam_prev)$type == "blank", 
                              ps_decontam_prev))[taxa]/nrow(subset(sample_data(ps_decontam_prev), type=="blank"))
  )
  
  plot_contam_prev_sample.vs.blank <- ggplot(data = contam_prev_data, 
                                             aes(x = prev.blank, y = prev.sample)) +
    geom_point() +
    scale_x_continuous(limits=c(0,1), 
                       breaks=seq(0,1, by=0.2), minor_breaks = seq(0,1, by=0.1)) +
    scale_y_continuous(limits=c(0,1), 
                       breaks=seq(0,1, by=0.2), minor_breaks = seq(0,1, by=0.1)) +
    xlab("Prevalence - blanks") + ylab("Prevalence - Samples")
  
  plot_contam_freq.conc <- plot_freq_function(c(taxa))
  
  ps_long_decontam <- ps_raw_decontam_full %>% 
    to_RA() %>% ps_to_df() %>%
    pivot_longer(names_to = "ASV", values_to = "RA", 
                 head(taxa_names(ps$raw), n=1):tail(taxa_names(ps$raw), n = 1)) %>% 
    mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
             fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
    mutate(visit = fct_explicit_na(visit, "blank") %>% 
             fct_relevel(., "V1", "V3", "V4", "V5", "V6", "blank")) %>% 
    mutate(participant = fct_explicit_na(participant, "blank"))
  
  plot_contam_RA.extrun <- ggplot(subset(ps_long_decontam, ASV == taxa), 
                                  aes(x = ext_run, y = RA)) +
    geom_jitter(width = 0.2, size=2.5, alpha=0.6) +
    geom_boxplot(alpha=0.5, outlier.colour = NA) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.y = ggtext::element_markdown()) +
    ggforce::facet_row(~ factor(ASV), space = "free") +
    labs(y = "Relative abundance", x="DNA extraction run") 
  
  plot_contam_RA.sampletype <- ggplot(subset(ps_long_decontam, ASV ==taxa), 
                                      aes(x = sample_type, y = RA)) +
    geom_jitter(width = 0.2, size=2.5, alpha=0.6) +
    geom_boxplot(alpha=0.5, outlier.colour = NA) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.y = ggtext::element_markdown()) +
    ggforce::facet_row(~ factor(ASV), space = "free") +
    labs(y = "Relative abundance", x="Sample type")
  
  plot_contam_RA.participant <- ggplot(subset(ps_long_decontam, ASV ==taxa), 
                                       aes(x = participant, y = RA)) +
    geom_jitter(width = 0.2, size=2.5, alpha=0.6) +
    geom_boxplot(alpha=0.5, outlier.colour = NA) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.y = ggtext::element_markdown()) +
    ggforce::facet_row(~ factor(ASV), space = "free") +
    labs(y = "Relative abundance", x="Participant")
  
  plot_contam_RA.visit <- ggplot(subset(ps_long_decontam, ASV==taxa), 
                                 aes(x = visit, y = RA)) +
    geom_jitter(width = 0.2, size=2.5, alpha=0.6) +
    geom_boxplot(alpha=0.5, outlier.colour = NA) +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.y = ggtext::element_markdown()) +
    ggforce::facet_row(~ factor(ASV), space = "free") +
    labs(y = "Relative abundance", x="Visit")
  
  plotgrid_decontam_bytaxa <- plot_grid(plot_contam_prev_sample.vs.blank, plot_contam_freq.conc, 
                                        plot_contam_RA.extrun, plot_contam_RA.sampletype,
                                        plot_contam_RA.visit, plot_contam_RA.participant,
                                        ncol=2, nrow=3)
  plotgrid_title <- ggdraw() + draw_label(taxa)
  plotgrid_decontam_bytaxa <- plot_grid(plotgrid_title, plotgrid_decontam_bytaxa,
                                        ncol=1, rel_heights = c(0.05,1))
}

plot_contam_grid <- plot_contam_grid_function("Pseudomonas_250") #Replace with any ASV name
plot_contam_grid

# Inspect remaining ASVs, any ?missed contam ---------------------------------

taxa_names(subset_taxa(ps_raw_decontam_full, !taxa_names(ps_raw_decontam_full) %in% taxa_contam.byfreq_TRUE))
#These are the ASVs NOT flagged as contam. Inspect top 200

plot_freq_function(taxa_names(
  subset_taxa(ps_raw_decontam_full, !taxa_names(ps_raw_decontam_full) %in% 
                taxa_contam.byfreq_TRUE))[1:64])
plot_freq_function(taxa_names(
  subset_taxa(ps_raw_decontam_full, !taxa_names(ps_raw_decontam_full) %in% 
                taxa_contam.byfreq_TRUE))[65:128])
plot_freq_function(taxa_names(
  subset_taxa(ps_raw_decontam_full, !taxa_names(ps_raw_decontam_full) %in% 
                taxa_contam.byfreq_TRUE))[129:200])

# Remove contaminants (comb model) ----------------------------------------

length(taxa_contam.byfreq_TRUE) #List ?contam ASV (by comb model)

not_contam <- c("x") #Define any ASVs manually identified as not contaminants / false +ve
is_contam <- c("y") #Define any ASVs manually identified as  contaminants / false -ve

taxa_contam_toremove <- c(
  taxa_contam.byfreq_TRUE[!taxa_contam.byfreq_TRUE %in% not_contam], is_contam) 

ps_raw_after_decontam <- ps_raw_decontam_full %>% 
  subset_taxa(!taxa_names(ps_raw_decontam_full) %in% taxa_contam_toremove) %>% 
  prune_samples(!sample_sums(.) == 0, .)

#How many reads removed by decontam?

sum(sample_sums(ps_raw_decontam_full)) #total number of reads before decontam
sum(sample_sums(ps_raw_after_decontam)) #reads left after decontam 
sum(sample_sums(ps_raw_decontam_full) - sample_sums(ps_raw_after_decontam)) #reads removed by decontam 
sample_data(ps_raw_after_decontam)$decontam_reads <- sample_sums(ps_raw_after_decontam)
(sum(sample_sums(ps_raw_decontam_full) - sample_sums(ps_raw_after_decontam)))/sum(sample_sums(ps_raw_decontam_full)) 
#Proportion reads removed by decontam

#Calculate % of reads removed by decontam per sample

sample_data(ps_raw_after_decontam)$prop_decontam <-
  1-(sample_data(ps_raw_after_decontam)$decontam_reads/sample_data(ps_raw_after_decontam)$raw_reads)

#Make table showing mean, min and max reads remaining after decontam per visit/sample type

table_decontam_reads_by_visit.type <- data.frame(sample_data(ps_raw_after_decontam)) %>% 
  group_by(visit.type) %>%
  summarise(meanreads=mean(decontam_reads),
            minreads=min(decontam_reads),maxreads=max(decontam_reads),
            propdecontam=mean(prop_decontam)) %>% 
  arrange(desc(propdecontam)) %>% kable(.)

table_decontam_reads_by_type <- data.frame(sample_data(ps_raw_after_decontam)) %>% 
  group_by(sample_type) %>%
  summarise(meanreads=mean(decontam_reads),
            minreads=min(decontam_reads),maxreads=max(decontam_reads),
            propdecontam=mean(prop_decontam)) %>% 
  arrange(desc(propdecontam)) %>% kable(.)

#What about mocks? How much of mocks ended up contaminating samples?

ps$RA %>% subset_samples(mock_type == "iso_mock") %>%
  prune_taxa(taxa_contam_toremove, .) %>% sample_sums()

#Plot ?contam (comb) ASVs in iso mocks

bar_mock_RA_contam <- ps$RA %>% subset_samples(mock_type == "iso_mock") %>%
  prune_taxa(taxa_contam_toremove, .) %>%
  create_bar(., n=6, ncol_legend=1) + 
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) + 
  labs(title = "Contaminants (combined model) in mocks")

bar_mock_RA_not.contam <- ps$RA %>% 
  subset_taxa(!taxa_names(ps$RA) %in% taxa_contam_toremove) %>% 
  prune_samples(!sample_sums(.) == 0, .) %>% subset_samples(mock_type == "iso_mock") %>% 
  create_bar(., n=6, ncol_legend=1) + 
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) + 
  labs(title = "Mocks after removal of contaminants") 

plot_grid(bar_mock_RA_contam, bar_mock_RA_not.contam)

# Inspect before & after decontam -----------------------------------------

#Check %reads removed/remaining vs raw reads per samples (colour by sample type)

decontam_removed_vs_readcount_type <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x = raw_reads, y = prop_decontam)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.8) +
  scale_x_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points +
  labs(y="Proportion reads removed", x="Raw read count") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  geom_vline(xintercept = 5000, linetype="dashed") + 
  theme(legend.position = "bottom")

decontam_remaining_vs_readcount_type <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x = raw_reads, y = decontam_reads)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.8) +
  scale_x_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points +
  labs(y="Reads remaining", x="Raw read count") +
  geom_hline(yintercept = 5000, linetype="dashed") +
  geom_vline(xintercept = 5000, linetype="dashed") + 
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500))

#Plot % contam reads against 16S qPCR DNA conc (by sample type)

decontam_removed_vs_conc_type <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.8) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(y="Proportion reads removed", x="DNA concentration (pg/ul)") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  theme(legend.position = "bottom")

#Plot remaining reads against 16S qPCR DNA conc (by sample type)

decontam_remaining_vs_conc_type <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.8) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(y = "Reads remaining", x="DNA concentration (pg/ul)") +
  #To add y-intercept, include: geom_hline(yintercept = x000, linetype="dashed") +
  geom_vline(xintercept = 0.3, linetype="dashed") + theme(legend.position = "bottom") +
  geom_hline(yintercept = 5000, linetype="dashed") +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500))

plotgrid_decontam_removed.vs.remaining.reads_type <-
  plot_grid(decontam_removed_vs_conc_type, decontam_remaining_vs_conc_type, labels="auto")

#Plot % contam reads against 16S qPCR DNA conc

decontam_removed_vs_conc_visit_IN <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("IN")) %>% 
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y="Proportion reads removed (IN)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0, 1.0, by = 0.25), 
                     minor_breaks = seq(0, 1.0, by = 0.125))

decontam_removed_vs_conc_visit_BM <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("BM")) %>% 
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y="Proportion reads removed (BM)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0, 1.0, by = 0.25), 
                     minor_breaks = seq(0, 1.0, by = 0.125))

decontam_removed_vs_conc_visit_MN <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V1", "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("MN")) %>% 
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y="Proportion reads removed (MN)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0, 1.0, by = 0.25), 
                     minor_breaks = seq(0, 1.0, by = 0.125))

decontam_removed_vs_conc_visit_IS <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("IS")) %>% 
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y="Proportion reads removed (IS)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0.0, 0.3), breaks = seq(0, 0.3, by = 0.1), 
                     minor_breaks = seq(0, 0.3, by = 0.05))

#Plot remaining reads against 16S qPCR DNA conc (by study visit for INMNBM only)

decontam_remaining_vs_conc_visit_IN <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("IN")) %>% 
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y = "Reads remaining (IN)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 50000), breaks = seq(0, 50000, by = 10000), 
                     minor_breaks = seq(0, 50000, by = 5000))

decontam_remaining_vs_conc_visit_IS <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("IS")) %>% 
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y = "Reads remaining (IS)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500))

decontam_remaining_vs_conc_visit_BM <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("BM")) %>% 
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y = "Reads remaining (BM)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 50000), breaks = seq(0, 50000, by = 10000), 
                     minor_breaks = seq(0, 50000, by = 5000))

decontam_remaining_vs_conc_visit_MN <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(visit = fct_relevel(visit, "V1","V3", "V4", "V5", "V6")) %>%
  subset(., sample_type %in% c("MN")) %>% 
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=visit), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = visit_colors) + # Color for data points
  labs(y = "Reads remaining (MN)", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 50000), breaks = seq(0, 50000, by = 10000), 
                     minor_breaks = seq(0, 50000, by = 5000))

plotgrid_decontam_removed.vs.remaining.reads_visit_MNBMINIS <- 
  plot_grid(decontam_removed_vs_conc_visit_MN, decontam_removed_vs_conc_visit_BM, decontam_removed_vs_conc_visit_IN, decontam_removed_vs_conc_visit_IS,
            decontam_remaining_vs_conc_visit_MN,  decontam_remaining_vs_conc_visit_BM, decontam_remaining_vs_conc_visit_IN, decontam_remaining_vs_conc_visit_IS,
            ncol = 4, nrow = 2)

# i.e. IN V3 most affected - only low conc low read count sample remaining

#Plot % contam reads against 16S qPCR DNA conc for each sample (by DNA ext run)

decontam_removed_vs_conc_extrun <- ps_raw_after_decontam %>% 
  meta_to_df() %>%
  ggplot(aes(x = qpcr_pgul, y = prop_decontam)) +
  geom_point(aes(color=ext_run), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  labs(y="Proportion reads removed", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom")

#Plot remaining reads against 16S qPCR DNA conc for each sample (by DNA ext run)

decontam_remaining_vs_conc_extrun <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  ggplot(aes(x = qpcr_pgul, y = decontam_reads)) +
  geom_point(aes(color=ext_run), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  labs(y = "Reads remaining", x="DNA concentration (pg/ul)") +
  geom_vline(xintercept = 0.3, linetype="dashed") + 
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500))

#i.e. later ext runs had more reads removed... but this is sample type confounding!
#Because lower biomass sample types extracted in later runs, & higher biomass in earlier runs
#Explore this in more detail: plot DNA ext run against reads removed & remaining by sample type

box_removedreads_extrun_type <- ps_raw_after_decontam %>%
  meta_to_df() %>% 
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x=ext_run, y=prop_decontam, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(x="DNA extraction run", y="Proportion reads removed") +
  geom_hline(yintercept = 0.8, linetype="dashed") +
  theme(legend.position = "bottom")

box_remainingreads_extrun_type <- ps_raw_after_decontam %>%
  meta_to_df() %>%
  mutate(sample_type = fct_explicit_na(sample_type, "blank") %>%
           fct_relevel(., "MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")) %>%
  ggplot(aes(x=ext_run, y=decontam_reads, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(x="DNA extraction run", y="Reads remaining") +
  geom_hline(yintercept = 5000, linetype="dashed") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500))

plotgrid_decontam_removed.vs.remaining.extrun_type <- 
  plot_grid(decontam_removed_vs_conc_extrun, decontam_remaining_vs_conc_extrun,
            box_removedreads_extrun_type, box_remainingreads_extrun_type,
            ncol=2, nrow=2)

#Plot blanks again after decontam

ps_after_decontam_RA <- ps_raw_after_decontam %>% to_RA

bar_decontam_iso_blanks_RA <- subset_samples(ps_after_decontam_RA, blank_type == "iso_blank") %>% 
  create_bar(., n = 20, ncol_legend = 1) +
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) +
  labs(title = "Taxonomic composition of isolation blanks (after decontam)")

bar_decontam_iso_blanks_raw_count <- subset_samples(ps_raw_after_decontam, blank_type == "iso_blank") %>% 
  meta_to_df() %>% rownames_to_column(var="sample_id") %>% 
  ggplot(aes(x = sample_id, y = decontam_reads)) +
  geom_col() +
  geom_text(aes(label = scales::comma(decontam_reads, accuracy = 1)),  angle = 90, 
            color = "white", hjust = 1, size = 2.5, nudge_y = -0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "", y = "Reads remaining")

cow_decontam_iso_blanks <- cowplot::plot_grid(bar_decontam_iso_blanks_RA, 
                                              bar_decontam_iso_blanks_raw_count,align = "v", axis = "lr", nrow = 2, 
                                              rel_heights = c(0.75, 0.25))

plotgrid_isoblank_before.vs.after.decontam <- plot_grid(cow_iso_blanks, 
                                                        cow_decontam_iso_blanks, 
                                                        ncol = 2, nrow = 1)

# Repeat decontam separately for V3.IS and V3.IN (+blanks) ----------

#Can INV3 be salvaged by performing decontam separately on V3 infant samples?

sl(contam.V3INIS, 
   isContaminant(subset_samples(ps_raw_decontam_full, 
                                visit.type %in% c("V3.IS","V3.IN","NA.NA")), 
                 conc = "qpcr_pgul", neg = "is_neg", 
                 method = "combined", threshold = 0.1), 
   dir_path = here(), overwrite = TRUE)

#Make lists of contaminants identified for IN/IS V3 samples

contam_TRUE_V3INIS <- contam.V3INIS %>% filter(contaminant == TRUE)

taxa_contam.byfreq_TRUE_V3INIS <- arrange(contam_TRUE_V3INIS, desc(freq)) %>% 
  rownames() #Arrange contamdf.comb output by desc freq & save taxa as vector

ps_raw_after_decontam_V3INIS <- ps_raw_decontam_full %>% 
  subset_samples(visit.type %in% c("V3.IS","V3.IN", "NA.NA")) %>%
  subset_taxa(!taxa_names(ps_raw_decontam_full) %in% taxa_contam.byfreq_TRUE_V3INIS) %>%
  prune_samples(!sample_sums(.) == 0, .) #phylo containing only decontam V3IS and V3IN

#Plot contam in low biomass datasets 

plot_contam_freq.vs.conc.V3INIS <- plot_frequency(subset_samples(ps_raw_decontam_full, 
                                                                 visit.type %in% c("V3.IS","V3.IN","NA.NA")), 
                                                  taxa_contam.byfreq_TRUE_V3INIS, 
                                                  conc="qpcr_pgul") + xlab("DNA concentration (pg/ul)") +
  ggtitle("Contaminants identified in INV3 and ISV3 samples (combined decontam)")
#Use this to plot ALL contam in V3 infant samples on one plot

#How many reads removed by decontam?

sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                               visit.type %in% c("V3.IS","V3.IN", "NA.NA")))) #total number of reads before decontam
sum(sample_sums(ps_raw_after_decontam_V3INIS)) #reads left after decontam
sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                               visit.type %in% c("V3.IS","V3.IN", "NA.NA")))
) - sum(sample_sums(ps_raw_after_decontam_V3INIS)) #reads removed by decontam
sample_data(ps_raw_after_decontam_V3INIS)$decontam_reads <- sample_sums(ps_raw_after_decontam_V3INIS)
(sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                                visit.type %in% c("V3.IS","V3.IN","NA.NA")))
) - sum(sample_sums(ps_raw_after_decontam_V3INIS))
)/sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                                 visit.type %in% c("V3.IS","V3.IN","NA.NA"))))  #Proportion reads removed by decontam

#Calculate % of reads removed by decontam per sample, add as column to ps

sample_data(ps_raw_after_decontam_V3INIS)$prop_decontam <-
  1-(sample_data(ps_raw_after_decontam_V3INIS)$decontam_reads/sample_data(ps_raw_after_decontam_V3INIS)$raw_reads)

#Make table showing mean, min and max reads remaining after decontam per visit/sample type

table_decontam_reads_by_visit.type_V3INIS <- 
  data.frame(sample_data(ps_raw_after_decontam_V3INIS)) %>% 
  group_by(visit.type) %>%
  summarise(meanreads=mean(decontam_reads),
            minreads=min(decontam_reads),maxreads=max(decontam_reads),
            propdecontam=mean(prop_decontam)) %>% 
  arrange(desc(propdecontam)) %>% kable(.)

#i.e. decontam V3.IN & V3.IS separately - 78% of INV3 reads still removed 
#(as opposed to 88% when done with whole dataset), 
#but only 67% of blanks removed. So INV3 really not looking good!

# Repeat decontam one last time - IN, MN and BM +blanks (i.e. all low biomass sample types) ----------

sl(contam.INMNBM, 
   isContaminant(subset_samples(ps_raw_decontam_full, 
                                !(sample_type %in% c("MO","MS","IS","CS"))), 
                 conc = "qpcr_pgul", neg = "is_neg", 
                 method = "combined", threshold = 0.1), 
   dir_path = here(), overwrite = TRUE)

contam_TRUE_INMNBM <- contam.INMNBM %>% filter(contaminant == TRUE)

taxa_contam.byfreq_TRUE_INMNBM <- arrange(contam_TRUE_INMNBM, desc(freq)) %>% 
  rownames() #Arrange contamdf.comb output by desc freq & save taxa as vector

ps_raw_after_decontam_INMNBM <- ps_raw_decontam_full %>% 
  subset_samples(!(sample_type %in% c("MO","MS","IS","CS"))) %>%
  subset_taxa(!taxa_names(ps_raw_decontam_full) %in% taxa_contam.byfreq_TRUE_INMNBM) %>%
  prune_samples(!sample_sums(.) == 0, .) #phylo containing only decontam MN,IN,BM

sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                               !(sample_type %in% c("MO","MS","IS","CS"))))) #total number of reads before decontam
sum(sample_sums(ps_raw_after_decontam_INMNBM)) #reads left after decontam
sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                               !(sample_type %in% c("MO","MS","IS","CS"))))
) - sum(sample_sums(ps_raw_after_decontam_INMNBM)) #reads removed by decontam
sample_data(ps_raw_after_decontam_INMNBM)$decontam_reads <- sample_sums(ps_raw_after_decontam_INMNBM)
(sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                                !(sample_type %in% c("MO","MS","IS","CS"))))
) - sum(sample_sums(ps_raw_after_decontam_INMNBM))
)/sum(sample_sums(subset_samples(ps_raw_decontam_full, 
                                 !(sample_type %in% c("MO","MS","IS","CS")))))  #Proportion reads removed by decontam
sample_data(ps_raw_after_decontam_INMNBM)$prop_decontam <-
  1-(sample_data(ps_raw_after_decontam_INMNBM)$decontam_reads/sample_data(ps_raw_after_decontam_INMNBM)$raw_reads)

table_decontam_reads_by_visit.type_INMNBM <- 
  data.frame(sample_data(ps_raw_after_decontam_INMNBM)) %>% 
  group_by(visit.type) %>%
  summarise(meanreads=mean(decontam_reads),
            minreads=min(decontam_reads),maxreads=max(decontam_reads),
            propdecontam=mean(prop_decontam)) %>% 
  arrange(desc(propdecontam)) %>% kable(.)

#Still decontam >80% V3IN, and only 71% blanks! i.e. decision to remove ALL V3IN

# Filter rare taxa & low read samples --------------------------------------------------------

#Subset samples with >80% decontam, <5000 reads remaining, and all of V3IN

samples_to_filter <- union(union(union(
  sample_names(subset_samples(ps_raw_after_decontam, prop_decontam>0.8)), 
  sample_names(subset_samples(ps_raw_after_decontam, decontam_reads<5000))), 
  sample_names(subset_samples(ps_raw_after_decontam, visit.type=="V3.IN"))), 
  sample_names(subset_samples(ps_raw_after_decontam, visit.type=="NA.NA"))) 

#Now remove these samples from phylo: 

ps_raw_after_decontam_low <- ps_raw_after_decontam %>% 
  subset_samples(!sample_names(ps_raw_after_decontam) %in% samples_to_filter) %>% 
  prune_samples(!sample_sums(.) == 0, .)

filtered_samples <- merge(as.data.frame(
  table(sample_data(ps_raw_after_decontam)$visit.type)), as.data.frame(
    table(sample_data(ps_raw_after_decontam_low)$visit.type)), 
  by = "Var1", all = TRUE)

colnames(filtered_samples) <- c("visit.type", "before.filter", "after.filter")
filtered_samples[is.na(filtered_samples)] <- 0
filtered_samples$prop.filter <- (filtered_samples$before.filter 
                                 - filtered_samples$after.filter) / filtered_samples$before.filter
filtered_samples <- filtered_samples[order(-filtered_samples$prop.filter), ]
print(filtered_samples)

#Next remove rare taxa:

taxa_notrare <- ps_raw_after_decontam_low %>% prune_taxa(taxa_sums(.)>0, .) %>% 
  to_RA %>%
  filter_taxa(function(x) sum(x > 0.001) >= 1, TRUE) %>%
  taxa_names #These are the taxa to KEEP (i.e. >0.1% RA & present in min 1 sample)
taxa_rare <- setdiff(taxa_names(ps_raw_after_decontam_low), taxa_notrare)

ps_raw_after_decontam_filt <- prune_taxa(taxa_notrare, ps_raw_after_decontam_low) %>% 
  prune_samples(!sample_sums(.) == 0, .)

ntaxa(ps_raw_after_decontam_low) # taxa remaining after decontam
ntaxa(ps_raw_after_decontam_filt) # taxa remaining after filtering rare taxa

#Add column to filtered ps: read count after filtering & prop reads removed by filtering (relative to decontam_reads):

sample_data(ps_raw_after_decontam_filt)$filt_reads <- sample_sums(ps_raw_after_decontam_filt)

sample_data(ps_raw_after_decontam_filt)$prop_filt <-
  1-(sample_data(ps_raw_after_decontam_filt)$filt_reads/sample_data(ps_raw_after_decontam_filt)$decontam_reads)

sum(sample_sums(ps_raw_after_decontam_low)) 
#  reads remaining after decontam & removing low read / high decontam samples
sum(sample_sums(ps_raw_after_decontam_filt)) 
#  reads remaining after filtering rare taxa
sum(sample_sums(ps_raw_after_decontam_low)) - sum(sample_sums(ps_raw_after_decontam_filt)) 
#  reads filtered
(sum(sample_sums(ps_raw_after_decontam_low)) - 
    sum(sample_sums(ps_raw_after_decontam_filt))) / sum(sample_sums(ps_raw_after_decontam_low)) 

min(sample_sums(ps_raw_after_decontam_filt)) #min reads per sample
max(sample_sums(ps_raw_after_decontam_filt)) #max reads per sample
median(sample_sums(ps_raw_after_decontam_filt)) #median reads per sample

#Now inspect effect of filtering by DNA conc, ext run, sample type & visit

#First add columns to show number of ASVs after decontam & filtering

sample_data(ps_raw_after_decontam_filt)$ASV_decontam <- rowSums(t(otu_table(ps_raw_after_decontam_low)) >0)
sample_data(ps_raw_after_decontam_filt)$ASV_filt <- rowSums(t(otu_table(ps_raw_after_decontam_filt)) >0)
sample_data(ps_raw_after_decontam_filt)$prop_filt_ASV <- 
  1-(sample_data(ps_raw_after_decontam_filt)$ASV_filt/sample_data(ps_raw_after_decontam_filt)$ASV_decontam)
colnames(sample_data(ps_raw_after_decontam_filt)) #Check all columns added correctly

#Next plot ASV prop removed vs ASV remaining after filtering, by DNA conc, ext run, type & visit

filt_removed_vs_readcount_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x = filt_reads, y = prop_filt_ASV)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.6) +
  scale_x_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500)) +
  scale_color_manual(values = sample_type_colors) +
  labs(y="Proportion ASV removed", x="Filtered read count") +
  geom_vline(xintercept = 5000, linetype="dashed") + theme(legend.position = "none")

filt_remaining_vs_readcount_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x = filt_reads, y = ASV_filt)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.6) +
  scale_x_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 15000), 
                     minor_breaks = seq(0, 60000, by = 7500)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points +
  labs(y="ASVs remaining", x="Filtered read count") +
  geom_vline(xintercept = 5000, linetype="dashed") + theme(legend.position = "none")

filt_removed_vs_conc_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x = qpcr_pgul, y = prop_filt_ASV)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(y="Proportion ASVs removed", x="DNA concentration (pg/ul)") +
  #To add y-intercept, include: geom_hline(yintercept = 0.xx, linetype="dashed") +
  geom_vline(xintercept = 0.3, linetype="dashed") + theme(legend.position = "none")

filt_remaining_vs_conc_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>%
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x = qpcr_pgul, y = ASV_filt)) +
  geom_point(aes(color=sample_type), fill="white", size=2.5, alpha=0.6) +
  scale_x_log10(labels = function(x) scales::comma(x, accuracy = 0.1)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(y = "ASVs remaining", x="DNA concentration (pg/ul)") +
  #To add y-intercept, include: geom_hline(yintercept = x000, linetype="dashed") +
  geom_vline(xintercept = 0.3, linetype="dashed") + theme(legend.position = "none")

plotgrid_filt_ASVs_reads.conc.type <- plot_grid(
  filt_removed_vs_readcount_type, filt_remaining_vs_readcount_type,
  filt_removed_vs_conc_type, filt_remaining_vs_conc_type)

#Now look at ASVs filtered (removed & remaining) by sample & visit

box_removedASVs_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>% 
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x=sample_type, y=prop_filt_ASV, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Sample type", y="Proportion ASVs removed") +
  theme(legend.position = "bottom")

box_remainingASVs_type <- ps_raw_after_decontam_filt %>%
  meta_to_df() %>% 
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")) %>%
  ggplot(aes(x=sample_type, y=ASV_filt, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Sample type", y="ASVs remaining") +
  theme(legend.position = "none")

box_removedASVs_visit <- 
  subset_samples(ps_raw_after_decontam_filt, sample_type != "CS") %>%
  meta_to_df() %>% 
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN")) %>%
  ggplot(aes(x=visit, y=prop_filt_ASV, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Proportion ASVs removed") +
  theme(legend.position = "none")

box_remainingASVs_visit <- 
  subset_samples(ps_raw_after_decontam_filt, sample_type != "CS") %>%
  meta_to_df() %>% 
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN")) %>%
  ggplot(aes(x=visit, y=ASV_filt, color=sample_type)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(alpha=0.8, 
              position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = sample_type_colors) + # Color for data points
  labs(x="Visit", y="ASVs remaining") +
  theme(legend.position = "none")

plotgrid_filt_ASVs_type_visit <- plot_grid(
  box_removedASVs_type, box_remainingASVs_type,
  box_removedASVs_visit, box_remainingASVs_visit)

# Before finalising phylo, check top 256 remaining taxa  --------

plot_freq_function_final <- function(taxa) {
  plot_frequency(ps_raw_after_decontam_filt, taxa, conc="qpcr_pgul") + 
    xlab("DNA concentration (pg/ul)")} 

plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[1:64]) 
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[65:128]) 
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[129:192]) 
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[193:256]) 

png(file=here("plots", "RemainingASV1.png"), width = 15, height = 10, units = "in",res=300)
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[1:64]) 
dev.off()

png(file=here("plots", "RemainingASV2.png"), width = 15, height = 10, units = "in",res=300)
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[65:128]) 
dev.off()

png(file=here("plots", "RemainingASV3.png"), width = 15, height = 10, units = "in",res=300)
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[129:192]) 
dev.off()

png(file=here("plots", "RemainingASV4.png"), width = 15, height = 10, units = "in",res=300)
plot_freq_function_final(taxa_names(ps_raw_after_decontam_filt)[193:256]) 
dev.off()

#Manually check individual ASVs as necessary. Top 256 all look fine!

# Compare with Salter, any missed contam? ----------------------------------------------------------

Salter <- read_tsv(here("notes/Salter_blank_contam.tsv")) 
Salter_genus <- Salter$Genus #93 contam genera from Salter paper
filt_genus <- sub("_.*","", taxa_names(ps_raw_after_decontam_filt)) %>% unique()

Salter.vs.filt_genus <- intersect(Salter_genus, filt_genus) 
#Remove some genera, as not to be treated as contam in this dataset:
Salter.vs.filt_genus <- Salter.vs.filt_genus[!(Salter.vs.filt_genus %in% 
                                                 c("Streptococcus", "Corynebacterium"))]

#Which of these are genera that have already been removed by decontam?
contam_removed_genus <- sub("_.*","", taxa_contam_toremove) %>% unique()

Salter.vs.filt.vs.removed_genus <- intersect(Salter.vs.filt_genus,contam_removed_genus)

Salter.vs.filt_ASV <- taxa_names(ps_raw_after_decontam_filt)[grep(paste(Salter.vs.filt_genus, collapse = "|"), 
                                                                  taxa_names(ps_raw_after_decontam_filt))]

Salter.vs.filt.vs.removed_ASV <- taxa_names(ps_raw_after_decontam_filt)[grep(paste(Salter.vs.filt.vs.removed_genus, collapse = "|"), 
                                                                             taxa_names(ps_raw_after_decontam_filt))]

plot_freq_function_final(Salter.vs.filt.vs.removed_ASV[1:20]) 

# Finalise phylo for stats analysis ---------------------------------------

#Save final phylo as ps_final & convert to RA, then write to results folder:

sample_data(ps_raw_after_decontam_filt)$visit <- fct_relevel(sample_data(ps_raw_after_decontam_filt)$visit, "V1", "V3", "V4", "V5", "V6")
sample_data(ps_raw_after_decontam_filt)$sample_type <- fct_relevel(sample_data(ps_raw_after_decontam_filt)$sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS")
sample_data(ps_raw_after_decontam_filt)$participant_type <- fct_relevel(sample_data(ps_raw_after_decontam_filt)$participant_type, "mother", "infant", "sibling")

ps_final <- ps_raw_after_decontam_filt
ps_final_RA <- ps_final %>% to_RA

sl(ps_final, ps_final, overwrite = T, dir_path = here("results"))
sl(ps_final_RA, ps_final_RA, overwrite = T, dir_path = here("results"))

table(sample_data(ps_final)$visit.type, sample_data(ps_final)$participant)

as.data.frame(rownames(sample_data(ps_final)))
nrow(sample_data(ps_final))

#Record session info:

sessionInfo()