# Load packages & functions -----------------------------------------------

library(tidyverse)
library(phyloseq)
library(here)
library(decontam)
library(knitr)
library(ggpubr)
library(cowplot)
library(data.table)
library(vegan)
library(RColorBrewer)

#Load additional functions (created by Dr Wouter de Steenhuijsen Piters - https://gitlab.com/wsteenhu/muis_trx)

source(here("src/load_dada.R"))
source(here("src/utils.R"))

# set paths
knitr::opts_knit$set(root.dir=".", aliases=c(h = "fig.height", w = "fig.width", ow = "out.width"))

#Load phylo

phylo <- readRDS(here("results/ps_final.Rds"))
phylo_RA <- readRDS(here("results/ps_final_RA.Rds"))
phylo_meta <- data.frame(sample_data(phylo))

# Beta diversity NMDS ----------------------------------------------------------

ASV_raw <- as.data.frame(t(as(otu_table(phylo), "matrix")))
ASV_RA <- as.data.frame(t(as(otu_table(phylo_RA), "matrix")))
braycurtis <- vegdist(ASV_RA, method="bray")
braycurtis.matrix <- as.matrix(braycurtis)

set.seed(27) #Set random seed as starting point. If don't set seed, wont' be reproducible!
ordBC <- metaMDS(ASV_RA, "bray", trymax=1000, k=2, trace = TRUE, autotransform = F)
print(ordBC)

#Aiming for stress 0-0.1 (low) or 0.1-0.2 (moderate). Poor if >0.2.
#Note best result not repeated, despite trymax=1000. But stress/dimensions fine.

# plot the samples as points in a 2 dimensional space, 
NMDS <- left_join(phylo_meta %>% rownames_to_column("sample_id") %>%
                    mutate(sample_id=as.character(sample_id)), 
                  ordBC$points %>%
                    data.frame() %>%
                    rownames_to_column("sample_id") %>%
                    select(sample_id, NMDS1=MDS1, NMDS2=MDS2), 
                  by="sample_id")

NMDS_infant <- NMDS %>% filter(participant_type == "infant")
NMDS_mother <- NMDS %>% filter(participant_type == "mother") %>% filter(visit %in% c("V3","V4","V5","V6"))
NMDS_mother_V1 <- NMDS %>% filter(participant_type == "mother") %>% filter(visit %in% c("V1"))

#Now plot & draw ellipses around groups of samples (participant type)

df_ellipse <- NMDS %>%
  group_by(participant_type) %>%
  do(vegan:::veganCovEllipse(
    cov.wt(
      data.frame(NMDS1=.$NMDS1, NMDS2=.$NMDS2))$cov %>% data.frame,
    center=c(mean(.$NMDS1), mean(.$NMDS2))) %>% data.frame(check.names=F)) %>% 
  ungroup

NMDS_participant <- ggplot(data=NMDS, aes(NMDS1, NMDS2)) + 
  geom_polygon(data=df_ellipse, aes(x=NMDS1, y=NMDS2, colour=participant_type, fill=participant_type), alpha=0.4, show.legend = F) +
  geom_point(aes(colour = participant_type), size=2, alpha=1) + 
  geom_text(aes(label = sprintf("Stress = %.3f", ordBC$stress)), x = Inf, y = -Inf, hjust = 1.1, vjust = -1, size = 3) +
  scale_color_manual(values = participant_type_colors) +
  scale_fill_manual(values = participant_type_colors) +
  theme_bw() +
  theme(legend.position = "bottom", aspect.ratio = 1)

#Function to plot NMDS for other factors

plot_NMDS <- function(NMDS, group, cols, ...)  {
  df_ellipse <- NMDS %>%
    group_by({{group}}) %>% # change; also possible to use two groups
    do(vegan:::veganCovEllipse(
      cov.wt(
        data.frame(NMDS1=.$NMDS1, NMDS2=.$NMDS2))$cov %>% data.frame,
      center=c(mean(.$NMDS1), mean(.$NMDS2))) %>% data.frame(check.names=F)) %>% 
    ungroup
  
  p1 <- ggplot(data=NMDS, aes(NMDS1, NMDS2)) + 
    geom_polygon(data=df_ellipse, aes(x=NMDS1, y=NMDS2, colour={{group}}, fill={{group}}), show.legend = F) +
    geom_point(aes(colour = {{group}}), size=2, alpha=1) +
    scale_color_manual(..., values=alpha(cols, 1)) +
    scale_fill_manual(..., values=alpha(cols, 0.4)) +
    theme_bw() +
    theme(legend.position = "bottom", aspect.ratio = 1)
  p1
}

NMDS_sampletype <- plot_NMDS(NMDS, sample_type, cols=sample_type_colors, name="Sample type") +
  ggtitle("NMDS by sample type")

#Function to plot NMDS by 2 factors

plot_NMDS2 <- function(NMDS, group1, group2, cols, ...)  {
  df_ellipse <- NMDS %>%
    group_by({{group1}}, {{group2}}) %>% 
    do(vegan:::veganCovEllipse(
      cov.wt(
        data.frame(NMDS1=.$NMDS1, NMDS2=.$NMDS2))$cov %>% data.frame,
      center=c(mean(.$NMDS1), mean(.$NMDS2))) %>% data.frame(check.names=F)) %>% 
    ungroup()
  
  p1 <- ggplot(data=NMDS, aes(NMDS1, NMDS2)) + 
    geom_polygon(data=df_ellipse, aes(x=NMDS1, y=NMDS2, 
                                      color = interaction({{group1}}, {{group2}}), 
                                      fill = interaction({{group1}}, {{group2}})), show.legend = F) +
    geom_point(aes(color = interaction({{group1}}, {{group2}})), size=2, alpha=0.8) +
    scale_color_manual(..., values=alpha(cols, 1)) +
    scale_fill_manual(..., values=alpha(cols, 0.5)) +
    theme_bw() +
    theme(legend.position = "bottom", aspect.ratio = 1)
  
  return(p1)
}

plot_NMDS_visit.type_infant <- plot_NMDS2(NMDS_infant, 
                                          group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) +
  ggtitle("Infants by sample type and visit")

plot_NMDS_visit.type_mother <- plot_NMDS2(NMDS_mother,  
                                          group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) +
  ggtitle("Mothers by sample type and visit")

plot_NMDS_visit.type_all <- plot_NMDS2(filter(NMDS, visit != "V1"),  
                                       group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) +
  ggtitle("All participants by sample type and visit")
plot_NMDS_visit.type_all + guides(color = guide_legend(nrow = 4, ncol = 6))

#Now do a separate plot for each niche to compare timpeoints:

NMDS_visit_MO <- plot_NMDS2(filter(NMDS, sample_type == "MO", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("MO (V3-V6)") + theme(legend.position = "none")
NMDS_visit_MS <- plot_NMDS2(filter(NMDS, sample_type == "MS", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("MS (V3-V6)") + theme(legend.position = "none")
NMDS_visit_MN <- plot_NMDS2(filter(NMDS, sample_type == "MN", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("MN (V3-V6)") + theme(legend.position = "none")
NMDS_visit_BM <- plot_NMDS2(filter(NMDS, sample_type == "BM", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("BM (V3-V6)") + theme(legend.position = "none")
NMDS_visit_IS <- plot_NMDS2(filter(NMDS, sample_type == "IS", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("IS (V3-V6)") + theme(legend.position = "none")
NMDS_visit_IN <- plot_NMDS2(filter(NMDS, sample_type == "IN", visit!="V1"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("IN (V4-V6)") + theme(legend.position = "none")

#And now a plot for each timepoint, to compare niches:

NMDS_visit_V3 <- plot_NMDS2(filter(NMDS, visit=="V3"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("Visit 3") + theme(legend.position = "none")
NMDS_visit_V4 <- plot_NMDS2(filter(NMDS, visit=="V4"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("Visit 4") + theme(legend.position = "none")
NMDS_visit_V5 <- plot_NMDS2(filter(NMDS, visit=="V5"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("Visit 5") + theme(legend.position = "none")
NMDS_visit_V6 <- plot_NMDS2(filter(NMDS, visit=="V6"), 
                            group1 = visit, group2 = sample_type, cols=sampletype_visit_colors) + 
  ggtitle("Visit 6") + theme(legend.position = "none")

plotgrid_NMDS_each.visit.niche1 <- plot_grid(NMDS_visit_MO, NMDS_visit_MS, NMDS_visit_MN, 
                                             NMDS_visit_BM, NMDS_visit_IS, NMDS_visit_IN, ncol=3, labels="auto")
plotgrid_NMDS_each.visit.niche2 <- plot_grid(NMDS_visit_V3, NMDS_visit_V4, 
                                             NMDS_visit_V5, NMDS_visit_V6, ncol=4, labels=c("g","h","i","j"))

#Now compare based on inoculated, maternal abx & sibling status (V3-V6)

plot_NMDS_inoc <- plot_NMDS2(filter(NMDS, sample_type!="CS"), 
                             group1 = inoculated, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Inoculated TRUE/FALSE") 

plot_NMDS_mother_inoc <- plot_NMDS2(NMDS_mother, 
                                    group1 = inoculated, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Maternal samples, inoculated TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_mother_nlac <- plot_NMDS2(NMDS_mother, 
                                    group1 = m_nlac_col, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Maternal samples, Nlac-colonised TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_mother_abx <- plot_NMDS2(NMDS_mother, 
                                   group1 = abx_mum_todate, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Maternal samples, antibiotics to date TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_mother_sib <- plot_NMDS2(NMDS_mother, 
                                   group1 = sib_5yr, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Maternal samples, child <5y TRUE/FALSE") +  
  theme(legend.position = "none")

plot_NMDS_infant_inoc <- plot_NMDS2(NMDS_infant, 
                                    group1 = inoculated, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Infant samples, maternal inoculation TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_infant_abx <- plot_NMDS2(NMDS_infant, 
                                   group1 = abx_mum_todate, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Infant samples, maternal antibiotics to date TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_infant_sib <- plot_NMDS2(NMDS_infant, 
                                   group1 = sib_5yr, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Infant samples, sibling <5y TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_infant_formula <- plot_NMDS2(NMDS_infant, 
                                       group1 = formula_todate, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Infant samples, formula to date TRUE/FALSE") +  
  theme(legend.position = "none")
plot_NMDS_infant_mat_nlac <- plot_NMDS2(NMDS_infant, 
                                        group1 = m_nlac_col, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Infant samples, Nlac-colonised mother TRUE/FALSE") +  
  theme(legend.position = "none")

plotgrid_NMDS_inoc_abx_sib <- 
  plot_grid(plot_NMDS_mother_inoc, plot_NMDS_mother_sib,
            plot_NMDS_mother_nlac, plot_NMDS_mother_abx,
            plot_NMDS_infant_inoc, plot_NMDS_infant_sib,
            plot_NMDS_infant_mat_nlac, plot_NMDS_infant_abx,
            ncol=4)

#And look for baseline differences between inoculated vs uninoculated mums, & those with vs without siblings
plot_NMDS_V1.inoc_mother <- plot_NMDS2(NMDS_mother_V1,  
                                       group1 = inoculated, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Visit 1 - inoculated vs baseline-colonised mothers")

plot_NMDS_V1.sibs_mother <- plot_NMDS2(NMDS_mother_V1,  
                                       group1 = sib_5yr, group2 = sample_type, cols=sampletype_TF_colors) +
  ggtitle("Visit 1 - siblings vs no siblings")

# Beta diversity PERMANOVA building ---------------------

meta_infant <- phylo_meta %>% filter(participant_type == "infant")
meta_mother <- phylo_meta %>% filter(participant_type == "mother")

braycurtis_infant <- braycurtis.matrix[
  rownames(braycurtis.matrix) %in% row.names(meta_infant),
  rownames(braycurtis.matrix) %in% row.names(meta_infant)]

braycurtis_mother <- braycurtis.matrix[
  rownames(braycurtis.matrix) %in% row.names(meta_mother),
  rownames(braycurtis.matrix) %in% row.names(meta_mother)]

box_braycurtis <- boxplot(
  list(all_samples = braycurtis.matrix,
       maternal_samples = braycurtis_mother,
       infant_samples = braycurtis_infant),
  names = c("All samples", "Maternal samples", "Infant samples"),
  main = "Bray-Curtis dissimilarity",
  ylab = "Distance",
  col = c("grey", "lightpink1", "darkseagreen2"))

#PERMANOVA

#Univariable to look at overall contribution to variance of niche (whole dataset):

set.seed(123)
permanova_niche_overall_V3to6 <- 
  adonis2(braycurtis.matrix[
    rownames(braycurtis.matrix) %in% row.names(subset(phylo_meta, visit!="V1", sample_type!="CS")),
    rownames(braycurtis.matrix) %in% row.names(subset(phylo_meta, visit!="V1", sample_type!="CS"))]
    ~ sample_type, 
    data = subset(phylo_meta, visit!="V1", sample_type!="CS"), 
    strata= subset(phylo_meta, visit!="V1", sample_type!="CS")$participant, 
    permutations = 1000, method="marginal")

#Notes that covariates that introduce multicolinearity will result in -inf F / 0 R2
#e.g. participant_type & sample_type, or participant & inoc / sib_5y.

#Create function to test multiple visit or sample types at once!

#First by visit for each niche...

permanova_results_visit <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(permanova_results_visit) <- c("niche", "R2", "pvalue")

for (niche in unique(subset(phylo_meta, sample_type!="CS")$sample_type)) {
  # Filter data based on sample_type
  subset_data_meta <- filter(subset(phylo_meta, sample_type!="CS",visit!="V1"), sample_type==niche)
  subset_data_braycurtis <- braycurtis.matrix[
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta),
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta)]
  # Construct the formula
  formula <- as.formula(paste("subset_data_braycurtis ~ visit"))
  # Perform PERMANOVA
  set.seed(27)
  model <- adonis2(formula, subset_data_meta, 
                   strata=subset_data_meta$participant, permutations = 1000)
  
  # Store the results
  tmp <- data.frame(niche = paste(niche), R2 = model$R2[1], pvalue = model$Pr[1])
  permanova_results_visit <- rbind(permanova_results_visit, tmp)
}

permanova_results_visit$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_visit$pvalue, method="BH")
permanova_results_visit <- permanova_results_visit %>%
  arrange(adj_pval)
permanova_results_visit #print results

#Now by niche for each visit...

permanova_results_niche <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(permanova_results_niche) <- c("timepoint", "R2", "pvalue")

for (timepoint in unique(subset(phylo_meta, visit!="V1", sample_type!="CS")$visit)) {
  # Filter data based on sample_type
  subset_data_meta <- filter(phylo_meta, visit==timepoint)
  subset_data_braycurtis <- braycurtis.matrix[
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta),
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta)]
  formula <- as.formula(paste("subset_data_braycurtis ~ sample_type"))
  set.seed(27)
  model <- adonis2(formula, subset_data_meta, 
                   strata=subset_data_meta$participant, permutations = 1000)
  tmp <- data.frame(timepoint = paste(timepoint), R2 = model$R2[1], pvalue = model$Pr[1])
  permanova_results_niche <- rbind(permanova_results_niche, tmp)
}

permanova_results_niche$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_niche$pvalue, method="BH")
permanova_results_niche <- permanova_results_niche %>%
  arrange(timepoint)
permanova_results_niche #print results

#Repeat this for mother vs infant separately (as expect more visit variance for infant)

permanova_results_niche_infant <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(permanova_results_niche_infant) <- c("timepoint", "R2", "pvalue")

for (timepoint in unique(meta_infant$visit)) {
  subset_data_meta <- filter(meta_infant, visit==timepoint)
  subset_data_braycurtis <- braycurtis.matrix[
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta),
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta)]
  formula <- as.formula(paste("subset_data_braycurtis ~ sample_type"))
  set.seed(27)
  model <- adonis2(formula, subset_data_meta, 
                   strata=subset_data_meta$participant, permutations = 1000)
  tmp <- data.frame(timepoint = paste(timepoint), R2 = model$R2[1], pvalue = model$Pr[1])
  permanova_results_niche_infant <- rbind(permanova_results_niche_infant, tmp)
}

permanova_results_niche_infant$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_niche_infant$pvalue, method="BH")
permanova_results_niche_infant <- permanova_results_niche_infant %>%
  arrange(timepoint)
permanova_results_niche_infant #print results

#... And again for mother:

permanova_results_niche_mother <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(permanova_results_niche_mother) <- c("timepoint", "R2", "pvalue")

for (timepoint in unique(meta_mother$visit)) {
  subset_data_meta <- filter(meta_mother, visit==timepoint)
  subset_data_braycurtis <- braycurtis.matrix[
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta),
    rownames(braycurtis.matrix) %in% row.names(subset_data_meta)]
  formula <- as.formula(paste("subset_data_braycurtis ~ sample_type"))
  set.seed(27)
  model <- adonis2(formula, subset_data_meta, 
                   strata=subset_data_meta$participant, permutations = 1000)
  tmp <- data.frame(timepoint = paste(timepoint), R2 = model$R2[1], pvalue = model$Pr[1])
  permanova_results_niche_mother <- rbind(permanova_results_niche_mother, tmp)
}

permanova_results_niche_mother$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_niche_mother$pvalue, method="BH")
permanova_results_niche_mother <- permanova_results_niche_mother %>%
  arrange(timepoint)
permanova_results_niche_mother #print results

#Now by meta variables: first screen using niche and timepoint.niche

permanova_variables_infant <- 
  c("abx_mum_todate","formula_todate","bf_since_last", 
    "inoculated","m_nlac_col", "sib_5yr", "cov_todate", "RTI_infant_todate")

# Create an empty dataframe for results
permanova_results_infant_variables <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(permanova_results_infant_variables) <- c("variable", "timepoint.niche", "R2", "pvalue")

# Loop through each variable
for (var in permanova_variables_infant) {
  # Loop through different niches and time points
  for (timepoint.niche in unique(meta_infant$visit.type)) {
    # Filter data based on sample_type and visit
    subset_data_meta <- filter(meta_infant, visit.type == timepoint.niche)
    subset_data_braycurtis <- braycurtis_infant[
      rownames(braycurtis_infant) %in% row.names(subset_data_meta),
      rownames(braycurtis_infant) %in% row.names(subset_data_meta)]
    formula <- as.formula(paste("subset_data_braycurtis ~", var))
    set.seed(27)
    model <- adonis2(formula, subset_data_meta, permutations = 1000)
    tmp <- data.frame(variable = var, timepoint.niche = paste(timepoint.niche), 
                      R2 = model$R2[1], pvalue = model$Pr[1])
    permanova_results_infant_variables <- rbind(permanova_results_infant_variables, tmp)
  }
}

permanova_results_infant_variables$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_infant_variables$pvalue, method="BH")
permanova_results_infant_variables <- permanova_results_infant_variables %>%
  arrange(pvalue)
permanova_results_infant_variables #print results

#Now same again for mum:

permanova_variables_mother <- 
  c("abx_mum_todate", "inoculated","m_nlac_col", "sib_5yr", "cov_todate", "RTI_mum_todate")

# Create an empty dataframe for results
permanova_results_mother_variables <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(permanova_results_mother_variables) <- c("variable", "timepoint.niche", "R2", "pvalue")

# Loop through each variable
for (var in permanova_variables_mother) {
  # Loop through different niches and time points
  for (timepoint.niche in unique(meta_mother$visit.type)) {
    # Filter data based on sample_type and visit
    subset_data_meta <- filter(meta_mother, visit.type == timepoint.niche)
    subset_data_braycurtis <- braycurtis_mother[
      rownames(braycurtis_mother) %in% row.names(subset_data_meta),
      rownames(braycurtis_mother) %in% row.names(subset_data_meta)]
    formula <- as.formula(paste("subset_data_braycurtis ~", var))
    set.seed(27)
    model <- adonis2(formula, subset_data_meta, permutations = 1000)
    tmp <- data.frame(variable = var, timepoint.niche = paste(timepoint.niche), 
                      R2 = model$R2[1], pvalue = model$Pr[1])
    permanova_results_mother_variables <- rbind(permanova_results_mother_variables, tmp)
  }
}

permanova_results_mother_variables$adj_pval <- #posthoc adjustment of pvalues
  p.adjust(permanova_results_mother_variables$pvalue, method="BH")
permanova_results_mother_variables <- permanova_results_mother_variables %>%
  arrange(pvalue)
permanova_results_mother_variables #print results

#Adj p-val >0.05 for all, and unadjusted <0.1 for abx_mum_todate, m_lac_col, inoc, covid.

#Now try multivariable model: use niche, visit, inoculated, m_nlac_col & abx for mum (visit 3-6 only)
#NOTE: covid & bf_since_last removed, as not significant in any multivariable model

set.seed(123)
permanova_multi_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1")),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1"))]
    ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate, 
    data = subset(meta_mother, sample_type=="MO" & visit!="V1"), 
    strata=subset(meta_mother, sample_type=="MO" & visit!="V1")$participant, 
    permutations = 1000, method="marginal")

#And so on, for each other niche...

#Try to untangle inoculated, sib_5yr, abx_mum_todate & m_nlac_col

#Effect of sibs on inoculated=TRUE subset (repeat for each niche...)

set.seed(123)
permanova_sib_inocTRUE_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & inoculated=="TRUE")),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & inoculated=="TRUE"))]
    ~ visit + sib_5yr + m_nlac_col + abx_mum_todate, 
    data = subset(meta_mother, sample_type=="MO" & visit!="V1" & inoculated=="TRUE"), 
    strata=subset(meta_mother, sample_type=="MO" & visit!="V1" & inoculated=="TRUE")$participant, 
    permutations = 1000, method="marginal")

#Effect of inoculation on sibs=TRUE subset (repeat for each niche...)

set.seed(123)
permanova_inoc_sibTRUE_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & sib_5yr=="TRUE")),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & sib_5yr=="TRUE"))]
    ~ visit + inoculated + m_nlac_col + abx_mum_todate, 
    data = subset(meta_mother, sample_type=="MO" & visit!="V1" & sib_5yr=="TRUE"), 
    strata=subset(meta_mother, sample_type=="MO" & visit!="V1" & sib_5yr=="TRUE")$participant, 
    permutations = 1000, method="marginal")

#Now m_nlac col for abx TRUE/FALSE... (repeat for each niche...)

set.seed(123)
permanova_abx_nlacTRUE_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & m_nlac_col=="TRUE")),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit!="V1" & m_nlac_col=="TRUE"))]
    ~ visit + inoculated + sib_5yr + abx_mum_todate, 
    data = subset(meta_mother, sample_type=="MO" & visit!="V1" & m_nlac_col=="TRUE"), 
    strata=subset(meta_mother, sample_type=="MO" & visit!="V1" & m_nlac_col=="TRUE")$participant, 
    permutations = 1000, method="marginal")

#Before & after inoc (V1 vs V3/V4) maternal (NO abx subset) (repeat for each maternal niche)

set.seed(123)
permanova_V1V3_noabx_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V3"))),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V3")))]
    ~ visit, 
    data = subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V3")), 
    strata=subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V3"))$participant, 
    permutations = 1000, method="marginal")

set.seed(123)
permanova_V1V4_noabx_MO <- 
  adonis2(braycurtis_mother[
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V4"))),
    rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V4")))]
    ~ visit, 
    data = subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V4")), 
    strata=subset(meta_mother, sample_type=="MO" & abx_mum_todate=="FALSE" & visit %in% c("V1","V4"))$participant, 
    permutations = 1000, method="marginal")

# Adding PERMANOVA terms iteratively --------------------------------------

set.seed(123)
permanova_build_IN_v <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vi <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vs <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + sib_5yr, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vm <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_va <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vis <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vim <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_via <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vsm <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + sib_5yr + m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vsa <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + sib_5yr + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vma <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + m_nlac_col + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vism <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_visa <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr +  abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vima <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + m_nlac_col + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_vsma <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + sib_5yr + m_nlac_col + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_visma <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_vi <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + visit:inoculated, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_vs <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + visit:sib_5yr, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_vm <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + visit:m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_va <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + visit:abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_is <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + inoculated:sib_5yr, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_im <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + inoculated:m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_ia <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + inoculated:abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_sm <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + sib_5yr:m_nlac_col, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_sa <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + sib_5yr:abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

set.seed(123)
permanova_build_IN_int_ma <- adonis2(braycurtis_infant[
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN")),
  rownames(braycurtis_infant) %in% row.names(subset(meta_infant, sample_type=="IN"))]
  ~ visit + inoculated + sib_5yr + m_nlac_col + abx_mum_todate + m_nlac_col:abx_mum_todate, 
  data = subset(meta_infant, sample_type=="IN"), 
  strata=subset(meta_infant, sample_type=="IN")$participant, 
  permutations = 1000, method="marginal")

#Repeat for rest of niches...

# V1 inoc vs uninoc? (maternal only)  ---------------------------------------

set.seed(123)
permanova_V1_inoc_MO <- adonis2(braycurtis_mother[
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit=="V1")),
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MO" & visit=="V1"))]
  ~ inoculated,
  data = subset(meta_mother, sample_type=="MO" & visit=="V1"), 
  permutations = 1000)

set.seed(123)
permanova_V1_inoc_MS <- adonis2(braycurtis_mother[
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MS" & visit=="V1")),
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MS" & visit=="V1"))]
  ~ inoculated, 
  data = subset(meta_mother, sample_type=="MS" & visit=="V1"), 
  permutations = 1000)

set.seed(123)
permanova_V1_inoc_MN <- adonis2(braycurtis_mother[
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MN" & visit=="V1")),
  rownames(braycurtis_mother) %in% row.names(subset(meta_mother, sample_type=="MN" & visit=="V1"))]
  ~ inoculated, 
  data = subset(meta_mother, sample_type=="MN" & visit=="V1"), 
  permutations = 1000)

# Check homogeneity of multivariate dispersions (variance/covariance) --------

#betadisper is a multivariate analogue of Levene's test for homogeneity of variances. 
#Non-euclidean distances between objects and group centres (centroids or medians) 
#are handled by reducing the original distances to principal coordinates

par(mar = c(1, 1, 1, 1))  # Adjust margins as needed

betadisper_niche <- betadisper(braycurtis, phylo_meta$sample_type)
anova(betadisper_niche)
permutest(betadisper_niche, pairwise=TRUE, permutations=100)
betadisper_niche.HSD <- TukeyHSD(betadisper_niche)
plot(betadisper_niche, col = rep(sample_type_colors, alpha = 0.1), ellipse=TRUE, hull=FALSE)
plot(betadisper_niche.HSD)
boxplot(betadisper_niche)

betadisper_timepoint <- betadisper(braycurtis, phylo_meta$visit)
anova(betadisper_timepoint)
permutest(betadisper_timepoint, pairwise=TRUE, permutations=100)
betadisper_timepoint.HSD <- TukeyHSD(betadisper_timepoint)
plot(betadisper_timepoint, col = rep(visit_colors, alpha = 0.1), ellipse=TRUE, hull=FALSE)
plot(betadisper_timepoint)
boxplot(betadisper_timepoint)

betadisper_inoc <- betadisper(braycurtis, phylo_meta$inoculated)
anova(betadisper_inoc)
permutest(betadisper_inoc, pairwise=TRUE, permutations=100)
betadisper_inoc.HSD <- TukeyHSD(betadisper_inoc)
plot(betadisper_inoc, ellipse=TRUE, hull=FALSE)
plot(betadisper_inoc.HSD)
boxplot(betadisper_inoc)

betadisper_nlac <- betadisper(braycurtis, phylo_meta$m_nlac_col)
anova(betadisper_nlac)
permutest(betadisper_nlac, pairwise=TRUE, permutations=100)
betadisper_nlac.HSD <- TukeyHSD(betadisper_nlac)
plot(betadisper_nlac, ellipse=TRUE, hull=FALSE)
plot(betadisper_nlac.HSD)
boxplot(betadisper_nlac)

betadisper_sib <- betadisper(braycurtis, phylo_meta$sib_5y)
anova(betadisper_sib)
permutest(betadisper_sib, pairwise=TRUE, permutations=100)
betadisper_sib.HSD <- TukeyHSD(betadisper_sib)
plot(betadisper_sib, ellipse=TRUE, hull=FALSE)
plot(betadisper_sib.HSD)
boxplot(betadisper_sib)

betadisper_abx <- betadisper(braycurtis, phylo_meta$abx_mum_todate)
anova(betadisper_abx)
permutest(betadisper_abx, pairwise=TRUE, permutations=100)
betadisper_abx.HSD <- TukeyHSD(betadisper_abx)
plot(betadisper_abx, ellipse=TRUE, hull=FALSE)
plot(betadisper_abx.HSD)
boxplot(betadisper_abx)

par(mar = c(5, 4, 4, 2))  # Reset margins to default

# Related vs unrelated mothers - Bray Curtis --------------------------------------------

braycurtis_related <- as.data.frame(matrix(NA, nrow = nrow(braycurtis.matrix), ncol = ncol(braycurtis.matrix)))
rownames(braycurtis_related) <- rownames(braycurtis.matrix)
colnames(braycurtis_related) <- colnames(braycurtis.matrix)

# Create a data frame to store the unrelated subset
braycurtis_unrelated <- as.data.frame(matrix(NA, nrow = nrow(braycurtis.matrix), ncol = ncol(braycurtis.matrix)))
rownames(braycurtis_unrelated) <- rownames(braycurtis.matrix)
colnames(braycurtis_unrelated) <- colnames(braycurtis.matrix)

# Fill in values for the related subset
for (row_name in rownames(braycurtis.matrix)) {
  for (col_name in colnames(braycurtis.matrix)) {
    if (substr(row_name, 5, 6) == substr(col_name, 5, 6)) { #where 5th & 6th char = participant ID!
      braycurtis_related[row_name, col_name] <- braycurtis.matrix[row_name, col_name]
    }
  }
}

# Fill in values for the unrelated subset
for (row_name in rownames(braycurtis.matrix)) {
  for (col_name in colnames(braycurtis.matrix)) {
    if (substr(row_name, 5, 6) != substr(col_name, 5, 6)) {
      braycurtis_unrelated[row_name, col_name] <- braycurtis.matrix[row_name, col_name]
    }
  }
}

#Now subset by sample types:
braycurtis_related_V4INMN <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V4")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]
braycurtis_related_V5INMN <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V5")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]
braycurtis_related_V6INMN <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V6")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]

braycurtis_related_V3ISMS <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V3")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_related_V4ISMS <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V4")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_related_V5ISMS <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V5")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_related_V6ISMS <- braycurtis_related[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V6")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]

braycurtis_unrelated_V4INMN <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V4")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]
braycurtis_unrelated_V5INMN <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V5")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]
braycurtis_unrelated_V6INMN <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IN" & visit=="V6")),
  row.names(subset(phylo_meta, sample_type == "MN"))
]

braycurtis_unrelated_V3ISMS <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V3")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_unrelated_V4ISMS <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V4")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_unrelated_V5ISMS <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V5")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]
braycurtis_unrelated_V6ISMS <- braycurtis_unrelated[
  row.names(subset(phylo_meta, sample_type == "IS" & visit=="V6")),
  row.names(subset(phylo_meta, sample_type =="MS"))
]

# List of braycurtis matrices
braycurtis_matrices_related <- list(
  braycurtis_related_V4INMN,
  braycurtis_related_V5INMN,
  braycurtis_related_V6INMN,
  braycurtis_related_V3ISMS,
  braycurtis_related_V4ISMS,
  braycurtis_related_V5ISMS,
  braycurtis_related_V6ISMS
)

braycurtis_matrices_unrelated <- list(
  braycurtis_unrelated_V4INMN,
  braycurtis_unrelated_V5INMN,
  braycurtis_unrelated_V6INMN,
  braycurtis_unrelated_V3ISMS,
  braycurtis_unrelated_V4ISMS,
  braycurtis_unrelated_V5ISMS,
  braycurtis_unrelated_V6ISMS
)

#Now compare related & unrelated braycurtis:

calculate_mean_braycurtis <- 
  function(related_matrix_braycurtis, unrelated_matrix_braycurtis, subsample_size = NULL) {
    df <- data.frame(
      RowNames = row.names(related_matrix_braycurtis),
      MeanBrayCurtisRelated = numeric(nrow(related_matrix_braycurtis)),
      NumRelatedSamples = numeric(nrow(related_matrix_braycurtis)),
      MeanBrayCurtisUnrelated = numeric(nrow(related_matrix_braycurtis)),
      NumUnrelatedSamples = numeric(nrow(related_matrix_braycurtis))
    )
    
    for (i in seq_along(df$RowNames)) {
      row_name <- df$RowNames[i]
      
      # For related
      row_values_related <- as.numeric(related_matrix_braycurtis[row_name, , drop = FALSE])
      non_na_values_related <- row_values_related[!is.na(row_values_related)]
      df$NumRelatedSamples[i] <- length(non_na_values_related)
      df$MeanBrayCurtisRelated[i] <- mean(non_na_values_related, na.rm = TRUE)
      
      # For unrelated
      row_values_unrelated <- as.numeric(unrelated_matrix_braycurtis[row_name, , drop = FALSE])
      non_na_values_unrelated <- row_values_unrelated[!is.na(row_values_unrelated)]
      
      if (!is.null(subsample_size) && length(non_na_values_unrelated) > subsample_size) {
        set.seed(1988)  # Ensure reproducibility
        non_na_values_unrelated <- sample(non_na_values_unrelated, subsample_size)
      }
      
      df$NumUnrelatedSamples[i] <- length(non_na_values_unrelated)
      df$MeanBrayCurtisUnrelated[i] <- mean(non_na_values_unrelated, na.rm = TRUE)
    }
    
    return(df)
  }

subsample_size <- 20

braycurtis_related_unrelated <- Map(
  function(related_matrix_braycurtis, unrelated_matrix_braycurtis) 
    calculate_mean_braycurtis(related_matrix_braycurtis, unrelated_matrix_braycurtis, subsample_size),
  braycurtis_matrices_related, 
  braycurtis_matrices_unrelated
)

braycurtis_related_values <- unlist(lapply(braycurtis_related_unrelated, function(x) x$MeanBrayCurtisRelated))
braycurtis_unrelated_values <- unlist(lapply(braycurtis_related_unrelated, function(x) x$MeanBrayCurtisUnrelated))

wilcoxon_braycurtis_related_unrelated <- list()

for (i in seq_along(braycurtis_related_unrelated)) {
  related_values_braycurtis <- braycurtis_related_unrelated[[i]]$MeanBrayCurtisRelated
  unrelated_values_braycurtis <- braycurtis_related_unrelated[[i]]$MeanBrayCurtisUnrelated
  
  wilcoxon_braycurtis <- wilcox.test(related_values_braycurtis, unrelated_values_braycurtis, paired=TRUE)
  
  wilcoxon_braycurtis_related_unrelated[[i]] <- wilcoxon_braycurtis
}

print(wilcoxon_braycurtis_related_unrelated)

