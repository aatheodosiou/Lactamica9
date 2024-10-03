# Load packages & data ----------------------------------------------------

library(tidyverse)
library(phyloseq)
library(here)
library(decontam)
library(knitr)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(lmerTest)
library(emmeans)

#Load additional functions (created by Dr Wouter de Steenhuijsen Piters - https://gitlab.com/wsteenhu/muis_trx)

source(here("src/load_dada.R"))
source(here("src/utils.R"))

# set paths
knitr::opts_knit$set(root.dir=".", aliases=c(h = "fig.height", w = "fig.width", ow = "out.width"))

# Rarefaction -------------------------------------------------------------

#Make an OTU table & use in rarefaction plot 

otu_table_decontam <- ps_final %>%
  otu_table() %>% as(., "matrix") %>% t

par(mar = c(1, 1, 1, 1))  # Adjust margins as needed

rarefaction5033 <- map2_dfr(vegan::rarecurve(otu_table_decontam, step = 10), 
                            rownames(otu_table_decontam), ~enframe(.x) %>% mutate(sample_id = .y)) %>%
  mutate(name = str_remove(name, "N") %>% as.numeric) %>%
  ggplot(aes(x = name, y = value, group = sample_id)) +
  geom_line(alpha = 0.15) +
  scale_y_log10() +
  geom_vline(xintercept = 5033, linetype="dashed") +
  labs(x = "Reads per sample", y = "Number of ASVs per sample") 

par(mar = c(5, 4, 4, 2))  # Reset margins to default

# Rarefy to even depth:choose based on rarefaction plot & min readcount sample

rarefy_even_depth_5033 <- 
  rarefy_even_depth(ps_final, 
                    rngseed =27, sample.size=5033)

#Now do alpha diversity for rarefied & non-rarefied

alpha_div_rarefied <- estimate_richness(rarefy_even_depth_5033, measures=c("Shannon", "Observed"))
alpha_div_unrarefied <- estimate_richness(ps_final, 
                                          measures=c("Shannon", "Observed"))

#Change column names of alpha-div tables, & then merge together to add back to phylo

colnames(alpha_div_rarefied) <- paste(colnames(alpha_div_rarefied), "rarefied",sep="_")
colnames(alpha_div_unrarefied) <- paste(colnames(alpha_div_unrarefied), "unrarefied",sep="_")

alpha_div_rarefied <- cbind(sample_id = rownames(alpha_div_rarefied), alpha_div_rarefied)
alpha_div_unrarefied <- cbind(sample_id = rownames(alpha_div_unrarefied), alpha_div_unrarefied)

alpha_div_all <- merge(alpha_div_unrarefied, alpha_div_rarefied, by="sample_id", all=T)

#How similar are alpha diversities for rarefied, non-rarefied filtered & non-rarefied unfiltered?

cor.test(alpha_div_all$Shannon_rarefied, alpha_div_all$Shannon_unrarefied, method="spearman") 
cor.test(alpha_div_all$Observed_rarefied, alpha_div_all$Observed_unrarefied, method="spearman") 

#Based on this, decision to proceed with rarefied

# Check normality & homoscedasticity --------------------------------------

lme_Shannon_visit.type_overall <- lmerTest::lmer(
  Shannon_rarefied ~ visit.type + (1|participant), 
  data = meta_alpha) #An overall dataset to test model

par(mfrow = c(1, 4)) #Make 4 plots in one
plot(lme_Shannon_visit.type_overall)
par(mar = c(1, 1, 1, 1))
qqnorm(ranef(lme_Shannon_visit.type_overall)$participant[, 1], 
       main = "qqnorm: Random effects of participant")
qqnorm(resid(lme_Shannon_visit.type_overall), main = "qqnorm: Residuals") 
hist(residuals(lme_Shannon_visit.type_overall), breaks=15, main = "Histogram of Residuals")
boxplot(residuals(lme_Shannon_visit.type_overall) ~ 
          meta_alpha$visit.type, main = "Boxplot of residuals")
par(mfrow = c(1, 1)) # Reset margins to default

#And again for observed...

# Plot alpha diversity ----------------------------------------------------

#Merge alpha diversity measures with phylo (after decontam, filt & even depth)

meta_alpha <- 
  merge(ps_final@sam_data,
        subset(alpha_div_rarefied, select = -sample_id), 
        by="row.names", all=TRUE) %>% 
  mutate(sample_type = fct_relevel(sample_type, "MO", "MS", "MN", "BM", "IS", "IN", "CS"))

meta_alpha$visit.type <- factor(meta_alpha$visit.type) #Change visit.type to factor

meta_alpha_mother <- subset(meta_alpha, participant_type=="mother")
meta_alpha_infant <- subset(meta_alpha, participant_type=="infant")

box_Shannon_sample <- ggplot(meta_alpha, aes(x=sample_type, y=Shannon_rarefied, color=sample_type)) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.5, aes(color = sample_type, fill = sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.3) + 
  scale_fill_manual(values = sample_type_colors) +
  scale_color_manual(values = sample_type_colors) + 
  xlab("Sample type") + ylab("Shannon index") + theme(legend.position = "none")

box_Observed_sample <- ggplot(meta_alpha, aes(x=sample_type, y=Observed_rarefied, color=sample_type)) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.5, aes(color = sample_type, fill = sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.3) + 
  scale_fill_manual(values = sample_type_colors) +
  scale_color_manual(values = sample_type_colors) + 
  xlab("Sample type") + ylab("Observed taxa") + theme(legend.position = "none")

box_Shannon_visit_sample <- 
  ggplot(subset(meta_alpha, sample_type != "CS"), aes(x=visit, y=Shannon_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Shannon_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, sample_type != "CS"), FUN = median), aes(group = sample_type, 
                                                                                                y = Shannon_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Shannon index") + theme(legend.position = "none")

box_Shannon_visit_sample_infant <- 
  ggplot(subset(meta_alpha, participant_type=="infant"), aes(x=visit, y=Shannon_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Shannon_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, participant_type=="infant"), FUN = median), aes(group = sample_type, 
                                                                                                       y = Shannon_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Shannon index") + theme(legend.position = "none")

box_Shannon_visit_sample_mother <- 
  ggplot(subset(meta_alpha, participant_type=="mother"), aes(x=visit, y=Shannon_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Shannon_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, participant_type=="mother"), FUN = median), aes(group = sample_type, 
                                                                                                       y = Shannon_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Shannon index") + theme(legend.position = "none")

box_Observed_visit_sample <- 
  ggplot(subset(meta_alpha, sample_type != "CS"), aes(x=visit, y=Observed_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Observed_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, sample_type != "CS"), FUN = median), aes(group = sample_type, 
                                                                                                y = Observed_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Observed taxa") + theme(legend.position = "none")

box_Observed_visit_sample_infant <- 
  ggplot(subset(meta_alpha, participant_type=="infant"), aes(x=visit, y=Observed_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Observed_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, participant_type=="infant"), FUN = median), aes(group = sample_type, 
                                                                                                       y = Observed_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Observed taxa") + theme(legend.position = "none")

box_Observed_visit_sample_mother <- 
  ggplot(subset(meta_alpha, participant_type=="mother"), aes(x=visit, y=Observed_rarefied, color=sample_type)) +
  geom_boxplot(outlier.colour = NA, alpha=0.8) +
  geom_jitter(size=2.4, alpha=0.5, 
              position = position_jitterdodge(dodge.width = 0.8)) + 
  geom_line(data = aggregate(Observed_rarefied ~ visit + sample_type, 
                             data = subset(meta_alpha, participant_type=="mother"), FUN = median), aes(group = sample_type, 
                                                                                                       y = Observed_rarefied, color=sample_type), 
            size = 1, alpha = 0.6, position = position_dodge(0.5)) +
  scale_color_manual(values = sample_type_colors) + 
  labs(x="Visit", y="Observed taxa") + theme(legend.position = "none")

plotgrid_alpha_diversity_rarefied <- plot_grid(
  box_Shannon_sample, box_Observed_sample, 
  ncol=2,labels="auto")

plotgrid_alpha_diversity_rarefied_visits <- plot_grid(
  box_Shannon_visit_sample, box_Observed_visit_sample, ncol=2, labels=c("e","f"))

# Correlation diversity vs readcount --------------------------------------

plot_Shannon_sample_readcount <- 
  ggplot(subset(meta_alpha, sample_type!="CS"), aes(x = filt_reads, y = Shannon_rarefied, color = sample_type)) +
  geom_point(alpha = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, aes(group = sample_type, color = sample_type), linetype = "dashed") +
  scale_color_manual(values = sample_type_colors) +
  labs(y = "Shannon index", x= "Filtered read count") + theme(legend.position = "none")

plot_Observed_sample_readcount <- 
  ggplot(subset(meta_alpha, sample_type!="CS"), aes(x = decontam_reads, y = Observed_rarefied, color = sample_type)) +
  geom_point(alpha = 0.8) + 
  scale_color_manual(values = sample_type_colors) +
  geom_smooth(method = "lm", se = FALSE, aes(group = sample_type, color = sample_type), linetype = "dashed") +
  labs(y = "Observed taxa", x= "Filtered read count") + theme(legend.position = "none")

plotgrid_cor_alpha_readcount_biomass <- plot_grid(
  plot_Shannon_sample_readcount, plot_Observed_sample_readcount)

cor.test(subset(meta_alpha, sample_type!="CS")$filt_reads, subset(meta_alpha, sample_type!="CS")$Shannon_rarefied, method="spearman")
cor.test(subset(meta_alpha, sample_type!="CS")$filt_reads, subset(meta_alpha, sample_type!="CS")$Observed_rarefied, method="spearman")

cor_Shannon_readcount_bysample <- subset(meta_alpha, sample_type!="CS") %>%
  group_by(sample_type) %>%
  summarize(correlation = cor(filt_reads, Shannon_rarefied, method="spearman"),
            p_value = cor.test(filt_reads, Shannon_rarefied, method="spearman")$p.value) %>% 
  mutate(fdr_p_value = p.adjust(p_value, method = "fdr")) %>% 
  mutate(fdr_sig = cut(fdr_p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                       labels = c("***", "**", "*", "^", "NS"))) 

cor_Observed_readcount_bysample <- subset(meta_alpha, sample_type!="CS") %>%
  group_by(sample_type) %>%
  summarize(correlation = cor(filt_reads, Observed_rarefied, method="spearman"),
            p_value = cor.test(filt_reads, Observed_rarefied, method="spearman")$p.value) %>% 
  mutate(fdr_p_value = p.adjust(p_value, method = "fdr")) %>% 
  mutate(fdr_sig = cut(fdr_p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                       labels = c("***", "**", "*", "^", "NS"))) 

# Beware missing samples! -------------------------------------------------

nrow(filter(meta_alpha, sample_type == "MO"))
filter(meta_alpha, sample_type=="MO") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31)

nrow(filter(meta_alpha, sample_type == "MS")) 
filter(meta_alpha, sample_type=="MS") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31) 

nrow(filter(meta_alpha, sample_type == "MN")) 
filter(meta_alpha, sample_type=="MN") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31) 

nrow(filter(meta_alpha, sample_type == "BM")) 
filter(meta_alpha, sample_type=="BM") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31) 

nrow(filter(meta_alpha, sample_type == "IS")) 
filter(meta_alpha, sample_type=="IS") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31) 

nrow(filter(meta_alpha, sample_type == "IN"))
filter(meta_alpha, sample_type=="IN") %>%
  group_by(participant) %>%
  summarise(unique_visits = toString(unique(visit))) %>% print(n=31)

# Shannon (filt, rarefied) stats analysis ------------------------------------------

# Is Shannon stat diff across visits & sample types?

#As dataset involves repeated measures (across visits), appropriate to
#fit linear mixed effects model (lmer), to use as input for ANOVA:
#Also, LME more robust to non-normal & heteroscedastic data

#F-statistic is a ratio of variances, specifically the ratio of the mean square 
#for the factor of interest to the mean square for the residuals (error); 
#whereas p-value is probability of observing an F-statistic as extreme 
#as the one calculated from the sample data, assuming null hypothesis true (=no effect).
#t-value in fixed effects output shows  how many standard errors the estimated 
#coefficient is away from zero. 

#Emmeans post-hoc test: estimates marginal means after fitting lmem
#i.e. post hoc prediction of means for different combinations of factors 
#while accounting for random effects in model.
#fdr = false discovery rate adjustment, to control expected proportion of falsely
#rejected null hypotheses in context of large number of comparisons

#First compare sample types (overall, not broken down by timepoint)
#Stats based on V3-6 (so that maternal & infant datasets comparable), and NOT CS

lme_Shannon_type <- lmerTest::lmer(
  Shannon_rarefied ~ sample_type + (1|participant), 
  data = subset(meta_alpha, sample_type!="CS", visit!="V1"))
anova(lme_Shannon_type)
summary(lme_Shannon_type)
fixef(lme_Shannon_type) 
emmeans(lme_Shannon_type, pairwise ~ sample_type, 
        adjust = "fdr", lmer.df = "satterthwaite")

# Shannon by meta variables ---------------------------------------

#Do univariable llmer to test each of these (collapse all visits, as not much difference):
#(Beware missing values MN & BM, not for lme)

lme_Shannon_sib_MO <- lmerTest::lmer(
  Shannon_rarefied ~ sib_5yr + (1|participant), 
  data = subset(meta_alpha, sample_type == "MO"))
anova(lme_Shannon_sib_MO)
summary(lme_Shannon_sib_MO)
fixef(lme_Shannon_sib_MO)
emmeans(lme_Shannon_sib_MO, pairwise ~ sib_5yr, 
        adjust = "fdr", lmer.df = "satterthwaite")

#Now by inoculation status... And so on for any other meta variables

lme_Shannon_inoc_MO <- lmerTest::lmer(
  Shannon_rarefied ~ inoculated + (1|participant), 
  data = subset(meta_alpha, sample_type == "MO" & visit != "V1"))
anova(lme_Shannon_inoc_MO)
summary(lme_Shannon_inoc_MO)
fixef(lme_Shannon_inoc_MO)
emmeans(lme_Shannon_inoc_MO, pairwise ~ inoculated, 
        adjust = "fdr", lmer.df = "satterthwaite")

#And now using visit.type...
#First, here's an example of how to decide if additive or interaction model more appropriate
#when comparing 2 or more variables (e.g. abx_mum_todate & visit.type)
#Do anova for additive vs interaction (likelihood ratio test p<0.05 suggests interaction useful), 
#and compare AIC (lower AIC value is better fit)

lme_Shannon_abx_visit.type_additive <- lmerTest::lmer(
  Shannon_rarefied ~ abx_mum_todate + visit.type + (1|participant), 
  data = subset(meta_alpha, sample_type != "CS", visit != "V1"))
lme_Shannon_abx_visit.type_interaction <- lmerTest::lmer(
  Shannon_rarefied ~ abx_mum_todate * visit.type + (1|participant), 
  data = subset(meta_alpha, sample_type != "CS", visit != "V1"))
anova(lme_Shannon_abx_visit.type_additive,lme_Shannon_abx_visit.type_interaction)
AIC(lme_Shannon_abx_visit.type_additive)
AIC(lme_Shannon_abx_visit.type_interaction)

#As comparing 2 variables, avoid pairwise across all results (too many comparisons!)...
#Instead, do emmeans without pairwise first (but with interaction), then pairs()...

emmeans_abx_visit.type_Shannon <- emmeans(lme_Shannon_abx_visit.type, ~ abx_mum_todate | visit.type, 
                                          adjust = "fdr", lmer.df = "satterthwaite") 

pairs(emmeans_abx_visit.type_Shannon, simple = "abx_mum_todate", adjust = "fdr")

#Maternal Nlac colonisation status (& so on for any other niches / variables to test)

lme_Shannon_nlac_MO <- lmerTest::lmer(
  Shannon_rarefied ~ m_nlac_col + (1|participant), 
  data = subset(meta_alpha, sample_type == "MO"))
anova(lme_Shannon_nlac_MO)
summary(lme_Shannon_nlac_MO)
fixef(lme_Shannon_nlac_MO)
emmeans(lme_Shannon_nlac_MO, pairwise ~ m_nlac_col, 
        adjust = "fdr", lmer.df = "satterthwaite")

#Now check V1 inoc vs uninoc = no baseline differences (Wilcoxon unpaired)
#(to make sure they weren't already different before inoculation!)

Shannon_wilcox_V1_MO_inocTF <-
  wilcox.test(Shannon_rarefied ~ inoculated, paired = FALSE, 
              data = subset(meta_alpha, sample_type == "MO" &
                              visit=="V1"))

# Now look at maternal Shannon V1 vs V3 and V4 for inoc (before & after)
#For 17 inoculated women:

Shannon_wilcox_V1V3_inoc_MO <-
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" & 
                              visit %in% c("V1","V3") & inoculated=="TRUE" & participant != "02"))
#NOTE: 02V1MO exists, but 02V3MO doesn't! sample missing from V3! So exclude 02.

Shannon_wilcox_V1V4_inoc_MO <- 
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" & 
                              visit %in% c("V1","V4") & inoculated == "TRUE"))

#So, MN sig diff V1 vs V4... Is this seen for uninoculated?
Shannon_wilcox_V1V4_uninoc_MN <-
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MN" &
                              visit %in% c("V1","V4") & inoculated=="FALSE" &
                              !(participant %in% c("01","11")))) 
#Note excluded participants with missing values, so only 5/7 included

#What about V3 & V4 inoc vs uninoc?
Shannon_wilcox_V3_inoc.uninoc_MN <-
  wilcox.test(Shannon_rarefied ~ inoculated, paired = FALSE, 
              data = subset(meta_alpha, sample_type == "MN" & visit == "V3"))
Shannon_wilcox_V4_inoc.uninoc_MN <-
  wilcox.test(Shannon_rarefied ~ inoculated, paired = FALSE, 
              data = subset(meta_alpha, sample_type == "MN" & visit == "V4"))

#Could pre/post inoculation effect actually be abx confounding?

#For 10 participants that had antibiotics between V1 & V3:

Shannon_wilcox_V1V3_abx_MO <-
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" & 
                              visit %in% c("V1","V3") & 
                              participant %in% c("05","09","14","16","20","25","31","11","17","24")))

#And now for V1V4, same participants as V1V3
Shannon_wilcox_V1V4_abx_MO <- #sig!!
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" &
                              visit %in% c("V1","V4") & 
                              participant %in% c("05","09","14","16","20","25","31","11","17","24"))) 

#To look at this another way, what about V4 abx vs no abx?
Shannon_wilcox_V4_MN_abx.noabx <-
  wilcox.test(Shannon_rarefied ~ abx_mum_todate, paired = FALSE, 
              data = subset(meta_alpha, sample_type == "MN" & visit=="V4"))
#NOT sig, suggesting apparent inoc diff at V4 ISN'T just an antibiotic effect

#As wilcox_V1V4_abx_MO sig, check that V1V4 & V1V3 NO abx not sig!
Shannon_wilcox_V1V4_noabx_MO <-
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" &
                              visit %in% c("V1","V4") & 
                              !(participant %in% c("02","05","09","14","16","20","25","31","11","17","24"))))

Shannon_wilcox_V1V3_noabx_MO <- #NS
  wilcox.test(Shannon_rarefied ~ visit, paired = TRUE, 
              data = subset(meta_alpha, sample_type == "MO" &
                              visit %in% c("V1","V3") & 
                              !(participant %in% c("02","05","09","14","16","20","25","31","11","17","24"))))

# Observed (filt, rarefied) stats analysis ------------------------------------------

#First compare sample types (overall, not broken down by timepoint)
#Stats based on V3-6 (so that maternal & infant datasets comparable), and NOT CS

lme_Observed_type <- lmerTest::lmer(
  Observed_rarefied ~ sample_type + (1|participant), 
  data = subset(meta_alpha, sample_type!="CS", visit!="V1"))
anova(lme_Observed_type)
summary(lme_Observed_type)
fixef(lme_Observed_type) 
emmeans(lme_Observed_type, pairwise ~ sample_type, 
        adjust = "fdr", lmer.df = "satterthwaite")

#Now by visit for each sample type (except BM & MN)

lme_Observed_MO_visit <- lmerTest::lmer(
  Observed_rarefied ~ visit + (1|participant), 
  data = filter(meta_alpha, sample_type=="MO"))
anova(lme_Observed_MO_visit)
summary(lme_Observed_MO_visit)
fixef(lme_Observed_MO_visit)
emmeans(lme_Observed_MO_visit, pairwise ~ visit, 
        adjust = "fdr", lmer.df = "satterthwaite") 

#And so on, as for Shannon above...

# Now some plots - Shannon ----------------------------------------------------------

#Only V3-6 (as no infant V1), and no CS (as only V6)

box_Shannon_inoc <- ggplot(subset(meta_alpha, visit !="V1" & sample_type !="CS"), 
                           aes(x=sample_type, y=Shannon_rarefied, color=inoculated)) +
  geom_boxplot(alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = c("TRUE"="deepskyblue", "FALSE"="pink2")) +  # Scale color manual
  labs(x = "Maternal inoculation status", y = "Shannon index") +
  theme(legend.position = "none")

box_Shannon_sib <- ggplot(subset(meta_alpha, visit !="V1" & sample_type !="CS"), 
                          aes(x=sample_type, y=Shannon_rarefied, color=sib_5yr)) +
  geom_boxplot(alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = c("TRUE"="goldenrod2", "FALSE"="dimgrey")) +  # Scale color manual
  labs(x = "Co-habiting siblings under 5", y = "Shannon index") +
  theme(legend.position = "none")

box_Shannon_abx <- ggplot(subset(meta_alpha, visit !="V1" & sample_type !="CS"), 
                          aes(x=sample_type, y=Shannon_rarefied, color=abx_mum_todate)) +
  geom_boxplot(alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = c("TRUE"="red2", "FALSE"="forestgreen")) +  # Scale color manual
  labs(x = "Maternal antibiotics to date", y = "Shannon index") +
  theme(legend.position = "bottom")

box_Shannon_nlac <- ggplot(subset(meta_alpha, visit !="V1" & sample_type !="CS"), 
                           aes(x=sample_type, y=Shannon_rarefied, color=m_nlac_col)) +
  geom_boxplot(alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = c("TRUE"="darkblue", "FALSE"="grey")) +  # Scale color manual
  labs(x = "Maternal N. lacatmica colonisation (Y92 or natural)", y = "Shannon index") +
  theme(legend.position = "bottom")

#Can also plot each metavariable by visit & niche, e.g. 

box_Shannon_abx_visit_MO <- ggplot(subset(meta_alpha, sample_type == "MO"), 
                                   aes(x=visit, y=Shannon_rarefied, color=abx_mum_todate)) +
  geom_boxplot(outlier.colour = NA, alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  labs(x = "Maternal antibiotics to date", y = "Shannon index (MO)") +
  theme(legend.position = "bottom")

#Plot before/after maternal abx or inoculation for maternal niches

box_Shannon_abx_V1V3V4 <- ggplot(subset(meta_alpha, 
                                        visit %in% c("V1","V3","V4") & sample_type %in% c("MO","MS","MN") &
                                          participant %in% c("05","09","14","16","20","25","31","11","17","24")), 
                                 aes(x = sample_type, y = Shannon_rarefied, color=visit)) +
  geom_boxplot(outlier.colour = NA, alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = visit_colors) + 
  labs(x = "Before and after maternal antibiotics", y = "Shannon index (mother)") +
  theme(legend.position = "bottom")

box_Shannon_inoc_V1V3V4 <- ggplot(subset(meta_alpha, 
                                         visit %in% c("V1","V3","V4") & sample_type %in% c("MO","MS","MN") &
                                           inoculated == TRUE), 
                                  aes(x = sample_type, y = Shannon_rarefied, color=visit)) +
  geom_boxplot(outlier.colour = NA, alpha=1) +
  geom_jitter(size=2.4, alpha=0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_color_manual(values = visit_colors) + 
  labs(x = "Before and after maternal inoculation", y = "Shannon index (mother)") +
  theme(legend.position = "bottom")

#What about V1vV4 for uninoculated? And V1inoc vs V1uninoc? And V4inoc vs V4 uninoc?

box_Shannon_uninoc_vs_inocV1V4MN <- meta_alpha %>% 
  subset(., visit %in% c("V1","V4") & sample_type == "MN") %>% 
  ggplot(aes(x = visit.type, y = Shannon_rarefied, color = inoculated)) +
  geom_boxplot(outlier.colour = NA, alpha = 1) +
  geom_jitter(aes(x = visit.type, y = Shannon_rarefied), 
              size = 2.4, alpha = 0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_y_continuous(limits = c(0, 5)) +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("TRUE" = "deepskyblue", "FALSE" = "pink2")) +  
  labs(x = "Timepoint.niche", y = "Shannon index") +
  theme(legend.position = "none")

box_Shannon_abx_vs_noabxV4MN <- meta_alpha %>% 
  subset(., visit %in% c("V1","V4") & sample_type == "MN") %>% 
  ggplot(aes(x = visit.type, y = Shannon_rarefied, color = abx_mum_todate)) +
  geom_boxplot(outlier.colour = NA, alpha = 1) +
  geom_jitter(aes(x = visit.type, y = Shannon_rarefied), 
              size = 2.4, alpha = 0.7, 
              position = position_jitterdodge(dodge.width = 0.8)) +  
  scale_y_continuous(limits = c(0, 5)) +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("TRUE" = "red2", "FALSE" = "forestgreen")) +  
  labs(x = "Timepoint.niche", y = "Shannon index") +
  theme(legend.position = "none")

plotgrid_Shannon_meta1 <- plot_grid(box_Shannon_sib, box_Shannon_inoc, labels=c("a","b"))
plotgrid_shannon_meta2 <- plot_grid(box_Shannon_uninoc_vs_inocV1V4MN, box_Shannon_abx_vs_noabxV4MN,
                                    nrow=2, labels=c("c","d"))
plotgrid_Shannon_meta <- plot_grid(plotgrid_Shannon_meta1, plotgrid_shannon_meta2, rel_widths = c(0.75,0.25))

#And can do same for Observed...