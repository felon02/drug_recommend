# plot5_tcga
# 2023.9.1


library(survivalROC)

options(digits = 4)
options(scipen = 10) # cancel 科学计数表示

library(rstatix) # have function: adjust_pvalue, add_significance
library(ggrepel)

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(export)
library(ggsci)
# library(scales)
#
ls("package:ggsci")
colorJama = pal_jama()(7)
colorJco = pal_jco()(10)
# show_col(colorJama)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
specify_decimal_2 <- function(x, k=2) trimws(format(round(x, k), nsmall=k))
specify_decimal_3 <- function(x, k=3) trimws(format(round(x, k), nsmall=k))

library(scales)
# show_col(colorJco)

color_cancers
color_datasets <- structure(colorJco[c(1,2,4)], names = c("GDSC", "CCLE", "TCGA"))
color_datasets

# load data: -----------------------------------------------------

response_tcga_v2 = read_csv("data_temp/response_tcga_v2.csv")
pheno_tcga_v2 = read_csv("data_temp/pheno_tcga_v2.csv")
snv_tcga_v1 = read_csv("data_temp/snv_tcga_v1.csv")
tpm_tcga_v2 = read_csv("data_temp/tpm_tcga_v2.csv")

ytrue_ypred_tcga = read_csv("data_temp/ytrue_ypred_tcga_20230828.csv")
df_ypred_byDrug_tcga = read_csv("data_temp/df_ypred_byDrug_tcga_20230828.csv")


# data from gdsc
pheno_gdsc_v1 = read_csv("data_temp/pheno_gdsc_v1.csv")
snv_gdsc_wide_v1 = read_csv("data_temp/snv_gdsc_wide_v1.csv")
tpm_gdsc_v1 = read_csv("data_temp/tpm_gdsc_v1.csv")


# TCGA feature and data description ---------------------------------------

# clean pheno

col_select = c("Overall Survival (Months)", "Overall Survival Status",
               "Progress Free Survival (Months)", "Progression Free Status",
               "Disease Free (Months)", "Disease Free Status",
               "Patient ID", "TCGA PanCanAtlas Cancer Type Acronym",
               "Neoplasm Disease Stage American Joint Committee on Cancer Code",
               "Diagnosis Age", "Sample Type",
               "Radiation Therapy", "Sex", "Patient Weight", "Subtype", "MSIsensor Score")
col_rename = c("OS_Months", "OS_Status", "PFS_Months", "PFS_Status", "DFS_Months", "DFS_Status",
               "Patient_ID", "cancerType_TCGA", "tumor_stage_raw", "Age_at_Diagnosis", "Sample_Type",
               "Radiation_Therapy", "Sex", "Patient_Weight", "Subtype", "MSIsensor Score")

any(duplicated(pheno_tcga_v2$`Patient ID`)) # False

pheno_tcga_v3 <- pheno_tcga_v2 %>% select(all_of(col_select)) 
colnames(pheno_tcga_v3) = col_rename


### cancer type frequency
pheno_tcga_v3 %>% count(cancerType_TCGA)
df1 <- pheno_tcga_v3 %>% 
  count(cancerType_TCGA) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_TCGA = n / sum(n))
df1

setdiff(pheno_gdsc_v1$cancerType_TCGA, pheno_tcga_v3$cancerType_TCGA)

df2 <- pheno_gdsc_v1 %>% 
  count(cancerType_TCGA) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_GDSC = n / sum(n))
df2

df_pheno <- df1 %>% 
  full_join(df2, by = "cancerType_TCGA") %>% 
  mutate(prop_TCGA = ifelse(is.na(prop_TCGA), 0, prop_TCGA))
df_pheno

p1 = ggplot(data = df_pheno, aes(x = prop_GDSC, y = prop_TCGA, color = cancerType_TCGA)) +
  geom_point(position = "identity", size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = color_cancers) +
  scale_x_continuous(n.breaks = 6, limits = c(0, 0.2)) +
  scale_y_continuous(n.breaks = 6, limits = c(0, 0.5)) +
  labs(title="Cancer type proportion", x="GDSC", y = "TCGA", color = "Cancer") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "right",
        legend.key.size = unit(0.1,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p1

### age density
pheno_gdsc_v1 %>% head(2)
pheno_tcga_v3 %>% head(2)
df_age <- pheno_gdsc_v1 %>% select(sample_id = SANGER_MODEL_ID, Age = Age) %>% mutate(dataset = "GDSC") %>% 
  bind_rows(pheno_tcga_v3 %>% select(sample_id = Patient_ID, Age = Age_at_Diagnosis) %>%  mutate(dataset = "TCGA")) %>% 
  mutate(dataset = factor(dataset, levels = c("GDSC", "TCGA")))
df_age %>% head(2)
p2 <- ggplot(data = df_age, aes(x = Age, color = dataset)) +
  geom_histogram(aes(y=after_stat(density),  fill = dataset),
                 alpha=0.4, bins = 50, position="identity",
                 color = "white") +
  geom_density(alpha = 1, linewidth = 1) +
  scale_x_continuous(n.breaks = 10) +
  scale_color_manual(values = color_datasets) +
  scale_fill_manual(values = color_datasets) +
  labs(title = NULL, x = "Age", y = "Density") +
  theme_pubr(base_size = 6) +
  theme(legend.position = c(0.85, 0.9),
        legend.title = element_blank(),
        # legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p2

### sex bar

pheno_gdsc_v1 %>% head(2)
pheno_tcga_v3 %>% head(2)
df_sex <- pheno_gdsc_v1 %>% count(Sex) %>% mutate(prop = n/sum(n), dataset = "GDSC") %>% 
  bind_rows(pheno_tcga_v3 %>% count(Sex) %>% mutate(prop = n/sum(n), dataset = "TCGA")) %>% 
  mutate(dataset = factor(dataset, levels = c("GDSC", "TCGA")))
df_sex

p_fisher = fisher.test(matrix(df_sex$n, nrow = 2))
p_fisher$p.value
color_sex = structure(colorJco[c(6,9)], names = c("Male", "Female"))
color_sex

p3 <- ggplot(data = df_sex, aes(x = dataset, y = prop)) +
  geom_bar(aes(fill = Sex), stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = color_sex) +
  geom_bracket(xmin = "GDSC", xmax = "TCGA", y.position = 1.02, 
               label = round(p_fisher$p.value, 3), label.size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = NULL, x = "", y = "Proportion") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        # legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

p3

p123 <- ggarrange(p1, p2, p3, nrow = 1)
p123



### snv frequency
snv_tcga_v1 %>% head(2)
snv_gdsc_wide_v1 %>% head(2)

calP <- function(x){round(fisher.test(matrix(as.integer(x), nrow = 2))$p.value, 3)}
# merge and add P value
snv_gdsc_tcga <- snv_gdsc_wide_v1 %>% 
  mutate(dataset = "GDSC", N_sample = nrow(.)) %>% 
  bind_rows(snv_tcga_v1 %>% mutate(dataset = "TCGA", N_sample = nrow(.))) %>% 
  select(-`...1`) %>% 
  pivot_longer(cols = -c(dataset, N_sample), names_to = "gene_symbol", values_to = "mutation_status") %>% 
  group_by(dataset, gene_symbol) %>% 
  summarize(N_sample = mean(N_sample), N_wt = sum(mutation_status == 0), 
            N_mut = sum(mutation_status == 1), prop_mut = N_mut / N_sample) %>% 
  group_by(gene_symbol) %>% 
  mutate(fisher_pvalue = calP(c(N_wt, N_mut))) %>% 
  ungroup() 

snv_gdsc_tcga_v1 <- snv_gdsc_tcga %>% 
  inner_join(snv_gdsc_tcga %>% 
               distinct(gene_symbol, .keep_all = TRUE) %>% 
               mutate(adj_pvalue = p.adjust(fisher_pvalue, method = "fdr")) %>% 
               select(gene_symbol, adj_pvalue)
             )
  
snv_gdsc_tcga_v1 %>% head()
snv_gdsc_tcga_v1 %>% filter(gene_symbol == "ACSL3")


snv_gdsc_tcga_v1 %>% head(2)
PROP_gene = 15
PROP_gene_tcga = 10
df_snv <- snv_gdsc_tcga_v1 %>% 
  mutate(prop_mut = prop_mut * 100) %>% 
  pivot_wider(id_cols = gene_symbol, names_from = dataset, values_from = prop_mut) %>% 
  mutate(gene_label = ifelse((GDSC > PROP_gene) | (TCGA > PROP_gene_tcga), gene_symbol, NA))
df_snv %>% filter(!is.na(gene_label))

p4 = ggplot(data = df_snv, aes(x = GDSC, y = TCGA)) +
  geom_point(position = "identity", size = 0.5, color = colorJama[1]) +
  stat_cor(method = "pearson", cor.coef.name = "R", size = 1.5,
           label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label = gene_label), size = 1.5, nudge_x = 1, nudge_y = 1) +
  scale_x_continuous(n.breaks = 6, limits = c(0, 75)) +
  scale_y_continuous(n.breaks = 6, limits = c(0, 75)) +
  labs(title="Mutation frequency", x="GDSC", y = "TCGA") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p4

### gep: pca
tpm_gdsc_v1 %>% head(2)
tpm_tcga_v2 %>% head(2)
tpm_full <- tpm_gdsc_v1 %>% 
  mutate(id_sample = paste("GDSC", `...1`, sep = "_")) %>% 
  select(-`...1`) %>% 
  bind_rows(tpm_tcga_v2 %>% 
              mutate(id_sample = paste("TCGA", `...1`, sep = "_")) %>% 
              select(-`...1`))
tpm_full %>% colnames()
tpm_full_df <- tpm_full %>% select(-id_sample) %>% as.data.frame()
dim(tpm_full_df)
rownames(tpm_full_df) = tpm_full$id_sample
colnames(tpm_full_df)[1:3]

run_pca <- function(mat){
  #calculate principal components
  res <- prcomp(mat, scale = TRUE)
  #reverse the signs
  res$rotation <- -1*res$rotation
  # principal components scores for each state are stored in results$x. 
  # We will also multiply these scores by -1 to reverse the signs
  #reverse the signs of the scores
  res$x <- -1*res$x
  return(res$x)
}

# visualize PCA_TPM
pca_tpm <- run_pca(tpm_full_df)
pca_tpm_tb <- as_tibble(pca_tpm, rownames="id_sample") %>% 
  mutate(Database = str_sub(id_sample, start = 1, end = 4)) %>% 
  select(id_sample, Database, everything())
pca_tpm_tb %>% head(2)

p5 <- ggplot(aes(x = PC1, y = PC2, color = Database), data = pca_tpm_tb) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_datasets) +
  # geom_abline(slope = 0.6, intercept = -10) +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(0.9, 0.9))
p5
# conclusion: GDSC and TCGA are two independent datasets


### gep: highly variable genes
tpm_tcga_v2 %>% head(2)
tpm_gdsc_v1 %>% head(2)
tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(mean) %>% summary()

# check the minium
tpm_tcga_v2 %>% select(-`...1`) %>% map_dbl(function(x){abs(mean(x))}) %>% min() # 0.005039
tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(function(x){abs(mean(x))}) %>% min() # 0.04723
#
gene_cv_tcga <- tpm_tcga_v2 %>% select(-`...1`) %>% map_dbl(function(x){sd(x)/abs(mean(x)+10^(-5))})
gene_cv_tcga = sort(gene_cv_tcga, decreasing = T)
gene_cv_tcga

gene_cv_gdsc <- tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(function(x){sd(x)/abs(mean(x)+10^(-5))})
gene_cv_gdsc = sort(gene_cv_gdsc, decreasing = T)
gene_cv_gdsc

CV_label_cutoff = 10
df_gepCV <- tibble(gene = names(gene_cv_gdsc), GDSC = gene_cv_gdsc) %>% 
  inner_join(tibble(gene = names(gene_cv_tcga), TCGA = gene_cv_tcga)) %>% 
  mutate(gene_label = ifelse((GDSC > CV_label_cutoff) | (TCGA > CV_label_cutoff), gene, NA))
df_gepCV %>% head(2)

p6 = ggplot(data = df_gepCV, aes(x = GDSC, y = TCGA)) +
  geom_point(position = "identity", size = 0.5, color = colorJama[1]) +
  stat_cor(method = "pearson", cor.coef.name = "R", size = 1.5,
           label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_text(aes(label = gene_label), size = 1.5, nudge_x = 0.08, nudge_y = 0.08) +
  # scale_x_continuous(n.breaks = 6, limits = c(0, 75)) +
  # scale_y_continuous(n.breaks = 6, limits = c(0, 75)) +
  labs(title="GEP: coefficient of variance", x="GDSC", y = "TCGA") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p6

# combine
p456 <- ggarrange(p4, p5, p6, nrow = 1)
p456

# Response ----------------------------------------------------------------

### drug usage
library(RColorBrewer)

df_drug_usage <- ytrue_ypred_tcga %>% 
  filter(response_status != "No_DATA") %>% 
  count(drug_name) %>% 
  arrange(desc(n)) %>% 
  mutate(prop = round(n / sum(n) * 100, 1), 
         drug_label = paste(drug_name, paste("N=", as.character(n), sep=""), sep = "\n"))

color_drug_tcga = brewer.pal(n=12, name="Paired")
color_drug_tcga = c(color_drug_tcga, "magenta4")
names(color_drug_tcga) = df_drug_usage$drug_name
show_col(color_drug_tcga)

hsize = 2
?ggrepel::geom_text_repel
p7 = ggplot(df_drug_usage, aes(x = hsize, y = prop, fill = drug_name)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  ggrepel::geom_text_repel(aes(label = drug_label), position = position_stack(vjust = 0.5), 
                           color = "black", size = 2.5) +
  # geom_text(aes(label = drug_label), position = position_stack(vjust = 0.5), 
  #           color = "black", size = 2.5)+
  scale_fill_manual(values = color_drug_tcga) +
  xlim(0, hsize+0.5) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + # reduce legend area
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void(base_size = 6) +
  theme(legend.position = "none")
p7


### IC50 Predicted

library(pheatmap)
library(viridis)

### heatmap
ytrue_ypred_tcga %>% head(2)

df <- ytrue_ypred_tcga %>% 
  pivot_wider(id_cols = bcr_patient_barcode, names_from = drug_name, values_from = ypred)
df1 <- df %>% select(-bcr_patient_barcode) %>% as.matrix()
rownames(df1) = df$bcr_patient_barcode

# sample annotation
anno_row <- ytrue_ypred_tcga %>% 
  filter(response_status != "No_DATA") %>% 
  mutate(Best_OR = factor(response_status, levels = c('CR', 'PR', 'SD', 'PD'))) %>% 
  select(bcr_patient_barcode, Best_OR, cancerType_TCGA, drug_name, Age, Sex)
anno_row_v1 <- anno_row %>% 
  select(Cancer_type = cancerType_TCGA, Drug = drug_name, Best_OR = Best_OR) %>% as.data.frame()
rownames(anno_row_v1) = anno_row$bcr_patient_barcode
# double check
sort(rownames(anno_row_v1)) == sort(rownames(df1))

# set patient order by bestOR
pt_order <- anno_row %>% arrange(Best_OR) %>% pull(bcr_patient_barcode)
j1 = match(pt_order, rownames(df1))
j1
j2 = match(pt_order, rownames(anno_row_v1))
df2 <- df1[j1,]
anno_row_v2 <- anno_row_v1[j2,]
# double check
df2[1:3, 1:3]
pt_order[1:3]
anno_row_v2[1:3, 1:3]


# set color of annotation
color_cancers
anno_row_v1 %>% count(Best_OR) %>% arrange(desc(n))
colore_bestOR = structure(pal_lancet()(4)[c(1,4,3,2)], names = c('CR', 'PR', 'SD', 'PD'))

anno_colour = list(
  Drug = color_drug_tcga,
  Cancer_type = color_cancers,
  Best_OR = colore_bestOR)
anno_colour


pheatmap(t(df2),
         annotation_col = anno_row_v2,
         annotation_colors = anno_colour,
         # color = viridis(n=100),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cluster_cols = F,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         scale = "none",
         show_rownames = T,
         show_colnames = F,
         treeheight_row = 10,
         treeheight_col = 10,
         legend = TRUE,
         # legend_breaks = seq(-1, 1, 0.2),
         display_numbers = F,
         fontsize = 6,
         # fontsize_number = 0.8 * 6,
         main = "TCGA: predicted IC50"
)



### 
### gemcitabin, Cisplatin
ytrue_ypred_tcga %>% head(2)
ytrue_ypred_tcga %>% filter(response_status!="No_DATA") %>% count(drug_name) %>% arrange(desc(n))
target_drug = c("Gemcitabine", "Cisplatin")

df_theDrug <- ytrue_ypred_tcga %>% 
  filter(response_status != "No_DATA", drug_name %in% target_drug) %>% 
  mutate(Best_OR = factor(response_status, levels = c('CR', 'PR', 'SD', 'PD')),
         Response = ifelse(response_status %in% c('CR', 'PR'), "R", "NR"),
         Response = factor(Response, levels = c("R", "NR")),
         CR_status = ifelse(response_status %in% c('CR'), "CR", "non-CR"),
         CR_status = factor(CR_status, levels = c("non-CR", "CR")),
         PD_status = ifelse(response_status %in% c('CR', 'PR', 'SD'), "non-PD", "PD"),
         PD_status = factor(PD_status, levels = c("non-PD", "PD")))
df_theDrug
unique(df_theDrug$response_status)
df_theDrug %>% count(Best_OR)
df_theDrug %>% count(cancerType_TCGA, drug_name)

df_theDrug %>% 
  # dplyr::filter(cancerType_TCGA != "SARC") %>% 
  # dplyr::filter(cancerType_TCGA == "PAAD") %>% 
  group_by(drug_name) %>% 
  rstatix::wilcox_test(ypred ~ PD_status) %>% 
  ungroup()

df_theDrug %>% 
  dplyr::filter(cancerType_TCGA != "SARC") %>% 
  group_by(cancerType_TCGA, drug_name, PD_status) %>% 
  rstatix::get_summary_stats(ypred)

comparisons_mannul <- list(c("R", "NR"))
comparisons_mannul <- list(c("PD", "non-PD"))

color_pd <- structure(c("magenta4", colore_bestOR['PD']),
                      names = c("non-PD", "PD"))
color_pd

df_gemcitabin <- df_theDrug %>% filter(drug_name == "Gemcitabine")
unique(df_gemcitabin$PD_status)
p8 <- ggplot(data = df_gemcitabin, aes(x = PD_status, y = ypred, color = PD_status)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(outlier.shape=NA, width = 0.1) + # use outlier.shape to remove outliers
  geom_point(# aes(shape = cancerType_TCGA),
             position = position_jitter(width = 0.2), 
             size = 1
             ) +
  scale_color_manual(values = color_pd) +
  # scale_fill_manual(values = color_cancers) +
  stat_compare_means(size=2,
                     label.y = 0.5,
                     # method = 'wilcox.test',
                     method = 't.test',
                     comparisons = comparisons_mannul) +
  scale_y_continuous(labels = specify_decimal_2,
                     expand = expansion(mult = c(0, 0.1))) +
  guides(color = "none") +
  labs(title=NULL, x="Best Overall Response", y = "Predicted Score", shape = "Cancer Type") +
  theme_pubr() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 0, hjust = 0),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 8),
        legend.position = "right")
# p8 = p8 + facet_wrap(facets = "drug_name", nrow = 1, scales = "fixed")
p8


# survival analysis by median values in target drug -----------------------

library(survival)
library(survminer)

source("src_r/script_KMplot.R")


pheno_tcga_v3 %>% head(2)
df_theDrug %>% head(2)
pheno_targetDrug <- pheno_tcga_v3 %>% 
  mutate(PFS_status = as.integer(str_sub(PFS_Status, start = 1, end = 1)),
         OS_status = as.integer(str_sub(OS_Status, start = 1, end = 1))) %>% 
  select(-c(PFS_Status, OS_Status)) %>% 
  inner_join(df_theDrug %>% 
               pivot_wider(id_cols = bcr_patient_barcode, names_from = drug_name, values_from = ypred),
             by = c("Patient_ID" = "bcr_patient_barcode"))
pheno_targetDrug %>% count(PFS_status)
pheno_targetDrug %>% filter(!is.na(Gemcitabine)) %>% count(cancerType_TCGA)
pheno_targetDrug %>% filter(!is.na(Cisplatin)) %>% count(cancerType_TCGA)

pheno_gemcitabin <- pheno_targetDrug %>% 
  filter(!is.na(Gemcitabine), cancerType_TCGA == "PAAD") %>% 
  mutate(Predicted_Score = ifelse(Gemcitabine > median(Gemcitabine), "High", "Low"))
pheno_gemcitabin$Gemcitabine %>% summary()

#
pheno_gemcitabin %>% head(2)
pheno_gemcitabin %>% count(PFS_status)
fit1_pfs <- survfit(Surv(PFS_Months, PFS_status) ~ Predicted_Score, data = pheno_gemcitabin)
fit1_pfs_cox <- coxph(Surv(PFS_Months, PFS_status) ~ Predicted_Score, data = pheno_gemcitabin)
summary(fit1_pfs_cox)
surv_median(fit1_pfs)

fig_pfs_gemcitabin <- km_plot_arranged(fit1_pfs, 
                                       pheno_gemcitabin, 
                                       xlab='Months', ylab='PFS Probability',
                             pval_coord=c(3, 0.2), legend_coord=c(0.7, 0.835),
                             legend_title = 'Predicted Score',
                             title = "") +
  geom_vline(xintercept = 0.22, linetype = 2, color = colorJama[4]) +
  annotate("text", x = 0.2, y = 0.4, label = "HR: 0.572, 95%CI: 0.26-1.26", size = 2)
fig_pfs_gemcitabin

### survival difference at 0.5, 1 year
ComparisonSurv::Fixpoint.test(time=pheno_gemcitabin$PFS_Months, 
              status=pheno_gemcitabin$PFS_status, 
              group=pheno_gemcitabin$Predicted_Score == "High", 
              t0=6)
# 0.9474 (0.8470-1.0478) vs. 0.6111 (0.3859-0.8363)
# method statistic  pvalue
# 1         Naive     7.144 0.00752
# 2           Log     5.021 0.02504
# 3       Cloglog     4.259 0.03905
# 4 Arcsin-square     7.219 0.00721
# 5         Logit     4.612 0.03176

pheno_gemcitabin$PFS_Months
pheno_gemcitabin$PFS_status


list_roc = list()
list_bestCutoff = list()
for (t in c(6, 12)) {
  t_char = paste("t", t, sep = "")
  
  roc_surv <- survivalROC(Stime=pheno_gemcitabin$PFS_Months, 
                          status=pheno_gemcitabin$PFS_status, 
                          marker=pheno_gemcitabin$Gemcitabine,
                          entry = NULL, 
                          predict.time = t,
                          # span = 0.25*nobs^(-0.20))
                          cut.values = NULL, 
                          method = "KM",
                          lambda = NULL, 
                          span = NULL, 
                          window = "symmetric")
  
  df_roc <- tibble(FP = roc_surv$FP, 
                      TP = roc_surv$TP, 
                      sum_FP_TP = TP + 1 - FP,
                      cutoff = roc_surv$cut.values, survival_rate = roc_surv$Survival,
                      auc = roc_surv$AUC,
                      timepoint = t_char) 
  bestCutoff <- df_roc %>% slice_max(sum_FP_TP)
  
  list_roc[[t_char]] = df_roc
  list_bestCutoff[[t_char]] = bestCutoff
}

list_bestCutoff
df_roc_all = list_roc$t6 %>% bind_rows(list_roc$t12)

df_auc = df_roc_all %>% 
  distinct(auc, .keep_all = TRUE) %>% 
  mutate(x = c(0.6, 0.6),
         y = c(0.2, 0.1),
         auc_label = paste("AUC: ", round(auc, 3), sep = "")
         )
df_auc
df_bestCutoff = do.call("bind_rows", list_bestCutoff) %>% 
  mutate(cutoff_label = paste("(", round(FP, 2), ", ", round(TP, 2), ")", sep = ""))
df_bestCutoff


fig_roc_tcga
fig_roc_tcga = ggplot(aes(x = FP, y = TP, color = timepoint), data  = df_roc_all) +
  geom_path(linetype = 1, linewidth = 0.7) +
  geom_abline(slope=1,
              intercept = 0,
              colour = colorJco[3],
              linewidth = 0.5,
              linetype = 2)+
  geom_point(size = 1.5, data = df_bestCutoff) +
  geom_text(aes(label = cutoff_label), size = 3, data = df_bestCutoff) +
  geom_text(aes(x = x, y = y, label = auc_label, color = timepoint), size = 3, data = df_auc) +
  scale_color_jama() +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme_pubr(base_size = 8) +
  theme(legend.position = c(0.8, 0.5))

fig_roc_tcga