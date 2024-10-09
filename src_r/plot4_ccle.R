# plot4_ccle
# 2023.8.25

options(digits = 4)
options(scipen = 10) # cancel 科学计数表示

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


# data --------------------------------------------------------------------

# ccle: data from m4
drug_used_info = read_csv("data_temp/out_ccle/v20230812/drug_used_info.csv")
y_ccle = read_csv("data_temp/out_ccle/v20230812/y_ccle.csv")
y_ccle_anno = read_csv("data_temp/out_ccle/v20230812/y_ccle_anno.csv")
pheno_ccle_v2 = read_csv("data_temp/out_ccle/v20230812/pheno_ccle_v2.csv")
pheno_ccle_v2_dummy = read_csv("data_temp/out_ccle/v20230812/pheno_ccle_v2_dummy.csv")
snv_v1 = read_csv("data_temp/out_ccle/v20230812/snv_v1.csv")
tpm_v1 = read_csv("data_temp/out_ccle/v20230812/tpm_v1.csv")

y_pred_v1_anno = read_csv("data_temp/out_ccle/v20230812/y_pred_v1_anno.csv")
y_pred_v1_pivot = read_csv("data_temp/out_ccle/v20230812/y_pred_v1_pivot.csv")


# data from m2
pheno_gdsc_v1 = read_csv("data_temp/pheno_gdsc_v1.csv")
snv_gdsc_wide_v1 = read_csv("data_temp/snv_gdsc_wide_v1.csv")
tpm_gdsc_v1 = read_csv("data_temp/tpm_gdsc_v1.csv")

# data from plot2
snv_freq_gdscANDccle <- read_csv("data_temp/freq2_df_v1.csv")
snv_freq_gdscANDccle %>% head()
snv_freq_gdscANDccle %>% dim()

color_cancers
length(color_cancers)
color_datasets <- structure(colorJco[1:2], names = c("GDSC", "CCLE"))
color_datasets


# CCLE feature and data description ---------------------------------------

### cancer type frequency
pheno_gdsc_v1 %>% count(cancerType_TCGA)
df1 <- pheno_gdsc_v1 %>% 
  count(cancerType_TCGA) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_GDSC = n / sum(n))
df1
pheno_ccle_v2 %>% head(2)

setdiff(pheno_gdsc_v1$cancerType_TCGA, pheno_ccle_v2$cancerType_TCGA)
length(unique(pheno_gdsc_v1$cancerType_TCGA))
length(unique(pheno_ccle_v2$cancerType_TCGA))

df2 <- pheno_ccle_v2 %>% 
  count(cancerType_TCGA) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_CCLE = n / sum(n))
df_pheno <- df1 %>% 
  full_join(df2, by = "cancerType_TCGA") %>% 
  mutate(prop_CCLE = ifelse(is.na(prop_CCLE), 0, prop_CCLE))
df_pheno

p1 = ggplot(data = df_pheno, aes(x = prop_GDSC, y = prop_CCLE, color = cancerType_TCGA)) +
  geom_point(position = "identity", size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = color_cancers) +
  scale_x_continuous(n.breaks = 6, limits = c(0, 0.16)) +
  scale_y_continuous(n.breaks = 6, limits = c(0, 0.16)) +
  labs(title="Cancer type proportion", x="GDSC(N=542)", y = "CCLE(N=93)", color = "Cancer") +
  theme_pubr(base_size = 6) +
  theme(legend.position = "right",
        legend.key.size = unit(0.1,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p1

### age density
pheno_gdsc_v1 %>% head(2)
pheno_ccle_v2 %>% head(2)
df_age <- pheno_gdsc_v1 %>% mutate(dataset = "GDSC") %>% 
  bind_rows(pheno_ccle_v2 %>% mutate(dataset = "CCLE")) %>% 
  mutate(dataset = factor(dataset, levels = c("GDSC", "CCLE")))
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
library(scales)
show_col(colorJco)
pheno_gdsc_v1 %>% head(2)
pheno_ccle_v2 %>% head(2)
df_sex <- pheno_gdsc_v1 %>% count(Sex) %>% mutate(prop = n/sum(n), dataset = "GDSC") %>% 
  bind_rows(pheno_ccle_v2 %>% count(Sex) %>% mutate(prop = n/sum(n), dataset = "CCLE")) %>% 
  mutate(dataset = factor(dataset, levels = c("GDSC", "CCLE")))
df_sex

p_fisher = fisher.test(matrix(df_sex$n, nrow = 2))
p_fisher$p.value
color_sex = structure(colorJco[c(6,9)], names = c("Male", "Female"))
color_sex

p3 <- ggplot(data = df_sex, aes(x = dataset, y = prop)) +
  geom_bar(aes(fill = Sex), stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = color_sex) +
  geom_bracket(xmin = "GDSC", xmax = "CCLE", y.position = 1.02, 
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
snv_freq_gdscANDccle %>% head(2)
PROP_gene = 20
df_snv <- snv_freq_gdscANDccle %>% 
  mutate(gene_label = ifelse((freq_GDSC > PROP_gene) | (freq_CCLE > PROP_gene-5), gene_symbol, NA))
df_snv %>% filter(!is.na(gene_label))
p4 = ggplot(data = df_snv, aes(x = freq_GDSC, y = freq_CCLE)) +
  geom_point(position = "identity", size = 0.3, color = colorJama[1]) +
  stat_cor(method = "pearson", cor.coef.name = "R", size = 1.5,
           label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label = gene_label), size = 1.5, nudge_x = 1, nudge_y = 1) +
  scale_x_continuous(n.breaks = 6, limits = c(0, 75)) +
  scale_y_continuous(n.breaks = 6, limits = c(0, 75)) +
  labs(title="Mutation frequency", x="GDSC", y = "CCLE") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p4

### gep: pca
tpm_full <- tpm_gdsc_v1 %>% 
  mutate(id_sample = paste("GDSC", `...1`, sep = "_")) %>% 
  select(id_sample, everything(), -`...1`) %>% 
  bind_rows(tpm_v1 %>% 
              mutate(id_sample = paste("CCLE", `level_0`, sep = "_")) %>% 
              select(id_sample, everything(), -`level_0`))
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
p5 <- ggplot(aes(x = PC1, y = PC2, color = Database), data = pca_tpm_tb) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_datasets) +
  # geom_abline(slope = 0.6, intercept = -10) +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(0.9, 0.9))
p5


### gep: highly variable genes
tpm_v1 %>% head(2)
tpm_gdsc_v1 %>% head(2)
tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(mean) %>% summary()

# check the minium
tpm_v1 %>% select(-`level_0`) %>% map_dbl(function(x){abs(mean(x))}) %>% min()
tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(function(x){abs(mean(x))}) %>% min()
#
gene_cv_ccle <- tpm_v1 %>% select(-`level_0`) %>% map_dbl(function(x){sd(x)/abs(mean(x)+10^(-5))})
gene_cv_ccle = sort(gene_cv_ccle, decreasing = T)
gene_cv_ccle

gene_cv_gdsc <- tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(function(x){sd(x)/abs(mean(x)+10^(-5))})
gene_cv_gdsc = sort(gene_cv_gdsc, decreasing = T)
gene_cv_gdsc

CV_label_cutoff = 15
df_gepCV <- tibble(gene = names(gene_cv_gdsc), GDSC = gene_cv_gdsc) %>% 
 inner_join(tibble(gene = names(gene_cv_ccle), CCLE = gene_cv_ccle)) %>% 
  mutate(gene_label = ifelse((GDSC > CV_label_cutoff) | (CCLE > CV_label_cutoff), gene, NA))
df_gepCV %>% head(2)

p6 = ggplot(data = df_gepCV, aes(x = GDSC, y = CCLE)) +
  geom_point(position = "identity", size = 0.3, color = colorJama[1]) +
  stat_cor(method = "pearson", cor.coef.name = "R", size = 1.5,
           label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_text(aes(label = gene_label), size = 1.5, nudge_x = 0.08, nudge_y = 0.08) +
  # scale_x_continuous(n.breaks = 6, limits = c(0, 75)) +
  # scale_y_continuous(n.breaks = 6, limits = c(0, 75)) +
  labs(title="GEP: coefficient of variance", x="GDSC", y = "CCLE") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
p6

# combine
p456 <- ggarrange(p4, p5, p6, nrow = 1)
p456


### MFI log2FC
library(pheatmap)

# note: Docetaxel with LOWER log2FC of MFI, "BRD-K30577245-341-01-9"
# Docetaxel with higher log2FC of MFI, "BRD-K30577245-001-04-3"
y_ccle %>% head(2)
y_ccle_anno %>% head(2)
unique(y_ccle_anno$drug_name_revised)
y_ccle_anno %>% filter(drug_name_revised == "Docetaxel") %>% count(broad_id)
pheno_ccle_v2 %>% head(2)

df_hp <- y_ccle %>% 
  pivot_wider(id_cols = DepMap_ID, names_from = drug_name_revised, values_from = log2FC)
df_hp_v1 <- df_hp %>% 
  inner_join(pheno_ccle_v2, by = c("DepMap_ID"="ModelID")) %>% 
  mutate(cancerType_TCGA = factor(cancerType_TCGA, levels = names(color_cancers)))
df_hp_v1 %>% head(2)
df_hp_mat <- df_hp_v1 %>% 
  select(-DepMap_ID, -Age, -Sex, -cancerType_TCGA) %>% 
  as.data.frame()
rownames(df_hp_mat) = df_hp_v1$DepMap_ID

# sample annotation
df_hp_v1 %>% select(Sex, cancerType_TCGA) %>% as.data.frame() -> df_hp_v1_anno_row
rownames(df_hp_v1_anno_row) = df_hp_v1$DepMap_ID

# set color of annotation
anno_colour = list(
  cancerType_TCGA = color_cancers,
  Sex = color_sex)

# HEATMAP
pheatmap(t(df_hp_mat),
         annotation_col = df_hp_v1_anno_row,
         annotation_colors = anno_colour,
         cluster_cols = F,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         treeheight_row = 10,
         treeheight_col = 10,
         scale = "none",
         show_colnames = FALSE,
         fontsize = 5,
         main = "CCLE: logFC of MFI"
)


### mean IMF FC per cancertype
df_hp %>% head(2)
drugs_ccle = setdiff(colnames(df_hp), c("DepMap_ID"))
drugs_ccle

df_meanMFI_byCancers <- df_hp_v1 %>% 
  select(cancerType_TCGA, all_of(drugs_ccle)) %>% 
  pivot_longer(cols = all_of(drugs_ccle), names_to = "drugs", values_to = "MFI") %>% 
  group_by(cancerType_TCGA, drugs) %>% 
  summarize(meanMFI = mean(MFI))
df_meanMFI_byCancers %>% head(2)

drug_order_ccle <- df_meanMFI_byCancers %>% group_by(drugs) %>% summarize(meanMFI = mean(meanMFI)) %>% 
  arrange(meanMFI) %>% pull(drugs)
drug_order_ccle

df_meanMFI_byCancers$drugs = factor(df_meanMFI_byCancers$drugs, levels = rev(drug_order_ccle))

fig_meanMFI <- ggplot(data = df_meanMFI_byCancers, aes(x = drugs, y = meanMFI, color = cancerType_TCGA)) +
  geom_boxplot(color = "gray20", fill = "transparent", size = 0.3, outlier.size = 0) +
  geom_point(size = 0.5, position = position_jitter(width = 0.2)) +
  scale_color_manual(values = color_cancers) +
  labs(y = "Mean MFI", x = NULL, color = NULL) +
  coord_flip() +
  theme_bw(base_size = 6) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        # legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

fig_meanMFI



# CCLE data prediction ----------------------------------------------------

y_pred_v1_anno %>% head(2)
pheno_ccle_v2 %>% head(2)

ytrue_ypred_ccle <- y_pred_v1_anno %>% 
  rename(y_true = log2FC, y_pred = ypred, drugs = drug_name_revised)
ytrue_ypred_ccle %>% head(2)

### re-order drugs and cancers by correlation coef
cal_pval_of_cortest <- function(x,y, method = "pearson"){
  if (length(x) < 3) {
    return(NA)
  }
  z = cor.test(x, y, method = method)
  z = z$p.value
  return(z)
}

corr_by_cancer <- ytrue_ypred_ccle %>% 
  group_by(cancerType_TCGA) %>% 
  summarize(coef_corr_pearson = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr_pearson = cal_pval_of_cortest(y_true, y_pred),
            coef_corr_spearman = cor(y_true, y_pred, method = "spearman"),
            pvalue_corr_spearman = cal_pval_of_cortest(y_true, y_pred, method = "spearman")
            ) %>% 
  arrange(desc(coef_corr_pearson))
corr_by_cancer

corr_by_drug <- ytrue_ypred_ccle %>% 
  group_by(drugs) %>% 
  summarize(coef_corr_pearson = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr_pearson = cal_pval_of_cortest(y_true, y_pred),
            coef_corr_spearman = cor(y_true, y_pred, method = "spearman"),
            pvalue_corr_spearman = cal_pval_of_cortest(y_true, y_pred, method = "spearman")
  ) %>% 
  arrange(desc(coef_corr_pearson))
corr_by_drug
# correlation coef of spearman is generally lower than pearson, so i use pearson directly

###
j = match(corr_by_cancer$cancerType_TCGA, names(color_cancers))
j
color_cancers_ccle = color_cancers[j]
color_cancers_ccle
length(color_cancers_ccle)
drug_order_ccle = corr_by_drug$drugs
drug_order_ccle

df_ytrue_ypred_ccle_v1 <- ytrue_ypred_ccle %>% 
  mutate(cancer_type = factor(cancerType_TCGA, levels = names(color_cancers_ccle)),
         drugs = factor(drugs, levels = drug_order_ccle))
df_ytrue_ypred_ccle_v1 %>% head(2)
df_ytrue_ypred_ccle_v1$y_pred %>% min()
df_ytrue_ypred_ccle_v1$y_pred %>% max()
df_ytrue_ypred_ccle_v1$y_true %>% min()
df_ytrue_ypred_ccle_v1$y_true %>% max()

scatter_plot_ccle <- function(df, drugCol=NULL, cancertypeCol=NULL){
  p = ggplot(data = df, aes(x = y_true, y = y_pred)) +
    geom_point(position = "identity", color = colorJama[1], size = 0.4) +
    geom_smooth(method='lm', se = FALSE, color=colorJama[4], linetype=1, size=0.6, fullrange = FALSE) +
    stat_cor(method = "pearson", cor.coef.name = "R", size = 1.0,
             label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(n.breaks = 6) +
    labs(title=NULL, x="MFI", y = "Score predicted")
    
  p_revise = p + ylim(-1, 1) +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6))
  p_small =  p +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(size = 4.5),
      axis.text.y = element_text(size = 4.5))
  
  p_drug = facet(p_small, facet.by = drugCol, nrow = 3, scales = "free") 
  p_cancer = facet(p_small, facet.by = cancertypeCol, nrow = 2, scales = "free") + ylim(-1,1)
  return(list(p_revise, p_drug, p_cancer))
}

fig_mfi <- scatter_plot_ccle(df=df_ytrue_ypred_ccle_v1, drugCol="drugs", cancertypeCol="cancer_type")
fig_mfi[[1]]
fig_mfi[[2]]
fig_mfi[[3]]
# 
fig_comb1_ccle <- ggarrange(ggarrange(fig_mfi[[1]], fig_mfi[[3]], nrow = 1, widths = c(1, 2.3)),
                            ggarrange(fig_mfi[[2]], NULL, nrow = 1, widths = c(1.6, 1)),
                            nrow = 2, heights = c(2,3))
fig_comb1_ccle



### heatmap of correlation per cancer per drug
df_ytrue_ypred_ccle_v1 %>% head(2)
df_ytrue_ypred_ccle_v1 %>% filter(cancer_type == "BRCA", drugs == "Dasatinib")


cal_cor_coef <- function(x,y, method="pearson"){
  z = ifelse(length(x) < 3, NA, cor(x, y, method = method))
  return(z)
}
corr_by_cancerANDdrug <- df_ytrue_ypred_ccle_v1 %>% 
  group_by(cancer_type, drugs) %>% 
  summarize(coef_corr = cal_cor_coef(y_true, y_pred, method = "pearson"),
            pvalue_corr = cal_pval_of_cortest(y_true, y_pred),
            coef_corr_spearman = cal_cor_coef(y_true, y_pred, method = "spearman"),
            pvalue_corr_spearman = cal_pval_of_cortest(y_true, y_pred, method = "spearman"))
corr_by_cancerANDdrug %>% head(2)
corr_by_cancerANDdrug_v1 <- corr_by_cancerANDdrug %>% filter(!is.na(coef_corr))
corr_by_cancerANDdrug_v1 %>% head()

library(pheatmap)
library(RColorBrewer)

df1 <- corr_by_cancerANDdrug_v1 %>% 
  pivot_wider(id_cols = drugs, names_from = cancer_type, values_from = coef_corr)

mat1 <- df1 %>% select(-drugs) %>% as.matrix()
rownames(mat1) = df1$drugs
colnames(mat1) = colnames(df1)[-1]
min(mat1)
max(mat1)


pheatmap(mat1,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100),
         cluster_cols = T,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         scale = "none",
         show_rownames = T,
         show_colnames = T,
         treeheight_row = 10,
         treeheight_col = 10,
         legend = TRUE,
         legend_breaks = seq(-1, 1, 0.2),
         display_numbers = TRUE,
         fontsize = 6,
         fontsize_number = 0.8 * 6,
         main = NA
)
