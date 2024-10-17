# plot3_gdsc_performance
# 2023.8.25

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


# 1. ytrue and ypred in training and testing dataset -------------------------

pheno_gdsc_v1 = read_csv("data_temp/pheno_gdsc_v1.csv")
df_ytrue_ypred = read_csv("data_temp/out_model/df_ytrue_ypred_20230828.csv")
df_ytrue_ypred %>% head(2)
df_ytrue_ypred %>% count(dataset)
pheno_gdsc_v1

df_ytrue_ypred_v1 <- df_ytrue_ypred %>% 
  inner_join(pheno_gdsc_v1) %>% 
  rename(drugs = level_1)


### re-order drugs and cancers by correlation coef
df_test <- df_ytrue_ypred_v1 %>% filter(dataset == "Testing")
df_test %>% head(2)
cal_pval_of_cortest <- function(x,y, method = "pearson"){
  z = cor.test(x, y, method = method)
  z = z$p.value
  return(z)
}
corr_by_cancer <- df_test %>% 
  group_by(cancerType_TCGA) %>% 
  summarize(coef_corr = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr = cal_pval_of_cortest(y_true, y_pred)) %>% 
  arrange(desc(coef_corr))
corr_by_cancer

corr_by_drug <- df_test %>% 
  group_by(drugs) %>% 
  summarize(coef_corr = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr = cal_pval_of_cortest(y_true, y_pred)) %>% 
  arrange(desc(coef_corr))
corr_by_drug

corr_by_cancerANDdrug <- df_test %>% 
  group_by(cancerType_TCGA, drugs) %>% 
  summarize(coef_corr = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr = cal_pval_of_cortest(y_true, y_pred)) %>% 
  arrange(desc(coef_corr))
corr_by_cancerANDdrug %>% head(2)


j = match(corr_by_cancer$cancerType_TCGA, names(color_cancers))
color_cancers_v1 = color_cancers[j]
color_cancers_v1
drug_order_v1 = corr_by_drug$drugs


df_ytrue_ypred_v1 <- df_ytrue_ypred_v1 %>% 
  mutate(cancer_type = factor(cancerType_TCGA, levels = names(color_cancers_v1)),
         drugs = factor(drugs, levels = drug_order_v1),
         dataset = factor(dataset, levels = c("Training", "Testing")))
df_ytrue_ypred_v1 %>% head(2)


scatter_plot <- function(df, datasetName, title, drugCol=NULL, cancertypeCol=NULL){
  df1 <- df %>% filter(dataset == datasetName)
  p = ggplot(data = df1, aes(x = y_true, y = y_pred)) +
    geom_point(position = "identity", color = colorJama[1], size = 0.1) +
    geom_smooth(method='lm', se = F, color=colorJama[4], linetype=1, size=0.6) +
    stat_cor(method = "pearson", cor.coef.name = "R", size = 1.5,
             label.x.npc = 0.02, label.y.npc = 0.98, digits = 3) +
    scale_x_continuous(n.breaks = 6) + 
    scale_y_continuous(n.breaks = 7) +
    labs(title=paste("Correlation: ", title, sep = ""), x="IC50", y = "Score predicted") + 
    theme_pubr(base_size = 6) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          # legend.key.size = unit(0.3,"cm"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))

  p_drug = facet(p, facet.by = drugCol, nrow = 3, scales = "free")
  
  p_cancer = facet(p, facet.by = cancertypeCol, nrow = 3, scales = "free")
  return(list(p, p_drug, p_cancer))
}


df_ytrue_ypred_v1 %>% count(dataset)
fig_trainset <- scatter_plot(df=df_ytrue_ypred_v1, datasetName="Training", title="Training-set", 
                             drugCol="drugs", cancertypeCol="cancer_type")
fig_trainset[[2]]

fig_testset <- scatter_plot(df=df_ytrue_ypred_v1, datasetName="Testing", title="Testing-set", 
                            drugCol="drugs", cancertypeCol="cancer_type")
fig_testset

# 
fig_overall = ggarrange(fig_trainset[[1]], fig_testset[[1]], nrow = 1)
fig_drugs = ggarrange(fig_trainset[[2]], fig_testset[[2]], nrow = 2)
fig_cancers = ggarrange(fig_trainset[[3]], fig_testset[[3]], nrow = 2)

p_trainset_cancer_2rows <- facet(fig_trainset[[1]], facet.by = "cancer_type", nrow = 2, scales = "free") +
  labs(title = "Correlation by cancer: Training-set")
fig_trainset_cancerANDdrug <- ggarrange(p_trainset_cancer_2rows, fig_trainset[[2]], nrow = 2, heights = c(2,3))
fig_testset_cancerANDdrug <- ggarrange(fig_testset[[3]], fig_testset[[2]], nrow = 2)
fig_trainset_cancerANDdrug

fig_comb1 <- ggarrange(ggarrange(ggarrange(fig_trainset[[1]], fig_testset[[1]], nrow = 2),
                                 fig_testset[[3]],
                                 nrow = 1, widths = c(1,2)),
                       fig_testset[[2]], nrow = 2)
fig_comb1

p_drug_5row_testset = facet(fig_testset[[1]], facet.by = "drugs", nrow = 5, ncol=4, scales = "free") +
  theme_pubr(base_size = 5.5) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        # legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4))
p_drug_5row_testset

fig_comb1_b1 <- ggarrange(ggarrange(ggarrange(fig_trainset[[1]], fig_testset[[1]], nrow = 2),
                                    fig_testset[[3]],
                                    nrow = 1, widths = c(1,2)),
                          ggarrange(p_drug_5row_testset, p_drug_5row_testset, nrow = 1),
                          nrow = 2,
                          heights = c(1, 1.2)
)
fig_comb1_b1


### correlation heatmap by cancer and drugs
df_ytrue_ypred_v1 %>% head(2)
corr_by_cancerANDdrug <- df_ytrue_ypred_v1 %>% 
  group_by(dataset, cancer_type, drugs) %>% 
  summarize(coef_corr = cor(y_true, y_pred, method = "pearson"),
            pvalue_corr = cal_pval_of_cortest(y_true, y_pred),
            coef_corr_label = specify_decimal(coef_corr, k=2)) %>% 
  arrange(desc(coef_corr)) %>% 
  mutate(drugs = factor(drugs, levels = rev(drug_order_v1)))
corr_by_cancerANDdrug %>% head(2)

fig_hp_corr <- ggplot(data = corr_by_cancerANDdrug, aes(x = cancer_type, y = drugs)) +
  geom_tile(aes(fill = coef_corr)) +
  scale_fill_viridis_c() +
  geom_text(aes(label = coef_corr_label), color = "black", size = 1)
fig_hp_corr

library(pheatmap)
library(viridis)
library(RColorBrewer)

### trainging-set
df1 <- corr_by_cancerANDdrug %>% 
  filter(dataset == "Training") %>% 
  pivot_wider(id_cols = drugs, names_from = cancer_type, values_from = coef_corr)

mat1 <- df1 %>% select(-drugs) %>% as.matrix()
rownames(mat1) = df1$drugs
colnames(mat1) = colnames(df1)[-1]

pheatmap(mat1,
         # color = viridis(n=100),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100),
         cluster_cols = T,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         scale = "none",
         show_rownames = T,
         show_colnames = T,
         # treeheight_row = 1,
         # treeheight_col = 1,
         legend = TRUE,
         legend_breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
         display_numbers = TRUE,
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 5,
         fontsize_number = 0.8 * 5,
         main = NA
)


## testing set
df2 <- corr_by_cancerANDdrug %>% 
  filter(dataset == "Testing") %>% 
  pivot_wider(id_cols = drugs, names_from = cancer_type, values_from = coef_corr)

mat2 <- df2 %>% select(-drugs) %>% as.matrix()
rownames(mat2) = df2$drugs
colnames(mat2) = colnames(df2)[-1]

pheatmap(mat2,
         # color = viridis(n=100),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100),
         cluster_cols = T,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         scale = "none",
         show_rownames = T,
         show_colnames = T,
         # treeheight_row = 1,
         # treeheight_col = 1,
         legend = TRUE,
         legend_breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
         display_numbers = TRUE,
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize = 5.5,
         fontsize_number = 0.8 * 5,
         main = NA
)
dev.off()


