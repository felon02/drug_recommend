# plot2_v2_features_and_outcomes
# 2023.8.21

getwd()
setwd("upload")

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
color_20 = pal_d3(palette = "category20", alpha = 1)(20)
# show_col(colorJama)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
specify_decimal_2 <- function(x, k=2) trimws(format(round(x, k), nsmall=k))


# 1. data.load ---------------------------------------------------------------

### data from mai2_v2
ic50_gdsc_wide = read_csv("data_temp/ic50_gdsc_wide.csv")
pheno_gdsc_v1 = read_csv("data_temp/pheno_gdsc_v1.csv")
snv_gdsc_wide_v1 = read_csv("data_temp/snv_gdsc_wide_v1.csv")
tpm_gdsc_v1 = read_csv("data_temp/tpm_gdsc_v1.csv")

### ccle
y_ccle = read_csv("data_temp/out_ccle/v20230812/y_ccle.csv")
y_ccle_anno.csv = read_csv("data_temp/out_ccle/v20230812/y_ccle_anno.csv")
pheno_ccle_v2 = read_csv("data_temp/out_ccle/v20230812/pheno_ccle_v2.csv")
snv_ccle_v1 = read_csv("data_temp/out_ccle/v20230812/snv_v1.csv")
tpm_ccle_v1 = read_csv("data_temp/out_ccle/v20230812/tpm_v1.csv")


# 2. snv frequency -----------------------------------------------------------

### 2.1 snv frequency 
snv_gdsc_wide_v1 %>% head(2)
freq_gdsc <- snv_gdsc_wide_v1 %>% 
  select(-`...1`) %>% 
  map_dbl(.f = sum) %>% 
  sort(decreasing = TRUE)
freq_gdsc %>% head()

snv_ccle_v1 %>% head()
freq_ccle <- snv_ccle_v1 %>% 
  select(-`DepMap_ID`) %>% 
  map_dbl(.f = sum) %>% 
  sort(decreasing = TRUE)

freq2_df <- tibble(gene_symbol = names(freq_gdsc), mut_GDSC = freq_gdsc) %>% 
  inner_join(tibble(gene_symbol = names(freq_ccle), mut_CCLE = freq_ccle)) %>% 
  mutate(N_GDSC = nrow(snv_gdsc_wide_v1), N_CCLE = nrow(snv_ccle_v1),
         wt_GDSC = N_GDSC - mut_GDSC, wt_CCLE = N_CCLE - mut_CCLE)
freq2_df %>% head()

## add fisher.exact p value to the df
freq2_df_forP <- freq2_df %>% 
  dplyr::select(mut_GDSC, mut_CCLE, wt_GDSC, wt_CCLE) %>% 
  mutate(across(.fns = as.integer))
  # type_convert(col_types = cols(.default = col_integer()))

calP <- function(x){round(fisher.test(matrix(as.integer(x), nrow = 2))$p.value, 3)}
p_per_row = apply(freq2_df_forP, 1, calP)
p_per_row

freq2_df_v1 <- freq2_df %>% 
  bind_cols(tibble(p_fisher = p_per_row)) %>% 
  adjust_pvalue(p.col = "p_fisher", method = "fdr") %>% 
  mutate(freq_GDSC = mut_GDSC/N_GDSC * 100, freq_CCLE = mut_CCLE/N_CCLE * 100)
freq2_df_v1 %>% head(2)
# save
freq2_df_v1 %>% write_csv("data_temp/freq2_df_v1.csv")


### 2.2 density plot
fig_snv_density <- ggplot(data = freq2_df_v1, aes(x = freq_GDSC)) +
  geom_histogram(aes(y=after_stat(density)),
                 alpha=0.6, bins = 70,
                 color = "white", position="identity") +
  geom_density(alpha = 1, linewidth = 0.6, color = colorJama[1]) +
  labs(title =NULL, x = "Mutation frequency", y = "Density") +
  theme_pubr(base_size = 6) +
  theme(legend.position = c(0.85, 0.9),
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
fig_snv_density

### 2.3 plot top10 in both gdsc and CCLE
top10Gene_gdsc = freq2_df_v1 %>% slice_max(order_by = mut_GDSC, n = 10) %>% mutate(gene_origin = "GDSC")
top10Gene_ccle = freq2_df_v1 %>% slice_max(order_by = mut_CCLE, n = 10) %>% arrange(mut_CCLE) %>% mutate(gene_origin = "CCLE")
df_plot1 <- top10Gene_gdsc %>% 
  bind_rows(top10Gene_ccle) %>% 
  mutate(gene_origin = ifelse(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE), 
                              "GDSC&CCLE", gene_origin)) %>% 
  distinct(gene_symbol, .keep_all = TRUE) %>% 
  mutate(gene_color = case_when(gene_origin == "GDSC" ~ colorJama[1],
                                gene_origin == "CCLE" ~ colorJama[2],
                                gene_origin == "GDSC&CCLE" ~ colorJama[3]))
df_plot1 %>% head(2)

df_plot1_v1 <- df_plot1 %>% 
  pivot_longer(cols = c(freq_GDSC, freq_CCLE), names_to = "Database", values_to = "Frequency") %>% 
  mutate(gene_symbol = factor(gene_symbol, levels = rev(df_plot1$gene_symbol)),
         Database = factor(Database, levels = c("freq_CCLE","freq_GDSC"), labels = c("CCLE", "GDSC")))
df_plot1_v1 %>% head()

df_plot1 %>% head(2)
df_plot1_p <- df_plot1 %>% 
  transmute(gene_symbol = gene_symbol,
            group1 = "GDSC", 
            group2 = "CCLE", 
            p = p_fisher.adj,
            p3decimal = round(p_fisher.adj, 3),
            y.position = pmax(freq_GDSC, freq_CCLE) + 3) %>% 
  add_significance("p") 
df_plot1_p %>% head(2)
#
?geom_bar
fig_freq <- ggplot(data = df_plot1_v1, aes(x = gene_symbol, y = Frequency)) +
  geom_bar(aes(fill = Database), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(GDSC = colorJama[1], CCLE = colorJama[2])) +
  stat_pvalue_manual(df_plot1_p,  x = "gene_symbol", label = "p.signif", size=3, hide.ns = T) +
  labs(x = NULL, y = "Mutation frequency, %", fill = NULL,
       title = "Top 10 mutant genes in GDSC/CCLE") +
  coord_flip() +
  theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, color = rev(df_plot1$gene_color)),
        legend.position = c(0.8, 0.2))
fig_freq


### 2.4 figtop10 genes in gdsc
top10Gene_gdsc
top10Gene_gdsc$gene_symbol
df_plot_gdsc_1 <- top10Gene_gdsc %>% mutate(gene_symbol = factor(gene_symbol, levels = rev(.$gene_symbol)))
fig_freqTop10_gdsc <- ggplot(data = df_plot_gdsc_1, aes(x = gene_symbol, y = freq_GDSC)) +
  geom_bar(stat = "identity", fill = colorJama[1], width = 0.6) +
  labs(x = NULL, y = "Mutation frequency, %", fill = NULL,
       title = NULL) +
  coord_flip() +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = c(0.8, 0.2))
fig_freqTop10_gdsc


# 3. Highly variable genes in gdsc -------------------------------------------

### 3.1 std calculation
tpm_gdsc_v1 %>% head(2)
tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(mean) %>% summary()

gene_cv <- tpm_gdsc_v1 %>% select(-`...1`) %>% map_dbl(function(x){sd(x)/abs(mean(x)+10^(-5))})
gene_cv = sort(gene_cv, decreasing = T)
gene_cv

### 3.2 heatmap of the top10 genes
tb_hvg_gdsc <- tpm_gdsc_v1 %>% select(`...1`, all_of(hvg_top10)) %>% 
  mutate(across(all_of(hvg_top10), .fns = function(x){(x-mean(x)/sd(x))})) %>% 
  pivot_longer(cols = all_of(hvg_top10), names_to = "Gene", values_to = "z_score")
  
tb_hvg_gdsc %>% head(2)
fig_hp <- ggplot(data = tb_hvg_gdsc, aes(y = Gene, x = `...1`, fill = z_score)) +
  geom_tile() +
  scale_fill_viridis_c()
fig_hp
# conclusion: useless, i can not cluster samples

### 3.3. density of std
gene_cv
df_cv = tibble(cv = gene_cv, Gene = names(gene_cv))

fig_cv_density <- ggplot(data = df_cv, aes(x = cv)) +
  geom_histogram(aes(y=after_stat(density)),
                 alpha=0.6, bins = 70,
                 color = "white", position="identity") +
  geom_density(alpha = 1, linewidth = 0.6, color = colorJama[1]) +
  scale_x_log10() +
  labs(title = NULL, x = "log10(Coefficient of Variation)", y = "Density") +
  theme_pubr(base_size = 5) +
  theme(legend.position = c(0.85, 0.9),
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5))
fig_cv_density


df_cv %>% head(2)
df_plot_cv <- df_cv %>% slice_max(order_by = cv, n = 10) %>% 
  mutate(gene_symbol = factor(Gene, levels = rev(.$Gene)))
fig_cv_top10 <- ggplot(data = df_plot_cv, aes(x = gene_symbol, y = cv)) +
  geom_bar(stat = "identity", fill = colorJama[1], width = 0.6) +
  labs(x = "Coefficient of Variation", y = "Coeficient of Variation", fill = NULL,
       title = NULL) +
  coord_flip() +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = c(0.8, 0.2))
fig_cv_top10


# 4 snv mutation frequency and gep variance-----------------------------------------------------------------------
 
freq2_df_v1 %>% head(2)
df_cv %>% head(2)
tb_freq_cv <- freq2_df_v1 %>% select(Gene = gene_symbol, Mutation_frequency = freq_GDSC) %>% 
  inner_join(df_cv, by = "Gene") %>% 
  mutate(gene_to_anno = ifelse((Mutation_frequency > 20) | cv > 10, Gene, NA),
         gene_color = case_when(Mutation_frequency > 20 ~ colorJco[1],
                                cv > 10 ~ colorJco[2],
                                (Mutation_frequency > 20) & (cv>10) ~ colorJco[3]))
tb_freq_cv
fig_scatter_mut_cv <- ggplot(data = tb_freq_cv, aes(x = Mutation_frequency, y = cv)) +
  geom_point(size = 0.5, color = colorJama[1])+
  geom_text_repel(aes(label = gene_to_anno, color = gene_color), size = 2, show.legend = F) +
  scale_color_jco() +
  scale_y_log10() +
  labs(x = "Mutation Frequency", y = "log10(Coeficient of Variation)") +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5))
fig_scatter_mut_cv
# conclusion: snv and gep are good supplement



# 5. cancer type: donut plot; sex -----------------------------------------------

### 5.1. cancer
pheno_gdsc_v1 %>% colnames()
pheno_gdsc_v1 %>% count(cancerType_TCGA)
df_plot_donut_cancer <- pheno_gdsc_v1 %>% 
  count(cancerType_TCGA) %>% 
  arrange(desc(n)) %>% 
  mutate(prop = round(n / sum(n) * 100, 1) , 
         cancer_label = paste(cancerType_TCGA, as.character(prop), sep = "\n"))

color_cancers = color_20[1:length(df_plot_donut_cancer$cancerType_TCGA)]
names(color_cancers) = df_plot_donut_cancer$cancerType_TCGA

#
hsize = 2

fig_pie_cancers = ggplot(df_plot_donut_cancer, aes(x = hsize, y = prop, fill = cancerType_TCGA)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(label = cancer_label), position = position_stack(vjust = 0.5), color = "black", size = 2.5)+
  scale_fill_manual(values = color_cancers) +
  xlim(0, hsize+0.5) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + # reduce legend area
  theme_void(base_size = 6) +
  theme(legend.position = "none")
fig_pie_cancers


### 5.2. sex
color_sex = structure(colorJco[c(6,9)], names = c("Male", "Female"))
color_sex

pheno_gdsc_v1 %>% colnames()
df_fig_pie_sex <- pheno_gdsc_v1 %>% 
  count(Sex) %>% 
  arrange(desc(n)) %>% 
  mutate(prop = round(n / sum(n) * 100, 1) , 
         cancer_label = paste(Sex, as.character(prop), sep = "\n"))
df_fig_pie_sex
fig_pie_sex = ggplot(df_fig_pie_sex, aes(x = hsize, y = prop, fill = Sex)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(label = Sex), position = position_stack(vjust = 0.5), color = "black", size = 2.5)+
  scale_fill_manual(values = color_sex) +
  xlim(0, hsize+0.5) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + # reduce legend area
  theme_void(base_size = 6) +
  theme(legend.position = "none")
fig_pie_sex


### Age
pheno_gdsc_v1 %>% count(Age)
fig_density_age <- ggplot(data = pheno_gdsc_v1, aes(x = Age)) +
  geom_histogram(aes(y=after_stat(density)),
                 alpha=0.6, bins = 50,
                 color = "white", position="identity") +
  geom_density(alpha = 1, linewidth = 1, color = colorJama[2]) +
  scale_x_continuous(n.breaks = 10) +
  labs(title = NULL, x = "Age", y = "Density") +
  theme_pubr(base_size = 6) +
  theme(legend.position = c(0.85, 0.9),
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
fig_density_age



# 6.outcome ---------------------------------------------------------------

ic50_gdsc_wide %>% head(2)
pheno_gdsc_v1 %>% colnames()

### 6.1 heatmap
library(pheatmap)

df_fig_hp <- ic50_gdsc_wide %>% 
  inner_join(pheno_gdsc_v1, by = "SANGER_MODEL_ID") %>% 
  mutate(cancerType_TCGA = factor(cancerType_TCGA, levels = names(color_cancers))) %>% 
  arrange(cancerType_TCGA)
df_fig_hp_v1 <- df_fig_hp %>% 
  select(-all_of(colnames(pheno_gdsc_v1))) %>% 
  as.data.frame()
rownames(df_fig_hp_v1) = df_fig_hp$SANGER_MODEL_ID

# sample annotation
df_fig_hp %>% select(Sex, cancerType_TCGA) %>% as.data.frame() -> df_fig_hp_anno_row
rownames(df_fig_hp_anno_row) = df_fig_hp$SANGER_MODEL_ID

# set color of annotation
color_cancers
pheno_gdsc_v1 %>% count(Sex)
anno_colour = list(
  cancerType_TCGA = color_cancers,
  Sex = color_sex)

# plot
pheatmap(t(df_fig_hp_v1),
         annotation_col = df_fig_hp_anno_row,
         annotation_colors = anno_colour,
         cluster_cols = F,
         cluster_rows = T,
         # clustering_distance_cols = "correlation",
         scale = "none",
         show_colnames = FALSE,
         fontsize = 5,
         main = "GDSC: IC50"
)


### mean IC50 per cancer type
df_fig_hp %>% head(2)
drugs = setdiff(colnames(ic50_gdsc_wide), c("SANGER_MODEL_ID"))
drugs

df_meanIC50_byCancers <- df_fig_hp %>% 
  select(cancerType_TCGA, all_of(drugs)) %>% 
  pivot_longer(cols = all_of(drugs), names_to = "drugs", values_to = "ic50") %>% 
  group_by(cancerType_TCGA, drugs) %>% 
  summarize(meanIC50 = mean(ic50))

drug_order <- df_meanIC50_byCancers %>% group_by(drugs) %>% summarize(meanIC50 = mean(meanIC50)) %>% 
  arrange(meanIC50) %>% pull(drugs)
drug_order


df_meanIC50_byCancers$drugs = factor(df_meanIC50_byCancers$drugs, levels = rev(drug_order))

fig_meanIC50 <- ggplot(data = df_meanIC50_byCancers, aes(x = drugs, y = meanIC50, color = cancerType_TCGA)) +
  geom_boxplot(color = "gray20", fill = "transparent", size = 0.3, outlier.size = 0) +
  geom_point(size = 0.5, position = position_jitter(width = 0.2)) +
  scale_color_manual(values = color_cancers) +
  labs(y = "Mean IC50", x = NULL, color = NULL) +
  coord_flip() +
  theme_bw(base_size = 5) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4,"cm"),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5))
  
fig_meanIC50
