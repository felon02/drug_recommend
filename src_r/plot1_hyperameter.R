# plot1_hyperameter_v2
# 2023.8.28

getwd()
setwd("upload")

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

# hyperparameter ----------------------------------------------------------

df_stat_train = list()
df_stat_test = list()

for (k in c("dot", "concat_sample", "concat_gep")) {
  df_stat_train[[k]] = read_csv(paste("data_temp/out_model/df_stat_train_", k, ".csv", sep = ""))
  df_stat_test[[k]] = read_csv(paste("data_temp/out_model/df_stat_test_", k, ".csv", sep = ""))
}


#
df_stat_train[['dot']] %>% head() # stat_train is only mse
df_stat_test[['concat_sample']] %>% head() # stat test contains loss and mse

df_stat_train_v1 = list()
df_stat_test_v1 = list()
df_stat = list()
for (k in c("dot", "concat_sample", "concat_gep")) {
  df_stat_train_v1[[k]] <- df_stat_train[[k]] %>% 
    mutate(mean_perFold = `mean`) %>% 
    group_by(alpha) %>% 
    summarize(mean_perAlpha = mean(mean_perFold),
              std_perAlpha = sd(mean_perFold)) %>% 
    mutate(dataset = 'Training-set',
           alpha = factor(alpha, levels = c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)))
  
  df_stat_test_v1[[k]] <- df_stat_test[[k]] %>% 
    filter(`...1` == "mse") %>% 
    transmute(alpha = alpha, mean_perAlpha = `mean`, 
              std_perAlpha = `std`, dataset = "Testing-set") %>% 
    mutate(alpha = factor(alpha, levels = c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)))
  
  df_stat[[k]] <- df_stat_train_v1[[k]] %>% 
    bind_rows(df_stat_test_v1[[k]]) %>% 
    mutate(alpha = factor(alpha, levels = c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)))
}
df_stat_test_v1[['concat_sample']]
df_stat_test[['concat_sample']]
# line plot with errorbar


fig_errorbar = list()
for (k in c("dot", "concat_sample", "concat_gep")) {
  fig_errorbar[[k]] <- ggplot(df_stat_test_v1[[k]], aes(x=alpha, y=mean_perAlpha, group=1)) + 
    geom_line(size=.5, color = colorJco[2]) +
    geom_point(size = 2, color = colorJco[1])+
    geom_errorbar(aes(ymin=mean_perAlpha-std_perAlpha, 
                      ymax=mean_perAlpha+std_perAlpha), 
                  width=.2) +
    scale_color_jama() +
    scale_y_continuous(limits = c(-0.03, 1), n.breaks = 10, labels = specify_decimal_2) +
    labs(title=paste("Model:", k, sep = ""), 
         x="Alpha: the regular coefficient ", 
         y = "Mean MSE")+
    theme_bw() +
    theme(legend.position = c(0.85,0.85), legend.title = element_blank())
}

fig_errorbar2 <- ggarrange(fig_errorbar[['dot']], 
                           fig_errorbar[['concat_sample']], 
                           nrow = 1)
fig_errorbar2
