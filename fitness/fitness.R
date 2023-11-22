library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

df <- read.csv(file="fitness/2022_summer_fitness_assay.csv")

df %>%
  tidyr::gather(-Condition, -Date, -Rep, key='trait', value='value') %>%
  ggplot(.) +
  geom_jitter(width=0.2, aes(color=Date)) +
  geom_boxplot(alpha=0.5) +
  aes(x=Condition, y=value/10) +
  theme_bw() +
  facet_grid(~factor(trait, levels=c("N_egg", "N_pupa"), labels = c("Egg", "Pupa"))) +
  labs(y="Counts / female", x=NULL)

df_fitness <- read.csv(file="fitness/2022_summer_fitness_assay.csv") %>%
     #dplyr::mutate(frac_pup = N_pupa/N_egg) %>%
     tidyr::gather(-Condition, -Date, -Rep, key='trait', value="counts")

df_fitness$Condition <- factor(df_fitness$Condition, levels=c("Control", "Starch_DW", "Starch_noni"), labels=c("Control", "+ Blue food", "+ Blue food \n+ Noni juice"))

df_fitness %>%
     ggplot(.)+
     geom_jitter(width=0.2, alpha=0.8) +
     geom_boxplot(alpha=0.5, aes(fill=Condition)) +
     aes(x=Condition, y=counts/10) +
     theme_bw() +
     theme(axis.title=element_text(size=12, color='black'),
                     axis.text=element_text(size=11.2, color='black'),
                     strip.text=element_text(size=12, color='black'),
                     legend.position = 'none') +
     facet_wrap(~factor(trait, levels = c("N_egg","N_pupa"), labels=c("Egg", "Pupa"))) +
     labs(y="Counts / female", x=NULL) +
     scale_fill_manual(values=c("grey80","grey10", "darkgreen"))

ggsave(file="Plots/Manuscript/Fig6a_boxplot_fitness.pdf", width=6.5, height=4)
