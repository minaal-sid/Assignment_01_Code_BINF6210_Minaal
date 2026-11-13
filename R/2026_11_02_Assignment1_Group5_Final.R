##***************************
##
## Assignment 1
##
## Group 5
##
## 2026-11-02
##
##***************************

## _ Packages used -------
library(tidyverse)
library(ggplot2)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::rename())
library(viridis)
theme_set(theme_light())

#####################
##### Setup ----#####
#####################

#Import the csv files
df_spiralis_te = read.csv("../data/GCA_001447595.1.csv")
df_nelsoni_te = read.csv("../data/GCA_001447455.1.csv")

#Filter out the subsets of Retrotransposons from transposons
ls_retro_te = c("LINE/CR1", "LINE/RTE-BovB", "LTR/DIRS", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad", "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE")
ls_transpo_te = c("DNA/CMC-Chapaev", "DNA/hAT", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc4", "RC/Helitron")
ls_retro_subset = c("PLE","PLE/Naiad","LTR/Pao", "LINE/CR1", "LTR/Gypsy")

#############################
#### Data Wrangling ----#####
#############################

#Summary of percentages for fig 1
df_spiralis_summary <- df_spiralis_te |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_transpo_te ~ "DNA transposon",
      matching_class %in% ls_retro_te ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "spiralis"
  )

# Repeating above coding steps for T. nelsoni

df_nelsoni_summary <- df_nelsoni_te |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_transpo_te ~ "DNA transposon",
      matching_class %in% ls_retro_te ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "nelsoni"
  )
df_te_summary <- rbind(df_spiralis_summary, df_nelsoni_summary)

# df_te_summary

#Filter the retrotransposon data for Figure 2

df_spiralis_retro_te = df_spiralis_te |> 
  filter(matching_class %in% ls_retro_te) |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n))


df_nelsoni_retro_te = df_nelsoni_te |> 
  filter(matching_class %in% ls_retro_te) |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n)) 

df_retros_long =
  left_join(df_nelsoni_retro_te, df_spiralis_retro_te, by="matching_class") |>
  replace_na(list (n.y = 0)) |>
  rename("T. nelsoni" = "n.x", "T. spiralis" = "n.y") |> 
  pivot_longer(cols = 2:3, names_to = "Species", values_to = "Count")

#Figure 3 data wrangling

df_spiralis_div <- df_spiralis_te |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. spiralis"
  ) 
df_nelsoni_div <- df_nelsoni_te |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. nelsoni"
  ) 

df_trich_div <- rbind(df_spiralis_div, df_nelsoni_div)

count(df_trich_div, species, matching_class)

####################
##### Plots ----####
####################

#Figure 1 - DNA Transposon vs Retrotransposon makeup between T. nelsoni and T. spiralis

ggplot(df_te_summary, aes(x = species, y = pct_copies, fill = te_type)) +
  geom_col() +
  labs(x = NULL, y = "% of TE copies", fill = NULL,
       title = "DNA transposons vs. Retrotransposons in Trichinella species") +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal()
ggsave("../output/2026-11-01_Fig1_retro_transpo_comparison.png", width = 6, dpi=600) 

#Figure 2 2 - Retrotransposon class composition between T. nelsoni and T. spiralis
ggplot(df_retros_long, aes(x = matching_class, y = Count, fill = `Species`)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Retrotransposon composition between T. spiralis and T. nelsoni", x = "Retrotransposon class", y = "Number of Retrotransposons") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave("../output/2026-11-01_Fig2_retro_comparison.png", width = 8, dpi=600) 

#Figure 3 - Percent divergence between certain retrotransposons between T. nelsoni and T. spiralis
ggplot(df_trich_div, aes(x = perc_div, fill = species)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~matching_class) +
  labs(x = "% divergence from consensus", y = "No. of copies", fill = NULL,
       title = "Repeat landscape for select retrotransposons in Trichinella species") +
  scale_fill_viridis_d(option = "viridis") +
  theme_dark()
ggsave("../output/2026-11-01_Fig3_retro_percdel_comparison.png", width = 8, dpi=600) 

#Mac users had to change permissions to Read and Write instead of Read only. Some of our group members had this error, but you probably already know this, but wanted to mention this just in case!

#####################
##### Stats-----#####
#####################

# Chi squared test for Fig 1
counts <- matrix(df_te_summary$n_copies, nrow = 2, byrow = TRUE)
#counts
chisq.test(counts)

# X-squared = 0.00018247, df = 1,
# p-value = 0.9892


#KS Test for Fig 3
spiralis_cr1 <- df_spiralis_div |>
  filter(matching_class == "LINE/CR1")

nelsoni_cr1 <- df_nelsoni_div |>
  filter(matching_class == "LINE/CR1")

ks_cr1 <- ks.test(spiralis_cr1$perc_div, nelsoni_cr1$perc_div)
rm(spiralis_cr1, nelsoni_cr1)

# LTR/Pao class

spiralis_pao <- df_spiralis_div |>
  filter(matching_class == "LTR/Pao")

nelsoni_pao <- df_nelsoni_div |>
  filter(matching_class == "LTR/Pao")

ks_pao <- ks.test(spiralis_pao$perc_div, nelsoni_pao$perc_div)
rm(spiralis_pao, nelsoni_pao)

# PLE class

spiralis_ple <- df_spiralis_div |>
  filter(matching_class == "PLE")

nelsoni_ple <- df_nelsoni_div |>
  filter(matching_class == "PLE")

ks_ple <- ks.test(spiralis_ple$perc_div, nelsoni_ple$perc_div)
rm(spiralis_ple, nelsoni_ple)

# PLE/Naiad class

spiralis_nai <- df_spiralis_div |>
  filter(matching_class == "PLE/Naiad")

nelsoni_nai <- df_nelsoni_div |>
  filter(matching_class == "PLE/Naiad")

ks_nai <- ks.test(spiralis_nai$perc_div, nelsoni_nai$perc_div)
rm(spiralis_nai, nelsoni_nai)

# LTR/Gypsy class

spiralis_gypsy <- df_spiralis_div |>
  filter(matching_class == "LTR/Gypsy")

nelsoni_gypsy <- df_nelsoni_div |>
  filter(matching_class == "LTR/Gypsy")

ks_gypsy <- ks.test(spiralis_gypsy$perc_div, nelsoni_gypsy$perc_div)
ks_gypsy
rm(spiralis_gypsy, nelsoni_gypsy)

# Adjusting p-values to ensure there are no false positives by multiplying by the number of KS tests performed (5)

ls_ks_pvalues <- c("LINE/CR1" = ks_cr1$p.value, "LTR/Pao" = ks_pao$p.value, "PLE" = ks_ple$p.value, "PLE/Naiad" = ks_nai$p.value, "LTR/Gypsy" = ks_gypsy$p.value)

p.adjust(ls_ks_pvalues, method = "bonferroni")

# P values:
# LINE/CR1   LTR/Pao       PLE PLE/Naiad 
# 1.0000000 0.2088786 1.0000000 1.0000000 
# LTR/Gypsy 
# 1.0000000 

rm(ks_cr1, ks_nai, ks_pao, ks_ple, ks_gypsy)