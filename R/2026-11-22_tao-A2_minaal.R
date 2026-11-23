### Assignment 2 
# Original Author: Minaal Siddiqui 
# - created original script file and figures along with Adhithya Dhivakar, Tomas Fischer, Purushoth Puvanendran, Bhumika Aparanji 

# Editor: William Tao
# - edits marked by #edit, suggestions marked by #suggestion
# 1st edit = #edit1 2nd edit = #edit2 ...




##***************************
##
## Assignment 1
##
## Group 5
##
## 2026-11-02
##
##***************************

#suggestion - consider making chunk titles consistent and removing all the "(untitled)" title labels which are caused when 5 # are put in a row. makes it a lot easier to read and follow in my opinion. instead of # maybe add something like ##****** instead if that format is wanted

##Packages used -------
library(tidyverse)
library(ggplot2)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::rename())
library(viridis)
theme_set(theme_light())

#edit - changing working directory for project 
setwd("./R/")

## *****
##### Setup ----#####
## *****

#Import the csv files
df_spiralis_te = read.csv("../data/GCA_001447595.1.csv")
df_nelsoni_te = read.csv("../data/GCA_001447455.1.csv")

#edit2: second edit pt 1 starts here (has 2 parts)
#command + f "edit2" to find them

#combine the data frames into 1 and then add a "species" column so that the code doesn't have to repeat for each step
df_combined_te = bind_rows(
  df_spiralis_te |> 
    mutate(species = "T. spiralis"),
  df_nelsoni_te |> 
    mutate(species = "T. nelsoni")
)
#edit2: pt 1 of second edit ends here

#Filter out the subsets of Retrotransposons from transposons
ls_retro_te = c("LINE/CR1", "LINE/RTE-BovB", "LTR/DIRS", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad", "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE")
ls_transpo_te = c("DNA/CMC-Chapaev", "DNA/hAT", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc4", "RC/Helitron")
ls_retro_subset = c("PLE","PLE/Naiad","LTR/Pao", "LINE/CR1", "LTR/Gypsy")

#edit1: first edit starts here
#code is missing some classes that exist in either species (DNA/TcMar-Tc2, LINE, DNA/TcMar-ISRm11)
#check to see how many copies were missed from original code
bind_rows(
  df_spiralis_te |> 
    mutate(species = "T. spiralis"),
  df_nelsoni_te |> 
    mutate(species = "T. nelsoni")) |> 
  filter(matching_class %in% c("DNA/TcMar-ISRm11", "DNA/TcMar-Tc2", "LINE")) |> 
  count(species, matching_class)

#only 10 total count missing, so would probably not impact the results much but still important to include in case it were more

#fixed filter to include missing retrotransposons and DNA transposons
#also removed classes that don't exist in either species (RC/Helitron, DNA/CMC-Chapaev)
ls_retro_te = c("LINE/CR1", "LINE/RTE-BovB", "LTR/DIRS", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad", "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE", "LINE")
ls_transpo_te = c("DNA/TcMar-ISRm11", "DNA/hAT", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc4", "DNA/TcMar-Tc2")
ls_retro_subset = c("PLE","PLE/Naiad","LTR/Pao", "LINE/CR1", "LTR/Gypsy")

#another potential way of doing the classes is to organize them by prefixes (similar to what we do in the programming class) Eg. all DNA prefix will be labelled DNA transposons so that it's not as likely to accidentally miss classes
#combined all unique classes 
all_classes = c(df_spiralis_te$matching_class, df_nelsoni_te$matching_class) |> 
  unique()
#retrotransposons
ls_retro_te = str_subset(all_classes, "^(LINE|LTR|SINE|PLE)")
ls_retro_te

#dna transposons
ls_transpo_te = str_subset(all_classes, "^(DNA)")
ls_transpo_te
#edit1: first edit ends here


## *****
#### Data Wrangling ----#####
## *****

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

#edit2: second edit pt2 starts here
#after combining the data into 1 df_combined_te, the code for the data wrangling of fig 1. changes

#essentially just uses the same code but with the df_combined_te data frame. one difference to note is that we have to add another group_by near the end so that the percentages add up to 100 within each respective species and not between them
df_te_summary = df_combined_te |> 
  mutate(te_type = case_when(
    matching_class %in% ls_transpo_te ~ "DNA transposon",
    matching_class %in% ls_retro_te ~ "Retrotransposon",
    TRUE ~ "Other"
    )
  ) |> 
  filter(te_type != "Other") |>
  group_by(species, te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  group_by(species) |> 
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp)
  )
df_te_summary
#edit2: second edit pt2 ends here

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

## *****
#### Plots ----####
## *****

#suggestion - in my opinion it's easier to follow and code when the same/related items are grouped together and back to back. for example having the code chunk for data wrangling of figure 1 and then right after that having the code chunk for the plotting of figure 1, instead of having all the data wrangling for all the figures and then having all the plotting. This way if something is wrong or there's an error you don't have to constantly scroll up and down because the code you need to change for the figure is right above it rather than in a different section entirely. (but also that's my preference)

#Figure 1 - DNA Transposon vs Retrotransposon makeup between T. nelsoni and T. spiralis

ggplot(df_te_summary, aes(x = species, y = pct_copies, fill = te_type)) +
  geom_col() +
  labs(x = NULL, y = "% of TE copies", fill = NULL,
       title = "DNA transposons vs. Retrotransposons in Trichinella species") +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal()
ggsave("../output/2026-11-01_Fig1_retro_transpo_comparison.png", width = 6, dpi=600) 

#Figure 2 - Retrotransposon class composition between T. nelsoni and T. spiralis
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

## *****
##### Stats-----#####
## *****

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

#edit3: third edit starts here 
#adding wilcox test could help support the hypothesis because it suggests looking into the te insertion ages.the ks is also good but is more broad and looks at all types of potential effects that could affect the distribution, so it might not show exactly what's wanted. I think the wilcox looks directly at the shift but also has its limitations so doing it with the KS test is probably worth.

#wilcox test
wilcox_results <- df_combined_te |>
  filter(
    matching_class %in% ls_retro_subset,
    is.finite(perc_div)
  ) |>
  group_by(matching_class) |>
  group_modify(~ {
    x <- .x$perc_div[.x$species == "T. spiralis"]
    y <- .x$perc_div[.x$species == "T. nelsoni"]
    
    tibble(
      n_spiralis = length(x),
      n_nelsoni = length(y),
      median_spiralis = median(x),
      median_nelsoni = median(y),
      p_value = wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
    )
  }) |>
  ungroup() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "bonferroni"),
    significant = p_adjusted < 0.05
  ) |>
  arrange(p_adjusted)

wilcox_results

#edit3: third edit ends here
