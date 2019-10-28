

### Project: Polio Mixing Experiment 2.0
### Purpose: Collapse output variant CSV by sequencing duplicate.

# ======================== Import packages, read in data =============================

library(tidyverse)

sampleInfo <- read.csv("./data/reference/MixingStudySamples2.csv")

variants_senslocal_binomial_100 <- read.csv("./data/raw/shiny2.senslocal.binomial.100.variants.csv")
variants_senslocal_onesided_100 <- read.csv("./data/raw/shiny2.senslocal.onesided.100.variants.csv")

variants_senslocal_binomial_50 <- read.csv("./data/raw/shiny2.senslocal.binomial.50.variants.csv")
variants_senslocal_onesided_50 <- read.csv("./data/raw/shiny2.senslocal.onesided.50.variants.csv")

variants_senslocal_binomial_25 <- read.csv("./data/raw/shiny2.senslocal.binomial.25.variants.csv")
variants_senslocal_onesided_25 <- read.csv("./data/raw/shiny2.senslocal.onesided.25.variants.csv")

variants_senslocal_binomial_10 <- read.csv("./data/raw/shiny2.senslocal.binomial.10.variants.csv")
variants_senslocal_onesided_10 <- read.csv("./data/raw/shiny2.senslocal.onesided.10.variants.csv")

variants_verysens_binomial_100 <- read.csv("./data/raw/shiny2.verysens.binomial.100.variants.csv")
variants_verysens_onesided_100 <- read.csv("./data/raw/shiny2.verysens.onesided.100.variants.csv")

variants_verysens_binomial_50 <- read.csv("./data/raw/shiny2.verysens.binomial.50.variants.csv")
variants_verysens_onesided_50 <- read.csv("./data/raw/shiny2.verysens.onesided.50.variants.csv")

variants_verysens_binomial_25 <- read.csv("./data/raw/shiny2.verysens.binomial.25.variants.csv")
variants_verysens_onesided_25 <- read.csv("./data/raw/shiny2.verysens.onesided.25.variants.csv")

variants_verysens_binomial_10 <- read.csv("./data/raw/shiny2.verysens.binomial.10.variants.csv")
variants_verysens_onesided_10 <- read.csv("./data/raw/shiny2.verysens.onesided.10.variants.csv")


# ======================== Collapse functions ============================

sift_dups <- function(df)
{
  if(nrow(df) > 2) stop("Too many mutations here")
  
  df <- dplyr::mutate(df, coverage = cov.tst.bw + cov.tst.fw)
  higher_qual <- subset(df, coverage == max(df$coverage))
  if(nrow(higher_qual) > 1)
  { 
    higher_qual <- higher_qual[1,]
  }
  return(higher_qual)
}

# Adapted for polio data.
collapse_localcov <- function(df)
{
  stopifnot(length(unique(df$group)) == 1)
  
  df %>% dplyr::group_by(mutation) %>% dplyr::summarize(found = length(mutation)) -> count_mutations
    
  stopifnot(max(count_mutations$found) < 3)
    
  count_mutations %>% dplyr::filter(found == 2) -> good_mut
    
  df %>% dplyr::filter(mutation %in% good_mut$mutation) -> good_var
    
  good_var %>% dplyr::group_by(mutation) %>% dplyr::do(sift_dups(.)) -> dups_good
    
  return(dplyr::ungroup(dups_good))
}

# ========================== Main ==============================

expectedVariants <- read.csv("./data/reference/MixingStudyExpectedVariants.csv")
sampleInfo <- mutate(sampleInfo, Id = SampleNumber)

### Binomial Data

# 100%
variants_with_meta <- left_join(variants_senslocal_binomial_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_100, "./data/processed/shiny2.senslocal.binomial.100.notcollapsed.variants.csv")
write_csv(bin_collapsed_100, "./data/processed/shiny2.senslocal.binomial.100.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_binomial_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_100, "./data/processed/shiny2.verysens.binomial.100.notcollapsed.variants.csv")
write_csv(bin_collapsed_100, "./data/processed/shiny2.verysens.binomial.100.collapsed.variants.csv")

# 50%
variants_with_meta <- left_join(variants_senslocal_binomial_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_50, "./data/processed/shiny2.senslocal.binomial.50.notcollapsed.variants.csv")
write_csv(bin_collapsed_50, "./data/processed/shiny2.senslocal.binomial.50.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_binomial_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_50, "./data/processed/shiny2.verysens.binomial.50.notcollapsed.variants.csv")
write_csv(bin_collapsed_50, "./data/processed/shiny2.verysens.binomial.50.collapsed.variants.csv")

# 25%
variants_with_meta <- left_join(variants_senslocal_binomial_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_25, "./data/processed/shiny2.senslocal.binomial.25.notcollapsed.variants.csv")
write_csv(bin_collapsed_25, "./data/processed/shiny2.senslocal.binomial.25.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_binomial_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_25, "./data/processed/shiny2.verysens.binomial.25.notcollapsed.variants.csv")
write_csv(bin_collapsed_25, "./data/processed/shiny2.verysens.binomial.25.collapsed.variants.csv")

# 10%
variants_with_meta <- left_join(variants_senslocal_binomial_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_10, "./data/processed/shiny2.senslocal.binomial.10.notcollapsed.variants.csv")
write_csv(bin_collapsed_10, "./data/processed/shiny2.senslocal.binomial.10.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_binomial_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_10, "./data/processed/shiny2.verysens.binomial.10.notcollapsed.variants.csv")
write_csv(bin_collapsed_10, "./data/processed/shiny2.verysens.binomial.10.collapsed.variants.csv")

### One-Sided Data

# 100%
variants_with_meta <- left_join(variants_senslocal_onesided_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_100, "./data/processed/shiny2.senslocal.onesided.100.notcollapsed.variants.csv")
write_csv(one_collapsed_100, "./data/processed/shiny2.senslocal.onesided.100.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_onesided_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_100, "./data/processed/shiny2.verysens.onesided.100.notcollapsed.variants.csv")
write_csv(one_collapsed_100, "./data/processed/shiny2.verysens.onesided.100.collapsed.variants.csv")

# 50%
variants_with_meta <- left_join(variants_senslocal_onesided_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_50, "./data/processed/shiny2.senslocal.onesided.50.notcollapsed.variants.csv")
write_csv(one_collapsed_50, "./data/processed/shiny2.senslocal.onesided.50.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_onesided_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_50, "./data/processed/shiny2.verysens.onesided.50.notcollapsed.variants.csv")
write_csv(one_collapsed_50, "./data/processed/shiny2.verysens.onesided.50.collapsed.variants.csv")

# 25%
variants_with_meta <- left_join(variants_senslocal_onesided_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_25, "./data/processed/shiny2.senslocal.onesided.25.notcollapsed.variants.csv")
write_csv(one_collapsed_25, "./data/processed/shiny2.senslocal.onesided.25.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_onesided_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_25, "./data/processed/shiny2.verysens.onesided.25.notcollapsed.variants.csv")
write_csv(one_collapsed_25, "./data/processed/shiny2.verysens.onesided.25.collapsed.variants.csv")

# 10%
variants_with_meta <- left_join(variants_senslocal_onesided_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_10, "./data/processed/shiny2.senslocal.onesided.10.notcollapsed.variants.csv")
write_csv(one_collapsed_10, "./data/processed/shiny2.senslocal.onesided.10.collapsed.variants.csv")

variants_with_meta <- left_join(variants_verysens_onesided_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_10, "./data/processed/shiny2.verysens.onesided.10.notcollapsed.variants.csv")
write_csv(one_collapsed_10, "./data/processed/shiny2.verysens.onesided.10.collapsed.variants.csv")

### Prepare coverage data ###

# sens local
coverage_sl_binomial_100 <- read.csv("./data/raw/shiny2.senslocal.binomial.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_sl_onesided_100 <- read.csv("./data/raw/shiny2.senslocal.onesided.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_sl_binomial_50 <- read.csv("./data/raw/shiny2.senslocal.binomial.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_sl_onesided_50 <- read.csv("./data/raw/shiny2.senslocal.onesided.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_sl_binomial_25 <- read.csv("./data/raw/shiny2.senslocal.binomial.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_sl_onesided_25 <- read.csv("./data/raw/shiny2.senslocal.onesided.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_sl_binomial_10 <- read.csv("./data/raw/shiny2.senslocal.binomial.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_sl_onesided_10 <- read.csv("./data/raw/shiny2.senslocal.onesided.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

sampleInfo <- mutate(sampleInfo, Id = as.character(Id))
coverage_sl_binomial_100_wmeta <- left_join(coverage_sl_binomial_100, sampleInfo, by = "Id")
coverage_sl_onesided_100_wmeta <- left_join(coverage_sl_onesided_100, sampleInfo, by = "Id")
coverage_sl_binomial_50_wmeta <- left_join(coverage_sl_binomial_50, sampleInfo, by = "Id")
coverage_sl_onesided_50_wmeta <- left_join(coverage_sl_onesided_50, sampleInfo, by = "Id")
coverage_sl_binomial_25_wmeta <- left_join(coverage_sl_binomial_25, sampleInfo, by = "Id")
coverage_sl_onesided_25_wmeta <- left_join(coverage_sl_onesided_25, sampleInfo, by = "Id")
coverage_sl_binomial_10_wmeta <- left_join(coverage_sl_binomial_10, sampleInfo, by = "Id")
coverage_sl_onesided_10_wmeta <- left_join(coverage_sl_onesided_10, sampleInfo, by = "Id")

write_csv(coverage_sl_binomial_100_wmeta, "./data/processed/shiny2.senslocal.binomial.100.coverage.csv")
write_csv(coverage_sl_onesided_100_wmeta, "./data/processed/shiny2.senslocal.onesided.100.coverage.csv")
write_csv(coverage_sl_binomial_50_wmeta, "./data/processed/shiny2.senslocal.binomial.50.coverage.csv")
write_csv(coverage_sl_onesided_50_wmeta, "./data/processed/shiny2.senslocal.onesided.50.coverage.csv")
write_csv(coverage_sl_binomial_25_wmeta, "./data/processed/shiny2.senslocal.binomial.25.coverage.csv")
write_csv(coverage_sl_onesided_25_wmeta, "./data/processed/shiny2.senslocal.onesided.25.coverage.csv")
write_csv(coverage_sl_binomial_10_wmeta, "./data/processed/shiny2.senslocal.binomial.10.coverage.csv")
write_csv(coverage_sl_onesided_10_wmeta, "./data/processed/shiny2.senslocal.onesided.10.coverage.csv")

# very sens
coverage_vs_binomial_100 <- read.csv("./data/raw/shiny2.verysens.binomial.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_vs_onesided_100 <- read.csv("./data/raw/shiny2.verysens.onesided.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_vs_binomial_50 <- read.csv("./data/raw/shiny2.verysens.binomial.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_vs_onesided_50 <- read.csv("./data/raw/shiny2.verysens.onesided.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_vs_binomial_25 <- read.csv("./data/raw/shiny2.verysens.binomial.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_vs_onesided_25 <- read.csv("./data/raw/shiny2.verysens.onesided.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_vs_binomial_10 <- read.csv("./data/raw/shiny2.verysens.binomial.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_vs_onesided_10 <- read.csv("./data/raw/shiny2.verysens.onesided.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_vs_binomial_100_wmeta <- left_join(coverage_vs_binomial_100, sampleInfo, by = "Id")
coverage_vs_onesided_100_wmeta <- left_join(coverage_vs_onesided_100, sampleInfo, by = "Id")
coverage_vs_binomial_50_wmeta <- left_join(coverage_vs_binomial_50, sampleInfo, by = "Id")
coverage_vs_onesided_50_wmeta <- left_join(coverage_vs_onesided_50, sampleInfo, by = "Id")
coverage_vs_binomial_25_wmeta <- left_join(coverage_vs_binomial_25, sampleInfo, by = "Id")
coverage_vs_onesided_25_wmeta <- left_join(coverage_vs_onesided_25, sampleInfo, by = "Id")
coverage_vs_binomial_10_wmeta <- left_join(coverage_vs_binomial_10, sampleInfo, by = "Id")
coverage_vs_onesided_10_wmeta <- left_join(coverage_vs_onesided_10, sampleInfo, by = "Id")

write_csv(coverage_vs_binomial_100_wmeta, "./data/processed/shiny2.verysens.binomial.100.coverage.csv")
write_csv(coverage_vs_onesided_100_wmeta, "./data/processed/shiny2.verysens.onesided.100.coverage.csv")
write_csv(coverage_vs_binomial_50_wmeta, "./data/processed/shiny2.verysens.binomial.50.coverage.csv")
write_csv(coverage_vs_onesided_50_wmeta, "./data/processed/shiny2.verysens.onesided.50.coverage.csv")
write_csv(coverage_vs_binomial_25_wmeta, "./data/processed/shiny2.verysens.binomial.25.coverage.csv")
write_csv(coverage_vs_onesided_25_wmeta, "./data/processed/shiny2.verysens.onesided.25.coverage.csv")
write_csv(coverage_vs_binomial_10_wmeta, "./data/processed/shiny2.verysens.binomial.10.coverage.csv")
write_csv(coverage_vs_onesided_10_wmeta, "./data/processed/shiny2.verysens.onesided.10.coverage.csv")
