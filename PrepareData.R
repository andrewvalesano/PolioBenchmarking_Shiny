

### Project: Polio Mixing Experiment 2.0
### Purpose: Collapse output variant CSV by sequencing duplicate.

# ======================== Import packages, read in data =============================

library(tidyverse)

sampleInfo <- read.csv("./data/reference/MixingStudySamples2.csv")

variants_binomial_100 <- read.csv("./data/raw/shiny2.binomial.100.variants.csv")
variants_onesided_100 <- read.csv("./data/raw/shiny2.onesided.100.variants.csv")

variants_binomial_50 <- read.csv("./data/raw/shiny2.binomial.50.variants.csv")
variants_onesided_50 <- read.csv("./data/raw/shiny2.onesided.50.variants.csv")

variants_binomial_25 <- read.csv("./data/raw/shiny2.binomial.25.variants.csv")
variants_onesided_25 <- read.csv("./data/raw/shiny2.onesided.25.variants.csv")

variants_binomial_10 <- read.csv("./data/raw/shiny2.binomial.10.variants.csv")
variants_onesided_10 <- read.csv("./data/raw/shiny2.onesided.10.variants.csv")

#variants_twosided <- read.csv("./data/raw/shiny2.twosided.variants.csv") # data file seems to be broken. Some columns shifted leftward. Ugh.
#filter(variants_twosided, chr != "pT7_S1F") %>% View()
#filter(variants_twosided, chr == "pT7_S1F") %>% View()

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
variants_with_meta <- left_join(variants_binomial_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_100, "./data/processed/shiny2.binomial.100.notcollapsed.variants.csv")
write_csv(bin_collapsed_100, "./data/processed/shiny2.binomial.100.collapsed.variants.csv")

# 50%
variants_with_meta <- left_join(variants_binomial_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_50, "./data/processed/shiny2.binomial.50.notcollapsed.variants.csv")
write_csv(bin_collapsed_50, "./data/processed/shiny2.binomial.50.collapsed.variants.csv")

# 25%
variants_with_meta <- left_join(variants_binomial_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_25, "./data/processed/shiny2.binomial.25.notcollapsed.variants.csv")
write_csv(bin_collapsed_25, "./data/processed/shiny2.binomial.25.collapsed.variants.csv")

# 10%
variants_with_meta <- left_join(variants_binomial_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
bin_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
bin_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(bin_notcollapsed_10, "./data/processed/shiny2.binomial.10.notcollapsed.variants.csv")
write_csv(bin_collapsed_10, "./data/processed/shiny2.binomial.10.collapsed.variants.csv")

### One-Sided Data

# 100%
variants_with_meta <- left_join(variants_onesided_100, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_100 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_100 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_100, "./data/processed/shiny2.onesided.100.notcollapsed.variants.csv")
write_csv(one_collapsed_100, "./data/processed/shiny2.onesided.100.collapsed.variants.csv")

# 50%
variants_with_meta <- left_join(variants_onesided_50, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_50 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_50 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_50, "./data/processed/shiny2.onesided.50.notcollapsed.variants.csv")
write_csv(one_collapsed_50, "./data/processed/shiny2.onesided.50.collapsed.variants.csv")

# 25%
variants_with_meta <- left_join(variants_onesided_25, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_25 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_25 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_25, "./data/processed/shiny2.onesided.25.notcollapsed.variants.csv")
write_csv(one_collapsed_25, "./data/processed/shiny2.onesided.25.collapsed.variants.csv")

# 10%
variants_with_meta <- left_join(variants_onesided_10, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
one_notcollapsed_10 <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
one_collapsed_10 <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
write_csv(one_notcollapsed_10, "./data/processed/shiny2.onesided.10.notcollapsed.variants.csv")
write_csv(one_collapsed_10, "./data/processed/shiny2.onesided.10.collapsed.variants.csv")

### Two-Sided Data
variants_with_meta <- left_join(variants_twosided, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed

two_notcollapsed <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
two_collapsed <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))

# These don't have MapQ, Phred, Read_pos, or p.val. enter here just to see data. Why is the two-sided run so terrible for the benchmarking data?!
two_notcollapsed <- mutate(two_notcollapsed, MapQ = 44, Phred = 39, p.val = 0, Read_pos = 100)
two_collapsed <- mutate(two_collapsed, MapQ = 44, Phred = 39, p.val = 0, Read_pos = 100)

write_csv(two_notcollapsed, "./data/processed/shiny2.twosided.notcollapsed.variants.csv")
write_csv(two_collapsed, "./data/processed/shiny2.twosided.collapsed.variants.csv")

### Prepare coverage data ###

coverage_binomial_100 <- read.csv("./data/raw/shiny2.binomial.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_onesided_100 <- read.csv("./data/raw/shiny2.onesided.100.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_binomial_50 <- read.csv("./data/raw/shiny2.binomial.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_onesided_50 <- read.csv("./data/raw/shiny2.onesided.50.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_binomial_25 <- read.csv("./data/raw/shiny2.binomial.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_onesided_25 <- read.csv("./data/raw/shiny2.onesided.25.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

coverage_binomial_10 <- read.csv("./data/raw/shiny2.binomial.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))
coverage_onesided_10 <- read.csv("./data/raw/shiny2.onesided.10.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

#coverage_twosided <- read.csv("./data/raw/shiny2.twosided.coverage.csv", colClasses = c("character", "integer", "integer", "integer", "character"))

sampleInfo <- mutate(sampleInfo, Id = as.character(Id))
coverage_binomial_100_wmeta <- left_join(coverage_binomial_100, sampleInfo, by = "Id")
coverage_onesided_100_wmeta <- left_join(coverage_onesided_100, sampleInfo, by = "Id")
coverage_binomial_50_wmeta <- left_join(coverage_binomial_50, sampleInfo, by = "Id")
coverage_onesided_50_wmeta <- left_join(coverage_onesided_50, sampleInfo, by = "Id")
coverage_binomial_25_wmeta <- left_join(coverage_binomial_25, sampleInfo, by = "Id")
coverage_onesided_25_wmeta <- left_join(coverage_onesided_25, sampleInfo, by = "Id")
coverage_binomial_10_wmeta <- left_join(coverage_binomial_10, sampleInfo, by = "Id")
coverage_onesided_10_wmeta <- left_join(coverage_onesided_10, sampleInfo, by = "Id")
#coverage_twosided_wmeta <- left_join(coverage_twosided, sampleInfo, by = "Id")

write_csv(coverage_binomial_100_wmeta, "./data/processed/shiny2.binomial.100.coverage.csv")
write_csv(coverage_onesided_100_wmeta, "./data/processed/shiny2.onesided.100.coverage.csv")
write_csv(coverage_binomial_50_wmeta, "./data/processed/shiny2.binomial.50.coverage.csv")
write_csv(coverage_onesided_50_wmeta, "./data/processed/shiny2.onesided.50.coverage.csv")
write_csv(coverage_binomial_25_wmeta, "./data/processed/shiny2.binomial.25.coverage.csv")
write_csv(coverage_onesided_25_wmeta, "./data/processed/shiny2.onesided.25.coverage.csv")
write_csv(coverage_binomial_10_wmeta, "./data/processed/shiny2.binomial.10.coverage.csv")
write_csv(coverage_onesided_10_wmeta, "./data/processed/shiny2.onesided.10.coverage.csv")
#write_csv(coverage_twosided_wmeta, "./data/processed/shiny2.twosided.coverage.csv")

