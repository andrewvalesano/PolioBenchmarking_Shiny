

### Project: Polio Mixing Experiment 2.0
### Purpose: Collapse output variant CSV by sequencing duplicate.

# ======================== Import packages, read in metadata =============================

library(tidyverse)
source("./supportFunctions.R")

sampleInfo <- read.csv("./data/reference/MixingStudySamples2.csv")
primer_fwd_outer <- 95 # 95 to 115
primer_rev_outer <- 7440 # 7415 to 7440
expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants.csv")
expectedTruePos <- nrow(expectedVariants) # Expect to see 66 variants.
possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.


# ======================= Options for reading in variant data ====================

dups <- "both" # options: both, first, second
subset <- "100" # options: 100, 50, 25, 10
disp <- "onesided" # options: onesided or binomial
inputLevel <- "9_E3" # options: 4_E4, 9_E3, 9_E2, 9_E1
tna <- "yes" # options: yes, no

p.val.cutoff <- 0.9
MapQ.cutoff <- 0
Phred.cutoff <- 0
freq.var.cutoff <- 0
Read_pos_cutoff_1 <- 0
Read_pos_cutoff_2 <- 250

# ============================== Read in variant data ===============================

replicate = ifelse(dups == "first" | dups == "second", "notcollapsed", "collapsed")
file <- paste0("./data/processed/shiny2.", disp, ".", subset, ".", replicate, ".variants.csv")
variants_raw <- read.csv(file)

if(dups == "first") {
  data <- filter(variants_raw, Rep == 1)
} else if(dups == "second") {
  data <- filter(variants_raw, Rep == 2)
} else {
  data <- variants_raw
}

data <- filter(data, InputLevel == inputLevel, inTNA == tna)
data <- mutate(data, SampleNumber = Id)
data <- mutate(data, exp.freq = PercWT/100)
data <- mutate(data, Id = as.factor(as.character(PercWT)))
data <- filter(data, pos > primer_fwd_outer & pos < primer_rev_outer)
variant_data <- filter(data, p.val <= p.val.cutoff & 
                              MapQ >= MapQ.cutoff & 
                              Phred >= Phred.cutoff & 
                              freq.var >= freq.var.cutoff & 
                              Read_pos <= Read_pos_cutoff_2 & 
                              Read_pos >= Read_pos_cutoff_1)
variant_data <- mutate(variant_data, level = ifelse(PercWT == 1, "1_percent", 
                                                    ifelse(PercWT == 2, "2_percent", 
                                                           ifelse(PercWT == 5, "5_percent", 
                                                                  ifelse(PercWT == 10, "10_percent", 
                                                                         ifelse(PercWT == 100, "100_percent", NA))))))
variant_data <- mutate(variant_data, mutation_level = paste0(mutation, "_", level))

# ============== Read in test single data ============
# Does running the pipeline with different parameters make a difference? More false positives to weed out?

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


variants_single <- read.csv("data/raw/test.verysens.onesided.100.variants.csv")

expectedVariants <- read.csv("./data/reference/MixingStudyExpectedVariants.csv")
sampleInfo <- mutate(sampleInfo, Id = SampleNumber)
variants_with_meta <- left_join(variants_single, sampleInfo, by = "Id")
variants_with_meta <- mutate(variants_with_meta, group = paste0(inTNA, "_", InputLevel, "_", PercWT))
variants_with_meta %>% group_by(group) %>% do(collapse_localcov(.)) -> variants_collapsed
notcollapsed <- mutate(variants_with_meta, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))
collapsed <- mutate(variants_collapsed, category = ifelse(mutation %in% expectedVariants$mutation, TRUE, FALSE))

variant_data <- notcollapsed # data to use below

variant_data <- mutate(variant_data, SampleNumber = Id)
variant_data <- mutate(variant_data, exp.freq = PercWT/100)
variant_data <- mutate(variant_data, Id = as.factor(as.character(PercWT)))
variant_data <- filter(variant_data, pos > primer_fwd_outer & pos < primer_rev_outer)
variant_data <- mutate(variant_data, level = ifelse(PercWT == 1, "1_percent", 
                                                    ifelse(PercWT == 2, "2_percent", 
                                                           ifelse(PercWT == 5, "5_percent", 
                                                                  ifelse(PercWT == 10, "10_percent", 
                                                                         ifelse(PercWT == 100, "100_percent", NA))))))
variant_data <- mutate(variant_data, mutation_level = paste0(mutation, "_", level))

variant_data <- filter(variant_data, freq.var > 0.5 & inTNA == "yes", InputLevel == "9_E2" & !(ref == var) & PercWT != 100 & category == FALSE)

# ===================================== Read in coverage data ================================

file <- paste0("./data/processed/shiny2.", disp, ".", subset, ".coverage.csv")
coverage_raw <- read.csv(file)

# ============================ Plot coverage ==============================

coverage <- filter(coverage_raw, Id %in% variant_data$SampleNumber) # get the samples currently being analyzed in the variant data

base = 12000
interval = 300
linesize = 1.2
ggplot(coverage, aes(x = chr.pos, y = coverage, group = Id)) + geom_line(aes(color = Id)) + 
  theme_bw() + 
  xlab("Genome Position") + 
  ylab("Coverage (Read Depth)") + 
  theme(legend.position = "none") +
  xlim(c(95, 7440)) +
  geom_segment(data = coverage, size = linesize, color = "orange1", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
  geom_segment(data = coverage, size = linesize, color = "mediumseagreen", aes(x = 1472, y = base + interval*1, xend = 4151, yend = base + interval*1)) +
  geom_segment(data = coverage, size = linesize, color = "indianred4", aes(x = 3214, y = base + interval*2, xend = 5861, yend = base + interval*2)) +
  geom_segment(data = coverage, size = linesize, color = "slateblue3", aes(x = 4966, y = base + interval*3, xend = 7400, yend = base + interval*3)) +
  geom_abline(intercept = 200, slope = 0, linetype = 3, size = 0.5, color = "black") + 
  geom_abline(intercept = 1000, slope = 0, linetype = 3, size = 0.5, color = "black")


# ================================== Make table ===========================

dd = 4
m.roc.table <- miseq.roc.table(variant_data, 1, expectedTruePos, possible_vars, ">")
m.roc.table <- rename(m.roc.table, c("exp.freq"="Frequency","adj.sensitivity"="Sensitivity","TP"="True\nPositives","adj.specificity"="Specificity","FP"="False\nPositives"))
m.roc.table$Frequency <- c("100%", "10%", "5%", "2%", "1%")
m.roc.table$Sensitivity <- round(m.roc.table$Sensitivity, digits = dd)
m.roc.table$Specificity <- round(m.roc.table$Specificity, digits = dd)
m.roc.table

# ============================ Frequency by Position =======================


variant_data_TP <- filter(variant_data, category == TRUE & !is.na(level))

levels <- c("100_percent", "10_percent", "5_percent", "2_percent", "1_percent")
chrs <- data.frame("level" = levels)
chrs <- mutate(chrs, start = 0, stop = 7440)
palette <- wes_palette("Darjeeling1")
variant_data_TP$level <- factor(variant_data_TP$level, levels = rev(c("100_percent","10_percent","5_percent","2_percent","1_percent")))
chrs$level <- factor(chrs$level, levels = levels(variant_data_TP$level))

base = 0.6
interval = 0.02
linesize = 1.2

ggplot(variant_data_TP)  +
  geom_point(aes(x = pos, y = freq.var, color = level)) +
  ylab("Frequency") +
  xlab("Genome Position") +
  theme_minimal() +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  theme(panel.grid.major = element_line(colour = "gray96"), panel.grid.minor = element_line(colour = "white")) + 
  theme(legend.position = "right") +
  ylim(c(0, 0.7)) +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 7500)) +
  scale_color_manual(values = palette, name = "Expected Frequency") +
  geom_abline(intercept = 0.1, slope = 0, linetype = 3, size = 0.5, color = palette[4]) + 
  geom_abline(intercept = 0.05, slope = 0, linetype = 3, size = 0.5, color = palette[3]) + 
  geom_abline(intercept = 0.02, slope = 0, linetype = 3, size = 0.5, color = palette[2]) + 
  geom_abline(intercept = 0.01, slope = 0, linetype = 3, size = 0.5, color = palette[1]) +
  geom_segment(data = data, size = linesize, color = "orange1", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
  geom_segment(data = data, size = linesize, color = "mediumseagreen", aes(x = 1472, y = base + interval*1, xend = 4151, yend = base + interval*1)) +
  geom_segment(data = data, size = linesize, color = "indianred4", aes(x = 3214, y = base + interval*2, xend = 5861, yend = base + interval*2)) +
  geom_segment(data = data, size = linesize, color = "slateblue3", aes(x = 4966, y = base + interval*3, xend = 7400, yend = base + interval*3))

# =================================== Plot true positives and false negatives by position =====================

expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants_ByLevel.csv")
expectedVariants <- mutate(expectedVariants, mutation_level = paste0(mutation, "_", level))
expectedVariants <- mutate(expectedVariants, found = ifelse(mutation_level %in% variant_data$mutation_level, "True Positive", "False Negative"))

levels <- c("100_percent", "10_percent", "5_percent", "2_percent", "1_percent")
chrs <- data.frame("level" = levels)
chrs <- mutate(chrs, start = 0, stop = 7440)
palette <- wes_palette("Darjeeling1")
expectedVariants$level <- factor(expectedVariants$level, levels = rev(c("100_percent","10_percent","5_percent","2_percent","1_percent")))
chrs$level <- factor(chrs$level, levels = levels(expectedVariants$level))

ggplot(expectedVariants, aes(x = Position, y = level)) +
  geom_point(aes(color = found), size = 5, shape = 108) +
  geom_segment(data = chrs, aes(x = start, y = level, xend = stop, yend = level)) +
  ylab("Expected Frequency Group") +
  xlab("") +
  theme_minimal() +
  scale_color_manual(name = "", values = palette[c(1,2)]) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = c()) +
  theme(panel.grid.major = element_line(colour = "gray96"), panel.grid.minor = element_line(colour = "white")) + theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 7500))
