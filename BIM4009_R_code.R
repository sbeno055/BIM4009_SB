##### Packages #####

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(readr)
library(ggrepel)

##### Data Loading #####

# read in excel file
path <- "/Users/sarahcerine/Desktop/BIM_BLAST_GenExp.xlsx" 

#  chromosome sheets to load
sheets <- c(
  "Chr01_FPKM","Chr02_FPKM","Chr03_FPKM","Chr04_FPKM","Chr05_FPKM","Chr06_FPKM","Chr07_FPKM",
  "Chr08_FPKM","Chr09_FPKM","Chr10_FPKM","Chr11_FPKM","Chr12_FPKM","Chr13_FPKM","Chr14_FPKM",
  "Chr15_FPKM","Chr16_FPKM","Chr17_FPKM","Chr18_FPKM","Chr19_FPKM","Chr20_FPKM","Chr21_FPKM",
  "Chr22_FPKM", "ChrX_FPKM", "ChrY_FPKM"
)

# maps each SRR accession to its experimental condition 
sample_key <- tibble::tribble(
  ~sample,           ~hur_status, ~drug,
  "SRR28288466_1",      "Hur",       "Eltro",
  "SRR28288467_1",      "Hur",       "Eltro",
  "SRR31393158_1",      "Hur",       "Eltro", # repl3
  
  "SRR28288468_1",      "Hur",       "Ctrl",
  "SRR28288469_1",      "Hur",       "Ctrl",
  "SRR31393159_1",      "Hur",       "Ctrl",# repl3
  
  "SRR28288470_1",      "HurKO",     "Eltro",
  "SRR28288471_1",      "HurKO",     "Eltro",
  "SRR31393160_1",      "HurKO",     "Eltro", #repl3
  
  "SRR28288472_1",      "HurKO",     "Ctrl",
  "SRR28288473_1",      "HurKO",     "Ctrl",
  "SRR31393161_1",      "HurKO",     "Ctrl" # repl3
)

# replicate 3 samples needed scaling
rep3_samples <- c(
  "SRR31393158_1", "SRR31393159_1",
  "SRR31393160_1", "SRR31393161_1"
)


##### Load FPKM values #####

# reads one chromosome sheet and returns a table with gene, chromosome, and one column per sample
read_fpkm_sheet <- function(sheet_name) {
  raw <- read_excel(path, sheet = sheet_name, col_names = FALSE) %>%
    as_tibble()
  srr_row_idx <- which( # find the row that contains SRR accession IDs - header
    apply(raw, 1, function(row) {
      sum(str_detect(as.character(row), "^(SRR|ERR|DRR)\\d+"), na.rm = TRUE) >= 2
    })
  )[1]
  
  # build column names from the SRR header row
  col_names      <- as.character(unlist(raw[srr_row_idx, ]))
  col_names[1]   <- "gene"
  blank          <- is.na(col_names) | col_names == ""
  col_names[blank] <- paste0("col_", which(blank))  # fill any blank columns
  col_names      <- make.unique(col_names)
  
  # find the row labelled "Gene" since data starts on the next row
  data_start_idx <- which(tolower(as.character(raw[[1]])) == "gene")[1]
  
  data <- raw[(data_start_idx + 1):nrow(raw), ]
  names(data) <- col_names
  
  data %>%
    mutate(
      chromosome = sheet_name,
      across(-c(gene, chromosome), as.numeric)
    ) %>%
    select(gene, chromosome, everything()) %>%
    filter(!is.na(gene), gene != "")
}

# apply the function to all 24 chromosome sheets and stack into one table
raw_all <- map_dfr(sheets, read_fpkm_sheet)


##### Label #####

# pivot to long format so each row is one gene × one sample
fpkm_long <- raw_all %>%
  pivot_longer(
    cols      = -c(gene, chromosome),
    names_to  = "sample",
    values_to = "fpkm"
  ) %>%
  left_join(sample_key, by = "sample")  # attach HuR status and drug condition


##### Scale replicate 3 ######

# replicate 3 was sequenced separately, so its total library size differs from
# replicates 1 and 2 so we scaled it so all three replicates are comparable

# compute total FPKM per sample
sample_totals <- fpkm_long %>%
  group_by(hur_status, drug, sample) %>%
  summarise(total_fpkm = sum(fpkm, na.rm = TRUE), .groups = "drop") %>%
  mutate(is_rep3 = sample %in% rep3_samples)

# scale factor = mean of rep1 & rep2 total / rep3 total
# per condition
scale_factors <- sample_totals %>%
  group_by(hur_status, drug) %>%
  summarise(
    target_total = mean(total_fpkm[!is_rep3], na.rm = TRUE),  # average of rep1+rep2
    rep3_total   = total_fpkm[is_rep3][1],
    rep3_sample  = sample[is_rep3][1],
    scale_factor = target_total / rep3_total,
    .groups = "drop"
  )

# apply the scale factor only to rep3
fpkm_scaled <- fpkm_long %>%
  left_join(
    scale_factors %>% select(hur_status, drug, rep3_sample, scale_factor),
    by = c("hur_status", "drug")
  ) %>%
  mutate(
    fpkm_scaled = if_else(sample == rep3_sample, fpkm * scale_factor, fpkm)
  )


##### Average replicates #####

# average the three replicates into one mean FPKM value per gene per condition
mean_fpkm <- fpkm_scaled %>%
  group_by(chromosome, gene, hur_status, drug) %>%
  summarise(
    mean_fpkm = mean(fpkm_scaled, na.rm = TRUE),
    .groups   = "drop"
  )

# convert to wide format: one row per gene, one column per condition
fpkm_wide <- mean_fpkm %>%
  mutate(
    condition  = paste(hur_status, drug, sep = "_"),
    gene_locus = paste(gene, chromosome, sep = "_")
  ) %>%
  select(gene_locus, gene, chromosome, condition, mean_fpkm) %>%
  pivot_wider(names_from = condition, values_from = mean_fpkm)


##### Filter low-expression genes ######

# constants
PSEUDOCOUNT      <- 1    # added to all raw fpkm to avoid log(0) errors
MIN_FPKM         <- 5    # minimum mean FPKM required in all four conditions
OUTPUT_DIR       <- "/Users/sarahcerine/Desktop/"

# add pseudocount
# keep mean FPKM >= 5 in every condition so FC are reliable
fpkm_wide <- fpkm_wide %>%
  mutate(across(-c(gene_locus, gene, chromosome), ~ . + PSEUDOCOUNT)) %>%
  filter(if_all(-c(gene_locus, gene, chromosome), ~ . >= MIN_FPKM + PSEUDOCOUNT))

# 12,249 genes remain after filtering

##### Compute log2 fold changes #####

# four contrasts and interaction calculated 
fpkm_contrasts <- fpkm_wide %>%
  mutate(
    log2FC_Eltro_Hur    = log2(Hur_Eltro    / Hur_Ctrl), #HuR KO,  untreated
    log2FC_Eltro_HurKO  = log2(HurKO_Eltro  / HurKO_Ctrl), #HuR KO,  Eltrombopag
    log2FC_Hur_Ctrl     = log2(Hur_Ctrl     / HurKO_Ctrl), #HuR WT, untreated
    log2FC_Hur_Eltro    = log2(Hur_Eltro    / HurKO_Eltro), #HuR WT, Eltrombopag
    # interaction: difference in drug effect depending on HuR status
    Interaction         = log2FC_Eltro_Hur - log2FC_Eltro_HurKO
  )

##### Classify ARE status ######

# ARE+ genes contain AU-rich elements (any of 5 clusters) from the ARED database
ared_raw  <- read_excel("/Users/sarahcerine/Desktop/Results.xlsx")
are_genes_db <- ared_raw %>% pull("Gene Symbol") %>% unique()

# classification
fpkm_contrasts <- fpkm_contrasts %>%
  mutate(
    ARE_status = case_when(
      gene %in% are_genes_db ~ "ARE+",
      TRUE                   ~ "ARE-"
    )
  )

##### Statistical tests ######

# Kolmogorov test: testing normality 
# ARE+ under control
ks_are_plus_ctrl <- ks.test(
  fpkm_contrasts %>% filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Ctrl),
  "pnorm",
  mean = mean(fpkm_contrasts %>% filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Ctrl), na.rm = TRUE),
  sd   = sd(fpkm_contrasts %>%   filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Ctrl), na.rm = TRUE)
)
cat("ARE+ log2FC_Hur_Ctrl:  D =", round(ks_are_plus_ctrl$statistic, 4),
    "  p =", signif(ks_are_plus_ctrl$p.value, 3), "\n")

# ARE- under control
ks_are_minus_ctrl <- ks.test(
  fpkm_contrasts %>% filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Ctrl),
  "pnorm",
  mean = mean(fpkm_contrasts %>% filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Ctrl), na.rm = TRUE),
  sd   = sd(fpkm_contrasts %>%   filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Ctrl), na.rm = TRUE)
)
cat("ARE- log2FC_Hur_Ctrl:  D =", round(ks_are_minus_ctrl$statistic, 4),
    "  p =", signif(ks_are_minus_ctrl$p.value, 3), "\n")

# ARE+ under Eltrombopag
ks_are_plus_eltro <- ks.test(
  fpkm_contrasts %>% filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Eltro),
  "pnorm",
  mean = mean(fpkm_contrasts %>% filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Eltro), na.rm = TRUE),
  sd   = sd(fpkm_contrasts %>%   filter(ARE_status == "ARE+") %>% pull(log2FC_Hur_Eltro), na.rm = TRUE)
)
cat("ARE+ log2FC_Hur_Eltro: D =", round(ks_are_plus_eltro$statistic, 4),
    "  p =", signif(ks_are_plus_eltro$p.value, 3), "\n")

# ARE- under Eltrombopag
ks_are_minus_eltro <- ks.test(
  fpkm_contrasts %>% filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Eltro),
  "pnorm",
  mean = mean(fpkm_contrasts %>% filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Eltro), na.rm = TRUE),
  sd   = sd(fpkm_contrasts %>%   filter(ARE_status == "ARE-") %>% pull(log2FC_Hur_Eltro), na.rm = TRUE)
)
cat("ARE- log2FC_Hur_Eltro: D =", round(ks_are_minus_eltro$statistic, 4),
    "  p =", signif(ks_are_minus_eltro$p.value, 3), "\n")

# Wilcoxon rank-sum test: drug effect distribution in HuR WT vs HuR KO (two independent groups, non-parametric)
mwu_drug <- wilcox.test(
  fpkm_contrasts$log2FC_Eltro_Hur,
  fpkm_contrasts$log2FC_Eltro_HurKO,
  paired      = FALSE,
  alternative = "two.sided"
)

# paired Wilcoxon test: HuR effect under control vs Eltrombopag (same set of loci measured under two conditions)
paired_hur <- wilcox.test(
  fpkm_contrasts$log2FC_Hur_Ctrl,
  fpkm_contrasts$log2FC_Hur_Eltro,
  paired      = TRUE,
  alternative = "two.sided"
)

# Wilcoxon rank-sum test: HuR effect on ARE+ vs ARE- under each drug condition
are_test_ctrl  <- wilcox.test(log2FC_Hur_Ctrl  ~ ARE_status, data = fpkm_contrasts)
are_test_eltro <- wilcox.test(log2FC_Hur_Eltro ~ ARE_status, data = fpkm_contrasts)

p_ctrl  <- signif(are_test_ctrl$p.value,  3)
p_eltro <- signif(are_test_eltro$p.value, 3)


##### Figure 3 — ARE+ vs ARE– boxplots #####

# side-by-side plotting
fig3_data <- fpkm_contrasts %>%
  select(gene_locus, ARE_status, log2FC_Hur_Ctrl, log2FC_Hur_Eltro) %>%
  pivot_longer(
    cols      = c(log2FC_Hur_Ctrl, log2FC_Hur_Eltro),
    names_to  = "condition",
    values_to = "log2FC_Hur"
  ) %>%
  mutate(
    condition = recode(condition,
                       "log2FC_Hur_Ctrl"  = "Control",
                       "log2FC_Hur_Eltro" = "Eltrombopag"),
    condition = factor(condition, levels = c("Control", "Eltrombopag"))
  )

# 150 random points per group for jittering
set.seed(42)
fig3_jitter <- fig3_data %>%
  group_by(condition, ARE_status) %>%
  slice_sample(n = 150) %>%
  ungroup()

# panel A — side-by-side boxplots with p-values annotated above each condition
p3A <- fig3_data %>%
  ggplot(aes(x = condition, y = log2FC_Hur, fill = ARE_status)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.7, width = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    data     = fig3_jitter,
    aes(x = condition, y = log2FC_Hur, fill = ARE_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 0.6, alpha = 0.35, color = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(ylim = c(-4, 4), clip = "off") +
  annotate("text", x = 1, y = 3.7, label = paste0("p = ", p_ctrl),  size = 3.5, color = "grey30") +
  annotate("text", x = 2, y = 3.7, label = paste0("p = ", p_eltro), size = 3.5, color = "grey30") +
  scale_y_continuous(breaks = seq(-4, 4, by = 1)) +
  scale_fill_manual(values = c("ARE+" = "#F4A899", "ARE-" = "#90C9A0")) +
  labs(
    title = "HuR-dependent transcript abundance: ARE+ vs ARE\u2212",
    x     = "Drug Condition",
    y     = "log2 Fold Change (HuR WT / HuR KO)",
    fill  = "ARE Status"
  ) +
  theme_bw()

# figure b — summary statistics table
table_fig3B <- fpkm_contrasts %>%
  group_by(ARE_status) %>%
  summarise(
    `Mean (Control)`       = round(mean(log2FC_Hur_Ctrl,    na.rm = TRUE), 3),
    `Median (Control)`     = round(median(log2FC_Hur_Ctrl,  na.rm = TRUE), 3),
    `SD (Control)`         = round(sd(log2FC_Hur_Ctrl,      na.rm = TRUE), 3),
    `Mean (Eltrombopag)`   = round(mean(log2FC_Hur_Eltro,   na.rm = TRUE), 3),
    `Median (Eltrombopag)` = round(median(log2FC_Hur_Eltro, na.rm = TRUE), 3),
    `SD (Eltrombopag)`     = round(sd(log2FC_Hur_Eltro,     na.rm = TRUE), 3),
    `n loci`               = n()
  ) 

##### Figure 4 - distribution plots #####

# A — density plot of the Eltrombopag drug effect in HuR WT vs HuR KO
fig4A_data <- fpkm_contrasts %>%
  select(gene_locus, log2FC_Eltro_Hur, log2FC_Eltro_HurKO) %>%
  pivot_longer(-gene_locus, names_to = "condition", values_to = "log2FC") %>%
  mutate(condition = recode(condition,
                            "log2FC_Eltro_Hur"   = "HuR WT",
                            "log2FC_Eltro_HurKO" = "HuR KO"))

p4A <- fig4A_data %>%
  ggplot(aes(x = log2FC, fill = condition, color = condition)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  coord_cartesian(xlim = c(-6, 6)) +
  annotate("text", x = 5, y = 0.65, label = "p < 0.05", size = 3.5, color = "grey30") +
  scale_fill_manual(values  = c("HuR WT" = "#808080", "HuR KO" = "#F94449")) +
  scale_color_manual(values = c("HuR WT" = "#808080", "HuR KO" = "#F94449")) +
  labs(
    title = "Distribution of Eltrombopag drug effect",
    x     = "log2 Fold Change (Eltrombopag / Control)",
    y     = "Density",
    fill  = "HuR Status",
    color = "HuR Status"
  ) +
  theme_classic(base_size = 13)

# B — density plot of the HuR effect under control vs Eltrombopag
fig4B_data <- fpkm_contrasts %>%
  select(gene_locus, log2FC_Hur_Ctrl, log2FC_Hur_Eltro) %>%
  pivot_longer(-gene_locus, names_to = "condition", values_to = "log2FC") %>%
  mutate(condition = recode(condition,
                            "log2FC_Hur_Ctrl"  = "Control",
                            "log2FC_Hur_Eltro" = "Eltrombopag"))

p4B <- fig4B_data %>%
  ggplot(aes(x = log2FC, fill = condition, color = condition)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  coord_cartesian(xlim = c(-6, 6)) +
  annotate("text", x = 5, y = 0.9, label = "p < 0.05", size = 3.5, color = "grey30") +
  scale_fill_manual(values  = c("Control" = "#F4A899", "Eltrombopag" = "#5A9E6F")) +
  scale_color_manual(values = c("Control" = "#F4A899", "Eltrombopag" = "#5A9E6F")) +
  labs(
    title = "Distribution of HuR effect",
    x     = "log2 Fold Change (HuR WT / HuR KO)",
    y     = "Density",
    fill  = "Drug Condition",
    color = "Drug Condition"
  ) +
  theme_classic(base_size = 13)

# C — summary statistics table for Figure 4 distributions
table_fig4C <- tibble(
  Comparison   = c("Eltrombopag effect", "", "HuR effect", ""),
  Distribution = c("HuR WT", "HuR KO", "Control", "Eltrombopag"),
  `Mean (x̄)`  = c(
    round(mean(fpkm_contrasts$log2FC_Eltro_Hur,   na.rm = TRUE), 3),
    round(mean(fpkm_contrasts$log2FC_Eltro_HurKO, na.rm = TRUE), 3),
    round(mean(fpkm_contrasts$log2FC_Hur_Ctrl,    na.rm = TRUE), 3),
    round(mean(fpkm_contrasts$log2FC_Hur_Eltro,   na.rm = TRUE), 3)
  ),
  `Median (m)` = c(
    round(median(fpkm_contrasts$log2FC_Eltro_Hur,   na.rm = TRUE), 3),
    round(median(fpkm_contrasts$log2FC_Eltro_HurKO, na.rm = TRUE), 3),
    round(median(fpkm_contrasts$log2FC_Hur_Ctrl,    na.rm = TRUE), 3),
    round(median(fpkm_contrasts$log2FC_Hur_Eltro,   na.rm = TRUE), 3)
  ),
  SD = c(
    round(sd(fpkm_contrasts$log2FC_Eltro_Hur,   na.rm = TRUE), 3),
    round(sd(fpkm_contrasts$log2FC_Eltro_HurKO, na.rm = TRUE), 3),
    round(sd(fpkm_contrasts$log2FC_Hur_Ctrl,    na.rm = TRUE), 3),
    round(sd(fpkm_contrasts$log2FC_Hur_Eltro,   na.rm = TRUE), 3)
  ),
  Test = c("Wilcoxon rank-sum", "", "Paired Wilcoxon", "")
)

# D — gene counts (upregulation/downregulation) for each contrast
table_fig4D <- fpkm_contrasts %>%
  summarise(
    `Total loci`                    = n(),
    `Drug effect HuR WT — up`       = sum(log2FC_Eltro_Hur   > 0, na.rm = TRUE),
    `Drug effect HuR WT — down`     = sum(log2FC_Eltro_Hur   < 0, na.rm = TRUE),
    `Drug effect HuR KO — up`       = sum(log2FC_Eltro_HurKO > 0, na.rm = TRUE),
    `Drug effect HuR KO — down`     = sum(log2FC_Eltro_HurKO < 0, na.rm = TRUE),
    `HuR effect Control — up`       = sum(log2FC_Hur_Ctrl    > 0, na.rm = TRUE),
    `HuR effect Control — down`     = sum(log2FC_Hur_Ctrl    < 0, na.rm = TRUE),
    `HuR effect Eltrombopag — up`   = sum(log2FC_Hur_Eltro   > 0, na.rm = TRUE),
    `HuR effect Eltrombopag — down` = sum(log2FC_Hur_Eltro   < 0, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "Contrast", values_to = "Number of Loci")

##### Figure 5 — dyad interaction  #####

# rank all genes by the absolute size of their interaction score
interaction_ranked <- fpkm_contrasts %>%
  arrange(desc(abs(Interaction)))

# labelling each gene by the pattern of its drug response in WT and KO cells
interaction_ranked <- interaction_ranked %>%
  mutate(
    Direction = case_when(
      log2FC_Eltro_Hur > 0  & log2FC_Eltro_HurKO < 0  ~ "(+/-)",  # up in WT, down in KO
      log2FC_Eltro_Hur < 0  & log2FC_Eltro_HurKO > 0  ~ "(-/+)",  # down in WT, up in KO
      log2FC_Eltro_Hur > 0  & log2FC_Eltro_HurKO >= 0 ~ "(+/+) WT", # stronger in WT
      log2FC_Eltro_Hur < 0  & log2FC_Eltro_HurKO <= 0 ~ "(-/-) WT", # stronger in WT
      log2FC_Eltro_Hur >= 0 & log2FC_Eltro_HurKO > 0  ~ "(+/+) KO", # stronger in KO
      log2FC_Eltro_Hur <= 0 & log2FC_Eltro_HurKO < 0  ~ "(-/-) KO", # stronger in KO
      TRUE ~ "Other"
    )
  )

top20 <- interaction_ranked %>% slice_max(abs(Interaction), n = 20)

# A — scatter plot of drug effect in HuR WT vs HuR KO for the top 20 genes
p5A <- top20 %>%
  ggplot(aes(x = log2FC_Eltro_HurKO, y = log2FC_Eltro_Hur, color = ARE_status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_label_repel(
    aes(label = gene),
    size = 3.2, fontface = "italic",
    box.padding = 0.5, max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("ARE+" = "#F4A899", "ARE-" = "#90C9A0")) +
  labs(
    title = "HuR-dependent Eltrombopag response for top 20 interaction genes",
    x     = "log2 Eltro effect in HuR KO cells",
    y     = "log2 Eltro effect in HuR WT cells",
    color = "ARE Status"
  ) +
  theme_bw(base_size = 13)

# B — table of top 20 interaction genes with all contrast values
table_fig5B <- top20 %>%
  select(
    Gene              = gene,
    Chromosome        = chromosome,
    `ARE Status`      = ARE_status,
    `Drug_HuR WT`     = log2FC_Eltro_Hur,
    `Drug_HuR KO`     = log2FC_Eltro_HurKO,
    `HuR_Ctrl`        = log2FC_Hur_Ctrl,
    `HuR_drug`        = log2FC_Hur_Eltro,
    `Interaction Score` = Interaction,
    Direction
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

