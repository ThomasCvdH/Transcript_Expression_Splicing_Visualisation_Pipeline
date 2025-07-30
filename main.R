# --- R Packages ---
# Install required packages if not already installed
required_packages <- c(
  "dplyr", "tidyr", "stringr", "ggplot2", "tibble",
  "ggtranscript", "rtracklayer", "patchwork", "svglite"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(vegan)
remotes::install_github("dzhang32/ggtranscript")
library(ggtranscript)
library(rtracklayer)
library(patchwork)

# --- Load Custom Functions ---
source("functions.R")

# --- 1. Set Dataset Context ---
dataset <- "timeseries" # Options: "adam", "coolLL2", "timeseries"

# --- 2. Set up output directory for figures ---
figures_dir <- file.path("figures", dataset)
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir)

# --- 3. Load Data Based on Context ---
data_list <- switch(dataset,
  "adam" = parse_data("Transcript_TPM_adam_TRIMMED.csv"),
  "coolLL2" = parse_data("Transcript TPMcoolLL.csv"),
  "timeseries" = parse_data("Transcript_TPM_timeSeries_trimmed.csv")
)
df <- data_list$df
averaged_df <- data_list$averaged_df
print(averaged_df)
# --- 4. GTF Import (shared) ---
gtf_path <- file.path("ALL_combined_annotations.gtf")
gtf_df <- rtracklayer::import(gtf_path) %>% as_tibble()


# --- 5. Usage Examples ---

# Generate and print a dashboard for a single gene
rve2_dashboard <- plot_expression_dashboard(
  gene_id = "AT5G37260",
  gene_name = "RVE2",
  gtf_data = gtf_df,
  raw_df = df,
  avg_df = averaged_df,
  dataset = dataset
)
print(rve2_dashboard)

rve8_dashboard <- plot_expression_dashboard(
  gene_id = "AT3G09600",
  gene_name = "RVE8",
  gtf_data = gtf_df,
  raw_df = df,
  avg_df = averaged_df,
  dataset = dataset
)
print(rve8_dashboard)

# Generate and print a shannon diversity comparison for a single gene
comparison_plot <- plot_diversity_comparison(
  tair_code = "AT3G09600",
  gene_name = "RVE8",
  time_h = 72,
  raw_df = df,
  dataset = dataset # Does not work with timeseries dataset
)
print(comparison_plot)

# Batch process list of genes for dashboards, grids, family grids.

genes_to_plot <- tibble::tribble(
  ~gene_name, ~tair_code,
  "RVE1", "AT5G17300 ",
  "RVE2", "AT5G37260",
  "RVE3", "AT1G01520",
  "RVE4", "AT5G02840",
  "RVE5", "AT4G01280",
  "RVE6", "AT5G52660",
  "RVE7", "AT1G18330",
  "RVE8", "AT3G09600",
  "LHY", "AT1G01060",
  "CCA1", "AT2G46830",
  "PHYB", "AT2G18790",
  "PRR9", "AT2G46790",
  "PRR7", "AT5G02810",
  "PRR5", "AT5G24470",
  "TOC1", "AT5G61380",
  "ZTL", "AT5G57360",
  "GI", "AT1G22770",
  "BOA", "AT5G59570",
  "LUX", "AT3G46640",
  "ELF3", "AT2G25930",
  "ELF4", "AT2G40080",
  "COP1", "AT2G32950",
  "LNK1", "AT5G64170",
  "LNK2", "AT3G54500",
  "LNK3", "AT3G12320",
  "LNK4", "AT5G06980",
  "LWD1", "AT1G12910",
  "LWD2", "AT3G26640",
  "TTG1", "AT5G24520",
  "COR27", "AT5G42900",
  "COR28", "AT4G33980",
  "PIF1", "AT2G20180",
  "PIF3", "AT1G09530",
  "PIF4", "AT2G43010",
  "PIF5", "AT3G59060",
  "TCP2", "AT3G27010",
  "TCP21", "AT5G08330",
  "TCP22", "AT1G72010",
  "PHYB", "AT2G18790"
)
# Dashboard batch
for (i in 1:nrow(genes_to_plot)) {
  current_gene_name <- genes_to_plot$gene_name[i]
  current_gene_id <- str_trim(genes_to_plot$tair_code[i])

  tryCatch(
    {
      cat(paste("--- Generating dashboard for", current_gene_name, "---\n"))
      dashboard_plot <- plot_expression_dashboard(
        gene_id = current_gene_id,
        gene_name = current_gene_name,
        gtf_data = gtf_df,
        raw_df = df,
        avg_df = averaged_df,
        dataset = dataset # Pass dataset
      )
    },
    error = function(e) {
      cat(paste(
        "Could not generate dashboard for",
        current_gene_name, ":",
        e$message, "\n\n"
      ))
    }
  )
}

# SD compare batch plot
for (i in 1:nrow(genes_to_plot)) {
  current_gene_name <- genes_to_plot$gene_name[i]
  current_gene_id <- str_trim(genes_to_plot$tair_code[i])

  tryCatch(
    {
      cat(paste("--- Generating SD comparison for", current_gene_name, "---\n"))
      comparison_plot <- plot_diversity_comparison(
        tair_code = current_gene_id,
        gene_name = current_gene_name,
        time_h = 16,
        raw_df = df,
        dataset = dataset # Does not work with timeseries dataset
      )
    },
    error = function(e) {
      cat(paste(
        "Could not generate SD comparison plot for",
        current_gene_name, ":",
        e$message, "\n\n"
      ))
    }
  )
}
# Grid in batch by family
# --- RVE Gene Family ---
rve_genes <- genes_to_plot %>% filter(str_starts(gene_name, "RVE"))
tryCatch(
  {
    cat("--- Generating RVE Family Grid Plots ---")
    total_expression_plot_rve <- plot_total_expression_grid(
      df,
      averaged_df,
      rve_genes,
      dataset
    )
    transcript_plot_rve <- plot_individual_transcripts_grid(
      averaged_df,
      rve_genes,
      dataset
    )
    ggsave(
      file.path(figures_dir, paste0(
        "total_expression_grid_RVE_family_",
        toupper(dataset),
        ".svg"
      )),
      plot = total_expression_plot_rve,
      width = 24,
      height = 10
    )
    ggsave(
      file.path(figures_dir, paste0(
        "individual_transcripts_grid_RVE_family_",
        toupper(dataset),
        ".svg"
      )),
      plot = transcript_plot_rve,
      width = 24,
      height = 12
    )
    message("RVE family grid plots saved successfully.")
  },
  error = function(e) {
    cat("Error generating RVE grid plots:", e$message, "\n")
  }
)

# --- PRR Gene Family ---
prr_genes <- genes_to_plot %>% filter(str_starts(gene_name, "PRR"))
tryCatch(
  {
    cat("--- Generating PRR Family Grid Plots ---")
    total_expression_plot_prr <- plot_total_expression_grid(
      df,
      averaged_df,
      prr_genes,
      dataset
    )
    transcript_plot_prr <- plot_individual_transcripts_grid(
      averaged_df,
      prr_genes,
      dataset
    )
    ggsave(
      file.path(figures_dir, paste0(
        "total_expression_grid_PRR_family_",
        toupper(dataset),
        ".svg"
      )),
      plot = total_expression_plot_prr,
      width = 24,
      height = 10
    )
    ggsave(
      file.path(figures_dir, paste0(
        "individual_transcripts_grid_PRR_family_",
        toupper(dataset),
        ".svg"
      )),
      plot = transcript_plot_prr,
      width = 15,
      height = 5
    )
    message("PRR family grid plots saved successfully.")
  },
  error = function(e) {
    cat("Error generating PRR grid plots:", e$message, "\n")
  }
)
