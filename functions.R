# --- Utility: Dataset-Specific Configurations ---

# Get Regex Patterns for Parsing Sample Names and Gene IDs
get_regex_patterns <- function(dataset) {
  # Define base patterns that are common unless specified otherwise
  base_patterns <- list(
    tair_code = "AT.G\\d+",
    transcript_suffix = "\\.\\d+$"
  )

  if (dataset == "adam") {
    c(list(
      condition_group = "^X\\d+",
      time_h = "(?<=\\.X)(\\d+)(?=H)",
      rep = "\\.X\\d+$",
      pattern_type = "direct"
    ), base_patterns)
  } else if (dataset == "coolLL2") {
    c(list(
      condition_group = "^X\\d+C",
      time_h = "(?<=\\.X)(\\d+)(?=h)",
      rep = "\\.rep\\d+$",
      pattern_type = "direct"
    ), base_patterns)
  } else if (dataset == "timeseries") {
    c(list(
      condition_group = "X\\d+",
      timepoint = "(?<=T)\\d+",
      rep = "\\.rep\\d+$",
      pattern_type = "mapped"
    ), base_patterns)
  } else {
    stop("Unknown dataset for regex patterns")
  }
}

# Get Time Mapping for "timeseries" Dataset
get_time_mapping <- function() {
  timepoints <- 1:26
  hours <- numeric(length(timepoints))

  hours[1] <- 3
  for (t in 2:length(timepoints)) {
    if (t <= 11) {
      hours[t] <- hours[t - 1] + 3
    } else if (t %in% 12:13) {
      hours[t] <- hours[t - 1] + 1.5
    } else if (t > 13 && t <= 18) {
      hours[t] <- hours[t - 1] + 3
    } else if (t == 19) {
      hours[t] <- hours[t - 1] + 48
    } else {
      hours[t] <- hours[t - 1] + 3
    }
  }

  tibble::tibble(
    timepoint = as.character(timepoints),
    time_h = hours
  )
}

# Get Plotting Aesthetics (Shapes, Colors, Linetypes)
get_plot_aesthetics <- function(dataset) {
  if (dataset == "adam") {
    list(
      labels = c("15°C", "20°C"),
      values = c("X15" = "#377EB8", "X20" = "#E41A1C"),
      shapes = c("X15" = 21, "X20" = 24),
      linetypes = c("X15" = "solid", "X20" = "dashed")
    )
  } else if (dataset == "coolLL2") {
    list(
      labels = c("12°C", "20°C"),
      values = c("X12C" = "#377EB8", "X20C" = "#E41A1C"),
      shapes = c("X12C" = 21, "X20C" = 24),
      linetypes = c("X12C" = "solid", "X20C" = "dashed")
    )
  } else if (dataset == "timeseries") {
    list(
      labels = c("4°C", "20°C"),
      values = c("X4" = "#377EB8", "X20" = "#E41A1C"),
      shapes = c("X4" = 21, "X20" = 24),
      linetypes = c("X4" = "solid", "X20" = "dashed")
    )
  } else {
    stop("Unknown dataset for plot aesthetics")
  }
}

# --- Data Preparation Functions ---
# Generic parsing function to reduce redundancy
parse_data <- function(file_path) {
  regex <- get_regex_patterns(dataset)
  df <- read.csv(file_path, row.names = 1, check.names = FALSE)
  averaged_df <- df %>%
    tibble::rownames_to_column(var = "gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "value") %>%
    mutate(condition = str_remove(sample, regex$rep)) %>%
    group_by(gene, condition) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = "condition", values_from = "mean_value")
  list(df = df, averaged_df = averaged_df)
}

# --- Plotting Functions ---

# Utility function to extract time (hours) from condition/sample string
extract_time <- function(x, regex_patterns, time_map = NULL) {
  sapply(x, function(val) {
    if (!is.null(time_map)) { # Mapped time (e.g., timeseries)
      tp <- str_extract(val, regex_patterns$timepoint)
      if (is.na(tp)) {
        return(NA_real_)
      }
      time_row <- time_map %>% filter(timepoint == tp)
      if (nrow(time_row) == 0) {
        return(NA_real_)
      }
      return(time_row$time_h)
    } else { # Direct time extraction
      as.numeric(str_extract(val, regex_patterns$time_h))
    }
  })
}


# Plot Total Gene Expression with Error Bars
plot_gene_expression <- function(tair_code, gene_name, raw_df, avg_df, dataset) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  if (!("gene" %in% colnames(raw_df))) raw_df <- tibble::rownames_to_column(raw_df, var = "gene")

  plot_data <- raw_df %>%
    dplyr::filter(str_starts(gene, tair_code)) %>%
    pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "abundance") %>%
    mutate(
      condition_group = str_extract(sample, regex$condition_group),
      time_h = extract_time(sample, regex, time_map)
    ) %>%
    group_by(sample, condition_group, time_h) %>%
    summarise(total_abundance_per_rep = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(condition_group, time_h) %>%
    summarise(
      mean_abundance = mean(total_abundance_per_rep),
      sd_abundance = sd(total_abundance_per_rep),
      .groups = "drop"
    ) %>%
    filter(!is.na(time_h))
  print(plot_data)
  ggplot(plot_data, aes(x = time_h, y = mean_abundance, group = condition_group, color = condition_group)) +
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.4) +
    geom_line(aes(linetype = condition_group), linewidth = 1.2) +
    geom_point(aes(shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
    scale_color_manual(name = "Condition", values = aesthetics$values, labels = aesthetics$labels) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Total Expression of", gene_name),
      subtitle = paste("TAIR ID:", tair_code),
      x = "Time (Hours)", y = "Mean Expression (TPM)"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, color = "gray40"), legend.position = "bottom")
}

# Plot Top 3 Individual Gene Transcripts
plot_individual_transcripts <- function(tair_code, gene_name, avg_df, dataset) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  plot_data_full <- avg_df %>%
    dplyr::filter(str_starts(gene, tair_code)) %>%
    pivot_longer(cols = where(is.numeric), names_to = "condition", values_to = "mean_value") %>%
    mutate(
      condition_group = str_extract(condition, regex$condition_group),
      time_h = extract_time(condition, regex, time_map),
      transcript_id = str_extract(gene, regex$transcript_suffix)
    ) %>%
    filter(!is.na(time_h) & !is.na(transcript_id))

  top_transcripts <- plot_data_full %>%
    group_by(transcript_id) %>%
    summarise(overall_mean = mean(mean_value, na.rm = TRUE)) %>%
    slice_max(order_by = overall_mean, n = 3) %>%
    pull(transcript_id)

  plot_data <- plot_data_full %>% dplyr::filter(transcript_id %in% top_transcripts)

  if (nrow(plot_data) == 0) {
    warning(paste("No transcript data to plot for", gene_name))
    return(ggplot() +
      theme_void() +
      labs(title = paste(gene_name, "\n(No transcript data)")))
  }

  ggplot(plot_data, aes(x = time_h, y = mean_value, group = interaction(transcript_id, condition_group))) +
    geom_line(aes(color = transcript_id, linetype = condition_group), linewidth = 1) +
    geom_point(aes(color = transcript_id, shape = condition_group), fill = "white", size = 3.5, alpha = 0.8) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Top 3", gene_name, "Transcripts"),
      subtitle = paste("TAIR ID:", tair_code),
      x = "Time (Hours)", y = "Mean Expression (TPM)", color = "Transcript ID"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, color = "gray40"), legend.position = "bottom", legend.box = "vertical")
}

# Plot Transcript Structures and a Comparison of the Top Two
plot_gene_structure_and_comparison <- function(gene_id, gene_name, gtf_data, avg_df, min_tpm_threshold = 1) {
  expressed_transcripts_df <- avg_df %>%
    dplyr::filter(str_starts(gene, gene_id)) %>%
    rowwise() %>%
    mutate(overall_mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::filter(overall_mean > min_tpm_threshold)

  expressed_transcripts_list <- expressed_transcripts_df %>% pull(gene)

  if (length(expressed_transcripts_list) == 0) {
    warning(paste("No transcripts with expression >", min_tpm_threshold, "TPM found for", gene_name))
    return(NULL)
  }

  expressed_gene_annotation <- gtf_data %>%
    dplyr::filter(!is.na(.data$gene_id) & .data$gene_id == !!gene_id) %>%
    dplyr::filter(transcript_id %in% expressed_transcripts_list)

  if (nrow(expressed_gene_annotation) == 0) {
    warning(paste("Expressed transcripts for", gene_name, "not found in GTF file."))
    return(NULL)
  }

  gene_exons <- expressed_gene_annotation %>% dplyr::filter(type == "exon")
  gene_cds <- expressed_gene_annotation %>% dplyr::filter(type == "CDS")
  gene_myb <- expressed_gene_annotation %>% dplyr::filter(type == "mybdomain_part")

  p1 <- ggplot(gene_exons, aes(xstart = start, xend = end, y = transcript_id)) +
    geom_intron(data = to_intron(gene_exons, "transcript_id"), aes(strand = strand)) +
    geom_range(aes(fill = "Exon"), height = 0.25) +
    geom_range(data = gene_cds, aes(fill = "CDS"), height = 0.5) +
    geom_range(data = gene_myb, aes(fill = "MYB Domain"), height = 0.4, color = "black", linewidth = 0.001) +
    scale_fill_manual(
      name = "Feature",
      values = c("Exon" = "white", "CDS" = "royalblue", "MYB Domain" = "coral"),
      breaks = c("Exon", "CDS", "MYB Domain")
    ) +
    labs(subtitle = "Expressed transcript structures", y = "Transcript ID") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  top_two_transcripts <- expressed_transcripts_df %>%
    slice_max(order_by = overall_mean, n = 2) %>%
    pull(gene)

  if (length(top_two_transcripts) < 2) {
    warning(paste("Fewer than two expressed transcripts found for", gene_name, ". Cannot create comparison plot."))
    return(
      p1 + plot_annotation(
        title = paste("Transcript Structure for", gene_name),
        subtitle = paste("TAIR ID:", gene_id),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12)
        )
      )
    )
  }

  transcript1_id <- top_two_transcripts[1]
  transcript2_id <- top_two_transcripts[2]

  t1_exons <- gene_exons %>% dplyr::filter(transcript_id == transcript1_id)
  t1_cds <- gene_cds %>% dplyr::filter(transcript_id == transcript1_id)
  t1_myb <- gene_myb %>% dplyr::filter(transcript_id == transcript1_id)

  t2_exons <- gene_exons %>% dplyr::filter(transcript_id == transcript2_id)
  t2_cds <- gene_cds %>% dplyr::filter(transcript_id == transcript2_id)
  t2_myb <- gene_myb %>% dplyr::filter(transcript_id == transcript2_id)

  y_axis_label <- paste(transcript1_id, "\n(Top)", "/", transcript2_id, "\n(Bottom)")

  p2 <- ggplot(mapping = aes(y = y_axis_label)) +
    geom_half_range(data = t2_exons, aes(xstart = start, xend = end), fill = "white", height = 0.2, color = "black") +
    geom_intron(data = to_intron(t2_exons, "transcript_id"), aes(xstart = start, xend = end, strand = strand), arrow.min.intron.length = 300) +
    geom_half_range(data = t2_cds, aes(xstart = start, xend = end), fill = "grey40", height = 0.2) +
    geom_half_range(data = t2_myb, aes(xstart = start, xend = end), fill = "coral", height = 0.19, color = "black", linewidth = 0.001) +
    geom_half_range(data = t1_exons, aes(xstart = start, xend = end), range.orientation = "top", fill = "white", height = 0.2, color = "black") +
    geom_intron(data = to_intron(t1_exons, "transcript_id"), aes(xstart = start, xend = end, strand = strand), arrow.min.intron.length = 300) +
    geom_half_range(data = t1_cds, aes(xstart = start, xend = end), fill = "purple", range.orientation = "top", height = 0.2) +
    geom_half_range(data = t1_myb, aes(xstart = start, xend = end), fill = "coral", range.orientation = "top", height = 0.19, color = "black", linewidth = 0.001) +
    labs(subtitle = "Comparison of two most abundant transcripts", x = "Genomic Coordinates", y = "") +
    theme_minimal(base_size = 12)

  combined_plot <- p1 / p2 +
    plot_layout(heights = c(length(unique(gene_exons$transcript_id)), 4)) +
    plot_annotation(
      title = paste("Transcript Structures for", gene_name),
      subtitle = paste("TAIR ID:", gene_id),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12)
      )
    )
  return(combined_plot)
}

# Create a Full Dashboard for a Gene
plot_expression_dashboard <- function(gene_id, gene_name, gtf_data, raw_df, avg_df, dataset) {
  p_expr <- plot_gene_expression(tair_code = gene_id, gene_name = gene_name, raw_df = raw_df, avg_df = avg_df, dataset = dataset)
  p_ind_transcripts <- plot_individual_transcripts(tair_code = gene_id, gene_name = gene_name, avg_df = avg_df, dataset = dataset)
  p_structure <- plot_gene_structure_and_comparison(gene_id = gene_id, gene_name = gene_name, gtf_data = gtf_data, avg_df = avg_df)

  combined_plot <- if (is.null(p_structure)) {
    p_expr + p_ind_transcripts
  } else {
    (p_expr + p_ind_transcripts) / p_structure + plot_layout(heights = c(1, 1.2))
  }

  final_dashboard <- combined_plot + plot_annotation(
    title = paste("Expression and Structure Dashboard for", gene_name, "-", toupper(dataset)),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
  )

  output_filename <- file.path(figures_dir, paste0("Dashboard_", gene_name, "_", toupper(dataset), ".svg"))
  ggsave(plot = final_dashboard, filename = output_filename, width = 12, height = 12)
  cat(paste("Dashboard saved as", output_filename, "\n"))

  return(final_dashboard)
}

# Plot a Grid of Total Gene Expression for a List of Genes
plot_total_expression_grid <- function(raw_df, avg_df, gene_list, dataset) {
  if (!("gene" %in% colnames(raw_df))) {
    raw_df <- tibble::rownames_to_column(raw_df, var = "gene")
  }

  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  all_tair_codes <- gene_list$tair_code

  plot_data <- raw_df %>%
    mutate(tair_code = str_extract(gene, regex$tair_code)) %>%
    filter(tair_code %in% all_tair_codes) %>%
    pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "abundance") %>%
    mutate(
      condition_group = str_extract(sample, regex$condition_group),
      time_h = extract_time(sample, regex, time_map)
    ) %>%
    group_by(tair_code, sample, condition_group, time_h) %>%
    summarise(total_abundance_per_rep = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(tair_code, condition_group, time_h) %>%
    summarise(mean_abundance = mean(total_abundance_per_rep), sd_abundance = sd(total_abundance_per_rep), .groups = "drop") %>%
    left_join(gene_list, by = "tair_code") %>%
    filter(!is.na(time_h) & !is.na(condition_group))

  # Stop if there is no data to plot after filtering
  if (nrow(plot_data) == 0) {
    warning("No data available for plotting in the specified gene list.")
    return(NULL)
  }

  ggplot(plot_data, aes(x = time_h, y = mean_abundance, group = condition_group, color = condition_group)) +
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax = mean_abundance + sd_abundance), width = 1.5, alpha = 0.4) +
    geom_line(aes(linetype = condition_group), linewidth = 1) +
    geom_point(aes(shape = condition_group), fill = "white", size = 2.5) +
    facet_wrap(~gene_name, scales = "free_y", ncol = 5) +
    scale_color_manual(name = "Condition", values = aesthetics$values, labels = aesthetics$labels) +
    scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
    scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
    labs(
      title = paste("Total Gene Expression Across All Genes -", toupper(dataset)),
      x = "Time (Hours)", y = "Mean Expression (TPM)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 10)
    )
}

# Plot a Grid of Top 3 Individual Transcripts for a List of Genes
plot_individual_transcripts_grid <- function(avg_df, gene_list, dataset) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

  if (!("gene" %in% colnames(avg_df))) {
    avg_df <- tibble::rownames_to_column(avg_df, var = "gene")
  }

  plot_list <- lapply(1:nrow(gene_list), function(i) {
    current_tair <- gene_list$tair_code[i]
    current_name <- gene_list$gene_name[i]

    plot_data_full <- avg_df %>%
      filter(str_starts(gene, current_tair)) %>%
      pivot_longer(cols = where(is.numeric), names_to = "condition", values_to = "mean_value") %>%
      mutate(
        condition_group = str_extract(condition, regex$condition_group),
        time_h = extract_time(condition, regex, time_map),
        transcript_id = str_extract(gene, regex$transcript_suffix)
      ) %>%
      filter(!is.na(time_h), !is.na(transcript_id))

    if (nrow(plot_data_full) == 0) {
      return(ggplot() +
        theme_void() +
        labs(title = paste(current_name, "\n(No Data)")))
    }

    top_transcripts <- plot_data_full %>%
      group_by(transcript_id) %>%
      summarise(overall_mean = mean(mean_value, na.rm = TRUE)) %>%
      slice_max(order_by = overall_mean, n = 3) %>%
      pull(transcript_id)

    plot_data <- plot_data_full %>% filter(transcript_id %in% top_transcripts)

    if (nrow(plot_data) == 0) {
      return(ggplot() +
        theme_void() +
        labs(title = paste(current_name, "\n(No top transcripts to plot)")))
    }

    ggplot(plot_data, aes(x = time_h, y = mean_value, group = interaction(transcript_id, condition_group))) +
      geom_line(aes(color = transcript_id, linetype = condition_group), linewidth = 0.8) +
      geom_point(aes(color = transcript_id, shape = condition_group), size = 2, alpha = 0.8) +
      scale_shape_manual(name = "Condition", values = aesthetics$shapes, labels = aesthetics$labels) +
      scale_linetype_manual(name = "Condition", values = aesthetics$linetypes, labels = aesthetics$labels) +
      guides(color = guide_legend(override.aes = list(linewidth = 1, shape = 16))) +
      labs(title = current_name, x = "Time (Hours)", y = "Mean Expression (TPM)", color = "Transcript ID") +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.box = "vertical"
      )
  })

  patchwork::wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = paste("Top 3 Transcripts for Each Gene -", toupper(dataset)),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22))
    )
}
# Plot diversity comparison of relevant gene
plot_diversity_comparison <- function(tair_code, gene_name, time_h, raw_df, dataset) {
  regex <- get_regex_patterns(dataset)
  aesthetics <- get_plot_aesthetics(dataset)
  time_map <- if (regex$pattern_type == "mapped") get_time_mapping() else NULL

diversity_summary_df <- raw_df %>%
  tibble::rownames_to_column(var = "transcript_id") %>%
  dplyr::filter(str_starts(transcript_id, tair_code)) %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    condition_group = str_extract(sample, regex$condition_group),
    time_h_calc = extract_time(sample, regex, time_map) # Assumes extract_time() is a custom function
  ) %>%
  dplyr::filter(time_h_calc == !!time_h) %>%
  # Pool Replicates 
  # Group by condition and transcript to average the replicates
  group_by(condition_group, transcript_id) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_abundance > 0) %>%
  group_by(condition_group) %>%
  # Calculate values needed for proportions (P)
  mutate(
    # N is the total abundance for the entire condition group
    N = sum(mean_abundance),
    # P is the proportion of each transcript within its group
    P = mean_abundance / N
  ) %>%
  # Calculate the components for the Shannon Index (H) and its variance
  mutate(
    P_lnP = P * log(P),
    P_lnP_sq = P * (log(P))^2
  ) %>%
  # Summarise Each Group for Final Test Statistics
  summarise(
    # Richness (S) is the count of unique transcripts in the group
    Richness = n(),
    # N is constant within the group, so we can take the mean (or first value)
    N = mean(N, na.rm = TRUE),
    # Calculate the Shannon Index (H)
    H = -sum(P_lnP),
    # Sum the squared term needed for variance calculation
    sum_P_lnP_sq = sum(P_lnP_sq),
    .groups = "drop"
  ) %>%
  # Calculate the Variance of H 
  mutate(
    Var_H = ((sum_P_lnP_sq - (H^2)) / N) + ((Richness - 1) / (2 * N^2))
  ) %>%
  select(condition_group, Richness, N, H, Var_H)

# --- Hutchinson's t-test Implementation ---

if (nrow(diversity_summary_df) == 2) {
  # Extract values for each group
  group1_H <- diversity_summary_df$H[1]
  group1_Var_H <- diversity_summary_df$Var_H[1]
  group1_N <- diversity_summary_df$N[1]
  
  group2_H <- diversity_summary_df$H[2]
  group2_Var_H <- diversity_summary_df$Var_H[2]
  group2_N <- diversity_summary_df$N[2]
  
  # Calculate the t-statistic
  t_statistic <- (group1_H - group2_H) / sqrt(group1_Var_H + group2_Var_H)
  
  # Calculate the degrees of freedom
  numerator_df <- (group1_Var_H + group2_Var_H)^2
  denominator_df <- ((group1_Var_H^2) / group1_N) + ((group2_Var_H^2) / group2_N)
  degrees_of_freedom <- numerator_df / denominator_df
  
  # Calculate the two-tailed p-value
  p_value <- 2 * pt(-abs(t_statistic), df = degrees_of_freedom)
} else {
  warning("Skipping t-test and plot generation because the number of groups is not equal to 2.")
  p_value <- NA # Set p_value to NA if test is not run
}

# Proceed with plotting only if the t-test was performed
if (!is.na(p_value)) {
  
  
  plot_data <- diversity_summary_df %>%
    mutate(
      CI_width = sqrt(Var_H),
      ymin = H - CI_width,
      ymax = H + CI_width
    )
  
  # Prepare labels and positions for the plot
  p_label <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", round(p_value, 3)))
  y_position <- max(plot_data$ymax) * 1.1

  # Define the significance bracket coordinates
  significance_bracket <- tibble(
    x = c(1, 1, 2, 2),
    y = c(y_position, y_position * 1.02, y_position * 1.02, y_position)
  )

  # General plot theme and labels
  plot_labs <- labs(
    title = paste("Shannon Diversity of", gene_name, "at", time_h, "h"),
    subtitle = paste("Hutchinson's t-test:", p_label),
    x = "Condition",
    y = "Shannon Index, H"
  )
  plot_theme <- theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )

  # Create the bar plot
  p_bar <- ggplot(plot_data, aes(x = condition_group, y = H, fill = condition_group)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.6) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linewidth = 0.7) +
    geom_line(data = significance_bracket, aes(x = x, y = y), inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = y_position * 1.05, label = p_label, size = 4.5) +
    scale_x_discrete(labels = aesthetics$labels) +
    scale_fill_manual(values = aesthetics$values, guide = "none") +
    plot_labs +
    plot_theme

  # Create the point plot
  p_point <- ggplot(plot_data, aes(x = condition_group, y = H, color = condition_group)) +
    geom_point(size = 4, shape = 18) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1, linewidth = 0.7) +
    geom_line(data = significance_bracket, aes(x = x, y = y), inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = y_position * 1.05, label = p_label, size = 4.5) +
    scale_x_discrete(labels = aesthetics$labels) +
    scale_color_manual(values = aesthetics$values, guide = "none") +
    plot_labs +
    plot_theme

  # Combine the plots using patchwork
  combined_plot <- p_bar + p_point
  
  # Create directory if it doesn't exist
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }

  # Save the combined plot
  output_filename <- file.path(figures_dir, paste0("SdComparison_", gene_name, "_", time_h, "_", toupper(dataset), ".svg"))
  ggsave(plot = combined_plot, filename = output_filename, width = 8, height = 8)
  cat(paste("SdComparison saved as", output_filename, "\n"))
  
  
}
return(combined_plot)
}

