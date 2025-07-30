======================================================
Project: name still TODO
======================================================

--- TODO ---

- Implement timeseries overlap line plot.
- Fix SDI comparison for timeseries.
- ?Error bars for top 10% transcripts plot?

--- Description ---

This R pipeline provides a comprehensive framework for processing, analyzing, and visualizing time-course gene expression data (in Transcripts Per Million, TPM) alongside gene structural information. It is designed to be modular and can handle multiple different datasets through a simple configuration switch at the top of the script.


--- Requirements ---

1. Software
- R (version 4.0 or newer is recommended).
- An R environment like RStudio or VSCode.

2. R Packages
The script will automatically attempt to install any missing packages. The required packages are:
- dplyr
- tidyr
- stringr
- ggplot2
- tibble
- ggtranscript
- rtracklayer
- patchwork
- svglite

3. Input Files
The script expects the following files to be present in the same directory as the script itself:

- Transcript TPM Data (CSV format): One CSV file for each dataset. The script is configured for:
    - Transcript_TPM_adam_TRIMMED.csv
    - Transcript TPMcoolLL.csv
    - Transcript_TPM_timeSeries_trimmed.csv
- Gene Annotation File (GTF format): A single annotation file.
    - ALL_combined_annotations.gtf


--- Setup & Usage ---

1. File Structure
Ensure your project directory is set up as follows:

/project_folder/
|-- functions.R
|-- main.R
|-- Transcript_TPM_adam_TRIMMED.csv
|-- Transcript TPMcoolLL.csv
|-- Transcript_TPM_timeSeries_trimmed.csv
|-- ALL_combined_annotations.gtf

2. Configure the Script
Open the main R script. The only line you need to change for standard operation is at the very top:

# --- 1. Set Dataset Context ---
dataset <- "coolLL2" # Options: "adam", "coolLL2", "timeseries"

Change "coolLL2" to "adam" or "timeseries" to select the dataset you wish to analyse.

3. Run the Script
Execute the script from within your R environment 


4. Check the Output
The script will create a "figures" directory (if it doesn't exist) and a sub-directory named after the chosen dataset. All output plots will be saved there as .svg files.

For example, if dataset <- "coolLL2", all plots will be saved in:
/your_project_folder/figures/coolLL2/


--- Customization ---

Analysing Different Genes
Change the gene_name and tair_code in the objects that call the plotting functions
To change the list of genes processed in the batch loop or in the gene family grids, change the tibble objects (genes_to_plot, rve_genes, prr_genes).
