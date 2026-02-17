
# Project:  DXY in Sliding Window Analysis

*Author: Maria Niki Chatzantoni*

*Class: BINP28*

*Date: February 2026*


This repository contains the scripts used to analyze absolute genetic divergence (dxy) between three populations: **8N**, **K**, and **Lesina**.

The pipeline takes raw VCF data, cleans it up, calculates divergence in sliding windows (using `pixy`), and generates the final plots.

## The Scripts

There are three main scripts that work together. Here is what each one does:

### 1. `project_dxy.sh` (The Main Pipeline)

This is the main script. It handles the entire workflow from raw data to final statistics:

* **Preprocessing:** Filters the input VCF to keep only biallelic SNPs (removes indels and multi allelic sites) and removes the `Naxos2` sample, which is not needed for the downstream analysis.
* **Stats Generation:** Calculates per site and per individual stats (Depth, Quality, Missingness, MAF) using `vcftools` so we can decide on filter cutoffs.
* **Filtering:** Creates two clean datasets based on Genotype Quality (GQ) thresholds:
* **Strict Set:** GQ $\ge$ 15
* **Relaxed Set:** GQ $\ge$ 10
* *Standard Filters:* Min Depth: 5, Max Depth: 27, Max Missing: 10%, MAF: 0.033.


*  Calculation: Runs `pixy` on the cleaned data using 10kb, 50kb, and 100kb sliding windows.

### 2. `filtering.R` (Quality Control)

We use this script *before* finalizing the parameters in the bash script.

* It visualizes the raw distributions of Quality, Depth, and Missing Data.
* **Goal:** It helps us spot outliers (e.g., finding the max depth cutoff to avoid paralogs) and confirm that our filter thresholds make sense.

### 3. `pixy.R` (Visualization)

Once `pixy` is done, this script makes the plots.

* It reads the output tables (`_dxy.txt`) and plots divergence across the genome.
* It adds a smoothing line to show trends.
* **Output:** Plots in a png format, comparing **Chr5 vs ChrZ** for all population pairs.
---
## Workflow Overview

### 1. Data Preparation & Filtering (`project_dxy.sh`)
This script is the main script of the analysis. It performs the following steps:
1.  **VCF Cleaning:** Keeps only SNPs (removes indels and multi allelic sites) and removes unwanted samples (e.g., `Naxos2`).
2.  **Statistics Generation:** Uses `vcftools` to calculate per site and per individual statistics (Depth, Quality, Missingness, MAF).
3.  **Hard Filtering:** Creates two filtered datasets based on different Genotype Quality (GQ) thresholds:
    * **Set 1:** GQ $\ge$ 15
    * **Set 2:** GQ $\ge$ 10
    * *Common Filters:* Max missing 10%, Min Depth 5, Max Depth 27, MAF 0.033.
4.  **$D_{xy}$ Calculation:** Runs `pixy` on both datasets using multiple window sizes (10kb, 50kb, and 100kb).
---

## How to Run It

### Step 1: Prepare the Data

Make sure the input VCF (`ProjTaxa.vcf`) is in the working directory. We'll also need the following tools installed: `bcftools`, `vcftools`, `htslib` (for tabix/bgzip), and `pixy`.

### Step 2: Run the Pipeline

Execute the main bash script. This will generate the filtered VCFs and run the `pixy` calculations.

```bash
chmod +x project_dxy.sh
./project_dxy.sh

```

### Step 3: Check Diagnostics (Optional)

If we want to see why we chose specific filter cutoffs, run the QC visualization script:

```bash
Rscript filtering.R

```

### Step 4: Generate Plots

Finally, visualize the  results:

```bash
Rscript pixy.R

```

---

## Populations DXY Analyzed

The analysis groups samples into three populations (defined in `populations.txt`):

* **8N** (5 samples)
* **K** (5 samples)
* **Lesina** (5 samples)

## Outputs

The pipeline produces:

* **Cleaned VCFs:** `filtered_gq15.vcf.gz` and `filtered_gq10.vcf.gz`
* **Data Tables:** `results_*_dxy.txt` (The raw numbers from pixy)
* **Figures:** `dxyplot_*.png` (The final plots for 10kb, 50kb, and 100kb windows)