# ascent: Functional Centroid Displacement Analysis in R 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20388186.svg)](https://doi.org/10.5281/zenodo.20388186)
The **`ascent`** (ASC-FCD) package provides a robust geometric and statistical framework to evaluate shifts in the functional structure of biological communities across environmental gradients, disturbances, or time series. 

Unlike traditional functional diversity metrics that merely quantify the magnitude of change, `ascent` identifies the statistical direction of functional displacement, assesses ecosystem resilience through null models, and taxonomically isolates the "winning" and "losing" species driving the turnover.

## 🚀 Key Features

* **Universal Trait Support:** Seamlessly handles binary, categorical, and continuous traits using Gower distances and Community Weighted Means (CWM) under the hood.
* **Paired and Network Designs:** Contrast specific temporal/spatial pairs (e.g., Before/After impact) or compute complex pairwise networks across multiple communities.
* **Null Model Integration:** Distinguish between neutral ecological drift and environmental filtering by assessing the structural resilience of the community.
* **Taxonomic Transitions:** Decompose critical axes of change to pinpoint the exact Functional Entities (FE) and species driving the ecological shift.
* **Functional Redundancy:** Quantify redundancy and identify unique Functional Entities to assess system vulnerability to local extinctions.

## 📦 Installation

You can install the development version of `ascent` from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("V3ndetta96/ascent")

💡 Quick Start (Usage)
Here is a basic workflow to evaluate the functional impact of an environmental disturbance on an estuarine community:

library(ascent)

# 1. Load your data (Traits matrix and Abundance matrix)
traits <- read.csv("traits_matrix.csv", row.names = 1)
abund <- read.csv("abundance_matrix.csv", row.names = 1)

# 2. Define experimental design vectors
time_vec <- ifelse(grepl("_before$", rownames(abund)), "Before", "After")
site_vec <- gsub("_before$|_after$", "", rownames(abund))

# 3. Compute paired functional centroid shifts
res_paired <- asc_paired(traits = traits, 
                         abund = abund, 
                         sites = site_vec, 
                         time = time_vec, 
                         ref_time = "Before")

# 4. Run null models (resilience) and extract species transitions
res_paired <- asc_null(res_paired, n_perm = 999)
transitions <- asc_transitions(res_paired, top_n = 3)

# 5. Review results and plot
summary(res_paired)
summary(transitions)
plot(res_paired, contrast = "Site_A")

📖 Documentation
For a detailed guide on the mathematical background and ecological interpretation of the ASC-FCD metric, please read the package vignette:

vignette("intro-to-ascent", package = "ascent")

## ✒️ Citation
While the official manuscript is in preparation, if you use `ascent` in your research, please cite the package directly:

Muñoz-Li, R. R. & F. Alvarez-Denis (2026). ascent: Functional Centroid Displacement Analysis in R. R package version 0.1.0. https://github.com/V3ndetta96/ascent
