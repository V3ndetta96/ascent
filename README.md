# ascent: A multi-layer framework for decomposing functional community restructuring into positional, dispersive and boundary components under hierarchical null models.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20480767.svg)](https://doi.org/10.5281/zenodo.20480767)
The **`ascent`** package implements the **ASC-FCD** (Asymmetric Shift in Centroid - Functional Centroid Displacement) framework.

Traditional functional ecology metrics often compress community dynamics into isolated indices, obscuring the underlying mechanisms of change. `ascent` provides a multi-layer topological approach to **decompose functional community restructuring into positional, dispersive, and boundary components under hierarchical null models.**

Whether you are evaluating temporal perturbations (e.g., deforestation, climate change) or spatial beta-diversity gradients, `ascent` allows you to track the exact geometric trajectory of ecosystems and identify the specific taxa driving these shifts.

## 🧠 Core Architecture

The framework evaluates ecological dynamics by deconstructing them into three orthogonal geometric layers within the functional hyperspace:

1. **Positional Component ($\Delta C$):** The net directional displacement of the community centroid. Measures the fundamental shift in the ecosystem's functional equilibrium.
2. **Dispersive Component ($\Delta FDis$):** The internal demographic reorganization. Quantifies whether biomass is concentrating into central redundant strategies or expanding towards peripheral phenotypes.
3. **Boundary Component ($\Delta FRic$):** The multidimensional Convex Hull volume. Captures the strict colonization of novel functional space or the extinction of extreme phenotypes.

### The Triad of Null Models
To isolate deterministic environmental filtering from stochastic noise, shifts are evaluated against a hierarchy of null models:
* **Structural Filter (Incidence/Curveball):** Does the taxonomic turnover alter the topology more than expected by chance?
* **Quantitative Filter (Demographic/SAD):** Is the internal biomass redirection deterministic?
* **Identity Filter (Trait Shuffle):** Given the observed richness and abundances, are the selected functional strategies significantly different from a random draw of the regional pool?

## 🚀 Installation

You can install the development version of `ascent` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("tu-usuario/ascent")

💻 Quick Start
Here is a basic workflow evaluating the functional impact of deforestation on a bird community using mixed trait data (quantitative and binary).
library(ascent)

library(ascent)

# 1. Define Traits and Abundances
# Note: Binary traits are set as factors to use the Gower metric.
traits <- data.frame(
  Body_Mass = c(15, 20, 30, 45, 60, 90, 150, 250, 400, 600),
  Frugivore = factor(c(0, 0, 1, 1, 1, 0, 1, 1, 1, 0))
)
rownames(traits) <- paste0("Sp", 1:10)

abund_time <- rbind(
  Reference = c( 0,  0, 10, 15, 25, 20, 10, 8, 4, 2), # Old-growth forest
  Impacted  = c(40, 35, 15,  5,  0,  0,  0, 0, 0, 0)  # Deforested
)

# 2. Execute the Paired Multi-Layer Framework
res <- asc_paired(
  traits = traits, abund = abund_time, 
  sites = c("Site1", "Site1"), time = c("Reference", "Impacted"),
  ref_time = "Reference", dist_method = "gower"
)

# 3. Evaluate Hierarchical Null Models
res <- asc_null(res, n_perm = 999, seed = 123)

# 4. View Automated Ecological Diagnostic
summary(res)

Visualizing Dynamic Topology and Functional Leverage
ascent includes built-in ggplot2/patchwork methods to visualize the PCoA trajectory, the 2D Convex Hull boundaries, and the Functional Leverage (the specific taxa pulling or releasing the centroid).

plot(res, contrast = "Site1", type = "both", n_sp = 5)

🛠️ Package API Overview
asc_entities(): Calculates the absolute baseline functional topology (CWM, FDis, FRic) for isolated communities.

asc_paired(): Evaluates temporal or paired experimental contrasts.

asc_pairwise(): Computes bidirectional spatial divergence across all combinations in a regional network (Beta-diversity).

asc_null(): Triggers the triad of null models (Structural, Quantitative, Identity).

asc_transitions(): Extracts the exact multidimensional leverage score for every species.

📖 Documentation & Vignettes
For a comprehensive guide, including spatial network analyses and detailed ecological interpretations, please read the official vignette:

vignette("ascent_tutorial")

## ✒️ Citation
While the official manuscript is in preparation, if you use `ascent` in your research, please cite the package directly:

Muñoz-Li, R. R. & F. Alvarez-Denis (2026). ascent: A multi-layer framework for decomposing functional community restructuring into positional, dispersive and boundary components under hierarchical null models. https://github.com/V3ndetta96/ascent
