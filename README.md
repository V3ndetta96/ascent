# ascent: A multi-layer framework for decomposing functional community restructuring into positional, dispersive and boundary components under hierarchical null models.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20388186.svg)](https://doi.org/10.5281/zenodo.20388186)

The **`ascent`** package implements the **ASC-CFD** (Assemblage Shift Characterization — Community Functional Dynamics) framework.

Most functional diversity metrics quantify how much communities differ,
but provide little information about how those differences are organized
within functional space.

Communities may exhibit identical levels of functional turnover while
following completely different ecological trajectories, including
directional shifts, internal reorganization, or expansion and contraction
of functional boundaries.

The `ascent` package implements the ASC-CFD framework to explicitly
decompose these processes. 

`ascent` provides a multi-layer topological approach to **decompose functional community restructuring into positional, dispersive, and boundary components under hierarchical null models.**

Whether you are evaluating temporal perturbations (e.g., deforestation, climate change) or spatial beta-diversity gradients, `ascent` allows you to reconstruct the geometric trajectory of ecosystems and identify the specific taxa driving these shifts.

## 🧠 Core Architecture

The framework evaluates ecological dynamics by deconstructing them into three analytically distinct geometric layers within the retained functional space:

1. **Positional Component ($\Delta C$):** The net directional displacement of the community centroid. Measures the displacement of the community functional centroid within the retained trait space.
2. **Dispersive Component ($\Delta FDis$):** The internal demographic reorganization. Quantifies whether biomass is concentrating into central redundant strategies or expanding towards peripheral phenotypes.
3. **Boundary Component ($\Delta FRic$):** The multidimensional Convex Hull volume. Captures changes in the outer functional boundaries of communities through Convex Hull geometry. Because this metric depends exclusively on extreme species, it should be interpreted as a topological descriptor rather than a robust abundance-weighted estimate of functional change.

What does each layer measure?

ΔC     → Where is the community moving?
ΔFDis  → How is biomass being redistributed?
ΔFRic  → Are functional boundaries expanding or collapsing?

### The Triad of Null Models
Rather than relying on a single null expectation, ascent implements a hierarchical inferential framework that sequentially isolates:
* **Structural Turnover (Incidence/Curveball):** Does the taxonomic turnover alter the topology more than expected by chance?
* **Demographic reorganization (Quantitative Filter - Demographic/SAD):** Is the internal biomass redirection deterministic?
* **Functional identity effects (Identity Filter - Trait Shuffle):** Given the observed richness and abundances, are the selected functional strategies significantly different from a random draw of the regional pool?

### Species-level Functional Leverage

The asc_transitions() module projects species-specific abundance changes onto the observed ecosystem trajectory, identifying the taxa responsible for driving functional restructuring.
This allows users to move beyond community-level metrics and directly quantify the contribution of individual species to directional ecosystem change.

### Why not just use Functional Beta Diversity?

Current frameworks quantify the magnitude of functional turnover between communities.
ascent instead decomposes restructuring into three distinct mechanisms:

| Question                                          | Classical beta diversity | ascent |
| ------------------------------------------------- | ------------------------ | ------ |
| How much changed?                                 | ✓                        | ✓      |
| In which direction did it change?                 | ✗                        | ✓      |
| Did the centroid move?                            | ✗                        | ✓      |
| Did internal functional dispersion change?        | ✗                        | ✓      |
| Did the functional boundaries expand or collapse? | ✗                        | ✓      |
| Which species drove the shift?                    | ✗                        | ✓      |



## 🚀 Installation

You can install the development version of `ascent` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("V3ndetta96/ascent")
```

## 💻 Quick Start

Here is a basic workflow evaluating the functional impact of deforestation on a bird community using mixed trait data (quantitative and binary).

```r
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
```

### Visualizing Dynamic Topology and Functional Leverage

`ascent` includes built-in `ggplot2`/`patchwork` methods to visualize the PCoA trajectory, the 2D Convex Hull boundaries, and the Functional Leverage (the specific taxa pulling or releasing the centroid).

```r
plot(res, contrast = "Site1", type = "both", n_sp = 5)
```

## 🛠️ Package API Overview

| Function | Purpose |
|---|---|
| `asc_baseline()` | Absolute baseline functional topology (CWM, FDis, FRic) for isolated communities |
| `asc_paired()` | Temporal or paired experimental contrasts |
| `asc_pairwise()` | Bidirectional spatial divergence across all community pairs |
| `asc_null()` | Triad of null models (Structural, Quantitative, Identity) |
| `asc_transitions()` | Multidimensional leverage score for every species |
| `asc_entities()` | Species clustering into discrete functional entities |
| `assess_functional_space()` | Diagnostic tool to evaluate PCoA quality and dimensionality |

## 📖 Documentation & Vignettes

For a comprehensive guide, including spatial network analyses and detailed ecological interpretations, please read the official vignette:

```r
vignette("ascent_tutorial")
```

## ⚠️ Methodological Notes

### FRic (Functional Richness)

FRic is computed as the volume of the convex hull in the retained PCoA space. Key considerations:

- **Depends exclusively on extreme species.** FRic captures the outer boundary of the functional space. It is insensitive to internal reorganization of biomass. A species that increases from 1% to 90% relative abundance will not alter FRic if it is not on the hull boundary.
- **Not abundance-weighted.** Unlike FDis or the centroid, FRic is a purely incidence-based metric within the framework. It should be interpreted as a **topological descriptor** of the functional boundary, not as a robust estimate of functional change.
- **Requires N > k.** When a community has fewer species than retained PCoA axes, the convex hull is geometrically undefined and FRic is reported as `NA`.

### ΔFDis (Functional Dispersion Change)

- **Not normalized.** ΔFDis is an absolute difference in the retained PCoA space. Its magnitude depends on the trait scales, the number of retained axes, and the Gower distance matrix. Direct comparisons across studies with different species pools or trait sets should be made with caution.
- **Normalized position:** Only the positional component (rΔC%) is normalized by the regional functional diameter (Dmax_regional), enabling cross-study comparisons.

### Null Model Scope

`asc_null()` operates on the full stacked community matrix, assuming a **shared regional species pool**. When analyzing biogeographically independent sites, contrasts should be evaluated separately by calling `asc_null()` on each `ascfcd` object independently.

### Functional Leverage: Additive Property

The Functional Leverage satisfies a strict mathematical decomposition:

$$\sum_{i=1}^{S} \text{Leverage}_i = \|\Delta C\|$$

where $\|\Delta C\|$ is the absolute Euclidean distance between centroids. Each species' leverage is the product of its demographic change ($\Delta p_i$) and its projection onto the unit direction of the centroid shift:

$$\text{Leverage}_i = \Delta p_i \cdot \text{proj}(\mathbf{f}_i, \hat{v})$$

This property ensures that the sum of all individual species contributions exactly reconstructs the total centroid displacement. Species with positive leverage **drive** the shift; species with negative leverage **resist** it.

## ✒️ Citation

While the official manuscript is in preparation, if you use `ascent` in your research, please cite the package directly:

Muñoz-Li, R. R. & F. Alvarez-Denis (2026). ascent: A multi-layer framework for decomposing functional community restructuring into positional, dispersive and boundary components under hierarchical null models. https://github.com/V3ndetta96/ascent
