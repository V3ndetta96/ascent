# Manuscript Outline: ASC-FCD Framework

**Target Journal Options (Free / Hybrid):** 
1. *Ecological Informatics* (Elsevier - Free subscription option)
2. *Ecological Modelling* (Elsevier - Free subscription option)
3. *Oikos* (Wiley - Hybrid, free to publish)
**Article Type:** Software / Methodological Paper
**Title Idea 1:** ascent: A multi-layer framework for decomposing functional community restructuring under hierarchical null models.
**Title Idea 2:** Deconstructing functional turnover: A multi-layer geometric approach to ecosystem restructuring.

## 1. Abstract (Max 250 words)
- **Background:** Current functional beta-diversity metrics collapse complex multi-dimensional changes into single scalars, masking divergent ecological trajectories.
- **Innovation:** We present the ASC-FCD framework (implemented in the R package `ascent`), which decomposes functional restructuring into three geometrically orthogonal layers: Position (centroid displacement), Dispersion (internal variance), and Boundary (convex hull volume).
- **Mechanism:** We introduce a novel "Functional Leverage" metric that attributes total ecosystem displacement to individual species' demographic changes.
- **Validation:** We propose a triad of hierarchical null models (Structural, Quantitative, Identity) to isolate specific assembly rules.
- **Conclusion:** `ascent` provides a unified, computationally efficient pipeline to track ecosystem trajectories under environmental change.

## 2. Introduction
- **The Problem:** Global change drives community restructuring. We use functional traits to understand this. However, traditional $\beta$-diversity indices (e.g., Jaccard, Bray-Curtis, or even functional turnover) cannot tell us *in which direction* an ecosystem is moving or *how* its internal structure is changing.
- **The Gap:** Existing frameworks (like `FD`, `betapart`, `TPD`) measure magnitude but not direction. They do not easily decompose the change into independent geometric components, nor do they map community-level shifts back to individual species drivers.
- **The Solution:** Introduce the ASC-FCD conceptual model. 
  - Layer 1: Directional Shift (Where is it going?)
  - Layer 2: Dispersive Reorganization (How is biomass rearranging?)
  - Layer 3: Boundary Dynamics (Are extremes being lost/gained?)
- **Objective:** Introduce the `ascent` R package, detail its mathematical foundations, and demonstrate its utility using a case study.

## 3. Mathematical Framework (Methodology)
*(This section will heavily reuse `ASCFCD_technical_document.md`)*
- **3.1 Functional Space Construction:** Gower distance + PCoA (handling mixed traits, quality metrics).
- **3.2 The Three Layers of Restructuring:** 
  - Formal definition of $\Delta C$ (Position).
  - Formal definition of $\Delta FDis$ (Dispersion).
  - Formal definition of $\Delta FRic$ (Boundary/Volume) and its caveats (incidence-based topological descriptor).
- **3.3 Species-Level Functional Leverage:** 
  - Mathematical definition (projection of species onto the directional vector).
  - Proof of additivity ($\sum Lev = \Delta C$).

## 4. Hierarchical Null Models for Inference
- Why a single null model is insufficient for functional ecology.
- **Model A (Structural / Curveball):** Isolating incidence and regional prevalence.
- **Model B (Quantitative / SAD):** Isolating demographic variance (explaining why FRic is NA here).
- **Model C (Identity / Trait Shuffle):** Isolating species-trait associations.

## 5. Software Implementation (The `ascent` Package)
- Brief overview of package architecture (S3 objects, speed, `ggplot2` integration).
- Core functions: `asc_paired`, `asc_pairwise`, `asc_transitions`, `assess_functional_space`.

## 6. Case Study / Application
*(We need a real or simulated ecological dataset here to show it in action)*
- E.g., The impact of deforestation on a bird community, or climate warming on alpine plants.
- Demonstrate how classical metrics would show "X% change", but `ascent` reveals that the centroid moved specifically towards "smaller, non-forest-dependent species", driven primarily by the collapse of Species A (highest leverage).

## 7. Discussion
- **Ecological Implications:** Why separating position, dispersion, and boundary matters for conservation (e.g., you can have zero boundary change but massive positional shift).
- **Limitations:** Assumptions of PCoA distortion, sample size requirements for null models, constraints of the convex hull.
- **Future Directions:** Integration with phylogenetic trees or continuous trait density functions.

## 8. Data Availability
- CRAN link, GitHub repository, and reproducible scripts for the case study.
