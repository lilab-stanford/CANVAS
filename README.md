# CANVAS

<table>
<tr>
<td width="60%" valign="top">

### A cross-modality AI framework for ecological habitat reconstruction from histopathology

**CANVAS** (***C**ellular **A**rchitecture and **N**eighborhood-informed **V**irtual **A**I-driven **S**patial profiling*) is a biologically grounded, cross-modality AI framework that integrates high-dimensional spatial proteomics with routine histopathology to reconstruct, quantify, and model tumor ecological habitats directly from H&E-stained tissue sections. By learning to transfer CODEX-defined cellular neighborhoods (CNs) onto native histology, CANVAS enables accurate and scalable inference of spatial habitats across whole-slide images. This cross-modality alignment is powered by vision–language foundation modeling, refined through AI-agent–guided biological interpretation, and optimized via efficient spatial feature selection and machine learning–based prognostic modeling. Spanning single-cell resolution to population-scale clinical inference, CANVAS builds a functional bridge between spatial biology and precision oncology, enabling interpretable and translatable spatial biomarker discovery.

</td>
<td width="60%" align="center" valign="middle">

<img src="https://github.com/lilab-stanford/CANVAS/blob/main/Abstruct_figure/CANVAS_image.png">

</td>
</tr>
</table>

### 🔧 Key Modules

#### (1) CN-to-habitat prediction via a vision–language foundation model

CANVAS leverages a vision–language foundation model trained on paired CODEX and H&E whole-slide images to infer CN-defined ecological habitats directly from unannotated histology. This module integrates:

- Use the palom package to co-register CODEX and H&E images, and to map CODEX-defined cellular neighborhoods onto histology
- Use cn_assignment to generate habitat labels for training and validation
- Train a vision–language model to predict habitat classes from H&E images based on the co-registered CNs
- Predict habitat classes on H&E images using the trained model

**Run:**

```bash
python cn_assignment.py
python habitat_training.py
python habitat_prediction.py
```

**Output:** Patch-level maps of predicted habitat classes across whole-slide H&E images.

#### (2) Habitat-level spatial feature generation

For each inferred habitat, CANVAS extracts a suite of biologically interpretable spatial features spanning six categories:

- **Composition:** Quantifies the relative abundance of each spatial habitat within the tissue section
- **Diversity:** Measures habitat heterogeneity using metrics such as Shannon index, Fisher’s alpha, and richness
- **Spatial dispersion:** Captures habitat-level spatial organization using Ripley’s K, L, and F functions, Clark–Evans index, and kernel density
- **Interaction:** Quantifies pairwise spatial relationships between habitats using permutation-based analysis
- **Distance:** Computes directional Euclidean distances between habitat pairs
- **Transition entropy:** Non-Euclidean metrics reflecting the internal complexity and interfacial dynamics of the spatial architecture

**Run:**

```bash
Rscript Feature_generation/Spatial_metrics.R
Python Feature_generation/Distance_calculation.py
```

**Output:** Multiscale spatial feature matrices per sample.

#### (3) Feature selection and prognostic modeling

To identify clinically meaningful spatial biomarkers, CANVAS performs:

- Bootstrap LASSO and RandomForestSRC for feature stability and importance
- Univariate Cox regression for interpretable survival association
- Multivariate modeling across cohorts for immunotherapy outcome prediction

**Run:**

```bash
Rscript Feature_selection_modeling/feature_selection_modeling.R
```

**Output:** Refined biomarker panels and risk prediction models for immunotherapy outcomes.

#### (4) AI-Agent module for spatial feature interpretation

CANVAS incorporates a large language model (LLM)-driven AI agent to facilitate the systematic interpretation of high-dimensional spatial features by linking them to underlying spatial biology and relevant clinical priors.

The agent generates structured, human-interpretable outputs across five key biological dimensions:

- (i) Feature category
- (ii) Characteristic cellular composition of the associated habitat
- (iii) Description of the encoded spatial property
- (iv) Tendency for topological coupling with other habitats
- (v) Potential biological and clinical relevance

**Run:**

```bash
python AI_Agent/run.py
```
**Note:** Replace the default API key with your own in spatial_agent.py.
**Output:** An interpretable annotation layer that contextualizes model-derived spatial features within tissue architecture, ecological topology, and clinical significance.

---

### 🔬 Applications

- Immunotherapy response stratification in NSCLC and other solid tumors using spatially resolved habitat features
- Habitat-informed prognostic modeling across large-scale cohorts including TCGA, PLCO, and NLST
- Patient subtyping independent of PD-L1 expression levels or canonical oncogenic mutations (e.g., EGFR, KRAS)
- Virtual spatial proteomic annotation directly inferred from standard H&E histology, enabling scalable tissue profiling without multiplex staining

---

### 📦 Installation

We recommend using `conda` to manage environments for CANVAS.

```bash
conda create -n canvas_env python=3.8
conda activate canvas_env
pip install -r requirements.txt
```

R packages used include:

```r
install.packages(c("survival", "glmnet", "randomForestSRC", "ggplot2", "vegan", "entropy", "spatstat"))
```

---

### 📂 Directory Structure

```
CANVAS/
├── README.md                          # Project documentation and usage instructions

├── Demo_data/                         # Example dataset
│   └── Spatial_feature_matrix.csv     # CANVAS-derived spatial feature matrix for each sample

├── Habitat_prediction/                # Module 1: CN-to-habitat prediction via a vision–language foundation model
│   ├── cn_assignment.py               # Co-registration of CODEX and H&E images at single-cell resolution
│   ├── habitat_prediction.py          # Predicts ecological habitats from CN annotations
│   └── habitat_training.py            # Trains vision–language model for habitat prediction

├── Feature_generation/                # Module 2: Habitat-level spatial feature generation
│   └── Spatial_metrics.R              # Calculates composition, diversity, interaction, and other spatial metrics

├── Feature_selection_modeling/        # Module 3: Feature selection and prognostic modeling
│   └── feature_selection_modeling.R   # Performs Bootstrap LASSO, random forest, and Cox regression modeling

├── AI_agent/                          # Module 4: AI-Agent module for spatial feature interpretation
│   └── AI_agent.py                    # Interactive agent for habitat-level biological and clinical annotation

├── Abstruct_figure/                   # Folder for figures used in the abstract or main text


```

### 📄 Citation

Li Y. et al. *Cellular Architecture and Neighborhood-informed Virtual Spatial Tumor Profiling from Histopathology* (2025)

---

### 📬 Contact

Project maintained by the Department of Radiation Oncology, Stanford University.\
📧 [ycli16@stanford.edu](mailto\:ycli16@stanford.edu)

---

### 🧠 Acknowledgements

We thank the developers of core tools and libraries including CODEX, Seurat, and PyTorch.\
This work was supported by funding from the National Institutes of Health (NIH).
