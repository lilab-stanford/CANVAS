# CANVAS

<table>
<tr>
<td width="60%" valign="top">

### A cross-modality AI framework for ecological habitat reconstruction from histopathology

**CANVAS** (***C**ellular **A**rchitecture and **N**eighborhood-informed **V**irtual **A**I-driven **S**patial profiling*) is a biologically grounded, cross-modality AI framework that integrates high-dimensional spatial proteomics with routine histopathology to reconstruct, quantify, and model tumor ecological habitats directly from H&E-stained tissue sections. By learning to transfer CODEX-defined cellular neighborhoods (CNs) onto native histology, CANVAS enables accurate and scalable inference of spatial habitats across whole-slide images. This cross-modality alignment is powered by vision–language foundation modeling, refined through AI-agent–guided biological interpretation, and optimized via efficient spatial feature selection and machine learning–based prognostic modeling. Spanning single-cell resolution to population-scale clinical inference, CANVAS builds a functional bridge between spatial biology and precision oncology, enabling interpretable and translatable spatial biomarker discovery.

</td>
<td width="60%" align="center" valign="middle">

<img src="https://github.com/lilab-stanford/CANVAS/blob/main/Abstract_figure/CANVAS_image.png">

</td>
</tr>
</table>


### 🔧 Key Modules

#### (1) CN-to-habitat prediction via a vision–language foundation model

CANVAS leverages a vision–language foundation model trained on paired CODEX and H&E whole-slide images to infer CN-defined ecological habitats directly from unannotated histology. This module integrates:

- Co-register CODEX and H&E images at single-cell resolution and map CODEX-defined cellular neighborhoods onto histology
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
- **Spatial dispersion:** Describes habitat spatial distribution based on Ripley’s K, L, and F functions, as well as related metrics
- **Interaction:** Quantifies pairwise spatial associations between habitats using permutation-based analysis
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

**Output:** Refined panels of spatial biomarkers and risk prediction models with translational relevance for immunotherapy.

#### (4) AI-Agent module for spatial feature interpretation

CANVAS incorporates a large language model (LLM)-driven AI agent to facilitate the systematic interpretation of high-dimensional spatial features by linking them to underlying spatial biology and relevant clinical priors.

The agent generates structured, human-interpretable outputs across five key biological dimensions:

- (i) Feature category
- (ii) Characteristic cellular composition of the associated habitat
- (iii) Description of the encoded spatial property
- (iv) Tendency for topological coupling with other habitats
- (v) Potential biological and clinical relevance

**Setup:**

```bash
export OPENAI_API_KEY="<your-openai-api-key>"
# Optional: only if you route through a custom endpoint
export OPENAI_BASE_URL="https://your-custom-endpoint"
```
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
- Virtual reconstruction of the tumor–immune ecosystem directly from routine H&E histology, enabling scalable tissue profiling without multiplex staining

---

### 📦 Installation and Basic Usage

We recommend using `conda` to manage environments for CANVAS.

```bash
#Please install MUSK first (required)
conda create -n canvas_env python=3.8
conda activate canvas_env
pip install -r requirements.txt
```

R packages used include:

```r
install.packages(c("survival", "glmnet", "randomForestSRC", "ggplot2", "vegan", "entropy", "spatstat"))
```

### Basic Usage
1. Load the model
```bash
from Habitat_prediction.api import load_model, predict_folder

model, device = load_model(
    #weights=str(WEIGHTS_PATH),    # optional: local weights only
    musk_source=MUSK_SOURCE,      # hf_hub or local path to MUSK backbone
)
print("Loaded model on:", device)
```
2. Run inference on all images under the folder
```bash
import pandas as pd
from pathlib import Path

df = predict_folder(
    model, device, DEMO_DIR,
    img_size=384, batch_size=64,
)

OUT_DIR = Path("Demo_data"); OUT_DIR.mkdir(exist_ok=True)
out_csv = OUT_DIR / "output.csv"
df.to_csv(out_csv, index=False)
print("Saved")
```


---

### 📂 Directory Structure

```
CANVAS/
│
│   demo.ipynb                 # Example Jupyter notebook demonstrating CANVAS workflow
│   README.md                  # Project documentation and usage instructions
│   requirements.txt           # Python dependencies for running CANVAS
│
├── Abstract_figure/           # Figures used for abstract or manuscript illustration
│   └── CANVAS_image.png       # Thumbnail of CANVAS
│
├── AI_Agent/                  # Module 4: AI-Agent for spatial feature interpretation
│   │   Instruction for AI agent.docx   # Documentation of AI-agent usage
│   │   README.md                        # Module-specific documentation
│   │   requirements.txt                 # Dependencies for AI-agent module
│   │   run.py                           # Main script to launch AI-agent
│   │   spatial_agent.py                 # Core implementation of AI-agent
│   │
│   └── data/                            # Supporting annotation files for AI-agent
│       ├── Feature_annotation.xlsx      # Definitions and biological annotation of spatial features
│       └── Habitat_annotation.docx      # Annotation of habitats and ecological interpretation
│
├── Demo_data/                 # Example datasets for demonstration
│   ├── 28000_56224.png        # Example tissue image (demo figure)
│   ├── 67424_15680.png        # Example tissue image (demo figure)
│   ├── output.csv             # Example model output (demo results)
│
├── Feature_generation/        # Module 2: Habitat-level spatial feature generation
│   ├── Distance_calculation.py  # Computes pairwise distances among habitats
│   ├── Habitat_freq.R           # Quantifies habitat frequency per sample
│   ├── Habitat_interaction.py   # Computes inter-habitat interactions
│   ├── Spatial_diversity.R      # Calculates diversity indices (Shannon, Simpson, etc.)
│   ├── Spatial_entropy.R        # Computes spatial transition entropy (STE)
│   └── Spatial_metrics.R        # Master script for habitat-level spatial dispersion
│
├── Feature_selection_modeling/   # Module 3: Feature selection and prognostic modeling
│   └── feature_selection_modeling.R  # Performs Bootstrap LASSO, random forest, Cox regression
│
└── Habitat_prediction/         # Module 1: CN-to-habitat prediction using vision–language model
    ├── api.py                  # API wrapper for model inference
    ├── cn_assignment.py        # Co-registers CODEX and H&E images at single-cell resolution
    ├── habitat_prediction.py   # Predicts ecological habitats from CN annotations
    ├── habitat_training.py     # Trains the vision–language model for habitat prediction
    ├── model.py                # Core model architecture definition
    └── reference_weight.pth    # Reference model weights for habitat prediction (available at Zenodo: [10.5281/zenodo.17220060](https://zenodo.org/records/17220060))
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
