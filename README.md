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
├── CANVAS/
│   │
│   │   demo.ipynb                 # Example Jupyter notebook demonstrating CANVAS workflow
│   │   README.md                  # Project documentation and usage instructions
│   │   requirements.txt           # Python dependencies for running CANVAS
│   │
│   ├── Abstract_figure/           # Figures used for abstract or manuscript illustration
│   │   └── CANVAS_image.png       # Thumbnail of CANVAS
│   │
│   ├── Demo_data/                 # Example datasets for demonstration
│   │   ├── 28000_56224.png        # Example tissue image (demo figure)
│   │   ├── 67424_15680.png        # Example tissue image (demo figure)
│   │   ├── output.csv             # Example model output (demo results)
│   │
│   ├── Habitat_prediction/        # Module 1: CN-to-habitat prediction using vision–language model
│   │   ├── api.py                 # API wrapper for model inference
│   │   ├── cn_assignment.py       # Co-registers CODEX and H&E images at single-cell resolution
│   │   ├── habitat_prediction.py  # Predicts ecological habitats from CN annotations
│   │   ├── habitat_training.py    # Trains the vision–language model for habitat prediction
│   │   └── model.py               # Core model architecture definition
│   │
│   ├── Feature_generation/        # Module 2: Habitat-level spatial feature generation
│   │   ├── Distance_calculation.py  # Computes pairwise distances among habitats
│   │   ├── Habitat_freq.R           # Quantifies habitat frequency per sample
│   │   ├── Habitat_interaction.py   # Computes inter-habitat interactions
│   │   ├── Spatial_diversity.R      # Calculates diversity indices (Shannon, Simpson, etc.)
│   │   ├── Spatial_entropy.R        # Computes spatial transition entropy (STE)
│   │   └── Spatial_metrics.R        # Master script for habitat-level spatial dispersion
│   │
│   └── Feature_selection_modeling/   # Module 3: Feature selection and prognostic modeling
│       └── feature_selection_modeling.R  # Performs Bootstrap LASSO, random forest, Cox regression
│
└── Downstream_Analysis/
    │
    ├── Figure1/                  # Cell composition, annotation, and compartment-level analyses
    │   ├── Cell_Type_Annotation_in_WSI_Cohort.R       # Cell type annotation in WSI cohort
    │   ├── Cell_Type_Annotation_in_TMA1.R             # Cell type annotation in TMA cohort 1
    │   ├── Cell_Type_Annotation_in_TMA2.R             # Cell type annotation in TMA cohort 2
    │   ├── Cell_Type_Fraction_Comparison.R            # Compares cell-type fractions across groups
    │   ├── Cellular_Compartment_Comparison.R          # Compares major cell compartments
    │   └── markers.csv                                # Marker reference table for annotation
    │
    ├── Figure2/                  # Cell interaction, co-occurrence, and network motif analyses
    │   ├── Cell_Cell_Co_Occurrence.R                  # Quantifies spatial co-occurrence across samples
    │   ├── Cell_Cell_Interaction.txt                  # Generates cell-cell interaction analysis
    │   ├── Cell_Type_Clinical_Parameters.R            # Associates cell abundance with clinical parameters
    │   └── GNN_Network_Triplet_Motif_Analysis.txt     # GNN-based triplet motif and cell network analysis
    │
    ├── Figure3/                  # Neighborhood, distance, ecological unit, and proximity analyses
    │   ├── Cell_Cell_Distance_Analysis.R              # Analyzes spatial distances between cell types
    │   ├── Cell_State_Ecological_Unit_Analysis.R      # Defines and analyzes cell-state ecological units
    │   ├── Cellular_Neighborhood_Comparison.R         # Compares cellular neighborhood patterns
    │   ├── Cellular_Neighborhood_Generation.txt       # Generates cellular neighborhood annotations
    │   ├── Distance_and_Proximity_Density.txt         # Generates distance and proximity density
    │   └── Proximity_Density_Analysis.R               # Analyzes proximity density of cells across CNs
    │
    ├── Figure4/                  # CANVAS model performance and molecular/clinical associations
    │   ├── CANVAS_model.txt                          # CANVAS model architecture
    │   ├── CANVAS_training.txt                       # CANVAS training procedure and settings
    │   ├── CN_Assignment.txt                         # CN assignment workflow
    │   ├── CNV_and_SNV_Associations.R                # Associates habitats with CNV and SNV
    │   ├── Habitat_Prediction.txt                    # Habitat prediction workflow summary
    │   ├── Pathway_and_IO_Signature_Analysis.R       # Pathway and immuno-oncology signature analysis
    │   ├── RNA_Seq_Cell_Fraction_Analysis.R          # RNA-seq-based cell fraction analysis
    │   └── Training_Validation_and_Testing_Metrics.R # Summarizes model performance metrics
    │
    └── Figure5/                  # Habitat-derived prognostic feature construction and modeling
        ├── Feature_Dispersion.txt                    # Habitat spatial dispersion feature generation
        ├── Feature_Distance.txt                      # Distance feature generation
        ├── Feature_Diversity.R                       # Diversity feature analysis
        ├── Feature_Entropy.R                         # Entropy feature analysis
        ├── Feature_Interaction.txt                   # Habitat interaction feature generation
        ├── Feature_Relative_Abundance.R              # Relative abundance feature analysis
        ├── Feature_Selection.R                       # Feature selection workflow
        └── Prognostic_Model_Construction.R           # Prognostic model construction and evaluation
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
