# STRATA: Spatial Transcription-factor Regulatory Architecture of Tissue Analysis

A differential-geometric framework for continuous regulon field analysis of spatial transcriptomics data.

## Overview

STRATA constructs continuous regulon activity fields from transcript coordinates and computes their coupling tensor to quantify local co-regulation between transcription factor programs. It derives a Regulon Stability Index (RSI) from the Jacobian singular value decomposition, identifying **coupling phase boundaries** — positions where the regulatory logic of tissue changes.

## Key Results

Applied to 10x Genomics Xenium human skin melanoma data (382 genes, 13.7M transcripts):

- **Conjecture 5.1 validated**: Phase boundaries track histological tissue architecture (r = 0.32 with DEJ, r = 0.51 with σ₁, P < 10⁻¹⁰)
- **Parameter robust**: Correlations significant across all 12 bandwidth × coupling window combinations tested
- **Biological discovery**: Melanoma does not abolish regulon coupling but **homogenizes** it — coupling variance ↓28%, phase boundary intensity ↓18% relative to epidermal zone

## Framework Layers

| Layer | Quantity | Biological meaning |
|-------|----------|-------------------|
| L1 | Regulon fields φᵢ(s) | Continuous TF program activity at each tissue position |
| L2 | Coupling tensor C(s) | Local covariance structure between regulon programs |
| L3 | RSI & σ₁ from J(s) SVD | Stability landscape — where regulatory configuration is rigid vs labile |
| Phase boundaries | ‖∇C‖ | Where the regulatory logic of tissue changes |

## Quick Start
```python
# Run in Google Colab (T4 GPU, high RAM)
# Full pipeline: ~20 seconds on 13.7M transcripts
```

Open `STRATA_Layer1_RealSkin_20260224_clean.ipynb` in Google Colab and run all cells.

## Data

10x Genomics Xenium Human Skin (Melanoma) — CC BY 4.0
- [Dataset page](https://www.10xgenomics.com/datasets/human-skin-preview-data-xenium-human-skin-gene-expression-panel-add-on-1-standard)
- 382 genes (282 pre-designed + 100 custom add-on)
- 87,499 cells, 13.7M transcripts
- FFPE melanoma skin section

## Repository Structure
```
STRATA/
├── README.md
├── LICENSE
├── STRATA_Layer1_RealSkin_20260224_clean.ipynb  # Complete pipeline
├── strata/
│   ├── __init__.py
│   ├── fields.py          # KDE, regulon field construction
│   ├── coupling.py        # Coupling tensor, phase boundaries
│   └── stability.py       # Jacobian SVD, RSI, σ₁
└── paper/
    └── figure_plan.md
```

## Citation

If you use STRATA in your research, please cite:

> Jiu JW. STRATA: Spatial Transcription-factor Regulatory Architecture of Tissue Analysis. 2026. GitHub: jengweitjiu/STRATA

## Related Work

- **DGSA** — Geometric stability decomposition for cell state feature ablation (Bioinformatics, under review)
- **SICAI** — Stromal-immune coupling analysis in psoriasis (Nature Communications, under review)
- **IPA** — Intercellular communication buffering genetic perturbations (Nature Communications, under review)

## License

MIT License
