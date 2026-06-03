# Trabajo Fin de Máster
## Caracterización multiόmica de la dormancia invernal en Prunus dulcis: identificación de biomarcadores moleculares 

**Autor: Juan José Belchí Navarro; Tutora: Raquel Sánchez Pérez; Cotutor: Ginés Almagro Hernández.**

Máster en Bioinformática, Universidad de Murcia, curso 2025-2026

### Resumen
Estudio multiómico de la dormancia invernal en yemas florales de almendro (_Prunus dulcis_ cv. Penta), integrando cuatro capas moleculares a lo largo de tres estados fisiológicos: Endodormancia (EndoD), Ruptura de endodormancia (EndoR) y Ecodormancia (EcoD), bajo tres tratamientos agroquímicos (Control, ácido ascórbico, ethephon).


### Estructura del repositorio
```
TFM/
├── RNA-seq/
│   ├── 01_alignment/          # Scripts HISAT2 v2.1.0
│   ├── 02_quantification/     # StringTie v2.1.3b
│   └── 03_DEG_analysis/       # pyDESeq2; modelo ~ treatment + dormancy_state
├── sRNA-seq/
│   ├── 01_mapping/            # Bowtie v1.1.2 + miRDeep2 v2.0.0.8
│   └── 02_target_prediction/  # psRNATarget v2 + validación Spearman
├── WGBS/
│   ├── 01_methylation_matrix/ # extract_methylation_matrix_v3.ipynb
│   └── 02_DMG_analysis/       # Welch t-test + corrección BH (padj<0.05, |Δ|>0.1)
├── metabolomics/
│   └── mummichog_enrichment/  # MetaboAnalyst; referencia Arabidopsis thaliana
├── integration/
│   ├── WGCNA/                 # Input: TPM de StringTie
│   └── MOFA2/                 # 5 capas: RNA-seq, miRNA, CG, CHG, CHH
├── figures/                   # Scripts Python/matplotlib; salida 600 DPI
└── data/
    └── Expression_Profile_Nonpareil_gene.xlsx
```




