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
│   ├── deseq2_analysis.R                  # análisis de los datos rnaseq por deseq2
│   ├── go_enrichment_rnaseq.R             # análisis de enriquecimiento GO ORA
├── sRNA-seq/
│   ├── PCA_analysis_sRNAseq.R             # pca y análisis de los datos de srnaseq
│   ├── extraer_fasta_psRNATarget.R        # extraer los miRNA para psRNAtarget
│   ├── mapeo_KAI-PRUDU.R                  # mapeo para cambio de anotación
│   └── spearman_analysis.ipynb            # psRNATarget v2 + validación Spearman
├── WGBS/
│   ├── methylation_context_analysis.ipynb # extract_methylation_matrix_v3.ipynb
│   └── DMR_analysis.R                     # Welch t-test corrección BH
├── treatment/
│   ├── GO_enrichment_treatments.ipynb     # análisis GO ORA enrichment x tratamiento 
│   └── treatments_DESeq2_analysis.ipynb   # análisis deseq2 x tratamiento
├── integration/
│   ├── wgcna_analysis.R                   # análisis redes coexpresión génica
│   ├── mofa2_train.ipynb                  # entrenamiento modelo para mofa2 
│   └── mofa2_analysis.R                   # multi-omic factor analysis 2
└── data/

```




