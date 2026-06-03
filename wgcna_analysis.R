# WGCNA — Análisis de co-expresión génica durante la dormancia
# Input:  Expression_Profile_Nonpareil_gene.xlsx (TPM, 39 muestras)
# Output: figuras + tabla de módulos + hub genes
 
# 0. Librerías
library(WGCNA)
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()   # usa todos los cores disponibles
 
#1. Cargar datos
cat("Cargando datos de expresión...\n")
 
expr_raw <- read_excel("Expression_Profile_Nonpareil_gene.xlsx")
 
# Extraer matriz TPM (columnas 1_TPM ... 39_TPM)
tpm_cols <- grep("_TPM$", colnames(expr_raw), value = TRUE)
tpm_mat  <- as.matrix(expr_raw[, tpm_cols])
rownames(tpm_mat) <- expr_raw$Gene_ID
colnames(tpm_mat) <- gsub("_TPM", "", tpm_cols)
 
cat(sprintf("Genes totales: %d | Muestras: %d\n", nrow(tpm_mat), ncol(tpm_mat)))
 
# 2. Metadatos de muestras 
# Basado en metadatos.xlsx
# Muestras 1-39 en orden del archivo de expresión
meta <- data.frame(
  sample = as.character(1:39),
  dormancy = c(
    # Control (muestras 1-18)
    rep("EndoD", 9),   # 1-9
    rep("EndoR",  3),   # 10-12
    rep("EcoD",  6),   # 13-18
    # Asc (muestras 19-27)
    rep("EndoD", 3),   # 19-21
    rep("EndoR",  3),   # 22-24
    rep("EcoD",  3),   # 25-27
    # Eth (muestras 28-39)
    rep("EndoD", 6),   # 28-33
    rep("EndoR",  3),   # 34-36
    rep("EcoD",  3)    # 37-39
  ),
  treatment = c(
    rep("Control", 18),   # 1-18
    rep("Asc",      9),   # 19-27
    rep("Eth",     12)    # 28-39
  ),
  stringsAsFactors = FALSE
)
# Codificar dormancia como numérico para correlación con módulos
meta$dormancy_num <- as.numeric(factor(meta$dormancy,
                      levels = c("EndoD","EndoR","EcoD")))
 
# 3. Filtrado de genes
cat("Filtrando genes...\n")
 
# Filtro 1: TPM medio > 1 en al menos 5 muestras
expressed <- rowSums(tpm_mat > 1) >= 5
tpm_filt  <- tpm_mat[expressed, ]
cat(sprintf("Genes tras filtro TPM>1 en >=5 muestras: %d\n", nrow(tpm_filt)))
 
# Filtro 2: usar solo DEGs para un análisis más informativo
var_genes <- apply(tpm_filt, 1, var)
top_n     <- 5000   # top 5000 genes más variables
tpm_top   <- tpm_filt[order(var_genes, decreasing = TRUE)[1:top_n], ]
cat(sprintf("Genes seleccionados (top varianza): %d\n", nrow(tpm_top)))
 
# Transformación log
datExpr <- t(log2(tpm_top + 1))   # muestras x genes
cat(sprintf("Matriz final: %d muestras x %d genes\n", nrow(datExpr), ncol(datExpr)))
 
# 4. QC de muestras
cat("Control de calidad de muestras...\n")
 
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("Algunas muestras/genes eliminados en QC\n")
}
 
# Clustering de muestras para detectar outliers
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("WGCNA_01_sample_clustering.pdf", width = 12, height = 6)
par(cex = 0.7, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering (outlier detection)",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
dev.off()
cat("Figura 1 guardada: WGCNA_01_sample_clustering.pdf\n")
 
# 5. Selección del poder de umbral (soft-thresholding)
cat("Calculando poder de umbral (puede tardar 5-10 min)...\n")
 
powers  <- c(1:10, seq(12, 20, 2))
sft     <- pickSoftThreshold(datExpr, powerVector = powers,
                              verbose = 5, networkType = "signed")
 
pdf("WGCNA_02_soft_threshold.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "red", lty = 2)
 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = 0.9, col = "red")
dev.off()
cat("Figura 2 guardada: WGCNA_02_soft_threshold.pdf\n")
 
# Seleccionar el poder mínimo con R² >= 0.80
power_sel <- sft$powerEstimate
if (is.na(power_sel)) power_sel <- 12   # fallback si no converge
cat(sprintf("Poder seleccionado: %d\n", power_sel))
 
# 6. Construcción de la red y detección de módulos
cat("Construyendo red (puede tardar 10-20 min)...\n")
 
net <- blockwiseModules(
  datExpr,
  power             = power_sel,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)
 
moduleColors <- net$colors
cat(sprintf("Módulos detectados: %d\n", length(unique(moduleColors)) - 1))
cat("Tamaño de módulos:\n")
print(sort(table(moduleColors), decreasing = TRUE))
 
# Dendrograma
pdf("WGCNA_03_dendrogram_modules.pdf", width = 14, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
cat("Figura 3 guardada: WGCNA_03_dendrogram_modules.pdf\n")
 
#7. Correlación módulo-rasgo 
cat("Calculando correlación módulo-rasgo...\n")
 
MEs       <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs_ord   <- orderMEs(MEs)
 
# Traits: dormancy (numérico) y treatment (dummy)
traits <- data.frame(
  row.names  = meta$sample,
  dormancy   = meta$dormancy_num,
  Asc        = as.numeric(meta$treatment == "Asc"),
  Eth        = as.numeric(meta$treatment == "Eth")
)
# Alinear orden de muestras
traits <- traits[rownames(datExpr), , drop = FALSE]
 
moduleTraitCor  <- cor(MEs_ord, traits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
 
# Heatmap módulo-rasgo
pdf("WGCNA_04_module_trait_heatmap.pdf", width = 8, height = 10)
textMatrix <- paste0(round(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPval, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix    = moduleTraitCor,
  xLabels   = colnames(traits),
  yLabels   = rownames(moduleTraitCor),
  ySymbols  = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors    = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text  = 0.6,
  zlim      = c(-1, 1),
  main      = "Module-trait relationships\nPrunus dulcis dormancy"
)
dev.off()
cat("Figura 4 guardada: WGCNA_04_module_trait_heatmap.pdf\n")
 
# 8. Módulos más correlacionados con dormancia
cat("\n=== Módulos más correlacionados con dormancia ===\n")
dorm_cor <- moduleTraitCor[, "dormancy"]
dorm_p   <- moduleTraitPval[, "dormancy"]
top_modules <- sort(abs(dorm_cor), decreasing = TRUE)[1:8]
for (m in names(top_modules)) {
  cat(sprintf("  %s: r=%.3f, p=%.3e\n", m, dorm_cor[m], dorm_p[m]))
}
 
# 9. Hub genes por módulo
cat("Calculando hub genes...\n")
 
# Gene significance y module membership
geneModMembership <- as.data.frame(cor(datExpr, MEs_ord, use = "p"))
MMpvalue          <- as.data.frame(corPvalueStudent(
                       as.matrix(geneModMembership), nrow(datExpr)))
 
geneSig_dorm <- as.numeric(cor(datExpr,
                  traits$dormancy, use = "p"))
GSdorm_p     <- as.numeric(corPvalueStudent(geneSig_dorm, nrow(datExpr)))
 
# Tabla completa de genes con módulo, MM y GS
gene_info <- data.frame(
  Gene_ID          = colnames(datExpr),
  module           = moduleColors,
  GS_dormancy      = round(geneSig_dorm, 3),
  GS_dormancy_pval = round(GSdorm_p, 4),
  stringsAsFactors = FALSE
)
 
# Añadir MM para cada módulo
for (col in colnames(geneModMembership)) {
  gene_info[[paste0("MM_", col)]] <- round(geneModMembership[[col]], 3)
}
 
# Añadir descripción
gene_info <- left_join(gene_info,
               expr_raw[, c("Gene_ID","Description")],
               by = "Gene_ID")
 
write.csv(gene_info, "WGCNA_gene_module_membership.csv",
          row.names = FALSE)
cat("Tabla guardada: WGCNA_gene_module_membership.csv\n")
 
# 10. Hub genes top por módulo
# Para módulos correlacionados con dormancia
sig_modules <- names(dorm_cor)[dorm_p < 0.05]
cat(sprintf("\nMódulos con correlación dormancia significativa (p<0.05): %d\n",
            length(sig_modules)))
 
hub_list <- list()
for (mod in sig_modules) {
  mod_name <- gsub("ME", "", mod)
  mod_genes <- gene_info[gene_info$module == mod_name, ]
  mm_col    <- paste0("MM_", mod)
  if (mm_col %in% colnames(mod_genes)) {
    hub <- mod_genes[order(abs(mod_genes[[mm_col]]),
                           decreasing = TRUE), ][1:min(20, nrow(mod_genes)), ]
    hub$hub_module <- mod_name
    hub_list[[mod_name]] <- hub[, c("Gene_ID","module","GS_dormancy",
                                     mm_col,"Description")]
  }
}
 
if (length(hub_list) > 0) {
  hub_table <- do.call(rbind, hub_list)
  write.csv(hub_table, "WGCNA_hub_genes.csv", row.names = FALSE)
  cat("Hub genes guardados: WGCNA_hub_genes.csv\n")
  cat("\nTop 5 hub genes por módulo significativo:\n")
  for (mod in names(hub_list)) {
    cat(sprintf("\n  Módulo %s (r_dormancy=%.3f):\n", mod, dorm_cor[paste0("ME",mod)]))
    mm_col <- paste0("MM_ME", mod)
    if (mm_col %in% colnames(hub_list[[mod]])) {
      print(hub_list[[mod]][1:5, c("Gene_ID","GS_dormancy",mm_col,"Description")])
    }
  }
}
 
# 11. Resumen final
cat("\n", rep("=",60), "\n", sep="")
cat("RESUMEN WGCNA\n")
cat(rep("=",60), "\n", sep="")
cat(sprintf("Genes analizados:   %d\n", ncol(datExpr)))
cat(sprintf("Muestras:           %d\n", nrow(datExpr)))
cat(sprintf("Poder seleccionado: %d\n", power_sel))
cat(sprintf("Módulos totales:    %d\n", length(unique(moduleColors))-1))
cat(sprintf("Módulos sig. dormancia: %d\n", length(sig_modules)))
cat("\nFiguras generadas:\n")
cat("  WGCNA_01_sample_clustering.pdf\n")
cat("  WGCNA_02_soft_threshold.pdf\n")
cat("  WGCNA_03_dendrogram_modules.pdf\n")
cat("  WGCNA_04_module_trait_heatmap.pdf\n")
cat("\nTablas generadas:\n")
cat("  WGCNA_gene_module_membership.csv\n")
cat("  WGCNA_hub_genes.csv\n")