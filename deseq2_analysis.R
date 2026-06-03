# ANÁLISIS DE EXPRESIÓN DIFERENCIAL - ESTADO DE DORMANCIA
# DESeq2; Modelo: treatment + dormancy_state

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(apeglm)   # shrinkage para coeficientes directos del modelo
library(ashr)     # shrinkage para contrastes no directos (EcoD vs EndoR)


# 1. CARGA DE DATOS
# ajustamos la ruta del archivo
expr_file <- "Expression_Profile_Nonpareil_gene.xlsx"

raw <- read_excel(expr_file)

# Extraemos las columnas de read count (muestras 1-39)
count_cols <- paste0(1:39, "_Read_Count")
count_matrix <- as.matrix(raw[, count_cols])
rownames(count_matrix) <- raw$Gene_ID
colnames(count_matrix) <- paste0("S", 1:39)

# Asegurarnos de que son enteros (requisito de DESeq2)
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

message(paste("Genes totales:", nrow(count_matrix)))
message(paste("Muestras:", ncol(count_matrix)))

# -------

# 2. METADATA (colData)
col_data <- data.frame(
  sample_id  = paste0("S", 1:39),
  replicate  = c(rep("PC4",3), rep("PC13",3), rep("PC14",3),
                 rep("PC16",3), rep("PC18",3), rep("PC20",3),  # Control
                 rep("PA1",3), rep("PA2",3), rep("PA4",3),     # Asc
                 rep("PE1",3), rep("PE13",3), rep("PE15",3), rep("PE17",3)), # Eth
  treatment  = factor(c(rep("Control",18),
                        rep("Asc",9),
                        rep("Eth",12)),
                      levels = c("Control","Asc","Eth")),
  dormancy   = factor(c(rep("EndoD",9),   # PC4+PC13+PC14 Control
                        rep("EndoR",3),     # PC16 Control
                        rep("EcoD",6),     # PC18+PC20 Control
                        rep("EndoD",3),    # PA1 Asc
                        rep("EndoR",3),     # PA2 Asc
                        rep("EcoD",3),     # PA4 Asc
                        rep("EndoD",6),    # PE1+PE13 Eth
                        rep("EndoR",3),     # PE15 Eth
                        rep("EcoD",3)),    # PE17 Eth
                      levels = c("EndoD","EndoR","EcoD")),
  row.names = paste0("S", 1:39)
)

# verificación
message("Distribución de muestras por grupo:")
print(table(col_data$treatment, col_data$dormancy))

#-----

# 3. FILTRADO PREVIO
# Eliminar genes con expresión muy baja (ruido)
# Criterio: al menos 10 reads en al menos 3 muestras

keep <- rowSums(count_matrix >= 10) >= 3
count_filtered <- count_matrix[keep, ]
message(paste("Genes tras filtrado:", nrow(count_filtered),
              "(eliminados:", sum(!keep), ")"))


# 3b. EXCLUSIÓN DE MUESTRA OUTLIER
# S10 (PC16, réplica 1, Control EndoR) se excluye por tener pocos reads

count_filtered <- count_filtered[, colnames(count_filtered) != "S10"]
col_data       <- col_data[rownames(col_data) != "S10", ]

message(paste("Muestras tras excluir S10:", ncol(count_filtered)))

# ----------

# 4. OBJETO DESeq2 Y MODELO
# Modelo: treatment + dormancy_state
# El tratamiento se controla como covariable; el interés principal es la dormancia


dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData   = col_data,
  design    = ~ treatment + dormancy
)

dds <- DESeq(dds)

message("Factores en el modelo:")
print(resultsNames(dds))

# -------

# 5. PCA DE EXPRESIÓN
vsd <- vst(dds, blind = FALSE)  # Variance Stabilizing Transformation

pca_data <- plotPCA(vsd, intgroup = c("dormancy","treatment"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                               color = dormancy, shape = treatment)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = c("EndoD" = "#2166AC",
                                "EndoR"  = "#F4A582",
                                "EcoD"  = "#D6604D")) +
  labs(title = "PCA - Expresión génica (VST)",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"),
       color = "Dormancy state", shape = "Treatment") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave("PCA_expression.pdf", p_pca, width = 7, height = 5)
ggsave("PCA_expression.png", p_pca, width = 7, height = 5, dpi = 300)
message("PCA guardado.")

# -----------
# 6. ANÁLISIS DE EXPRESIÓN DIFERENCIAL - COMPARACIONES DE DORMANCIA
# Referencia: EndoD (nivel base del factor)
#
# NOTA SOBRE LFC SHRINKAGE:
# DESeq2 sin shrinkage produce LFC inflados para genes poco expresados
# (ej: LFC=29 es biológicamente imposible). lfcShrink() modera estos valores
# hacia cero según su nivel de evidencia estadística.
# - apeglm: para coeficientes directos del modelo (más preciso)
# - ashr:   para contrastes no directos como Eco vs EndoR

message("Coeficientes disponibles en el modelo:")
print(resultsNames(dds))

# 6a. EndoR vs EndoD (coeficiente directo -> apeglm)
res_EndoR_vs_EndoD_raw <- results(dds,
                               name  = "dormancy_EndoR_vs_EndoD",
                               alpha = 0.05)
res_EndoR_vs_EndoD <- lfcShrink(dds,
                              coef = "dormancy_EndoR_vs_EndoD",
                              type = "apeglm",
                              res  = res_ER_vs_Endo_raw)
res_EndoR_vs_EndoD <- res_EndoR_vs_EndoD[order(res_ER_vs_Endo$padj), ]

#6b. EcoD vs EndoD (coeficiente directo -> apeglm)
res_EcoD_vs_EndoD_raw <- results(dds,
                                name  = "dormancy_EcoD_vs_EndoD",
                                alpha = 0.05)
res_EcoD_vs_EndoD <- lfcShrink(dds,
                               coef = "dormancy_EcoD_vs_EndoD",
                               type = "apeglm",
                               res  = res_EcoD_vs_EndoD_raw)
res_EcoD_vs_EndoD <- res_EcoD_vs_EndoD[order(res_EcoD_vs_EndoD$padj), ]

# 6c. EcoD vs EndoR (contraste no directo -> ashr)
# No es un coeficiente del modelo, hay que usar type="ashr" con contrast
res_EcoD_vs_EndoR_raw <- results(dds,
                              contrast = c("dormancy","EcoD","EndoR"),
                              alpha    = 0.05)
res_EcoD_vs_EndoR <- lfcShrink(dds,
                             contrast = c("dormancy","EcoD","EndoR"),
                             type     = "ashr",
                             res      = res_EcoD_vs_EndoR_raw)
res_EcoD_vs_EndoR <- res_EcoD_vs_EndoR[order(res_EcoD_vs_EndoR$padj), ]

# Resumen de resultados
message("\n=== EndoR vs EndoD (lfcShrink apeglm) ===")
summary(res_ER_vs_Endo)
message("\n=== EcoD vs EndoD (lfcShrink apeglm) ===")
summary(res_Eco_vs_Endo)
message("\n=== EcoD vs EndoR (lfcShrink ashr) ===")
summary(res_EcoD_vs_EndoR)

# Comparar LFC antes/después del shrinkage para verificar el efecto
message("\n=== Efecto del shrinkage - Top 5 genes más cambiados ===")
top5 <- head(rownames(res_EcoD_vs_EndoD), 5)
comparison <- data.frame(
  Gene    = top5,
  LFC_raw = round(res_EcoD_vs_EndoD_raw[top5, "log2FoldChange"], 2),
  LFC_shr = round(res_EcoD_vs_EndoD[top5,     "log2FoldChange"], 2),
  padj    = signif(res_EcoD_vs_EndoD[top5, "padj"], 3)
)
print(comparison)

# --------------

# 7. EXPORTAR RESULTADOS A CSV

# función para exportar los resultados con anotación
export_results <- function(res, raw_data, filename,
                           lfc_threshold = 1, basemean_min = 0) {
  df <- as.data.frame(res)
  df$Gene_ID <- rownames(df)
  
  # Añadir anotación del archivo excel original
  annot <- raw_data[, c("Gene_ID","Gene_Symbol","Description","gene_biotype")]
  df <- merge(df, annot, by = "Gene_ID", all.x = TRUE)
  
  # Clasificar como UP/DOWN/NS usando los umbrales de esta comparación
  df$regulation <- ifelse(is.na(df$padj) | is.na(df$log2FoldChange), "NS",
                   ifelse(df$padj >= 0.05, "NS",
                   ifelse(df$baseMean < basemean_min, "NS",
                   ifelse(df$log2FoldChange >  lfc_threshold, "UP",
                   ifelse(df$log2FoldChange < -lfc_threshold, "DOWN", "NS")))))
  
  df <- df[order(df$padj), ]
  write.csv(df, filename, row.names = FALSE)
  message(paste("Guardado:", filename))
  
  sig <- df[df$regulation %in% c("UP","DOWN"), ]
  up   <- sum(df$regulation == "UP")
  down <- sum(df$regulation == "DOWN")
  message(paste0("  DEGs (padj<0.05, |LFC|>", lfc_threshold,
                 ", baseMean>=", basemean_min, "): UP=", up, " DOWN=", down))
  return(df)
}

df_ER   <- export_results(res_ER_vs_Endo,  raw,
                           "DEG_EndoR_vs_EndoD.csv",
                           lfc_threshold = 1, basemean_min = 0)

df_Eco1 <- export_results(res_Eco_vs_Endo, raw,
                           "DEG_EcoD_vs_EndoD.csv",
                           lfc_threshold = 1, basemean_min = 0)

# EcoD vs EndoR: criterios más estrictos por n=2 en EndoR Control
# tras excluir S10. El shrinkage (ashr) no es suficiente para estabilizar
# los LFC con tan pocas réplicas. Se aplica baseMean >= 50 y |LFC| > 2.
message("\nNOTA: EcoD vs EndoR usa criterios más estrictos (baseMean>=50, |LFC|>2)")
message("Motivo: EndoR Control tiene n=2 réplicas tras excluir S10 (outlier)")
df_Eco2 <- export_results(res_Eco_vs_ER,   raw,
                           "DEG_EcoD_vs_EndoR.csv",
                           lfc_threshold = 2, basemean_min = 50)

# ----------
# 8. VOLCANO PLOTS

make_volcano <- function(df, title, lfc_threshold = 1,
                         color_up="#D6604D", color_down="#2166AC") {
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  up   <- sum(df$regulation == "UP",   na.rm=TRUE)
  down <- sum(df$regulation == "DOWN", na.rm=TRUE)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                 color = regulation)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("UP"=color_up, "DOWN"=color_down, "NS"="grey70"),
                       labels = c(paste0("UP (n=", up, ")"),
                                  paste0("DOWN (n=", down, ")"),
                                  "NS")) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype="dashed", color="black", linewidth=0.4) +
    geom_hline(yintercept = -log10(0.05),
               linetype="dashed", color="black", linewidth=0.4) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)", color = "") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

v1 <- make_volcano(df_ER,   "EndoR vs EndoD",      lfc_threshold = 1)
v2 <- make_volcano(df_Eco1, "EcoD vs EndoD",       lfc_threshold = 1)
v3 <- make_volcano(df_Eco2, "EcoD vs EndoR\n(baseMean≥50, |LFC|>2)", lfc_threshold = 2)

ggsave("Volcano_EndoR_vs_EndoD.pdf", v1, width=6, height=5)
ggsave("Volcano_EcoD_vs_EndoD.pdf",  v2, width=6, height=5)
ggsave("Volcano_EcoD_vs_EndoR.pdf",   v3, width=6, height=5)


# -------------

# 9. HEATMAP - TOP 50 DEGs MÁS VARIABLES

# Unir DEGs significativos de todas las comparaciones
sig_genes <- unique(c(
  rownames(res_EndoR_vs_EndoD[!is.na(res_EndoR_vs_EndoD$padj) & 
                           res_EndoR_vs_EndoD$padj < 0.05 & 
                           abs(res_EndoR_vs_EndoD$log2FoldChange) > 1, ]),
  rownames(res_EcoD_vs_EndoD[!is.na(res_EcoD_vs_EndoD$padj) & 
                            res_EcoD_vs_EndoD$padj < 0.05 & 
                            abs(res_EcoD_vs_EndoD$log2FoldChange) > 1, ]),
  rownames(res_EcoD_vs_EndoR[!is.na(res_EcoD_vs_EndoR$padj) & 
                          res_EcoD_vs_EndoR$padj < 0.05 & 
                          abs(res_EcoD_vs_EndoR$log2FoldChange) > 1, ])
))

message(paste("\nTotal DEGs únicos en alguna comparación:", length(sig_genes)))

# Seleccionar top 50 por varianza si hay muchos
if (length(sig_genes) > 50) {
  vsd_mat <- assay(vsd)[sig_genes, ]
  gene_var <- apply(vsd_mat, 1, var)
  top_genes <- names(sort(gene_var, decreasing=TRUE))[1:50]
} else {
  top_genes <- sig_genes
}

if (length(top_genes) > 1) {
  mat_heatmap <- assay(vsd)[top_genes, ]
  mat_scaled  <- t(scale(t(mat_heatmap)))  # Z-score por gen
  
  # Anotación de columnas
  annot_col <- data.frame(
    Dormancy  = col_data$dormancy,
    Treatment = col_data$treatment,
    row.names = col_data$sample_id
  )
  
  annot_colors <- list(
    Dormancy  = c("EndoD"="#2166AC", "EndoR"="#F4A582", "EcoD"="#D6604D"),
    Treatment = c("Control"="#1B7837",      "Asc"="#762A83",         "Eth"="#E08214")
  )
  
  pdf("Heatmap_top_DEGs.pdf", width=12, height=10)
  pheatmap(mat_scaled,
           annotation_col  = annot_col,
           annotation_colors = annot_colors,
           show_rownames   = (length(top_genes) <= 50),
           show_colnames   = FALSE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           main = paste0("Top DEGs por estado de dormición (n=", length(top_genes), ")"))
  dev.off()
}

message("\n=== ANÁLISIS COMPLETADO ===")
message("Archivos generados:")
message("  - PCA_expression.pdf/.png")
message("  - DEG_EndoR_vs_EndoD.csv")
message("  - DEG_EcoD_vs_EndoD.csv")
message("  - DEG_EcoD_vs_EndoR.csv")
message("  - Volcano_*.pdf (3 comparaciones)")
message("  - Heatmap_top_DEGs.pdf")
