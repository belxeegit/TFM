library(DESeq2)
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)

# 1. CARGAR DATOS

# Matriz de expresión: extraer solo Read_Counts
srna_raw <- read.table("1_MatureBased_expression.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Columnas de Read_Count: X1_Read_Count ... X39_Read_Count
count_cols <- grep("Read_Count", colnames(srna_raw), value = TRUE)
counts <- srna_raw[, count_cols]
rownames(counts) <- srna_raw$Mature_ID

# Renombrar columnas a IDs numéricos limpios
colnames(counts) <- gsub("X(\\d+)_Read_Count", "\\1", colnames(counts))


#2. METADATOS

meta_clean <- read_excel("metadatos.xlsx")

# Verificación crítica: deben ser exactamente 39 filas únicas
stopifnot(nrow(meta_clean) == 39)
stopifnot(!duplicated(meta_clean$ID))

# Ordenar metadatos igual que columnas de counts
meta_clean <- meta_clean %>% arrange(ID)
counts <- counts[, as.character(meta_clean$ID)]

# Limpiar nombres de columnas de factores
meta_clean$dormancy_state <- trimws(meta_clean$`Dormancy state`)
meta_clean$treatment       <- meta_clean$Treatment
meta_clean$treatment[is.na(meta_clean$treatment)] <- "Control"
meta_clean$treatment       <- factor(meta_clean$treatment, 
                                     levels = c("None", "Asc", "Eth"))
meta_clean$dormancy_state  <- factor(meta_clean$dormancy_state,
                                     levels = c("EndoD", 
                                                "EndoR", 
                                                "EcoD"))

# 3. FILTRADO DE BAJA EXPRESIÓN 

# Mantener miRNAs con >= 10 reads en al menos 3 muestras
keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, ]
cat("miRNAs antes del filtro:", nrow(counts), "\n")
cat("miRNAs tras el filtro:  ", nrow(counts_filt), "\n")

#4. DESEQ2 

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = meta_clean,
  design    = ~ dormancy_state + treatment
)

dds <- DESeq(dds)

#5. PCA DE CALIDAD 

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("dormancy_state", "treatment"), 
                    returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = dormancy_state, shape = treatment)) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  theme_bw(base_size = 13) +
  ggtitle("PCA - sRNA-seq (miRNAs conocidos, ppe)")

ggsave("PCA_sRNAseq.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# 6. ANÁLISIS DIFERENCIAL 
# Comparación principal: EcoD vs EndoD
res_eco_vs_endo <- results(dds, 
                           contrast = c("dormancy_state", 
                                        "EcoD", 
                                        "EndoD"),
                           alpha = 0.05)

res_eco_vs_endo <- res_eco_vs_endo[order(res_eco_vs_endo$padj), ]

# Resumen
summary(res_eco_vs_endo)

# Tabla de significativos
de_mirnas <- as.data.frame(res_eco_vs_endo) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

cat("DE-miRNAs significativos (padj<0.05, |LFC|>=1):", nrow(de_mirnas), "\n")
print(de_mirnas)

# Guardar resultados
write.csv(de_mirnas, "DE_miRNAs_EcoD_vs_EndoD.csv", quote = FALSE)
write.csv(as.data.frame(res_eco_vs_endo), "All_miRNAs_results.csv", quote = FALSE)

# Endodormancy release vs Endodormancy
res_release_vs_endo <- results(dds,
                               contrast = c("dormancy_state",
                                            "EndoR",
                                            "EndoD"),
                               alpha = 0.05)
res_release_vs_endo <- res_release_vs_endo[order(res_release_vs_endo$padj), ]

de_mirnas_release_vs_endo <- as.data.frame(res_release_vs_endo) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

cat("DE-miRNAs EndoR vs EndoD:", nrow(de_mirnas_release_vs_endo), "\n")

# Ecodormancy vs Endodormancy release
res_eco_vs_release <- results(dds,
                              contrast = c("dormancy_state",
                                           "EcoD",
                                           "EndoR"),
                              alpha = 0.05)
res_eco_vs_release <- res_eco_vs_release[order(res_eco_vs_release$padj), ]

de_mirnas_eco_vs_release <- as.data.frame(res_eco_vs_release) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

cat("DE-miRNAs EcoD vs EndoR:", nrow(de_mirnas_eco_vs_release), "\n")

# Guardar las tres comparaciones
write.csv(as.data.frame(res_release_vs_endo), 
          "All_miRNAs_EndoR_vs_EndoD.csv", quote = FALSE)
write.csv(as.data.frame(res_eco_vs_release),  
          "All_miRNAs_EcoD_vs_EndoR.csv", quote = FALSE)
write.csv(de_mirnas_release_vs_endo, 
          "DE_miRNAs_EndoR_vs_EndoD.csv", quote = FALSE)
write.csv(de_mirnas_eco_vs_release,  
          "DE_miRNAs_EcoD_vs_EndoR.csv", quote = FALSE)

library(ggVennDiagram)

# Listas de DE-miRNAs por comparación
venn_dormancy <- list(
  "EndoR vs EndoD"     = rownames(de_mirnas_release_vs_endo),
  "EcoD vs EndoD"              = rownames(de_mirnas),
  "EcoD vs EndoR"      = rownames(de_mirnas_eco_vs_release)
)

venn_plot <- ggVennDiagram(venn_dormancy, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "#f7f7f7", high = "#2171b5") +
  scale_color_manual(values = rep("grey40", 3)) +
  ggtitle("DE-miRNAs por transición de dormancia") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14))

ggsave("Venn_DE_miRNAs_dormancy.png", plot = venn_plot, 
       width = 8, height = 7, dpi = 300)




######LOS MAS ABUNDANTES
# Counts normalizados por DESeq2 (size-factor normalization)
counts_norm <- counts(dds, normalized = TRUE)  # dds, no dds_srna

# Media de expresión normalizada a través de todas las muestras
mean_expr <- rowMeans(counts_norm)

# Top 30 más expresados globalmente
top30 <- sort(mean_expr, decreasing = TRUE)[1:30]

# Construir tabla informativa
top30_df <- data.frame(
  miRNA        = names(top30),
  mean_norm_counts = round(top30, 1),
  is_DE_Eco_vs_Endo = names(top30) %in% rownames(de_mirnas),
  log2FC_Eco_vs_Endo = round(res_eco_vs_endo[names(top30), "log2FoldChange"], 2),
  padj_Eco_vs_Endo   = signif(res_eco_vs_endo[names(top30), "padj"], 3)
)

print(top30_df)
write.csv(top30_df, "Top30_expressed_miRNAs.csv", quote = FALSE, row.names = FALSE)

# Gráfico de barras horizontal - top 30
top30_df$miRNA <- factor(top30_df$miRNA, levels = rev(top30_df$miRNA))

ggplot(top30_df, aes(x = mean_norm_counts, y = miRNA, fill = is_DE_Eco_vs_Endo)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d62728"),
                    labels = c("No DE", "DE (Eco vs Endo)"),
                    name = "") +
  labs(x = "Counts normalizados (media)", y = NULL,
       title = "Top 30 miRNAs más expresados en yemas de almendro") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("Top30_expressed_miRNAs.png", width = 8, height = 8, dpi = 300)
