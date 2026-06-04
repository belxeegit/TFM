# MOFA2 - Análisis e interpretación de resultados
# Capas: RNA-seq ; miRNA ; CG_promoter ; CHG_promoter ; CHH_promoter
# Grupos: Control ; Asc (Ascórbico) ; Eth (Ethephon)

library(MOFA2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# 1. Cargar modelo
model <- load_model("mofa2_model.hdf5")

# Añadir metadataos de muestras al modelo
# Dormancy state por sample ID (IDs 1-39)
dormancy_map <- c(
  # Control
  "1"="Endodormancy","2"="Endodormancy","3"="Endodormancy",
  "4"="Endodormancy","5"="Endodormancy","6"="Endodormancy",
  "7"="Endodormancy","8"="Endodormancy","9"="Endodormancy",
  "10"="EndoRelease","11"="EndoRelease","12"="EndoRelease",
  "13"="Ecodormancy","14"="Ecodormancy","15"="Ecodormancy",
  "16"="Ecodormancy","17"="Ecodormancy","18"="Ecodormancy",
  # Asc
  "19"="Endodormancy","20"="Endodormancy","21"="Endodormancy",
  "22"="EndoRelease","23"="EndoRelease","24"="EndoRelease",
  "25"="Ecodormancy","26"="Ecodormancy","27"="Ecodormancy",
  # Eth
  "28"="Endodormancy","29"="Endodormancy","30"="Endodormancy",
  "31"="Endodormancy","32"="Endodormancy","33"="Endodormancy",
  "34"="EndoRelease","35"="EndoRelease","36"="EndoRelease",
  "37"="Ecodormancy","38"="Ecodormancy","39"="Ecodormancy"
)

samples_meta <- samples_metadata(model)
samples_meta$dormancy <- dormancy_map[samples_meta$sample]
samples_meta$dormancy <- factor(samples_meta$dormancy,
                                levels = c("EndoD","EndoR","EcoD"))
samples_metadata(model) <- samples_meta

col_dormancy <- c("EndoD"="#2166AC",
                  "EndoR"  ="#F4A736",
                  "EcoD"  ="#1A9641")
col_group    <- c("Control"="#555555","Asc"="#4DAF4A","Eth"="#E41A1C")

# 2. Varianza explicada

# 2a. Heatmap global: factores × vistas, por grupo
# Aquí no se usa plot_total = TRUE para que mantenga la columna 'factor'
p_var <- plot_variance_explained(
  model,
  x = "view", y = "factor" 
)
ggsave("mofa2_variance_explained.pdf", p_var, width=10, height=5)

# 2b. Varianza total por vista (barplot)
# Aquí usamos plot_total = TRUE pero ELIMINAMOS los argumentos x e y.
# MOFA2 ya sabe internamente que debe poner 'view' en x y 'variance' en y.
p_var_total <- plot_variance_explained(
  model,
  plot_total = TRUE
)
ggsave("mofa2_variance_total.pdf", p_var_total, width=8, height=4)

# 3. Scores de los factores
# 3a. Scatter Factor1 vs Factor2 coloreado por estado de dormancia
p_f1f2 <- plot_factors(
  model,
  factors = c(1, 2),
  color_by = "dormancy",
  shape_by  = "group",
  dot_size  = 3
) +
  scale_color_manual(values = col_dormancy, name = "Dormancy state") +
  scale_shape_manual(values = c(Control=16, Asc=17, Eth=15), name = "Treatment") +
  theme_classic(base_size = 12) +
  labs(title = "MOFA2 – Factor 1 vs Factor 2")
ggsave("mofa2_factors_F1F2_dormancy.pdf", p_f1f2, width=6, height=5)

# 3b. Scatter Factor1 vs Factor2 coloreado por grupo (tratamiento)
p_f1f2_group <- plot_factors(
  model,
  factors = c(1, 2),
  color_by = "group",
  dot_size  = 3
) +
  scale_color_manual(values = col_group, name = "Treatment") +
  theme_classic(base_size = 12) +
  labs(title = "MOFA2 – Factor 1 vs Factor 2 (by treatment)")
ggsave("mofa2_factors_F1F2_treatment.pdf", p_f1f2_group, width=6, height=5)

# 3c. Beeswarm de Factor 1 por estado de dormancia
p_bee1 <- plot_factor(
  model,
  factor = 1,
  color_by = "dormancy",
  dot_size = 3,
  stroke = 0.4
) +
  scale_color_manual(values = col_dormancy) +
  scale_fill_manual(values = col_dormancy) +
  theme_classic(base_size = 12) +
  labs(title = "Factor 1 – Dormancy gradient")
ggsave("mofa2_factor1_beeswarm.pdf", p_bee1, width=5, height=4)

# 3d. Panel con todos los factores vs dormancy
p_all_factors <- plot_factors(
  model,
  factors = 1:10,
  color_by = "dormancy"
) +
  scale_color_manual(values = col_dormancy)
ggsave("mofa2_all_factors_dormancy.pdf", p_all_factors, width=12, height=10)

# 4. Top pesos por factor y vista 
views_list <- c("RNA-seq","miRNA","CG_promoter","CHG_promoter","CHH_promoter")

# Función: top features para una vista y factor dados
plot_top_weights_custom <- function(model, view, factor, n=20) {
  p <- plot_top_weights(
    model,
    view   = view,
    factor = factor,
    nfeatures = n,
    scale  = TRUE,
    abs    = FALSE
  ) +
    theme_classic(base_size=10) +
    labs(title = paste0(view, " – Factor ", factor))
  return(p)
}

# Factor 1: top pesos en cada vista
for (v in views_list) {
  p <- plot_top_weights_custom(model, v, factor=1, n=20)
  fname <- paste0("mofa2_weights_F1_", gsub("[^a-zA-Z0-9]","_",v), ".pdf")
  ggsave(fname, p, width=7, height=5)
}

# Factor 2: top pesos en RNA-seq y miRNA
for (v in c("RNA-seq","miRNA")) {
  p <- plot_top_weights_custom(model, v, factor=2, n=20)
  fname <- paste0("mofa2_weights_F2_", gsub("[^a-zA-Z0-9]","_",v), ".pdf")
  ggsave(fname, p, width=7, height=5)
}

# 5. Limpiar feature names (quitar prefijo de vista) 
# Corrección del regex: ".*__" elimina todo hasta el doble guion bajo
clean_feature_name <- function(x) sub(".*__", "", x)

# 6. Tabla de pesos para genes candidatos estrella 
star_candidates <- c("PRUDU35267","PRUDU20869","PRUDU17011","PRUDU31204",
                     "PRUDU20878","PRUDU26110","PRUDU12379","PRUDU1176",
                     "PRUDU41546","PRUDU45936")

weights_rna <- get_weights(model, views="RNA-seq", as.data.frame=TRUE)
# Aplicamos la nueva función de limpieza
weights_rna$feature_clean <- clean_feature_name(weights_rna$feature)

# Ahora sí filtrará correctamente y pivotará sin errores
star_weights <- weights_rna %>%
  filter(feature_clean %in% star_candidates) %>%
  pivot_wider(names_from=factor, values_from=value) %>%
  arrange(desc(abs(Factor1)))

write.csv(star_weights, "mofa2_star_candidates_weights.csv", row.names=FALSE)
print("Top star candidate weights:")
print(star_weights)

# 7. Correlación factores con dormancy (test estadístico) 
# Extraer scores
Z <- get_factors(model, as.data.frame=TRUE)

# IMPORTANTE: Eliminamos la columna 'group' por defecto de Z para evitar 
# que colisione con la columna 'group' real de samples_meta al hacer el join.
if("group" %in% colnames(Z)) {
  Z <- Z %>% select(-group)
}

Z_meta <- left_join(Z, samples_meta, by="sample")

# ANOVA por factor - dormancy state
Z_ctrl <- Z_meta %>% filter(group == "Control")

anova_results <- Z_ctrl %>%
  group_by(factor) %>%
  summarise(
    p_anova = tryCatch(
      summary(aov(value ~ dormancy))[[1]][["Pr(>F)"]][1],
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(p_anova, method="BH"),
         significant = padj < 0.05)

print("ANOVA factor ~ dormancy (Control group):")
print(anova_results)
write.csv(anova_results, "mofa2_anova_factors_dormancy.csv", row.names=FALSE)

#8. Correlación entre factores
# 1. Abrimos el PDF
pdf("mofa2_factor_correlation.pdf", width=6, height=5)

# 2. Dibujamos el gráfico (esta función lo pinta directamente en el PDF)
plot_factor_cor(model) 

# 3. Cerramos y guardamos el archivo
dev.off()

#9. Heatmap de top features en RNA-seq para Factor 1
p_heat <- plot_data_heatmap(
  model,
  view   = "RNA-seq",
  factor = 1,
  features = 25,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_samples = "dormancy",
  annotation_colors = list(dormancy = col_dormancy),
  scale = "row"
)
ggsave("mofa2_heatmap_RNAseq_F1.pdf", p_heat, width=8, height=7)

# 10. Scatter de scores con gradiente de dormancia
# Figuras
dormancy_numeric <- c("Endodormancy"=1, "EndoRelease"=2, "Ecodormancy"=3)

Z_plot <- Z_meta %>%
  filter(factor %in% c("Factor1","Factor2")) %>%
  pivot_wider(names_from=factor, values_from=value) %>%
  mutate(dormancy_n = dormancy_numeric[as.character(dormancy)])

p_final <- ggplot(Z_plot, aes(x=Factor1, y=Factor2,
                               color=dormancy, shape=group)) +
  geom_point(size=3.5, stroke=0.5, alpha=0.9) +
  scale_color_manual(values=col_dormancy, name="Dormancy state") +
  scale_shape_manual(values=c(Control=16, Asc=17, Eth=15), name="Treatment") +
  stat_ellipse(aes(group=dormancy, color=dormancy), 
               level=0.80, linetype="dashed", linewidth=0.6) +
  theme_classic(base_size=13) +
  theme(legend.position="right",
        axis.title=element_text(size=12),
        plot.title=element_text(size=13, face="bold")) +
  labs(
    title = "MOFA2 multi-omic integration",
    subtitle = "RNA-seq | miRNA | CG/CHG/CHH promoter methylation",
    # CORRECCIÓN AQUÍ: usamos r2_per_factor en lugar de r2_total
    x = paste0("Factor 1 (", round(get_variance_explained(model)$r2_per_factor$Control["Factor1","RNA-seq"], 1), "% RNA-seq var.)"),
    y = "Factor 2"
  )
ggsave("mofa2_scatter_final.pdf", p_final, width=7, height=5.5)

cat("\nOK Análisis MOFA2 completado. Archivos generados:\n")
cat("  - mofa2_variance_explained.pdf\n")
cat("  - mofa2_factors_F1F2_dormancy.pdf\n")
cat("  - mofa2_factors_F1F2_treatment.pdf\n")
cat("  - mofa2_factor1_beeswarm.pdf\n")
cat("  - mofa2_weights_F1_[vista].pdf  (x5)\n")
cat("  - mofa2_heatmap_RNAseq_F1.pdf\n")
cat("  - mofa2_factor_correlation.pdf\n")
cat("  - mofa2_scatter_final.pdf\n")
cat("  - mofa2_star_candidates_weights.csv\n")
cat("  - mofa2_anova_factors_dormancy.csv\n")


############################################
# Inspección Factor 8
p_f8 <- plot_factor(model, factor=8, color_by="dormancy") +
  scale_color_manual(values=col_dormancy)
ggsave("mofa2_factor8_beeswarm.pdf", p_f8, width=5, height=4)

# Inspección Factor 2
samples_meta$collection_group <- ifelse(
  samples_meta$sample %in% c("1","2","3","28","29","30"),
  "Batch1_Nov21", "Batch2_Jan-Mar22"
)
samples_metadata(model) <- samples_meta
p_f2_batch <- plot_factor(model, factor=2, color_by="collection_group") +
  labs(title="Factor 2 – posible efecto de lote")
ggsave("mofa2_factor2_batch_check.pdf", p_f2_batch, width=5, height=4)
############################################


# 11. Extracción de Top Genes del Factor 1 para Análisis Funcional (GO/KEGG)

library(dplyr)
library(clusterProfiler)

# 1. Extraer TODOS los pesos del RNA-seq para el Factor 1
weights_F1 <- get_weights(model, views="RNA-seq", factors=1, as.data.frame=TRUE)
weights_F1$feature_clean <- clean_feature_name(weights_F1$feature)

# 2. Separar los drivers POSITIVOS (Promueven Endodormancia)
# Cogemos el Top 100 de genes con los valores más altos
top_positivos <- weights_F1 %>%
  filter(value > 0) %>%
  arrange(desc(value)) %>%
  slice_head(n = 100)

# 3. Separar los drivers NEGATIVOS (Promueven Ecodormancia)
# Cogemos el Top 100 de genes con los valores más bajos (más negativos)
top_negativos <- weights_F1 %>%
  filter(value < 0) %>%
  arrange(value) %>% # Ordenamos de más negativo a menos
  slice_head(n = 100)

# 4. Exportar listas limpias para usar en herramientas web (AgriGO, etc.) si lo prefieres
write.table(top_positivos$feature_clean, "mofa2_F1_Top100_Endodormancia.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(top_negativos$feature_clean, "mofa2_F1_Top100_Ecodormancia.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("\nOK Listas de genes Top 100 extraídas y guardadas en .txt\n")


