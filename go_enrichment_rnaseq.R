# GO ENRICHMENT ANALYSIS
# Anotaciones: P. dulcis Nonpareil v1.0 (GDR) via InterProScan
# Anotaciones nativas de almendro
# Referencia: D'Amico-Willman et al. 2022 G3 (Bethesda) 12(5):jkac065

library(clusterProfiler)
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ggpubr)

# 1. CONSTRUIR TABLA gene2GO DESDE ANOTACIÓN NATIVA DE ALMENDRO

go_raw <- read_excel("pdulcis_Nonpareil_v1_0_genes2Go.xlsx", skip = 2,
                     col_names = c("Query", "GO_ID", "Description"))

# Extraer PRUDU ID limpio: "rna-PRUDU1000.1_pno_v1" -> "PRUDU1000"
go_raw <- go_raw %>%
  filter(!is.na(Query), !is.na(GO_ID)) %>%
  filter(str_starts(GO_ID, "GO:")) %>%
  mutate(Gene_ID = str_extract(Query, "PRUDU[0-9]+")) %>%
  filter(!is.na(Gene_ID))

message(paste("Genes con GO terms:", length(unique(go_raw$Gene_ID))))
message(paste("Total anotaciones GO:", nrow(go_raw)))

# Tabla TERM2GENE: GO_ID -> Gene_ID (formato clusterProfiler)
term2gene <- go_raw %>%
  select(GO_ID, Gene_ID) %>%
  distinct()

# Tabla TERM2NAME: GO_ID -> descripción legible
term2name <- go_raw %>%
  select(GO_ID, Description) %>%
  mutate(Description = str_remove(Description, "^(Biological Process|Molecular Function|Cellular Component):")) %>%
  distinct(GO_ID, .keep_all = TRUE)

# Separar por ontología para análisis independiente
term2gene_BP <- term2gene %>%
  filter(GO_ID %in% go_raw$GO_ID[str_detect(go_raw$Description, "Biological Process")])
term2gene_MF <- term2gene %>%
  filter(GO_ID %in% go_raw$GO_ID[str_detect(go_raw$Description, "Molecular Function")])
term2gene_CC <- term2gene %>%
  filter(GO_ID %in% go_raw$GO_ID[str_detect(go_raw$Description, "Cellular Component")])

term2name_BP <- term2name %>% filter(GO_ID %in% term2gene_BP$GO_ID)
term2name_MF <- term2name %>% filter(GO_ID %in% term2gene_MF$GO_ID)
term2name_CC <- term2name %>% filter(GO_ID %in% term2gene_CC$GO_ID)

message(paste("  Términos BP:", length(unique(term2gene_BP$GO_ID))))
message(paste("  Términos MF:", length(unique(term2gene_MF$GO_ID))))
message(paste("  Términos CC:", length(unique(term2gene_CC$GO_ID))))


# -----

# 2. CARGAR DEGs
deg_ER  <- read.csv("DEG_EndoR_vs_EndoD.csv")
deg_Eco <- read.csv("DEG_EcoD_vs_EndoD.csv")

# Background: todos los genes testados en DESeq2
universe <- unique(c(deg_ER$Gene_ID, deg_Eco$Gene_ID))
message(paste("\nGenes en universe:", length(universe)))
message(paste("Universe con GO terms:", sum(universe %in% go_raw$Gene_ID),
              paste0("(", round(100*sum(universe %in% go_raw$Gene_ID)/length(universe),1), "%)")))

# -----------

# 3. FUNCIÓN GO ENRICHMENT CON ANOTACIONES NATIVAS

run_go_native <- function(gene_ids, universe_ids, t2g, t2n, label) {

  if (length(gene_ids) < 5) {
    message(paste("  Saltando", label, "- muy pocos genes"))
    return(NULL)
  }

  n_with_go <- sum(gene_ids %in% t2g$Gene_ID)
  message(paste0("  ", label, ": ", n_with_go, "/", length(gene_ids),
                 " genes con GO terms (", round(100*n_with_go/length(gene_ids),1), "%)"))

  tryCatch({
    result <- enricher(
      gene          = gene_ids,
      universe      = universe_ids,
      TERM2GENE     = t2g,
      TERM2NAME     = t2n,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2,
      minGSSize     = 5,
      maxGSSize     = 500
    )

    if (!is.null(result) && nrow(result@result) > 0) {
      n_sig <- sum(result@result$p.adjust < 0.05)
      message(paste0("    -> ", n_sig, " términos significativos (padj < 0.05)"))
      return(result)
    } else {
      message("    -> Sin términos significativos")
      return(NULL)
    }

  }, error = function(e) {
    message(paste("    ERROR:", e$message))
    return(NULL)
  })
}

# -----------------

# 4. EJECUTAR ENRICHMENT PARA AMBAS COMPARACIONES
comparisons <- list(
  list(name="Eco_vs_Endo",  deg=deg_Eco,
       label_up="Activated in Ecodormancy",
       label_down="Repressed in Ecodormancy"),
  list(name="ER_vs_Endo",   deg=deg_ER,
       label_up="Activated in EndoRelease",
       label_down="Repressed in EndoRelease")
)

all_results <- list()

for (comp in comparisons) {
  message(paste0("\n", rep("=",50)))
  message(paste("Procesando:", comp$name))

  up_genes   <- comp$deg$Gene_ID[comp$deg$regulation == "UP"]
  down_genes <- comp$deg$Gene_ID[comp$deg$regulation == "DOWN"]
  all_genes  <- comp$deg$Gene_ID[comp$deg$regulation %in% c("UP","DOWN")]

  message(paste("  UP:", length(up_genes), "| DOWN:", length(down_genes)))

  # Biological Process
  message("\n-- Biological Process --")
  all_results[[paste0(comp$name, "_BP_UP")]]   <- run_go_native(up_genes,   universe, term2gene_BP, term2name_BP, "UP")
  all_results[[paste0(comp$name, "_BP_DOWN")]] <- run_go_native(down_genes, universe, term2gene_BP, term2name_BP, "DOWN")
  all_results[[paste0(comp$name, "_BP_ALL")]]  <- run_go_native(all_genes,  universe, term2gene_BP, term2name_BP, "ALL")

  # Molecular Function
  message("\n-- Molecular Function --")
  all_results[[paste0(comp$name, "_MF_UP")]]   <- run_go_native(up_genes,   universe, term2gene_MF, term2name_MF, "UP")
  all_results[[paste0(comp$name, "_MF_DOWN")]] <- run_go_native(down_genes, universe, term2gene_MF, term2name_MF, "DOWN")

  # Cellular Component
  message("\n-- Cellular Component --")
  all_results[[paste0(comp$name, "_CC_UP")]]   <- run_go_native(up_genes,   universe, term2gene_CC, term2name_CC, "UP")
  all_results[[paste0(comp$name, "_CC_DOWN")]] <- run_go_native(down_genes, universe, term2gene_CC, term2name_CC, "DOWN")
}

# ------------------

# 5. FIGURAS - DOTPLOTS
make_dotplot <- function(result_obj, title_text, top_n = 20,
                          color_low = "#4393C3", color_high = "#D6604D") {
  if (is.null(result_obj)) return(NULL)

  df <- result_obj@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      Description = str_wrap(Description, 45),
      GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text=x)))
    )

  if (nrow(df) == 0) return(NULL)

  # Ordenar por GeneRatio
  df$Description <- factor(df$Description,
                            levels = df$Description[order(df$GeneRatio_num)])

  ggplot(df, aes(x = GeneRatio_num, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = color_low, high = color_high,
                         name = "p.adjust",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) +
    labs(title = title_text, x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      axis.text.y   = element_text(size = 9),
      legend.position = "right"
    )
}

# --- Figura principal: Eco vs Endo BP UP + DOWN combinado ---
p_bp_up <- make_dotplot(
  all_results[["EcoD_vs_EndoD_BP_UP"]],
  "GO:BP — Activated genes\n(Ecodormancy vs Endodormancy)",
  top_n = 20, color_low = "#FDBB84", color_high = "#D7301F"
)

p_bp_down <- make_dotplot(
  all_results[["EcoD_vs_EndoD_BP_DOWN"]],
  "GO:BP — Repressed genes\n(Ecodormancy vs Endodormancy)",
  top_n = 20, color_low = "#9ECAE1", color_high = "#08519C"
)

if (!is.null(p_bp_up) && !is.null(p_bp_down)) {
  combined <- ggarrange(p_bp_up, p_bp_down, ncol = 2, labels = c("A", "B"))
  ggsave("GO_EcoD_vs_EndoD_BP_combined.pdf", combined, width = 20, height = 9)
  ggsave("GO_EcoD_vs_EndoD_BP_combined.png", combined, width = 20, height = 9, dpi = 300)
  message("Figura combinada BP guardada.")
}

# Guardar dotplots individuales para todos los resultados
plot_config <- list(
  list(key="EcoD_vs_EndoD_BP_ALL",   w=11, h=9,
       title="GO:BP — All DEGs (Ecodormancy vs Endodormancy)"),
  list(key="EcoD_vs_EndoD_MF_UP",    w=11, h=8,
       title="GO:MF — Activated genes (Ecodormancy vs Endodormancy)"),
  list(key="EcoD_vs_EndoD_MF_DOWN",  w=11, h=8,
       title="GO:MF — Repressed genes (Ecodormancy vs Endodormancy)"),
  list(key="EcoD_vs_EndoD_CC_UP",    w=10, h=7,
       title="GO:CC — Activated genes (Ecodormancy vs Endodormancy)"),
  list(key="EndoR_vs_EndoD_BP_UP",     w=11, h=8,
       title="GO:BP — Activated genes (EndoRelease vs Endodormancy)"),
  list(key="EndoR_vs_EndoD_BP_DOWN",   w=11, h=8,
       title="GO:BP — Repressed genes (EndoRelease vs Endodormancy)")
)

for (cfg in plot_config) {
  p <- make_dotplot(all_results[[cfg$key]], cfg$title)
  if (!is.null(p)) {
    ggsave(paste0(cfg$key, ".pdf"), p, width = cfg$w, height = cfg$h)
    ggsave(paste0(cfg$key, ".png"), p, width = cfg$w, height = cfg$h, dpi = 300)
    message(paste("Guardado:", cfg$key))
  }
}

# -----------

# 6. EXPORTAR TABLAS COMPLETAS
export_results <- function(results_list, pattern, outfile) {
  keys <- names(results_list)[str_detect(names(results_list), pattern)]
  dfs  <- list()
  for (k in keys) {
    res <- results_list[[k]]
    if (!is.null(res)) {
      df <- res@result %>%
        filter(p.adjust < 0.05) %>%
        mutate(analysis = k)
      dfs[[k]] <- df
    }
  }
  if (length(dfs) > 0) {
    out <- bind_rows(dfs) %>% arrange(analysis, p.adjust)
    write.csv(out, outfile, row.names = FALSE)
    message(paste("Tabla guardada:", outfile, "-", nrow(out), "términos"))
    return(out)
  }
  return(NULL)
}

tbl_eco <- export_results(all_results, "EcoD_vs_EndoD", "GO_EcoD_vs_EndoD_results.csv")
tbl_er  <- export_results(all_results, "EndoR_vs_EndoD",  "GO_EndoR_vs_EndoD_results.csv")

# -----------

# 7. RESUMEN FINAL
message("\n", rep("=", 60))
message("RESUMEN GO ENRICHMENT")
message(rep("=", 60))
for (key in names(all_results)) {
  res <- all_results[[key]]
  if (!is.null(res)) {
    n <- sum(res@result$p.adjust < 0.05)
    message(sprintf("  %-35s  %3d términos significativos", key, n))
  }
}

message("\nArchivos generados:")
message("  GO_EcoD_vs_EndoD_BP_combined.pdf/.png")
message("  GO_EcoD_vs_EndoD_BP_results.csv")
message("  GO_EndoR_vs_EndoD_results.csv")

