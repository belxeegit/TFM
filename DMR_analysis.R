#!/usr/bin/env Rscript
# DMR — Análisis de metilación diferencial
# Comparación principal: Ecodormancy vs Endodormancy
# Contextos: CG, CHG, CHH | Regiones: promotor y gene body

library(BiocParallel)

# CONFIGURACIÓN — AJUSTAR ESTAS RUTAS
INPUT_DIR  <- "."          # carpeta con los CSV de metilación
OUTPUT_DIR <- "./DSS_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Parámetros DSS
DELTA     <- 0.1    # diferencia mínima de metilación para llamar DMR
P_CUTOFF  <- 0.05   # p-value para DML
MIN_CG    <- 3      # mínimo de CpGs en un DMR
MIN_LEN   <- 50     # longitud mínima de DMR en bp

# ------------
# METADATOS — muestras excluidas: 3 y 38


meta <- data.frame(
  ID = c(1,2,4,5,6,7,8,9,          # EndoD control (sin muestra 3)
         10,11,12,                   # EndoR control
         13,14,15,16,17,18,          # EcoD control
         19,20,21,                   # EndoD Asc
         22,23,24,                   # EndoR Asc
         25,26,27,                   # EcoD Asc
         28,29,30,31,32,33,          # EndoD Eth
         34,35,36,                   # EndoR Eth
         37,39),                     # EcoD Eth (sin muestra 38)
  state = c(rep("EndoD", 8),
            rep("EndoR", 3),
            rep("EcoD", 6),
            rep("EndoD", 3),
            rep("EndoR", 3),
            rep("EcoD", 3),
            rep("EndoD", 6),
            rep("EndoR", 3),
            rep("EcoD", 2)),
  stringsAsFactors = FALSE
)

# Índices por estado
idx_endo <- which(meta$state == "EndoD")   # n=17
idx_er   <- which(meta$state == "EndoR")    # n=9
idx_eco  <- which(meta$state == "EcoD")    # n=11

cat("Muestras EndoD:  ", length(idx_endo), "\n")
cat("Muestras EndoR:   ", length(idx_er),   "\n")
cat("Muestras EcoD:   ", length(idx_eco),  "\n")

# Comparaciones a realizar
COMPARISONS <- list(
  list(name = "EcoD_vs_EndoD", idx_A = idx_eco,  idx_B = idx_endo, label_A = "EcoD",  label_B = "EndoD"),
  list(name = "EndoR_vs_EndoD",  idx_A = idx_er,   idx_B = idx_endo, label_A = "EndoR",   label_B = "EndoD"),
  list(name = "EcoD_vs_EndoR",   idx_A = idx_eco,  idx_B = idx_er,   label_A = "EcoD",  label_B = "EndoR")
)

# ---------------
# FUNCIÓN PRINCIPAL


run_dss <- function(context, region, comp) {

  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("Contexto: %s | Región: %s | Comparación: %s\n",
              context, region, comp$name))
  cat(strrep("=", 60), "\n")

  # --- Leer matriz de metilación ---
  infile <- file.path(INPUT_DIR, sprintf("methyl_%s_%s.csv", context, region))
  if (!file.exists(infile)) {
    cat("AVISO: archivo no encontrado:", infile, "\n")
    return(NULL)
  }

  mat <- read.csv(infile, row.names = 1, check.names = FALSE)
  mat <- mat[as.character(meta$ID), ]
  cat(sprintf("Matriz cargada: %d muestras x %d genes\n", nrow(mat), ncol(mat)))

  genes  <- colnames(mat)
  n_genes <- length(genes)
  results <- vector("list", n_genes)

  for (i in seq_along(genes)) {
    g <- genes[i]
    vals_A <- as.numeric(mat[comp$idx_A, g])
    vals_B <- as.numeric(mat[comp$idx_B, g])

    vals_A <- vals_A[!is.na(vals_A)]
    vals_B <- vals_B[!is.na(vals_B)]

    if (length(vals_A) < 3 || length(vals_B) < 3) next

    tt <- tryCatch(
      t.test(vals_A, vals_B, var.equal = FALSE),
      error = function(e) NULL
    )
    if (is.null(tt)) next

    results[[i]] <- data.frame(
      gene_id          = g,
      mean_A           = mean(vals_A),
      mean_B           = mean(vals_B),
      delta            = mean(vals_A) - mean(vals_B),
      pvalue           = tt$p.value,
      n_A              = length(vals_A),
      n_B              = length(vals_B),
      stringsAsFactors = FALSE
    )

    if (i %% 5000 == 0) cat(sprintf("  Procesados: %d/%d\n", i, n_genes))
  }

  res <- do.call(rbind, results[!sapply(results, is.null)])
  cat(sprintf("Genes testados: %d\n", nrow(res)))

  # Renombrar columnas con etiquetas de la comparación
  names(res)[names(res) == "mean_A"] <- paste0("mean_", comp$label_A)
  names(res)[names(res) == "mean_B"] <- paste0("mean_", comp$label_B)

  # FDR
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$padj), ]

  # DMGs
  dmgs <- res[res$padj < 0.05 & abs(res$delta) > DELTA, ]
  cat(sprintf("DMGs (padj<0.05, |delta|>%.1f): %d\n", DELTA, nrow(dmgs)))
  cat(sprintf("  Hipermetilados en %s (delta>0): %d\n",
              comp$label_A, sum(dmgs$delta > 0)))
  cat(sprintf("  Hipometilados en %s (delta<0):  %d\n",
              comp$label_A, sum(dmgs$delta < 0)))

  # Guardar
  outfile_all <- file.path(OUTPUT_DIR,
    sprintf("DSS_%s_%s_%s_all.csv", context, region, comp$name))
  write.csv(res, outfile_all, row.names = FALSE)
  cat(sprintf("Guardado: %s\n", outfile_all))

  if (nrow(dmgs) > 0) {
    outfile_dmg <- file.path(OUTPUT_DIR,
      sprintf("DSS_%s_%s_%s_DMGs.csv", context, region, comp$name))
    write.csv(dmgs, outfile_dmg, row.names = FALSE)
    cat(sprintf("Guardado: %s\n", outfile_dmg))
  }

  return(list(all = res, dmgs = dmgs))
}

# ----------------
# EJECUTAR PARA LOS 6 CONTEXTOS × REGIONES

all_results <- list()

for (ctx in c("CG", "CHG", "CHH")) {
  for (reg in c("promoter", "genebody")) {
    for (comp in COMPARISONS) {
      key <- paste(ctx, reg, comp$name, sep = "_")
      all_results[[key]] <- run_dss(ctx, reg, comp)
    }
  }
}

# ------------
# RESUMEN FINAL

cat("\n", strrep("=", 60), "\n")
cat("RESUMEN FINAL\n")
cat(strrep("=", 60), "\n")
cat(sprintf("%-35s %8s %8s\n", "Análisis", "Testados", "DMGs"))
for (key in names(all_results)) {
  if (is.null(all_results[[key]])) next
  res  <- all_results[[key]]$all
  dmgs <- all_results[[key]]$dmgs
  cat(sprintf("%-35s %8d %8d\n", key, nrow(res), nrow(dmgs)))
}
