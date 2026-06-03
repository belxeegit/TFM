# Instalar en caso de ser necesario
install.packages("Biostrings") # o via BiocManager

library(Biostrings)

# Leer el FASTA limpio
mirnas <- readLines("mature.fa.txt")

# DE-miRNAs significativos (los de las 3 comparaciones)
de_names <- c("ppe-miR159", "ppe-miR319a", "ppe-miR167d",
              "ppe-miR166a", "ppe-miR166b", "ppe-miR166c",
              "ppe-miR166d", "ppe-miR166e", "ppe-miR164a",
              "ppe-miR164b", "ppe-miR164c", "ppe-miR156a",
              "ppe-miR156b", "ppe-miR156f", "ppe-miR156g",
              "ppe-miR156h", "ppe-miR156i", "ppe-miR398a-3p",
              "ppe-miR398b", "ppe-miR398a-5p", "ppe-miR160a",
              "ppe-miR160b", "ppe-miR482d-5p", "ppe-miR482d-3p",
              "ppe-miR393a", "ppe-miR393b", "ppe-miR395a-3p",
              "ppe-miR395b-3p", "ppe-miR395c", "ppe-miR3627-5p",
              "ppe-miR6285", "ppe-miR7122a-5p")

# Extraer headers e índices
header_lines <- which(startsWith(mirnas, ">"))
selected <- c()

for (name in de_names) {
  idx <- header_lines[grepl(name, mirnas[header_lines], fixed = TRUE)]
  if (length(idx) > 0) {
    selected <- c(selected, idx, idx + 1)
  }
}

selected <- sort(unique(selected))
writeLines(mirnas[selected], "DE_miRNAs_for_psRNATarget.fasta")
cat("Entradas guardadas:", length(selected)/2, "\n")
