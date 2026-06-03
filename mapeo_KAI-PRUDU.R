library(dplyr)

# 1. Leer el GFF y extraer el mapeo KAI -> PRUDU
gff <- read.table("GCA_021292205.2_OSU_Pdul_2.5_genomic.gff",
                  sep="\t", comment.char="#", quote="",
                  col.names=c("seqid","source","type","start","end","score","strand","phase","attributes"))

# Filtrar solo CDS o genes
gff_cds <- gff[gff$type == "CDS", ]

# Extraer protein_id (KAI) y gene ID (PRUDU)
extract_attr <- function(attrs, key) {
  sub(paste0(".*", key, "=([^;]+).*"), "\\1", attrs)
}

gff_cds$protein_id <- extract_attr(gff_cds$attributes, "protein_id")
gff_cds$gene_id    <- extract_attr(gff_cds$attributes, "gene")

# Guardar tabla de mapeo
mapping <- gff_cds[, c("protein_id","gene_id")] %>% distinct()
write.csv(mapping, "KAI_to_PRUDU_mapping.csv", row.names=FALSE)
