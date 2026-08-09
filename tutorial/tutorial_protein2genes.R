# ==============================================================================
# Protein marker annotation and feature harmonization
#
# Purpose:
#   Convert protein markers into corresponding gene symbols and Ensembl IDs for 
#   feature harmonization
#   with transcriptomic datasets.
#
# Important:
#   Protein markers are usually antibody names rather than official gene symbols.
#   Therefore, the recommended workflow is:
#
#   Protein marker
#          |
#          ↓
#   Manual curation of protein-to-gene mapping
#          |
#          ↓
#   biomaRt annotation (gene symbol -> Ensembl ID)
#
# For small numbers of markers (<100), manual annotation using UniProt/HGNC
# is recommended.
# For large-scale protein datasets, biomaRt, mygene, or UniProt REST API
# can be used.
#
# biomaRt documentation:
# https://github.com/Huber-group-EMBL/biomaRt
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Install and load packages
# ------------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}


library(biomaRt)
library(dplyr)



# ------------------------------------------------------------------------------
# 1. Input protein marker list
#
# Example:
# Protein markers from CODEX spatial proteomics
# ------------------------------------------------------------------------------

protein_marker_list <- c(
  "MUC2", "SOX9", "MUC1", "CD31", "Synapto",
  "CD49f", "CD15", "CHGA", "CDX2", "ITLN1",
  "CD4", "CD127", "Vimentin", "HLADR", "CD8",
  "CD11c", "CD44", "CD16", "BCL2", "CD3",
  "CD123", "CD38", "CD90", "aSMA", "CD21",
  "NKG2D", "CD66", "CD57", "CD206", "CD68",
  "CD34", "aDef5", "CD7", "CD36", "CD138",
  "CD45RO", "Cytokeratin", "CD117", "CD19",
  "Podoplanin", "CD45", "CD56", "CD69",
  "Ki67", "CD49a", "CD163", "CD161",
  "OLFM4", "FAP", "CD25", "CollIV", "CK7",
  "MUC6"
)



# ------------------------------------------------------------------------------
# 2. Manual curation:
#    Protein marker -> official gene symbol
#
# Note:
#   Many CODEX/IMC markers are antibody names:
#
#   CD31      -> PECAM1
#   Synapto   -> SYP
#   aSMA      -> ACTA2
#   HLADR     -> HLA-DRA
#   Cytokeratin -> KRT family
#
# For markers already equal to gene symbols,
# they can be directly retained.
# ------------------------------------------------------------------------------

manual_mapping <- data.frame(
  
  protein_marker = c(
    "CD31",
    "CD49f",
    "Synapto",
    "Vimentin",
    "HLADR",
    "aSMA",
    "CD206",
    "Cytokeratin",
    "Podoplanin",
    "Ki67",
    "CollIV",
    "CK7",
    "aDef5"
  ),
  
  gene_symbol = c(
    "PECAM1",
    "ITGA6",
    "SYP",
    "VIM",
    "HLA-DRA",
    "ACTA2",
    "MRC1",
    "KRT",
    "PDPN",
    "MKI67",
    "COL4A1",
    "KRT7",
    "DEFA5"
  ),
  
  stringsAsFactors = FALSE
)



# ------------------------------------------------------------------------------
# 3. Combine manually curated markers with directly matched markers
# ------------------------------------------------------------------------------

# Markers not manually annotated
direct_markers <- setdiff(
  protein_marker_list,
  manual_mapping$protein_marker
)


direct_mapping <- data.frame(
  protein_marker = direct_markers,
  gene_symbol = direct_markers,
  stringsAsFactors = FALSE
)


# Combine mapping table

protein_gene_mapping <- bind_rows(
  manual_mapping,
  direct_mapping
)


# Remove duplicated entries

protein_gene_mapping <- protein_gene_mapping %>%
  distinct()



# View mapping

print(protein_gene_mapping)



# ------------------------------------------------------------------------------
# 4. Connect to Ensembl database
#
# Human:
#     hsapiens_gene_ensembl
#
# Mouse:
#     mmusculus_gene_ensembl
#
# ------------------------------------------------------------------------------

options(timeout = 3000)


ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)



# ------------------------------------------------------------------------------
# 5. Query Ensembl annotation
#
# Convert:
#
# gene symbol
#      |
#      ↓
# Ensembl gene ID
#      +
# UniProt accession
#
# ------------------------------------------------------------------------------

annotation <- getBM(
  
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "uniprot_swissprot",
    "gene_biotype"
  ),
  
  filters = "hgnc_symbol",
  
  values = unique(protein_gene_mapping$gene_symbol),
  
  mart = ensembl
)



# Remove empty gene symbols

annotation <- annotation %>%
  filter(hgnc_symbol != "")



# ------------------------------------------------------------------------------
# 6. Merge protein marker information
# ------------------------------------------------------------------------------

final_map <- protein_gene_mapping %>%
  
  left_join(
    annotation,
    by = c(
      "gene_symbol" = "hgnc_symbol"
    )
  )



# Rename columns

final_map <- final_map %>%
  
  rename(
    ensembl_gene_id = ensembl_gene_id,
    uniprot_id = uniprot_swissprot
  )



# ------------------------------------------------------------------------------
# 7. Check unmapped markers
# ------------------------------------------------------------------------------

unmapped <- final_map %>%
  
  filter(
    is.na(ensembl_gene_id)
  )



cat(
  sprintf(
    "Total protein markers: %d\n",
    length(protein_marker_list)
  )
)


cat(
  sprintf(
    "Mapped markers: %d\n",
    sum(!is.na(final_map$ensembl_gene_id))
  )
)


cat(
  sprintf(
    "Unmapped markers: %d\n",
    nrow(unmapped)
  )
)



if (nrow(unmapped) > 0) {
  
  cat("\nMarkers requiring manual inspection:\n")
  
  print(
    unmapped$protein_marker
  )
  
}



# ------------------------------------------------------------------------------
# 8. Export annotation table
# ------------------------------------------------------------------------------

output_file <- "protein_marker_gene_annotation.csv"


write.csv(
  final_map,
  output_file,
  row.names = FALSE
)


cat(
  sprintf(
    "\nAnnotation table saved to: %s\n",
    output_file
  )
)



# ==============================================================================
# Output columns:
#
# protein_marker
#     Original protein/antibody marker name
#
# gene_symbol
#     Corresponding official gene symbol
#
# ensembl_gene_id
#     Ensembl gene identifier
#
# uniprot_id
#     UniProt accession (if available)
#
# gene_biotype
#     Gene annotation
#
# ==============================================================================