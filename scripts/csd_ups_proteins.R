## Charlotte Dawson 2023
## This script generates a better formatted UPS proteins fasta.

## Start script
library(here)
library(dplyr)
library(tidyr)
library(Biostrings)
library(stringr)
library(uniprotREST)

ups <- readAAStringSet(here("ups/ups1-ups2-sequences.fasta"))

ups_headers <- sub("ups", "", names(ups), ignore.case = FALSE)
ups_headers <- sub("UPS ", "UPS\\|", ups_headers)
ups_headers <- sub("HUMAN_UPS", "UPS", ups_headers, ignore.case = FALSE)
ups_headers <- sub("- Homo sapiens \\(Human\\)", "", ups_headers, ignore.case = FALSE)
ups_headers <- str_trim(ups_headers)
ups_headers <- data.frame(header = ups_headers) %>% 
  separate_wider_delim(header, delim = "|", names = c("Accession", "Entry.name", "Description"))
ups_headers

ups_accessions <- sub("^(.+?)\\|.+", "\\1", ups_headers$Accession)
ups_accessions

# Every UPS accession still maps to a valid entry in 2023-12-30, except P62988.
# Note that P62988 (UBIQ) was made obsolete in 2010 and was de-merged into 
# P0CG47, P0CG48, P62979, and P62987.
# Note that 

actual_sequences <- uniprot_map(
  ups_accessions,
  format = "tsv",
  fields = c("accession", "id", "gene_primary"),
  isoform = TRUE
)
actual_sequences

actual_sequences2 <- actual_sequences %>% 
  `colnames<-`(c("From", "Entry", "Entry.name", "Gene.name")) %>% 
  filter(From == Entry) %>% 
  select(-From) %>% 
  add_row(Entry = "P62988", Entry.name = "UBIQ_HUMAN", Gene.name = "UBIQ") %>% 
  mutate(Gene.name = sub(";.*", "", Gene.name))

ups_headers2 <- ups_headers %>% 
  left_join(actual_sequences2, by = c("Accession" = "Entry")) %>% 
  mutate(Entry.name.y = sub("HUMAN", "UPS", Entry.name.y))

# Reconstruct fasta headers to my liking
output_headers <- paste0(
  "sp|UPS_",
  ups_headers2$Accession,
  "|",
  ups_headers2$Entry.name.y,
  " ",
  ups_headers2$Description,
  " OS=Synthetic OX=0000 GN=",
  ups_headers2$Gene.name,
  " PE=1 SV=1"
)

output_headers

output_ups <- ups
names(output_ups) <- output_headers

# Order should be preserved
ups
output_ups

writeXStringSet(
  output_ups,
  filepath = here("2023_12_csd_ups_proteins.fasta.gz"),
  compress = TRUE
)
