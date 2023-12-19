## Charlotte Dawson 2023
## This script I used to create my own universal contaminants FASTA file with
## minimal sequence redundancy

## Pre-processing:
# I passed Hao Group Universal Contaminants fasta through cd-hit v4.8.1 with the
# following command to remove duplicate sequences
# ~/Downloads/cd-hit-v4.8.1-2019-0228/cd-hit -i hao_lab/contaminants_universal.fasta -o hao_lab/contaminants_universal_nr.fasta -c 1.00 -n 5 -d 0 -M 12000 -T 4

# Subsequently I passed this through cd-hit v4.8.1 again with the following 
# command to remove all sequences in hao_lab/contaminants_universal_nr.fasta
# which are identical to sequences in the UniProt human
# reference proteome (release 2023_05). Note, this will get rid of human keratins:
# ~/Downloads/cd-hit-v4.8.1-2019-0228/cd-hit-2d -i uniprot/2023_05_UP000005640_9606.fasta.gz -i2 hao_lab/contaminants_universal_nr.fasta -o hao_lab/contaminants_universal_nr_nh.fasta -c 1.00 -n 5 -d 0 -M 12000 -T 4

# Then I used Phil Wilmarth's fasta_utilities/FASTA_digest_unique.py script with
# default settings on hao_lab/contaminants_universal_nr.fasta
# python scripts/FASTA_digest_unique.py

## Script:
library(here)
library(dplyr)
library(tidyr)
library(Biostrings)
library(uniprotREST)

digest_log_filepath <- here("hao_lab/FASTA_digest_unique_log.txt")
contaminants_fasta_filepath <- here("hao_lab/contaminants_universal_nr.fasta")
contaminants_excel_filepath <- here("hao_lab/contaminants_universal.xlsx")
contaminants_nh_fasta_filepath <- here("hao_lab/contaminants_universal_nr_nh.fasta")

# Load contaminants excel sheet
xl <- readxl::read_xlsx(contaminants_excel_filepath, skip = 1) %>% 
  janitor::clean_names()

# Get results from cd-hit
nr_fasta <- readAAStringSet(contaminants_fasta_filepath)
length(nr_fasta)

nh_fasta <- readAAStringSet(contaminants_nh_fasta_filepath)
length(nh_fasta)
# Indicate proteins in nr_fasta
xl <- xl %>% 
  mutate(in_nr_fasta = uniprot_id %in% sub("(.*)\\|Cont_(.*)\\|(.*)", "\\2", names(nr_fasta)),
         in_nh_fasta = uniprot_id %in% sub("(.*)\\|Cont_(.*)\\|(.*)", "\\2", names(nh_fasta)),
         is_human = grepl("HUMAN", entry_name))

# Get results of FASTA_digest_unique.py
# (Tryptic peptides include those with up to 2 missed cleavages)
tryptic <- read.delim(
  digest_log_filepath,
  skip = 2
) %>% 
  group_by(Accession) %>% 
  mutate(Total_Peptides = n())

# Count number of unique and non-unique peptides per protein
tryptic_n_unique <- tryptic %>% 
  group_by(Accession, Unique) %>% 
  mutate(N_Peptides = n()) %>% 
  select(Accession, Total_Peptides, Unique, N_Peptides) %>% 
  distinct() %>% 
  pivot_wider(
    id_cols = c(Accession, Total_Peptides), 
    names_from = Unique,
    names_prefix = "N_Unique_",
    values_from = "N_Peptides"
  )

# Grab accessions of 'tryptically non-unique' proteins
trypsin_not_unique <- tryptic_n_unique %>% 
  filter(is.na(N_Unique_TRUE))

# Indicate 'tryptically unique' proteins
xl <- xl %>% 
  mutate(tryptic_unique = !uniprot_id %in% sub("(.*)\\|Cont_(.*)\\|(.*)", "\\2", trypsin_not_unique$Accession))

# Indicate 'Other' proteins
xl <- xl %>% 
  mutate(other = source_of_contamination == "Others")

# Grab UniProt accessions for contaminants of interest
poi <- xl %>% 
  filter(
    in_nr_fasta, # Non-redundant proteins, within Hao Lab contaminants
    in_nh_fasta | is_human, # Non-redundant proteins, compared with the human proteome (except human keratins)
    tryptic_unique, # Proteins with non-redundant tryptic peptides
    !other, # Protein not in Hao Lab contaminants "Other" category which seems to be a holdover from GPM cRAP?
    status != "manually added" # Will add back FLAG and HA tag later
  )

aoi <- poi %>% 
  pull(uniprot_id)

# Get FASTA sequences for these proteins
aoi_sequences <- uniprotREST::uniprot_map(
  ids = c(
    aoi,
    "P00651" # OOPS RNase T1
  ),
  format = "fasta"
)
length(aoi) + 1 == length(aoi_sequences)

# Get UniProt metadata for these proteins
aoi_metadata <- uniprotREST::uniprot_map(
  ids = c(
    aoi,
    "P00651" # OOPS RNase T1
  ),
  format = "tsv",
  fields = c("id", "accession", "gene_primary", "protein_name", "organism_name", "length", "uniparc_id", "reviewed", "sequence_version")
) %>% 
  select(-From)

aoi_metadata2 <- aoi_metadata %>% 
  `colnames<-`(c("entry_name", "uniprot_id", "gene_name", "protein_name", "organism", "length", "uniparc_id", "uniprot_status", "sequence_version")) %>% 
  left_join(
    select(xl, entry_name, source_of_contamination),
    by = "entry_name"
  ) %>% 
  select(-entry_name) %>% 
  mutate(
    source_of_contamination = case_when(
      uniprot_id == "P00761" ~ "Enzyme",
      uniprot_id == "P00651" ~ "Enzyme",
      TRUE ~ source_of_contamination
    )
  )

# Manually add FLAG, HA tags
aoi_metadata3 <- aoi_metadata2 %>% 
  bind_rows(
    data.frame(
      uniprot_id = c("AAAAA1", "AAAAA2"),
      gene_name = c("FLAG", "HA"),
      protein_name = c("FLAG Tag", "HA Tag"),
      organism = c("Synthetic", "Synthetic"),
      length = c(8, 9),
      uniparc_id = c("", ""),
      uniprot_status = c("", ""),
      sequence_version = c(1, 1),
      source_of_contamination = c("FLAG beads", "HA beads")
    )
  )

tags <- AAStringSet(
  c("sp|AAAAA1|FLAG Tag OS=Synthetic OX=0000 GN=FLAG PE=1 SV=1" = "DYKDDDDK",
    "sp|AAAAA2|HA Tag OS=Synthetic OX=0000 GN=HA PE=1 SV=1" = "YPYDVPDYA")
)

aoi_sequences2 <- c(aoi_sequences, tags)
aoi_sequences2
# Check lengths
nrow(aoi_metadata3) == length(aoi_sequences2)

# Add Cont_ prefix
names(aoi_sequences2)
names(aoi_sequences2) <- sub("^([a-z]{2}\\|)", "\\1Cont_", names(aoi_sequences2))
names(aoi_sequences2)

# Write outputs
write.table(
  x = aoi_metadata3,
  file = here("2023_12_csd_universal_contaminants.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

writeXStringSet(
  x = aoi_sequences2,
  filepath = here("2023_12_csd_universal_contaminants.fasta.gz"),
  compress = TRUE
)
