#### Load packages ####
library(here)
library(Biostrings)

#### Load standard contaminants ####
cont_fasta_path <- list.files(pattern = "contaminants.fasta")

cont_fasta_release <- unlist(regmatches(cont_fasta_path, regexec("[0-9]{4}_[0-9]{2}", cont_fasta_path)))

cont_fasta <- Biostrings::readAAStringSet(here(cont_fasta_path))

#### Append RNases used in OOPS ####
rnase_fasta_path <- list.files(pattern = "rnase_a_t1.fasta")

cont_rnases <- Biostrings::readAAStringSet(here(rnase_fasta_path))

cont_comb <- c(cont_fasta, cont_rnases)

#### Save fasta ####
out_fasta_path <- sub("contaminants", "contaminants_w_rnase_a_t1", cont_fasta_path)

Biostrings::writeXStringSet(cont_comb, out_fasta_path, append = FALSE)

#### Save metadata ####
out_metadata <- data.frame(
  "File" = out_fasta_path,
  "Number of Sequences" = length(cont_comb),
  "UniProt Release" = cont_fasta_release,
  "Date" = Sys.time()
)

write.table(out_metadata, 
            here("metadata", "contaminants_w_rnase_a_t1-metadata.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
