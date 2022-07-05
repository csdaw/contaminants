#### Load packages ####
library(here)
library(curl)
library(uniprotREST)
library(Biostrings)
library(magrittr)
library(httr2)
library(jsonlite)
library(dplyr)

#### Download source files ####
# Download Hao group universal contaminants excel file
curl::curl_download(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/blob/main/Universal%20protein%20contaminant%20FASTA/Contaminant%20protein%20information%20in%20the%20FASTA%20library%20and%20the%20potential%20source%20of%20contaminations.xlsx?raw=true",
  destfile = here::here(
    "hao_lab",
    "contaminants_universal.xlsx"
  )
)

# Download Hao group universal contaminants FASTA
curl::curl_download(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/main/Universal%20protein%20contaminant%20FASTA/Universal%20Protein%20Contaminants.fasta",
  destfile = here::here(
    "hao_lab",
    "contaminants_universal.fasta"
  )
)

# Download Hao Group cell culture contaminants FASTA
curl::curl_download(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/main/Sample-type%20specific%20contaminant%20FASTA/Cell%20Culture%20Contaminants.fasta",
  destfile = here::here(
    "hao_lab",
    "contaminants_cell_culture.fasta"
  )
)

# Download Hao Group mouse tissue contaminants FASTA
curl::curl_download(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/main/Sample-type%20specific%20contaminant%20FASTA/Mouse%20Tissue%20Contaminants.fasta",
  destfile = here::here(
    "hao_lab",
    "contaminants_mouse_tissue.fasta"
  )
)

# Download Hao Group rat tissue contaminants FASTA
curl::curl_download(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/main/Sample-type%20specific%20contaminant%20FASTA/Rat%20Tissue%20Contaminants.fasta",
  destfile = here::here(
    "hao_lab",
    "contaminants_rat_tissue.fasta"
  )
)

#### Check Universal FASTA sequences ####
cont_uni <- Biostrings::readAAStringSet(here("hao_lab", "contaminants_universal.fasta"))
cont_cell <- Biostrings::readAAStringSet(here("hao_lab", "contaminants_cell_culture.fasta"))
cont_mouse <- Biostrings::readAAStringSet(here("hao_lab", "contaminants_mouse_tissue.fasta"))
cont_rat <- Biostrings::readAAStringSet(here("hao_lab", "contaminants_rat_tissue.fasta"))

#### Extract UniProt accessions from FASTA headers ####
headers <- lapply(list(cont_uni, cont_cell, cont_mouse, cont_rat), names)
accessions_list <- lapply(headers, function(txt) {
  regmatches(txt, regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", txt))
})
accessions <- unique(unlist(accessions_list))

### Get data from UniProt ####
# Send UniProt accessions for mapping
rest_url <- "https://rest.uniprot.org/idmapping"

post_req <- httr2::request(rest_url) %>%
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url_path_append("run") %>%
  httr2::req_body_form(
    `from` = "UniProtKB_AC-ID",
    `to` = "UniProtKB",
    `ids` = paste(accessions, collapse = ",")
  ) %>%
  httr2::req_retry(max_tries = 5)

post_resp <- httr2::req_perform(post_req)
jobid <- httr2::resp_body_json(post_resp)
message(paste("Job ID:", jobid))
Sys.sleep(3)

# Get URL to get results
status_req <- httr2::request(rest_url) %>%
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url_path_append("status", jobid) %>%
  httr2::req_retry(is_transient = ~ httr2::resp_status(.x) %in% c("200", "429", "503")) %>%
  httr2::req_method("HEAD") %>% 
  httr2::req_options(followlocation = FALSE)

status_resp <- httr2::req_perform(status_req)

result_url <- status_resp %>% 
  httr2::resp_header("location") %>% 
  gsub("results/", "results/stream/", .)

cols_to_get <- c("accession", "id", "reviewed", "protein_name", "gene_names", 
                 "organism_name", "length", "sequence")

# Get results via stream (should only do this for small queries)
result_req <- httr2::request(rest_url) %>% 
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url(result_url) %>% 
  httr2::req_url_query(
    `format` = "tsv",
    `fields` = paste(cols_to_get, collapse = ",")
  )

result_resp <- httr2::req_perform(result_req)

result <- result_resp %>% 
  httr2::resp_body_string() %>% 
  read.delim(text = ., check.names = FALSE) %>% 
  select(-From)

#### Combine with existing data ####
cont_hao <- readxl::read_xlsx(here("hao_lab/contaminants_universal.xlsx"), skip = 1)

out <- left_join(
  result, 
  select(cont_hao, c("Uniprot ID", "Source of Contamination")), 
  by = c("Entry" = "Uniprot ID")
)

missing_rows <- cont_hao %>% 
  filter(`Uniprot ID` %in% c("AAAA1", "AAAA2")) %>% 
  mutate(Sequence = c("DYKDDDDK", "YPYDVPDYA"), .before = "Source of Contamination") %>% 
  mutate(Length = nchar(Sequence)) %>% 
  `colnames<-`(
    c(
      "Entry", 
      "Entry Name",
      "Reviewed",
      "Protein names",
      "Gene Names",
      "Organism",
      "Length",
      "Sequence",
      "Source of Contamination"
    )
  )

# Manually add data for FLAG and HA tags
# and also standardise column names
out2 <- bind_rows(out, missing_rows) %>% 
  janitor::clean_names()

#### Add FASTA column ####
# Indicate which FASTA should contain which proteins of the total set
out3 <- out2 %>% 
  mutate(universal = entry %in% c(accessions_list[[1]], "AAAA1", "AAAA2"),
         cell_culture = entry %in% accessions_list[[2]],
         mouse_tissue = entry %in% c(accessions_list[[3]], "AAAA1", "AAAA2"),
         rat_tissue = entry %in% c(accessions_list[[4]], "AAAA1", "AAAA2"))

#### Deal with duplicate sequences
# 3 duplicated sequences but only 2 are potentially problematic
# Also for now I'll only be using the 'universal' proteins
out4 <- out3 %>% 
  filter(!entry %in% c("Q2FZL2", "Q9BYR9")) %>% 
  filter(universal) %>% 
  select(-universal, -cell_culture, -mouse_tissue, -rat_tissue)

#### Save as CSV ####
readr::write_csv(out4, file = "contaminants.csv")

#### Generate and save FASTA ####
# Send UniProt accessions for mapping
post_req <- httr2::request(rest_url) %>%
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url_path_append("run") %>%
  httr2::req_body_form(
    `from` = "UniProtKB_AC-ID",
    `to` = "UniProtKB",
    `ids` = paste(out4$entry, collapse = ",")
  ) %>%
  httr2::req_retry(max_tries = 5)

post_resp <- httr2::req_perform(post_req)
jobid <- httr2::resp_body_json(post_resp)
message(paste("Job ID:", jobid))
Sys.sleep(3)

# Get URL to get results
status_req <- httr2::request(rest_url) %>%
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url_path_append("status", jobid) %>%
  httr2::req_retry(is_transient = ~ httr2::resp_status(.x) %in% c("200", "429", "503")) %>%
  httr2::req_method("HEAD") %>% 
  httr2::req_options(followlocation = FALSE)

status_resp <- httr2::req_perform(status_req)

result_url <- status_resp %>% 
  httr2::resp_header("location") %>% 
  gsub("results/", "results/stream/", .)

# Get results via stream (should only do this for small queries)
result_req <- httr2::request(rest_url) %>% 
  httr2::req_user_agent("Charlotte Dawson csd51@cam.ac.uk") %>%
  httr2::req_url(result_url) %>% 
  httr2::req_url_query(
    `format` = "fasta",
    `compressed` = FALSE
  )

fasta_path <- paste(uniprotREST::current_release(), "contaminants.fasta", sep = "_")
result_resp <- httr2::req_perform(
  result_req, 
  path = fasta_path
)

# Add 'Cont_' prefix to sequence headers
Biostrings::readAAStringSet(fasta_path) %>% 
  `names<-`(gsub("(^sp\\||^tr\\|)", "\\1\\Cont_", names(.))) %>% 
  Biostrings::writeXStringSet(fasta_path)

#### Manually add extra sequences ####
manual_seqs <- Biostrings::readAAStringSet(here("hao_lab", "contaminants_universal.fasta")) %>% 
  .[grep("AAAA[1-2]", names(.))]

Biostrings::writeXStringSet(
  manual_sequences,
  fasta_path,
  append = TRUE
)

#### Save metadata ####
out_metadata <- data.frame(
  "File" = fasta_path,
  "Number of Sequences" = nrow(out4),
  "UniProt Release" = uniprotREST::current_release(),
  "Date" = Sys.time()
)

write.table(out_metadata, 
            here("metadata", "contaminants-metadata.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
