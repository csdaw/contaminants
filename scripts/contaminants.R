#### Load packages ####
library(here)
library(curl)

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
