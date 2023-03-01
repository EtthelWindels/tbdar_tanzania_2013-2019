## ------------------------------------------------------------------------
## Replace names of sequences with complete names for 
## BEAST analyses (GNUMBER/date) 
## 2022-02-02 Etthel Windels
## ------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(plyr)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(tibble)))
suppressMessages(suppressWarnings(require(seqinr)))
suppressMessages(suppressWarnings(require(stringr)))


# Parser ------------------------------------------------------------------

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata file")
parser$add_argument("--alignment", type="character", help="Alignment file")
parser$add_argument("--output_alignment", type = "character",
                    help = "Output file for updated alignment")
parser$add_argument("--intro_name", type = "character",
                    help = "Name of introduction corresponding to lineage")
args <- parser$parse_args()


# Read arguments ----------------------------------------------------------

ALIGNMENT <- args$alignment
METADATA <- args$metadata
OUTPUT_ALIGNMENT <- args$output_alignment
INTRO_NAME <- args$intro_name

print(paste("metadata: ", METADATA))
print(paste("alignment: ", ALIGNMENT))
print(paste("output alignment: ", OUTPUT_ALIGNMENT))
print(paste("introduction: ", INTRO_NAME))


# Read files --------------------------------------------------------------

alignment <- seqinr::read.fasta(file = ALIGNMENT, seqtype="DNA", forceDNAtolower = F)
full_metadata <- read.table(file = METADATA, header=T,sep= '\t')


# Adjust metadata ---------------------------------------------------------

metadata <- 
  full_metadata %>%
# 1. Filter metadata for only sequences in the alignment
  filter(G_NUMBER %in% names(alignment)) %>%
# 2. Adjust date format
  mutate(ids_date_interview = str_split_fixed(ids_date_interview,' ',2)[,1]) %>% 
  mutate(ids_date_interview=paste(str_split_fixed(ids_date_interview,pattern="\\.",3)[,1],str_split_fixed(ids_date_interview,pattern="\\.",3)[,2],str_split_fixed(ids_date_interview,pattern="\\.",3)[,3], sep='-')) %>%
# 3. Create new_id column
  mutate(new_id = paste(G_NUMBER, ids_date_interview, sep = "/"))


# Adjust alignment --------------------------------------------------------

# 1. Filter alignment for introduction and remove outgroup G00157
alignment <- 
  alignment[names(alignment) %in% metadata$G_NUMBER[metadata$Introduction==INTRO_NAME]]

# 2. Adjust alignment names
gnumber <- names(alignment)
full_names <- metadata$new_id[match(gnumber, metadata$G_NUMBER)]
print(full_names)
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names


# Save output files -------------------------------------------------------

seqinr::write.fasta(alignment, names = names(alignment), file.out = OUTPUT_ALIGNMENT)
