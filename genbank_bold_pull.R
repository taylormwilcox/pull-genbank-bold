### REQUIRED PACKAGES (will need to install each one time on your machine) 
#install.packages("rentrez") # R package to call Entrez search tool (NCBI, including GenBank)
require("rentrez")           # Tutorial: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

#install.packages("bold")    # R package to pull sequences from BOLD
require("bold")              # Tutorial: https://github.com/ropensci/bold

### Input and output files - PUT YOUR FILE NAMES HERE
setwd("~/pull-genbank-bold") # this will set the working directory to the example folder
taxon_file  <- "example_taxon_list.txt"                           # must be in working directory or provide full path
out_table   <- "name_of_file_for_merged_genbank_bold_entries.csv" # these print to your working directory, or provide full path
out_genbank <- "name_of_file_for_genbank_sequences.fasta"
out_bold    <- "name_of_file_for_bold_sequences.fasta"
out_genbank_bold <- "name_of_file_for_merged_sequences.fasta"     # requested by Tommy

### NO NEED TO EDIT BELOW THIS LINE - RUN ENTIRE SCRIPT #####################################################################################

### Format search query for GenBank
taxon_list <- read.table(taxon_file, sep = '\n')$V1                             # this alone adequate for BOLD query

orgn_term <- vector(, length = 0)                                               # empty vector to fill with organism search terms
number_OR <- length(taxon_list) - 1                                             # number of organisms to paste 'OR' after
for(i in 1:length(number_OR)){
  orgn_term <- append(orgn_term, paste(taxon_list[i], "[ORGN] OR ", sep = ""))
}

last_org <- taxon_list[length(taxon_list)]                                      # final organism in list
genbank_term <- paste("mitochondrion[FILT] AND (",                      
                      orgn_term, 
                      last_org,
                      "[ORGN]",
                      ")",
                      sep = "", collapse = "")                                  # this is the search term for GenBank

### Get the GenBank sequence entries
r_search <- entrez_search(db = "nucleotide",                                    # database; run 'entrez_db_summary("nucleotide")' for details
                          term = genbank_term,
                          use_history = T) 

r_fetched <- entrez_summary(db = "nucleotide", 
                            web_history = r_search$web_history,
                            retmode = "xml")  

# list of lists with records
uid_genbank <- names(r_fetched) # integer ids (NCBI-specific)
orgn_genbank <- c()
accession_genbank <- c()
location_genbank <- c()
for(i in 1:length(r_fetched)){
  orgn_genbank[i] <- r_fetched[[i]]$Organism                                    # taxon name
  accession_genbank[i] <- r_fetched[[i]]$Caption                                # accession numbers
  ifelse(is.na(r_fetched[[i]]$SubName) == T,
         location_genbank[i] <- NA,
         location_genbank[i] <- r_fetched[[i]]$SubName)                         # sometimes contains lat/lon
}

records_genbank <- data.frame(cbind(organism = orgn_genbank,
                                    uid = uid_genbank,
                                    genbank_accession = accession_genbank, 
                                    loc = location_genbank))

### Get the BOLD sequence entries
records_bold <- bold_seqspec(taxon = taxon_list)[, c('species_name',   
                                                     'processid',               # BOLD identifier
                                                     'genbank_accession', 
                                                     'lat', 
                                                     'lon')]

### Merge GenBank and BOLD entries (and create output table)
joint_records <- merge(records_genbank, 
                       records_bold, 
                       by = "genbank_accession", all = T) 
in_genbank <- subset(joint_records, is.na(joint_records$uid) == F)              # there is a genbank record (Entrez indexing faster than BOLD)
bold_only <- subset(joint_records, is.na(joint_records$uid) == T)               # only avaiable at bold

### GenBank FASTA - for large sets of records (>300, have to do batches)
ids <- in_genbank$uid 
chunk_size <- 300
ids_chunked <- split(ids, ceiling(seq_along(ids)/chunk_size))

genbank_fasta <- c()
for(i in 1:length(ids_chunked)){
  genbank_fasta <- append(genbank_fasta,
                          entrez_fetch(db = "nucleotide",
                                       id = ids_chunked[[i]],
                                       rettype = "fasta"))
  print(paste("Chunk #",i,"of",length(ids_chunked),"done!"))                    # progress message
}

write(genbank_fasta, out_genbank) # EXPORT (fasta-format)

### BOLD FASTA - NOTE: Empty BOLD entries exist
if(nrow(bold_only) > 0) {
  bold_fasta <- bold_seq(ids = bold_only$processid)                             
  fasta_vector <- vector(, length = length(bold_fasta))                           # vector to put fasta sequences in
  for(i in 1:length(fasta_vector)){
    fasta_vector[i] <- paste('>', bold_fasta[[i]]$id, '|', bold_fasta[[i]]$name, "\r", 
                             bold_fasta[[i]]$sequence, sep = '')                  # reconstruct fasta with formatted headers
  }
    write(fasta_vector, out_bold) # EXPORT (fasta-format)
    
    ### Joint records
    out_joint_records <- subset(joint_records, joint_records$processid %in% sapply(bold_fasta, "[[", "id") |
                                  joint_records$genbank_accession %in% in_genbank$genbank_accession)
    write.csv(out_joint_records, out_table) # EXPORT (readable in Excel
    
    ### Joint fasta
    genbank_bold_fasta <- append(genbank_fasta, paste(fasta_vector, sep="", collapse=""))
    write(genbank_bold_fasta, file = out_genbank_bold)
  
} else{
  print("All BOLD entries also in GenBank")
  write.csv(joint_records, out_table)
}


