### REQUIRED PACKAGES (will need to install each one time on your machine) 
#install.packages("rentrez") # R package to call Entrez search tool (NCBI, including GenBank)
require("rentrez")           # Tutorial: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

#install.packages("bold")    # R package to pull sequences from BOLD
require("bold")              # Tutorial: https://github.com/ropensci/bold

### MODIFIED to make recursive over a list of taxa
### Prints out one fasta per taxon in the list 

### Input and output files - PUT YOUR FILE NAMES HERE
taxon_file  <- "example_taxon_list.txt"                           # must be in working directory or provide full path
out_table   <- "name_of_file_for_merged_genbank_bold_entries.csv" # these print to your working directory, or provide full path

### Give ending name of fasta files. E.g., "my_fun_project.fasta" will become "Oncorhynchus_mykiss_my_fun_project.fasta"
out_genbank_bold_end <- "end_name_file_fasta_merged.fasta"     

### NO NEED TO EDIT BELOW THIS LINE - RUN ENTIRE SCRIPT #####################################################################################
### Dataframe to append for out summary table
all_taxa_out <- data.frame(
  genbank_accession	= vector(, length = 0),
  organism	= vector(, length = 0),
  uid	= vector(, length = 0),
  loc	= vector(, length = 0),
  species_name	= vector(, length = 0),
  sampleid	= vector(, length = 0),
  lat	= vector(, length = 0),
  lon	= vector(, length = 0))

### taxa to loop through
taxa_to_loop <- read.table(taxon_file, sep = '\n')$V1                           # this alone adequate for BOLD query

### loop through taxa list
for(i in 1:length(taxa_to_loop)){
  taxon_list <- taxa_to_loop[i]                                                 # taxon list becomes entry 'i' 
  print(paste("Starting", taxon_list))                                          # print progress
  
  fasta_name <- paste(gsub('\\s', '_', taxon_list), '_', out_genbank_bold_end, sep = '') # needs an identifiable name for outputs

  ### GenkBank search query (simplified with length = 1)
  genbank_term <- paste("mitochondrion[FILT] AND ", 
                        taxon_list, 
                        "[ORGN]",
                        sep = "", collapse = "")                                # this is the search term for GenBank
  
  ### Get the GenBank sequence entries
  r_search <- entrez_search(db = "nucleotide",                                  # database; run 'entrez_db_summary("nucleotide")' for details
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
  for(j in 1:length(r_fetched)){
    orgn_genbank[j] <- r_fetched[[j]]$Organism                                  # taxon name
    accession_genbank[j] <- r_fetched[[j]]$Caption                              # accession numbers
    ifelse(is.na(r_fetched[[j]]$SubName) == T,
           location_genbank[j] <- NA,
           location_genbank[j] <- r_fetched[[j]]$SubName)                       # sometimes contains lat/lon
  }
  
  records_genbank <- data.frame(cbind(organism = orgn_genbank,
                                      uid = uid_genbank,
                                      genbank_accession = accession_genbank, 
                                      loc = location_genbank))
  
  ### Get the BOLD sequence entries
  records_bold <- bold_seqspec(taxon = taxon_list)[, c('species_name',   
                                                       'processid',             # BOLD identifier
                                                       'genbank_accession', 
                                                       'lat', 
                                                       'lon')]
  
  ### Merge GenBank and BOLD entries (and create output table)
  joint_records <- merge(records_genbank, 
                         records_bold, 
                         by = "genbank_accession", all = T) 
  in_genbank <- subset(joint_records, is.na(joint_records$uid) == F)            # there is a genbank record (Entrez indexing faster than BOLD)
  bold_only <- subset(joint_records, is.na(joint_records$uid) == T)             # only avaiable at bold

  ### GenBank FASTA - for large sets of records (>300, have to do batches)
  ids <- in_genbank$uid 
  chunk_size <- 300
  ids_chunked <- split(ids, ceiling(seq_along(ids)/chunk_size))
  
  genbank_fasta <- c()
  for(j in 1:length(ids_chunked)){
    genbank_fasta <- append(genbank_fasta,
                            entrez_fetch(db = "nucleotide",
                                         id = ids_chunked[[j]],
                                         rettype = "fasta"))
    print(paste("Chunk #",j,"of",length(ids_chunked),"done!"))                  # progress message
  }
  
  ### if BOLD-only entries, join, otherwise, use GenBank only
  if(nrow(bold_only) > 0) {
    bold_fasta <- bold_seq(ids = bold_only$processid)                             
    fasta_vector <- vector(, length = length(bold_fasta))                       # vector to put fasta sequences in
    for(j in 1:length(fasta_vector)){
      fasta_vector[j] <- paste('>', bold_fasta[[j]]$id, '|', bold_fasta[[j]]$name, "\r", 
                               bold_fasta[[j]]$sequence, sep = '')              # reconstruct fasta with formatted headers
    }
    
    #### Joint records
    out_joint_records <- subset(joint_records, joint_records$processid %in% sapply(bold_fasta, "[[", "id") |
                                  joint_records$genbank_accession %in% in_genbank$genbank_accession)
    all_taxa_out <- rbind(all_taxa_out, out_joint_records)
    
    ### Joint fasta
    genbank_bold_fasta <- append(genbank_fasta, paste(fasta_vector, sep="", collapse=""))
    write(genbank_bold_fasta, file = fasta_name)
    
  } else{
    print("All BOLD entries also in GenBank")
    all_taxa_out <- rbind(all_taxa_out, in_genbank)
    write(genbank_fasta, fasta_name)
  }
  
  print(paste("Finished", taxon_list))                                          # progress print
}

write.csv(all_taxa_out, out_table) # ALTOGETHER for summary table





