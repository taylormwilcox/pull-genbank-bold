README for genbank_bold_pull.R

- You need to install the packages rentrez and bold one time on your machine and load them each session (code in script)
- You only need to modify the first chunk of code to set your working directory, input file name, and output file names
- Input file is a simple text file where each line is a taxon name 
	- Should be fine to copy-paste from Excel column
	- If you get error 'incomplete final line...' simply go to last entry, then hit 'Enter' to create an empty line at the end of the file
- The output table (csv) is readable in Excel and has the following columns (ignore first column without header)
	- genbank_accession	
	- organism = species name from GenBank
	- uid = entry identifier for NCBI (can ignore)	
	- loc = GenBank field that sometimes contains location
	- species_name	= species name from BOLD
	- processid = BOLD identifier
    - lat = BOLD latitude	
	- lon = BOLD longitude
- There are two FASTA outputs: One for GenBank and one for BOLD 
	- GenBank FASTA have original headers (typically "Accesion_number taxon_name additional_info"
	- BOLD FASTA have header "bold_id|taxon_name"
