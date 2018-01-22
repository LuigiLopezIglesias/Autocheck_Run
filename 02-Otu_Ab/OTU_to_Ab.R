library(dplyr)
require(xlsx)
library(biomformat)
library(optparse)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-v', '--confidenceValue'), action='store',
              dest='confVal', type='numeric', default = 0.7),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(paste0("Run date to analyze is: ",opt$date))
print(paste0("Your confidence value is: ",opt$confVal))
print(paste0("The path where input file and output file an folder is: ",opt$Path))
#########################################################################################

for (chain in c("16S","ITS")) {
	# Loading OTUs File
	dat <- read_biom(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_mapped_reads_tax.biom"))
	otu_table <- as.data.frame(as.matrix(biom_data(dat)))
	otu_table$OTU.ID <- rownames(otu_table)
	
	MetaDat <- observation_metadata(dat)
	MetaDat$OTU.ID <- rownames(MetaDat)
	
	BiomFile <- merge(otu_table, MetaDat, by = "OTU.ID") 
	
	CleanBiomFile <- BiomFile %>%
	  filter(confidence > opt$confVal)
	
	SamplesFile <- CleanBiomFile %>%
	  select(-OTU.ID, -confidence, -taxonomy1, -taxonomy2, -taxonomy3,- taxonomy4, -taxonomy5, -taxonomy6)
	
	## Selection Species taxa
	SamplesFile$Species <-  gsub('.*s:','',SamplesFile$taxonomy7)
	
	## summarise reads of same Specie but different otus
	AbundanceFile <- SamplesFile %>%
	  select(-taxonomy7) %>%
	  group_by(Species) %>%
	  summarise_all(funs(sum))
	
	## Write abundance file to Anita
	write.table(AbundanceFile, file = paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance.csv"), 
	            col.names = TRUE, row.names = FALSE, sep = ",")
	
	## Individual Abundance by sample
	colnames(AbundanceFile) <- gsub('-','.',colnames(AbundanceFile))
	
	dir.create(paste0(opt$Path,opt$date,"/",chain,"_Individuales/"), showWarnings = FALSE, recursive = TRUE)
	
	for(i in colnames(AbundanceFile[-1])){
	  Samp <- AbundanceFile %>%
	    select_(.dots = "Species",i) %>%
	    filter_(.dots = paste0(i, " > 0")) %>% 
	    arrange(Species)
	  Samp[2] <- Samp[2] * 100/colSums(Samp[2])
	  colnames(Samp) <- gsub('\\.','-',colnames(Samp))
	  #print(colnames(Samp[2]))
	  #write.table(Samp, file = paste0(opt$Path,opt$date,"/",chain,"_Individuales/",colnames(Samp[2]),"_",chain,".csv"), 
	  #            col.names = TRUE, row.names = FALSE, sep = ",")
	  write.xlsx(x = as.data.frame(Samp), file = paste0(opt$Path,opt$date,"/",chain,"_Individuales/",colnames(Samp[2]),"_",chain,".xlsx"), 
	             row.names = FALSE)
	}
}

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |   STEP   | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    2     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")