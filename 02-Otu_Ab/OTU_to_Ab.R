library(dplyr)
require(xlsx)
library(biomformat)
library(optparse)
library(crayon)

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

cat(blue("Run date to analyze is: "%+%green$bold(opt$date)%+%"\n"))
cat(blue("Your confidence value is: "%+%green$bold(opt$confVal)%+%"\n"))
cat(blue("The path where input file and output file an folder is: "%+%green$bold(opt$Path)%+%"\n"))

#########################################################################################

for (chain in c("16S","ITS")) {
	# Loading OTUs File
	ifelse(file.exists(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_mapped_reads_tax.biom")),
	dat <- read_biom(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_mapped_reads_tax.biom")),
	next)
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
	    arrange_(.dots = paste0("desc(",i,")")) 
	    #arrange(Species)
	  Samp[2] <- Samp[2] * 100/colSums(Samp[2])
	  colnames(Samp) <- gsub('\\.','-',colnames(Samp))
	  #print(colnames(Samp[2]))
	  #write.table(Samp, file = paste0(opt$Path,opt$date,"/",chain,"_Individuales/",colnames(Samp[2]),"_",chain,".csv"), 
	  #            col.names = TRUE, row.names = FALSE, sep = ",")
	  write.xlsx(x = as.data.frame(Samp), file = paste0(opt$Path,opt$date,"/",chain,"_Individuales/",colnames(Samp[2]),"_",chain,".xlsx"), 
	             row.names = FALSE)
	}
	OTUSamplesFile <- CleanBiomFile %>% 
		select( -confidence)

	colnames(OTUSamplesFile) <- gsub("\\-",".",colnames(OTUSamplesFile))
        cat(yellow("creating abundances of "%+%green$bold(chain)%+%" samples of "%+%green$bold(opt$date)%+%" Run\n"))
	#OTUs Abundance by sample
	dir.create(paste0(opt$Path,opt$date,"/",chain,"_OTUs_Ab/"), showWarnings = FALSE, recursive = TRUE)

	for(OTUSamp in colnames(OTUSamplesFile[2:(length(colnames(OTUSamplesFile))-7)])){
	  #print(OTUSamp)
          Samp <- OTUSamplesFile %>%
            select_(.dots = "OTU.ID","taxonomy1","taxonomy2","taxonomy3","taxonomy4","taxonomy5","taxonomy6","taxonomy7",OTUSamp) %>%
          filter_(.dots = paste0(OTUSamp, " > 0")) %>% 
            arrange(taxonomy7)
          Samp[9] <- Samp[9] * 100/colSums(Samp[9])
          write.xlsx(x = as.data.frame(Samp), file = paste0(opt$Path,opt$date,"/",chain,"_OTUs_Ab/",OTUSamp,"_",chain,".xlsx"), row.names = FALSE)
	}
}

cat(red$bold("O       o O       o O       o          O       o O       o O       o\n"))
cat(red$bold("| O   o | | O   o | | O   o |   STEP   | O   o | | O   o | | O   o |\n"))
cat(red$bold("| | O | | | | O | | | | O | |    2     | | O | | | | O | | | | O | |\n"))
cat(red$bold("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |\n"))
cat(red$bold("o       O o       O o       O          o       O o       O o       O\n"))
