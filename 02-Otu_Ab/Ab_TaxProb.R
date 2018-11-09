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

for (chain in c("16S")) {
        # Loading OTUs File
        if (file.exists(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_mapped_reads_tax.biom"))==TRUE){
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
	  
   	  ## Loading Tax_Problem csv
	  TaxProblems <- read.csv(paste0("/home/yamishakka/Escritorio/github/Luigi/Autocheck_Run/02-Otu_Ab/DB_Analysis/",chain,"/RealibleNames_10_10.csv"))
          
  	  TaxProblems <- TaxProblems %>%
	    select(ReliableName,specie) %>%
	    unique

          TaxProblems$Species <- substr(TaxProblems$specie, 3, length(TaxProblems$specie))

	  UnionTaxProbAnalysis <- merge(AbundanceFile, TaxProblems, by = "Species", all.x = TRUE)
	
	  print(head(TaxProblems))
	  print(head(AbundanceFile$Species))

	  print(head(UnionTaxProbAnalysis))

## Write abundance file to Anita
          write.table(UnionTaxProbAnalysis, file = paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance_tax_problems.csv"),
                      col.names = TRUE, row.names = FALSE, sep = ",")
}
}
