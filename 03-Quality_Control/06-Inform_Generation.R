library(optparse)
library(data.table)
library(dplyr)
require(xlsx)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/")
)
### Load arguments
opt <- parse_args(OptionParser(option_list=option_list))

##### ANALYSIS
for (chain in c("16S","ITS")) {
 FERM_size <- file.size(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/Informs/Good_Fer_",opt$date,"_",chain,".csv"))
 ifelse(FERM_size < 5, 
        FERMENTED_GOOD <- data.frame(),
        FERMENTED_GOOD <- read.csv(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/Informs/Good_Fer_",opt$date,"_",chain,".csv")))
 
 LOW_size <- file.size(paste0(opt$Path,opt$date,"/Bad_Reads/",chain,"/Informs/Bad_Reads_",opt$date,"_",chain,".csv"))
 ifelse(LOW_size < 5,
        LOW_READS <- data.frame(),
        LOW_READS <- read.csv(paste0(opt$Path,opt$date,"/Bad_Reads/",chain,"/Informs/Bad_Reads_",opt$date,"_",chain,".csv")))
 
 BMK_size <- file.size(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/Informs/Good_BMK_",opt$date,"_",chain,".csv"))
 ifelse(BMK_size < 5, 
        BMK_GOOD <- data.frame(),
        BMK_GOOD <- read.csv(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/Informs/Good_BMK_",opt$date,"_",chain,".csv")))
 
 SG_size <- file.size(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/Good_SoilGrape_",opt$date,"_",chain,".csv"))
 ifelse(SG_size < 5,
        SOIL_GRAPE_GOOD <- data.frame(),
        SOIL_GRAPE_GOOD <- read.csv(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/Good_SoilGrape_",opt$date,"_",chain,".csv")))
 
 ANITA_INFO <- rbindlist(list(SOIL_GRAPE_GOOD, BMK_GOOD, FERMENTED_GOOD, LOW_READS), fill = TRUE) 
 ANITA_INFO[ANITA_INFO == "NoA"] <- NA
 dir.create(paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/"), showWarnings = FALSE, recursive = TRUE)
 write.table(ANITA_INFO, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 ANIITA_INFO <- read.csv(paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"))
 ANIITA_INFO <- ANIITA_INFO %>%
   arrange(Sample) %>%
   select(Sample, Chain_Type, Project, Sample_type,  Initial_reads, Final_reads, Contamination_sp, Sp_Num, Species, Sp_Max_perc,  Control_Num_Sp, Control_Max_Sp, Control_Samp_Type_SP, Diagnostic)
 
 write.table(ANIITA_INFO, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 write.xlsx(x = ANIITA_INFO, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".xlsx"), row.names = FALSE)
 print(paste0("Finished analysis of ",chain))
}

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    6     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")