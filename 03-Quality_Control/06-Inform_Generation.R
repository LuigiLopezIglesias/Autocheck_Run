library(optparse)
library(data.table)
library(dplyr)
require(xlsx)
library(crayon)
library(RPostgreSQL)
library(stringr)
#########################################################################################
###     Cargar argumentos//Load arguments      ###
#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/"),
  make_option(c('-U', '--user'), action='store',
              dest='user', type='character'),
  make_option(c('-P', '--password'), action='store',
              dest='password', type='character'),
  make_option(c('-H', '--host'), action='store',
              dest='host', type='character'),
  make_option(c('-g','--port'), action='store',
              dest='port', type='integer'),
  make_option(c('-n','--dbname'), action='store',
              dest='dbname', type='character')
)
### Load arguments
opt <- parse_args(OptionParser(option_list=option_list))

cat(blue("Run date to analyze is: "%+%green$bold(opt$date)%+%"\n"))
cat(blue("The path where input file and output file an folder is: "%+%green$bold(opt$Path)%+%"\n"))

##### ANALYSIS
for (chain in c("16S","ITS")) {
   if(file.exists(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/Informs/Good_Fer_",opt$date,"_",chain,".csv"))==TRUE){
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
     if (FERM_size < 5 && LOW_size < 5 && BMK_size < 5 && SG_size < 5) {
       cat(red$bold("Don't have information to make inform file\n"))
     } else {
       ANITA_INFO <- rbindlist(list(SOIL_GRAPE_GOOD, BMK_GOOD, FERMENTED_GOOD, LOW_READS), fill = TRUE)
         ANITA_INFO[ANITA_INFO == "NoA"] <- NA
         dir.create(paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/"), showWarnings = FALSE, recursive = TRUE)
         write.table(ANITA_INFO, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")

         ANIITA_INFO <- read.csv(paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"))	
       if (FERM_size < 5) {
         ANIITA_INFO <- ANIITA_INFO %>%
           arrange(Sample) %>%
           select(Sample, Chain_Type, Project, Sample_type, Initial_reads, Final_reads, Contamination_sp, Sp_Num, Species, Sp_Max_perc,  Control_Num_Sp, Control_Max_Sp, Control_Samp_Type_SP, Diagnostic)
       } else { 
         ANIITA_INFO <- ANIITA_INFO %>%
           arrange(Sample) %>%
           select(Sample, Chain_Type, Project, Sample_type, Stage, Substage, Initial_reads, Final_reads, Contamination_sp, Sp_Num, Species, Sp_Max_perc,  Control_Num_Sp, Control_Max_Sp, Control_Samp_Type_SP, Diagnostic)
       }
     #########################################################################################
     ### Conexion con la base de datos
     drv <- dbDriver("PostgreSQL")
     local_DB <- dbConnect(drv, user=opt$user,
                          password=opt$password,
                          host=opt$host,
                          port=opt$port, dbname=opt$dbname)

     Repetitions <- paste0("select m.c_muestra_wineseq, repeat_dna_extraction, repeat_pcr_16s, repeat_pcr_its
                            from muestra m
                            join muestra_internal mi on m.id = mi.id_muestra")
     Rep_set <- dbGetQuery(local_DB, Repetitions)
     colnames(Rep_set) <- c("Sample", "DNA_rep", "16S_PCR_Rep", "ITS_PCR_Rep")
     if (chain == "16S") {
       Rep_set$Sample <- paste0(Rep_set$Sample,"b")
     }
     Report <- merge(ANIITA_INFO, Rep_set)
    
     write.table(Report, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
     write.xlsx(x = Report, file = paste0(opt$Path,opt$date,"/RUN_Inform/",chain,"/Inform_",opt$date,"_",chain,".xlsx"), row.names = FALSE)
     cat(blue("Finished analysis of "%+%green$bold(chain)%+%"\n"))
   }
   } else {
    cat(red$bold(paste0("\n Dont have the files needed to analyze ",chain,"\n")))
   }
}

cat(magenta$bold("O       o O       o O       o          O       o O       o O       o\n"))
cat(magenta$bold("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |\n"))
cat(magenta$bold("| | O | | | | O | | | | O | |    6     | | O | | | | O | | | | O | |\n"))
cat(magenta$bold("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |\n"))
cat(magenta$bold("o       O o       O o       O          o       O o       O o       O\n"))
