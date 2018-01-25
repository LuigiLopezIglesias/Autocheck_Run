library(dplyr)
library(optparse)
library(RPostgreSQL)

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

opt <- parse_args(OptionParser(option_list=option_list))

##### FUNCTIONS
## Column Selection 
Sample_Row <- function(Muestra){
  Samp <- Bad_Reads_number %>%
    filter(Sample %in% Muestra)
}
## Text generation 
TEXTO <- function(Muestra, Sample_Type, Cadena, Projecto, Reads_F, Resultado)
{Texto <- data.frame(Sample=as.character(Muestra),
                     Sample_type=as.character(Sample_Type[1,2]),
                     Chain_Type=Cadena,
                     Project=Projecto,
                     Initial_reads=as.integer(Sample_Type[1,4]),
                     Final_reads=Reads_F,
                     Diagnostic = as.character(Resultado))
}

for (chain in c("16S","ITS")) {
 #########################################################################################
 ### Conexion con la base de datos
 drv <- dbDriver("PostgreSQL")
 local_DB <- dbConnect(drv, user=opt$user,
                      password=opt$password,
                      host=opt$host,
                      port=opt$port, dbname=opt$dbname)

 ##### FILE LOAD
 ### Runs del projecto
 File <- read.csv(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Reads_Raw.csv"))
 File$Sample <- gsub("\\-",".",File$Sample)
 File$ws_sample_name <- substring(File$Sample, 1, 6)
 ### Selection od samples in abundance file that match with Reads file
 OTUs <- read.csv(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance.csv"))
 MUESTRAS <- colnames(File)

 ##### ANALYSIS
 Samps_types <- paste0("select muestra.c_muestra_wineseq, tipo_muestra.d_tipo_muestra from muestra, tipo_muestra
                       where muestra.c_tipo_muestra = tipo_muestra.c_tipo_muestra
                       order by c_muestra_wineseq")
 
 Tipos_muestra <- dbGetQuery(local_DB, Samps_types)
 colnames(Tipos_muestra) <- c("ws_sample_name", "tipo_muestra")
 
 Bad_Reads_number <- File %>%
   filter(Reads < 20000)
 Bad_Reads_number <- merge(Tipos_muestra, Bad_Reads_number) %>%
   filter(!grepl("BMKCPF",ws_sample_name))
 
 Bad_Samples <- OTUs %>%
   select(Species, one_of(as.character(Bad_Reads_number$Sample)))
 
 BadReadstxt <- data.frame()
 for(i in as.character(Bad_Reads_number$Sample)) {
   Samp_data <- Sample_Row(i)
   Samp <- Bad_Samples %>%
     select_(.dots = i,"Species") %>%
     filter_(.dots = paste0(i, "> 0"))
   Reads_finales <- sum(Samp[i])
   TXT <- TEXTO(i, Samp_data, chain, opt$date, Reads_finales, "BAD")
   BadReadstxt <- rbind(BadReadstxt, TXT)
 }
 dir.create(paste0(opt$Path,opt$date,"/Bad_Reads/",chain,"/Informs/"), showWarnings = FALSE, recursive = TRUE)
 write.table(BadReadstxt, file = paste0(opt$Path,opt$date,"/Bad_Reads/",chain,"/Informs/Bad_Reads_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 print(paste0("Finished analysis of ",chain))
 dbDisconnect(local_DB)
} 

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    5     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")

