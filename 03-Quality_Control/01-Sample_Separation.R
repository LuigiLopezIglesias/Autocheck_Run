library(optparse)
library(dplyr)
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

#print(opt)
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
 MUESTRAS <- File$Sample
 File$ws_sample_name <- gsub('-.*', '', File$Sample)
 
 File$ws_sample_name <- gsub('b.*', '', File$ws_sample_name) 
 
 #### ANALYSIS
 WS_Samps <- paste0("select muestra.c_muestra_wineseq, tipo_muestra.d_tipo_muestra from muestra, tipo_muestra
                    where muestra.c_tipo_muestra = tipo_muestra.c_tipo_muestra
                    and (tipo_muestra.d_tipo_muestra = 'Soil'
                    or tipo_muestra.d_tipo_muestra = 'Grape'
                    or tipo_muestra.d_tipo_muestra = 'Fermented')
                    order by c_muestra_wineseq")
 
 WineSeq <- dbGetQuery(local_DB, WS_Samps)
 colnames(WineSeq) <- c("ws_sample_name", "tipo_muestra")
 
 NWS_Samps <- paste0("select muestra.c_muestra_wineseq, tipo_muestra.d_tipo_muestra from muestra, tipo_muestra
                     where muestra.c_tipo_muestra = tipo_muestra.c_tipo_muestra
                     and (tipo_muestra.d_tipo_muestra != 'Soil'
                     and tipo_muestra.d_tipo_muestra != 'Grape'
                     and tipo_muestra.d_tipo_muestra != 'Fermented')
                     order by c_muestra_wineseq")
 
 NoWineSeq <- dbGetQuery(local_DB, NWS_Samps)
 colnames(NoWineSeq) <- c("ws_sample_name", "tipo_muestra")
 
 ### Separacion muestras WS y NWS
 WS <- merge(File, WineSeq)
 NWS <- merge(File, NoWineSeq)
 
 
 ### Separacion de muestras con Bajo reads de buen numero de reads (profundidad)
 ###~~ Muestras Wineseq
 BR_WS <- WS %>%
   filter(Reads < 20000)
 dir.create(paste0(opt$Path,opt$date,"/Wineseq/Bad_Reads/",chain), showWarnings = FALSE, recursive = TRUE)
 write.table(BR_WS, file = paste0(opt$Path,opt$date,"/Wineseq/Bad_Reads/",chain,"/",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 GR_WS <- WS %>%
   filter(Reads >= 20000)
 
 ###-- Separacion por tipo de muestra 
 ## Suelo uva
 ST_WS <- GR_WS %>%
   filter(tipo_muestra %in% c('Soil','Grape'))
 dir.create(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain), showWarnings = FALSE, recursive = TRUE)
 write.table(ST_WS, file = paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 ## Fermentados
 FER_WS <- GR_WS %>%
   filter(tipo_muestra %in% 'Fermented')
 dir.create(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain), showWarnings = FALSE, recursive = TRUE)
 write.table(FER_WS, file = paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 ##~~ Muestras no Wineseq
 
 BR_NWS <- NWS %>%
   filter(Reads < 20000)
 dir.create(paste0(opt$Path,opt$date,"/NoWineseq/Bad_Reads/",chain), showWarnings = FALSE, recursive = TRUE)
 write.table(BR_NWS, file = paste0(opt$Path,opt$date,"/NoWineseq/Bad_Reads/",chain,"/",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 GR_NWS <- NWS %>%
   filter(Reads >= 20000)
 dir.create(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain), showWarnings = FALSE, recursive = TRUE)
 write.table(GR_NWS, file = paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 
 print(paste0("Finished analysis of ",chain))
 dbDisconnect(local_DB)
}

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    1     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")
