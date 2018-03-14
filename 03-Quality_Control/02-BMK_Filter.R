library(optparse)
library(dplyr)
library(crayon)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat(blue("Run date to analyze is: "%+%green$bold(opt$date)%+%"\n"))
cat(blue("The path where input file and output file an folder is: "%+%green$bold(opt$Path)%+%"\n"))

#### FUNCIONES
## Especies por muestra positiva
BMK_Gr <- function(Muestra){
  Samp <- BMK_GR %>%
    select(Muestra, Species) %>%
    filter_(.dots = paste0(Muestra, "> 0"))
}
## Seleccion de columna 
BMK_Row <- function(Muestra){
  Samp <- Reads %>%
    filter(Sample %in% Muestra)
}
### Text de Biomemakers
TEXTO_bmk <- function(Muestra, Sample_type, Cadena, Projecto, Reads_F, Sps_Control, Num_Sp, Top, Resultado)
{Texto <- data.frame(Sample=as.character(Muestra), 
                     Sample_type=as.character(Sample_type[1,5]),
                     Chain_Type=Cadena,
                     Project=Projecto,
                     Initial_reads=Sample_type[1,3],
                     Final_reads=Reads_F,
                     Contamination_sp=Sps_Control,
                     Sp_Num=nrow(Num_Sp),
                     Sp_Max=Top[1,2],
                     Sp_Max_perc=as.character(Top[1,1]),
                     Diagnostic=as.character(Resultado))
}
for (chain in c("16S","ITS")) {
 ### ARCHIVOS
 ### Muestras de BMK que pasan el corte de reads
 if(file.exists(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/",opt$date,"_",chain,".csv"))== TRUE) {
   Reads <- read.csv(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/",opt$date,"_",chain,".csv"))
   Reads$Sample <- gsub("\\-",".",Reads$Sample)
   
   ### seleccion de muestras que estan en Reads
   OTUs <- read.csv(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance.csv"))
   BMK_GR<- OTUs %>%
     select(Species, one_of(as.character(Reads$Sample)))
   
   if (chain == "ITS") {
     BMK <- data.frame()
     for(i in Reads$Sample) {
       Samp_data <- BMK_Row(i)
       Samp <- BMK_Gr(i)
       Reads_Finales <- sum(Samp[i])
       Samp[i] <- round((Samp[i]/sum(Samp[i]))*100, digits = 7)
       TOP_10 <- Samp %>% 
         arrange_(.dots = paste0("desc(",i,")")) %>%
         slice(1:10)
       Fun_control <- TOP_10 %>%
         filter(grepl("Malassezia|Epicoccum",Species))  ## <- Add contamination Sp
       Control <- ifelse(nrow(Fun_control) > 0 ,
                         ifelse(sum(Fun_control[1]) > 5,
                                "High",
                                "Low"),
                         "Low")
       Result <- ifelse(Control == "High",
                        "REVIEW",
                        "GOOD")
       TXT <- TEXTO_bmk(i, Samp_data, chain, opt$date, Reads_Finales, Control, Samp, TOP_10, Result)
       BMK <- rbind(BMK,TXT)
     }
   } else {
     BMK <- data.frame()
     for(i in Reads$Sample) {
       Samp_data <- BMK_Row(i)
       Samp <- BMK_Gr(i)
       Reads_Finales <- sum(Samp[i])
       Samp[i] <- round((Samp[i]/sum(Samp[i]))*100, digits = 7)
       TOP_10 <- Samp %>% 
         arrange_(.dots = paste0("desc(",i,")")) %>% 
         slice(1:10)
       ###### Si aparece alguno para usar como control volver a usar
       #Bac_control <- TOP_10 %>%
       #  filter(grepl("Nitrososphaera",Species))  ## <- Add contamination Sp
       #Control <- ifelse(nrow(Bac_control) > 0 ,
       #                  ifelse(sum(Bac_control[1]) > 5,
       #                         "High",
       #                         "Low"),
       #                  "Low")
       Control <- "Low"
       Result <- ifelse(Control == "High",
                        "REVIEW",
                        "GOOD")
       TXT <- TEXTO_bmk(i, Samp_data, chain, opt$date,  Reads_Finales, Control, Samp, TOP_10, Result)
       BMK <- rbind(BMK,TXT)
     }
   }
   dir.create(paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/Informs/"), showWarnings = FALSE, recursive = TRUE)
   write.table(BMK, file = paste0(opt$Path,opt$date,"/NoWineseq/Good_Reads/",chain,"/Informs/Good_BMK_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
   cat(blue("Finished analysis of "%+%green$bold(chain)%+%"\n"))
  } else {
    cat(red$bold(paste0("\n Dont have files needed to analyze ",chain,"\n")))
  }
}

cat(magenta$bold("O       o O       o O       o          O       o O       o O       o\n"))
cat(magenta$bold("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |\n"))
cat(magenta$bold("| | O | | | | O | | | | O | |    2     | | O | | | | O | | | | O | |\n"))
cat(magenta$bold("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |\n"))
cat(magenta$bold("o       O o       O o       O          o       O o       O o       O\n"))
