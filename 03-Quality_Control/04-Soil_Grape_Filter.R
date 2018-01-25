library(optparse)
library(dplyr)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/"),
  make_option(c('-f', '--filterdate'), action='store',
              dest='filDate', type='integer')
)
### Load arguments
opt <- parse_args(OptionParser(option_list=option_list))

##### FUNCTIONS
### Column Selection 
Sample_Row <- function(Muestra)
{Samp <- Reads %>%
  filter(Sample %in% Muestra)
}
## Positive samples by species
Sample_Gr <- function(Muestra)
{Samp <- Samp_GR %>%
  select(Muestra, Species) %>%
  filter_(.dots = paste0(Muestra, "> 0"))
}
### Text generation
## Text generation
TEXTO_Soil_Grape <- function(Muestra, Sample_type, Cadena, Projecto, Reads_F, Sps_Control, Num_Sp, Top, Cont_Num, Cont_perc, Cont_Sp, Resultado)
{Texto <- data.frame(Sample=as.character(Muestra),
                     Sample_type=as.character(Sample_type[1,5]),
                     Chain_Type=Cadena,
                     Project=Projecto,
                     Initial_reads=Sample_type[1,3],
                     Final_reads=Reads_F,
                     Contamination_sp=Sps_Control,
                     Sp_Num = nrow(Num_Sp),
                     Sp_Max = Top[1,2],
                     Sp_Max_perc = as.character(Top[1,1]),
                     Control_Num_Sp = as.character(Cont_Num),
                     Control_Max_Sp = Cont_perc,
                     Control_Samp_Type_SP = as.character(Cont_Sp),
                     Diagnostic = as.character(Resultado))
}
##### ANALYSIS
for (chain in c("16S","ITS")) {
 ##### FILE LOAD
 ### Fermentative samples that pass the reads number limit
 Reads <- read.csv(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/",opt$date,"_",chain,".csv"))
 Reads$Sample <- gsub("\\-",".",Reads$Sample)

 ### Selection od samples in abundance file that match with Reads file
 OTUs <- read.csv(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance.csv"))
 Samp_GR <- OTUs %>%
   select(Species, one_of(as.character(Reads$Sample)))

 if (chain == "ITS") {
   Fun_Soil_Grape <- data.frame()
   for(i in as.character(Reads$Sample)) {
     ### Sample selection
     Samp <- Sample_Gr(i)
     Reads_Finales <- sum(Samp[i])
     Samp[i] <- round((Samp[i]/sum(Samp[i]))*100, digits = 7)
     Samp_data <- Sample_Row(i)
     statistics <- Samp %>% summarise(num_sp = n(),
                                      max = max(Samp[i]))
      ### Filter collection selection
      Filt_Path <- paste0(getwd(),"/03-Quality_Control/Filter_Parameters/",chain,"/")
      if (!is.null(opt$filtDate)) {
        F_Date <- opt$filtDate
      } else {
        F_Date <- list.files(path = Filt_Path,
                   pattern = as.character(Samp_data$tipo_muestra), all.files = FALSE,
                   full.names = FALSE, recursive = FALSE,
                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
        F_Date <- gsub("^.*?_","",F_Date[length(F_Date)])
      }

     ### Contamination check
     TOP_10 <- Samp %>% 
       arrange_(.dots = paste0("desc(",i,")")) %>% 
       slice(1:10)
     Fun_control <- TOP_10 %>%
       filter(grepl("Malassezia|Epicoccum",Species))
     Control <- ifelse(nrow(Fun_control) > 0 ,
                       ifelse(sum(Fun_control[1]) > 5,
                              "High",
                              "Low"),
                       "Low")
     ### Limits files
     Limits <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/limit_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     Limits <- Limits %>%
       select(genus_specie, avg_PER, sd) %>%
       `colnames<-`(c("Species", "avg_PER", "sd"))
     Num_SP_control <- Limits %>% summarise(num_sp = n())
     mean <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/number_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     MAX <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/max_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     ### Species number analysis
     SpeciesNUM <- ifelse(statistics[,1] < (mean[2]-mean[3]), 
                          "Low", 
                          "GOOD")
     ### Abundance specie analysis
     SPMAX <- Samp %>%
       filter_(.dots = paste0(i,"==",statistics[,2]))
     SpeciesMAX <- ifelse(grepl(" sp",SPMAX[1,2]),
                          ifelse(statistics[,2] > (MAX[2]+(2*MAX[3])*100), #<- cambiar para diferentes sd con sp. 
                                 "ALERT", 
                                 "GOOD"),
                          ifelse(statistics[,2] > (MAX[2]+MAX[3]*100), 
                                 "ALERT", 
                                 "GOOD"))
     ### Analisis especies control
     limit_x <- Num_SP_control[,1]*10/100 #limite en eje para muy graves
     limit_y <- Num_SP_control[,1]*20/100 #limite en eje para normales
     
     Area_control <- ((limit_x*limit_y) - limit_x^2)
     
     ### Analisis con especies control de muestra
     comparaTION <- merge(Samp, Limits)
     Num_SP_CM <- comparaTION %>%
       summarise(number = n())
     
     if (Num_SP_CM[,1] < (Num_SP_control[,1]/2)) {
       Samp_data$tipo_muestra <- paste0("Isn't ",Samp_data$tipo_muestra)
       ControlSpecies <- "Have less than 50% of Sp"
       Diagnostic <- "REVIEW"
     }else {
       ### Filtro para Genero sp.
       Genero_sp <- comparaTION %>%
         filter(grepl(" sp", Species))
       
       Genero_sp <- Genero_sp %>% 
         mutate(Norm_limit = (avg_PER + 1.5*sd)*100,
                Alert_limit = (avg_PER + 2.5*sd)*100)
       
       ### Filtro para especies
       Specie_sp <- comparaTION %>%
         filter(!grepl(" sp", Species))
       
       Specie_sp <- Specie_sp %>% 
         mutate(Norm_limit = (avg_PER + sd)*100,
                Alert_limit = (avg_PER + 2*sd)*100)
       
       comparaTION <- merge(Specie_sp, Genero_sp, all = TRUE)
       
       comparaTION$reason <- ifelse(comparaTION[i] > comparaTION$Norm_limit, 
                                    ifelse(comparaTION[i] > comparaTION$Alert_limit, 
                                           "Great_problem", 
                                           "Normal_problem"),
                                    "No_problem")
       
       Problems <- comparaTION %>%
         group_by(reason) %>%
         summarise(number = n())
       
       GP <- Problems %>%
         filter(reason == "Great_problem") %>%
         summarise(n = sum(number))
       
       NP <- Problems %>%
         filter(reason == "Normal_problem") %>%
         summarise(n = sum(number))
       
       Area_problema <- (((-2*GP[,1]-NP[,1])/-2) * (NP[,1]-(-2*GP[,1])))-(((-2*GP[,1]-NP[,1])/-2)^2)
       
       ControlSpecies <- ifelse(Area_control < Area_problema[,1], 
                                "ALERT", 
                                "GOOD")
       
       ### hay que cambiarlo por todo el proceso de control de integrales y rectas
       Diagnostic <- ifelse(SpeciesMAX  == "GOOD", 
                            ifelse(ControlSpecies == "GOOD", 
                                   ifelse(SpeciesNUM == "GOOD", 
                                          "GOOD", 
                                          "REVIEW"),
                                   "REVIEW"),
                            "REVIEW")
     }
     
     TXT <- TEXTO_Soil_Grape(i, Samp_data, chain, opt$date, Reads_Finales, Control, Samp, TOP_10, SpeciesNUM, SpeciesMAX, ControlSpecies, Diagnostic)
     Fun_Soil_Grape <- rbind(Fun_Soil_Grape,TXT)
   }
   dir.create(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/"), showWarnings = FALSE, recursive = TRUE)
   write.table(Fun_Soil_Grape, file = paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/Good_SoilGrape_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 } else {
   Bac_Soil_Grape <- data.frame()
   for(i in as.character(Reads$Sample)) {
     ### Sample selection
     Samp <- Sample_Gr(i)
     Reads_Finales <- sum(Samp[i])
     Samp[i] <- round((Samp[i]/sum(Samp[i]))*100, digits = 7)
     Samp_data <- Sample_Row(i)
     statistics <- Samp %>% summarise(num_sp = n(),
                                      max = max(Samp[i]))

     ### Filter collection selection
     Filt_Path <- paste0(getwd(),"/03-Quality_Control/Filter_Parameters/",chain,"/")
     if (!is.null(opt$filtDate)) {
       F_Date <- opt$filtDate
     } else {
       F_Date <- list.files(path = Filt_Path,
                  pattern = as.character(Samp_data$tipo_muestra), all.files = FALSE,
                  full.names = FALSE, recursive = FALSE,
                  ignore.case = FALSE, include.dirs = FALSE, no.. = TRUE)
       F_Date <- gsub("^.*?_","",F_Date[length(F_Date)])
     }

     ### Contamination check
     TOP_10 <- Samp %>% 
       arrange_(.dots = paste0("desc(",i,")")) %>% 
       slice(1:10)
     ### Si aparece algun MO contaminate ponerlo
     #Bac_control <- TOP_10 %>%
     #  filter(grepl("Nitrososphaera",Species))   ## <- Add contamination Sp
     #Control <- ifelse(nrow(Bac_control) > 0 ,
     #                  ifelse(sum(Bac_control[1]) > 5,
     #                         "High",
     #                         "Low"),
     #                  "Low")
     Control <- "Low"
     ### Limits files
     Limits <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/limit_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     Limits <- Limits %>%
       select(genus_specie, avg_PER, sd) %>%
       `colnames<-`(c("Species", "avg_PER", "sd"))
     Num_SP_control <- Limits %>% summarise(num_sp = n())
     mean <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/number_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     MAX <- read.csv(paste0(Filt_Path,Samp_data$tipo_muestra,"_",F_Date,"/max_Species_",chain,"_",Samp_data$tipo_muestra,".csv"))
     ### Species number analysis
     SpeciesNUM <- ifelse(statistics[,1] < (mean[2]-mean[3]), 
                          "Low", 
                          "GOOD")
     ### Abundance specie analysis
     SPMAX <- Samp %>%
       filter_(.dots = paste0(i,"==",statistics[,2]))
     
     SpeciesMAX <- ifelse(grepl(" sp",SPMAX[1,2]),
                          ifelse(statistics[,2] > ((MAX[2]+(2*MAX[3]))*100), #<- cambiar para diferentes sd con sp. 
                                 "ALERT", 
                                 "GOOD"),
                          ifelse(statistics[,2] > (MAX[2]+MAX[3]*100), 
                                 "ALERT", 
                                 "GOOD"))
     ### Analisis especies control
     limit_x <- Num_SP_control[,1]*10/100 #limite en eje para muy graves
     limit_y <- Num_SP_control[,1]*20/100 #limite en eje para normales
     
     Area_control <- ((limit_x*limit_y) - limit_x^2)
     
     ### Analisis con especies control de muestra
     comparaTION <- merge(Samp, Limits)
     Num_SP_CM <- comparaTION %>%
       summarise(number = n())
     
     if (Num_SP_CM[,1] < (Num_SP_control[,1]/2)) {
       Samp_data$tipo_muestra <- paste0("Isn't ",Samp_data$tipo_muestra)
       ControlSpecies <- "Have less than 50% of Sp"
       Diagnostic <- "REVIEW"
     }else {
       ### Filtro para Genero sp.
       Genero_sp <- comparaTION %>%
         filter(grepl(" sp", Species))
       
       Genero_sp <- Genero_sp %>% 
         mutate(Norm_limit = (avg_PER + 1.5*sd)*100,
                Alert_limit = (avg_PER + 2.5*sd)*100)
       
       ### Filtro para especies
       Specie_sp <- comparaTION %>%
         filter(!grepl(" sp", Species))
       
       Specie_sp <- Specie_sp %>% 
         mutate(Norm_limit = (avg_PER + sd)*100,
                Alert_limit = (avg_PER + 2*sd)*100)
       
       comparaTION <- merge(Specie_sp, Genero_sp, all = TRUE)
       
       comparaTION$reason <- ifelse(comparaTION[i] > comparaTION$Norm_limit, 
                                    ifelse(comparaTION[i] > comparaTION$Alert_limit, 
                                           "Great_problem", 
                                           "Normal_problem"),
                                    "No_problem")
       
       Problems <- comparaTION %>%
         group_by(reason) %>%
         summarise(number = n())
       
       GP <- Problems %>%
         filter(reason == "Great_problem") %>%
         summarise(n = sum(number))
       
       NP <- Problems %>%
         filter(reason == "Normal_problem") %>%
         summarise(n = sum(number))
       
       Area_problema <- (((-2*GP[,1]-NP[,1])/-2) * (NP[,1]-(-2*GP[,1])))-(((-2*GP[,1]-NP[,1])/-2)^2)
       
       ControlSpecies <- ifelse(Area_control < Area_problema[,1], 
                                "ALERT", 
                                "GOOD")
       
       ### hay que cambiarlo por todo el proceso de control de integrales y rectas
       Diagnostic <- ifelse(SpeciesMAX  == "GOOD", 
                            ifelse(ControlSpecies == "GOOD", 
                                   ifelse(SpeciesNUM == "GOOD", 
                                          "GOOD", 
                                          "REVIEW"),
                                   "REVIEW"),
                            "REVIEW")
     }
     
     TXT <- TEXTO_Soil_Grape(i, Samp_data, chain, opt$date, Reads_Finales, Control, Samp, TOP_10, SpeciesNUM, SpeciesMAX, ControlSpecies, Diagnostic)
     Bac_Soil_Grape <- rbind(Bac_Soil_Grape,TXT)
   }
   dir.create(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/"), showWarnings = FALSE, recursive = TRUE)
   write.table(Bac_Soil_Grape, file = paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Soil_Grape/",chain,"/Informs/Good_SoilGrape_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",")
 }
 print(paste0("Finished analysis of ",chain))
}

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    4     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")
 
