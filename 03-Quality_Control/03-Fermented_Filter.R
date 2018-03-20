library(optparse)
library(dplyr)
library(RPostgreSQL)
library(crayon)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/"),
  make_option(c('-f', '--filterdate'), action='store',
              dest='filDate', type='integer'),
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
cat(blue("The Database info is: \n   User - "%+%green$bold(opt$user,"\n   ")
        %+%"Password - "%+%green$bold("***************\n   ")
        %+%"Host - "%+%green$bold(opt$host,"\n   ")
        %+%"Password - "%+%green$bold("****\n   ")
        %+%"DB name - "%+%green$bold(opt$dbname,"\n")))
PwD <- getwd()
##### FUNCTIONS
Stage_Substage <- function(Muestra,conexion) 
{SampTypeSamps <- paste0("select m.c_muestra_wineseq, m.c_tipo_muestra, fs.stage, fss.substage from  muestra m
                            left join  muestra_ferm mf on mf.id_muestra = m.id
                            left join ferm_stage fs on fs.id = mf.stage
                            left join ferm_substage fss on fss.id = mf.substage
                            where m.c_muestra_wineseq = '",Muestra,"'
                            order by m.c_muestra_wineseq")
CheckID <- dbGetQuery(conexion, SampTypeSamps)
}
### Column Selection 
Sample_Row <- function(Muestra)
{Samp <- Reads %>%
  filter(Sample %in% Muestra)
}
## Positive samples by species
Ferm_Gr <- function(Muestra)
{Samp <- Fer_GR %>%
  select(Muestra, Species) %>%
  filter_(.dots = paste0(Muestra, "> 0"))
}
### Text generation
TEXTO <- function(Muestra, Sample_Type, Estado, Cadena, Projecto, Reads_F, Sps_Control, Ferm_sp, Num_Sp, Top, Cont_Num, Cont_perc, Cont_Sp, Resultado)
{Texto <- data.frame(Sample=as.character(Muestra),
                     Sample_type=as.character(Sample_Type[1,5]),
                     Stage=as.character(Estado[3]),
                     Substage=as.character(Estado[4]),
                     Chain_Type=Cadena,
                     Project=Projecto,
                     Initial_reads=Sample_Type[1,3],
                     Final_reads=Reads_F,
                     Contamination_sp=Sps_Control,
                     Fermentativ_num=Ferm_sp,
                     Sp_Num = nrow(Num_Sp),
                     Sp_Max = Top[1,2],
                     Sp_Max_perc = as.character(Top[1,1]),
                     Control_Num_Sp = as.character(Cont_Num),
                     Control_Max_Sp = Cont_perc,
                     Control_Samp_Type_SP = as.character(Cont_Sp),
                     Diagnostic = as.character(Resultado))
}
### Conexion con la base de datos
 drv <- dbDriver("PostgreSQL")
 local_DB <- dbConnect(drv, user=opt$user,
                       password=opt$password,
                       host=opt$host,
                       port=opt$port, dbname=opt$dbname)

for (chain in c("ITS","16S")) {
  ##### FILE LOAD
  ### Fermentative samples that pass the reads number limit
  if(file.exists(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/",opt$date,"_",chain,".csv"))==TRUE) {
    Reads <- read.csv(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/",opt$date,"_",chain,".csv"))
    Reads$Sample <- gsub("\\-",".",Reads$Sample)
    ### Selection od samples in abundance file that match with Reads file
    OTUs <- read.csv(paste0(opt$Path,opt$date,"/",chain,"_",opt$date,"_Abundance.csv"))
    Fer_GR <- OTUs %>%
      select(Species, one_of(as.character(Reads$Sample)))
    
    ### Filter collection selection
    Filt_Path <- paste0(PwD,"/03-Quality_Control/Filter_Parameters/",chain,"/")
    if (!is.null(opt$filtDate)) {
      F_Date <- opt$filtDate
    } else {
      F_Date <- list.files(path = Filt_Path,
                 pattern = "Must", all.files = FALSE,
                 full.names = FALSE, recursive = FALSE,
                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
      F_Date <- gsub("^.*?_","",F_Date[length(F_Date)])
    }
    cat(blue("Control files used in "%+%chain%+%" are on date: "%+%green$bold(F_Date)%+%"\n"))
    ##### ANALISIS
    Analysis <- data.frame()
    for(i in as.character(Reads$Sample)) {
      j <- sub('\\..*', '', i)
      j <- sub('b', '', j)
      CheckID <- Stage_Substage(j,local_DB)
      ### Sample selection
      Samp <- Ferm_Gr(i)
      Reads_Finales <- sum(Samp[i])
      Samp[i] <- round((Samp[i]/sum(Samp[i]))*100, digits = 7)
      Samp_data <- Sample_Row(i)
      #Contamination and Fermentative check
      if (chain == "ITS") {
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
        
        ### Fermentative check 
        Ferm_sp <- Samp %>%
          filter(grepl("Saccharomyces|Hanseniaspora|Metschnikowia|Torulaspora",Species))
        S_P <- ifelse(nrow(Ferm_sp)>0, sum(Ferm_sp[i]), 0)
        
        Control_2 <- ifelse(S_P < 50,"Low", "High")
        
        Controles <-ifelse(Control == "Low",
                           ifelse(Control_2 == "High",
                                  "GOOD", 
                                  "REVIEW" ), 
                           "REVIEW")
      } else {
        ### Contamination check
        TOP_10 <- Samp %>% 
          arrange_(.dots = paste0("desc(",i,")")) %>% 
          slice(1:10)
        #### Si aparece algun MO que sea contaminacion ponerlo de neuvo
        #Bac_control <- TOP_10 %>%
        #  filter(grepl("Nitrososphaera",Species))   ## <- Add contamination Sp
        #Control <- ifelse(nrow(Bac_control) > 0 ,
        #                  ifelse(sum(Bac_control[1]) > 5,
        #                         "High",
        #                         "Low"),
        #                  "Low")
        Control <- "Low"
        ### Fermentative check 
        Ferm_sp <- Samp %>%
          filter(grepl("Erwinia|Lactobacillus|Lactococcus|Leuconostoc|Pediococcus|Bacillus|Gluconobacter|Gluconacetobacter|Acetobacter|Oenococcus",Species))
        S_P <- ifelse(nrow(Ferm_sp)>0, sum(Ferm_sp[i]), 0)
        
        TOP_20 <- Samp %>% 
          arrange_(.dots = paste0("desc(",i,")")) %>% ##eliminar si solo es presencia y poner abajo TOP_20 con samp
          slice(1:20)
        Ferm_Sp <- TOP_20 %>%                        
          filter(grepl('Erwinia|Lactobacillus|Lactococcus|Leuconostoc|Pediococcus|Bacillus|Gluconobacter|Gluconacetobacter|Acetobacter|Oenococcus', Species))
        Control_2 <- ifelse(nrow(Ferm_Sp) == 0,"Low", "High")
        Controles <-ifelse(Control == "Low",
                           ifelse(Control_2 == "High",
                                  "GOOD", 
                                  "REVIEW" ), 
                           "REVIEW")
      }
      if (is.na(CheckID[3]) == F) {
        if (CheckID[3] == "Must") {
          ### Limits files
          Limits <- read.csv(paste0(Filt_Path,CheckID[3],"_",F_Date,"/limit_Species_",chain,"_",CheckID[3],".csv"))
          Limits <- Limits %>%
            select(genus_specie, avg_PER, sd) %>%
            `colnames<-`(c("Species", "avg_PER", "sd"))
          Num_SP_control <- Limits %>% summarise(num_sp = n())
          mean <- read.csv(paste0(Filt_Path,CheckID[3],"_",F_Date,"/number_Species_",chain,"_",CheckID[3],".csv"))
          MAX <- read.csv(paste0(Filt_Path,CheckID[3],"_",F_Date,"/max_Species_",chain,"_",CheckID[3],".csv"))
          ### Species number analysis
          SpeciesNUM <- ifelse(nrow(Samp) < (mean[2]-mean[3]), 
                               "Low", 
                               "GOOD")
          ### Abundance specie analysis
          SPMAX <- Samp %>%
            filter_(.dots = paste0(i,"==",max(Samp[,1])))
          
          SpeciesMAX <- ifelse(grepl(" sp",SPMAX[1,2]),
                               ifelse(SPMAX[,1] > (MAX[2]+(2*MAX[3])), 
                                      "ALERT", 
                                      "GOOD"),
                               ifelse(SPMAX[,1] > (MAX[2]+MAX[3]), 
                                      "ALERT", 
                                      "GOOD"))
          
          ### Control species type analysis
          limit_x <- Num_SP_control[,1]*10/100 #limite en eje para muy graves
          limit_y <- Num_SP_control[,1]*20/100 #limite en eje para normales
          
          Area_control <- ((limit_x*limit_y) - limit_x^2)
          ### Analisis con especies control de muestra
          comparaTION <- merge(Samp, Limits)
          Num_SP_CM <- comparaTION %>%
            summarise(number = n())
          
          if (Num_SP_CM[,1] < (Num_SP_control[,1]/2)) {
            CheckID[,3] <- paste0("Isn't ",CheckID[3]," ",CheckID[3])
            ControlSpecies <- "Have less than 50% of Sp"
            Diagnostic <- "Visual analysis is needed"
          }else {
            ### Genero sp. filter
            Genero_sp <- comparaTION %>%
              filter(grepl(" sp", Species))
            
            Genero_sp <- Genero_sp %>% 
              mutate(Norm_limit = (avg_PER + 1.5*sd)*100,
                     Alert_limit = (avg_PER + 2.5*sd)*100)
            
            ### Species filter
            Specie_sp <- comparaTION %>%
              filter(!grepl(" sp", Species))
            
            Specie_sp <- Specie_sp %>% 
              mutate(Norm_limit = avg_PER + sd,
                     Alert_limit = avg_PER + 2*sd)
            
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
            Result <- ifelse(SpeciesMAX  == "GOOD", 
                             ifelse(ControlSpecies == "GOOD", 
                                    ifelse(SpeciesNUM == "REVIEW", 
                                           "REVIEW", 
                                           "GOOD"),
                                    "REVIEW"),
                            "REVIEW")
              
             TXT <- TEXTO(i, Samp_data, CheckID, chain, opt$date, Reads_Finales, Control, Control_2, Samp, TOP_10, SpeciesNUM, SpeciesMAX, ControlSpecies, Result)
             Analysis <- rbind(Analysis,TXT)
          }
        } else {
          #print("Fermentation with stage")
          TXT <- TEXTO(i, Samp_data, CheckID, chain, opt$date, Reads_Finales, Control, Control_2, Samp, TOP_10, "NA", "NA", "NA", Controles)
          Analysis <- rbind(Analysis,TXT)
        }
      } else {
        print(paste0(i," have unknown stage"))
        TXT <- TEXTO(i, Samp_data, "NA", chain, opt$date, Reads_Finales, Control, Control_2, Samp, TOP_10, "NA", "NA", "NA", Controles)
        Analysis <- rbind(Analysis,TXT)
      }
    }
    dir.create(paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/Informs/"), showWarnings = FALSE, recursive = TRUE)
    write.table(Analysis, file = paste0(opt$Path,opt$date,"/Wineseq/Good_Reads/Fermented/",chain,"/Informs/Good_Fer_",opt$date,"_",chain,".csv"), col.names = TRUE, row.names = FALSE, sep = ",") 
  } else {
    cat(red$bold(paste0("\n Dont have files needed to analyze ",chain,"\n")))
  }   
}
  
cat(magenta$bold("O       o O       o O       o          O       o O       o O       o\n"))
cat(magenta$bold("| O   o | | O   o | | O   o |  Q_STEP  | O   o | | O   o | | O   o |\n"))
cat(magenta$bold("| | O | | | | O | | | | O | |    3     | | O | | | | O | | | | O | |\n"))
cat(magenta$bold("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |\n"))
cat(magenta$bold("o       O o       O o       O          o       O o       O o       O\n"))
