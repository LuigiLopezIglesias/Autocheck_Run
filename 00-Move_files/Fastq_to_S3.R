require(optparse)
require(stringr)
library(crayon)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-w', '--warehouse'), action='store',
              dest='Warehouse', type='character',
              default="Hard Disk")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat(blue("Run date to upload is: "%+%green$bold(opt$date)%+%"\n"))
cat(blue("Location of run folder is on: "%+%green$bold(opt$Warehouse)%+%"\n"))
#########################################################################################

### Archivos descargados de sFTP
cat(silver("\n Making list of Run \n"))
PATH <- ifelse(opt$Warehouse == "Local",
       paste0("/home/yamishakka/Escritorio/Runs_rename/",opt$date,"/"),
       paste0("/media/yamishakka/Elements/Runs_Good/",opt$date,"/"))
comando <- paste0("ls ",PATH)
Samples <- system (comando, inter = TRUE)

### Move no Fastq files
#Path <- paste0(Path,Date)
Folder <- paste0("s3://raw-illumina-runs/RUN_",opt$date,"/FILES/")
TXTFiles <- paste0("aws s3 mv ",PATH,"/FastqSummaryF1L1.txt ",Folder)
  system(TXTFiles, inter = TRUE)

XMLFiles <- paste0("aws s3 mv ",PATH,"/config.xml ",Folder)
  system(XMLFiles, inter = TRUE)

RunInfo <- paste0("aws s3 mv ",PATH,"/RunInfo.xml ",Folder)
  system(RunInfo, inter = TRUE)

RunParameter <- paste0("aws s3 mv ",PATH,"/runParameters.xml ",Folder)
  system(RunParameter, inter = TRUE)

SampSheet <- paste0("aws s3 mv ",PATH,"/SampleSheet.csv ",Folder)
  system(SampSheet, inter = TRUE)

### Seleccion de nombre de muestra
cat(silver("\n Making Sample name \n"))
mysamps <- vector()
for (i in Samples) {
  Split_pos <- str_sub(i, 1, (str_locate(i, "_")-1))
  mysamps = rbind(mysamps, Split_pos[1])
}
mysamps <- mysamps[!is.na(mysamps)]
mysamps <- unique(mysamps, incomparables = FALSE)

#### Nombre de las muestras a guardar

for (SamP in mysamps) {
  cat("Uploadin sample "%+%yellow(SamP)%+%" to ")
  Folder <- paste0("s3://raw-illumina-runs/RUN_",opt$date,"/FASTQ/",opt$date,"_",SamP,"/")
  UpCommand <- paste0("aws s3 sync ",PATH," ",Folder," --exclude '*' --include '",SamP,"_*'")
  system(UpCommand, inter = TRUE)
  cat(blue(Folder)%+%"\n")
}

cat(red$bold("O       o O       o O       o          O       o O       o O       o\n"))
cat(red$bold("| O   o | | O   o | | O   o |  UPLOAD  | O   o | | O   o | | O   o |\n"))
cat(red$bold("| | O | | | | O | | | | O | |          | | O | | | | O | | | | O | |\n"))
cat(red$bold("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |\n"))
cat(red$bold("o       O o       O o       O          o       O o       O o       O\n"))
