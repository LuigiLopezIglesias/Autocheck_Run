require(optparse)
require(stringr)

#########################################################################################
###     Cargar argumentos//Load arguments      ###
option_list <- list(
  make_option(c('-d', '--date'), action='store',
              dest='date', type='character'),
  make_option(c('-w', '--warehouse'), action='store',
              dest='Warehouse', type='character',
              default="Hard Disk"),
  make_option(c('-p', '--path'), action='store',
              dest='Path', type='character',
              default="/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(paste0("Run date to upload is: ",opt$date))
print(paste0("Location of run folder is on: ",opt$Warehouse))
print(paste0("Location od pipeline storage files is on: ",opt$Path))
#########################################################################################

### Archivos descargados de sFTP
writeLines("\n Making list of Run \n")
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

### Seleccion de nombre de muestra
writeLines("\n Making Sample name \n")
mysamps <- vector()
for (i in Samples) {
Split_pos <- str_sub(i, 1, (str_locate(i, "_")-1))
mysamps = rbind(mysamps, Split_pos[1])
}
mysamps <- mysamps[!is.na(mysamps)]
mysamps <- unique(mysamps, incomparables = FALSE)

#### Nombre de las muestras a guardar

for (k in mysamps) {
print(k)
  Folder <- paste0("s3://raw-illumina-runs/RUN_",opt$date,"/FASTQ/",opt$date,"_",k,"/")
  UpCommand <- paste0("aws s3 sync ",PATH," ",Folder," --exclude '*' --include '",k,"_*'")
  system(UpCommand, inter = TRUE)
  print(Folder)
}

print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |  UPLOAD  | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |          | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")

