library(dplyr)
library(optparse)

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

print(paste0("Run date to analyze is: ",opt$date))
print(paste0("Location of run folder is on: ",opt$Warehouse))
print(paste0("Location od pipeline storage files is on: ",opt$Path))
#########################################################################################


PATH <- ifelse(opt$Warehouse == "Local",
       paste0("/home/yamishakka/Escritorio/Runs_rename/",opt$date,"/"),
       paste0("/media/yamishakka/Elements/Runs_Good/",opt$date,"/"))

Fastq <- list.files(path = PATH, 
                    pattern = "_R1_", 
                    all.files = FALSE, 
                    full.names = FALSE, 
                    recursive = FALSE, 
                    ignore.case = FALSE, 
                    include.dirs = FALSE, 
                    no.. = FALSE)

writeLines(c("","||||","Ready to anlyze reads number by sample in run","||||",""), con = stdout(), sep = "\n")
#Fastq<-gsub("Undetermined.*","",Fastq)

Reads <- character()
Samples <- character()
for (i in Fastq) {
  sample <- sub('_S.*', '',i)
  #print(sample)
  Samples <- rbind(Samples, sample)
  Count <- system(paste0("cp ",PATH,i," /tmp/fastatmp.gz && gunzip /tmp/fastatmp.gz && wc -l /tmp/fastatmp && rm /tmp/fastatmp"), intern = TRUE)
  Raw_Reads <- as.integer(sub(' /.*', '',Count))
  reads <- Raw_Reads/4
  Reads <- cbind(Reads, reads)
 # print(reads)
}

ReadsNumber = data.frame(Sample = Samples,
                         Reads = as.integer(Reads),
                         stringsAsFactors = FALSE)
ReadsNumber["Date"] <- opt$date

Reads_file <- unique(ReadsNumber) %>%
  filter(!grepl("Undeter",Sample))

Reads_16S <- Reads_file %>%
  filter(grepl("b", Sample))
print(Reads_16S)

Reads_ITS <- Reads_file %>%
  filter(!grepl("b", Sample))
print(Reads_ITS)
dir.create(paste0(opt$Path,opt$date,"/"), showWarnings = FALSE, recursive = TRUE)
write.table(Reads_ITS, file = paste0(opt$Path,opt$date,"/ITS_",opt$date,"_Reads_Raw.csv"), sep=",",row.names = FALSE, quote = FALSE)
write.table(Reads_16S, file = paste0(opt$Path,opt$date,"/16S_",opt$date,"_Reads_Raw.csv"), sep=",",row.names = FALSE, quote = FALSE)


print("O       o O       o O       o          O       o O       o O       o")
print("| O   o | | O   o | | O   o |   STEP   | O   o | | O   o | | O   o |")
print("| | O | | | | O | | | | O | |    1     | | O | | | | O | | | | O | |")
print("| o   O | | o   O | | o   O | FINISHED | o   O | | o   O | | o   O |")
print("o       O o       O o       O          o       O o       O o       O")

