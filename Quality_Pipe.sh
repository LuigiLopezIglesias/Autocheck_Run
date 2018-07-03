# To add parameters differents to standard you must to put in variable name -letter and info
#
# Conf="-v 0.8"
# PATHx="-p Path"
# QDate="-f date" 
Date=$1

DATE="-d ${Date}"

Conf=
PATHx=
QDate=

Rscript ./02-Otu_Ab/OTU_to_Ab.R $DATE $Conf $PATHx 

Rscript ./03-Quality_Control/01-Sample_Separation.R $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_NAME $PATHx

Rscript ./03-Quality_Control/02-BMK_Filter.R $DATE $PATHx

Rscript ./03-Quality_Control/03-Fermented_Filter.R $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_NAME $PATHx

Rscript ./03-Quality_Control/04-Soil_Grape_Filter.R $DATE $QDate $PATHx

Rscript ./03-Quality_Control/05-Bad_Reads_Inform.R $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_NAME $PATHx

Rscript ./03-Quality_Control/06-Inform_Generation.R $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_NAME $PATHx 
