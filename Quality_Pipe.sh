DATE=$1

Rscript ./02-Otu_Ab/OTU_to_Ab.R -d $DATE #-v 0.8 -p PATH

Rscript ./03-Quality_Control/01-Sample_Separation.R -d $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_DBNAME #-p PATH

Rscript ./03-Quality_Control/02-BMK_Filter.R -d $DATE #-p PATH

Rscript ./03-Quality_Control/03-Fermented_Filter.R -d $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_DBNAME #-p PATH

Rscript ./03-Quality_Control/04-Soil_Grape_Filter.R -d $DATE #-f 20171112 -p PATH

Rscript ./03-Quality_Control/05-Bad_Reads_Inform.R -d $DATE -U $DB_USER -P $DB_PASSWORD -H $DB_HOST -g $DB_PORT -n $DB_DBNAME #-p PATH

Rscript ./03-Quality_Control/06-Inform_Generation.R -d $DATE #-p PATH
