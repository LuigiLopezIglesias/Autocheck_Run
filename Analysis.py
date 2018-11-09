import optparse
import PyScripts.GeneralTools    as PyGT
import PyScripts.ExcelGeneration as PyEG
###############################################################################
##  Load arguments
###############################################################################

parser = optparse.OptionParser()

parser.add_option('-p', '--Project', action="store", dest="Project", help="Project name to work", default=False)
#parser.add_option('-X', '--FastqPath', action="store", dest="FastqPath", help="Path where was storaged the fastq downloaded", default=os.getcwd())
#parser.add_option('-G', '--GitPath', action="store", dest="GitPath", help="Path of folder where will be created git repository", default='/home/yamishakka/Escritorio/New_Pipeline/RunAnalysis/')
#parser.add_option('-D', '--HardDisk', action="store_true", dest="Storage", help="Store files in HardDisk")

options, args = parser.parse_args()

print('Project to analyse is \x1b[1;31;10m'+options.Project+'\x1b[0m')

###############################################################################
##  First step: creation of abundance (OTU/Percentaje) files
###############################################################################
PyGT.create_dir('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+options.Project)

PyEG.readsCount(options.Project)
#PyEG.fileDownloader(options.Project)
#PyEG.mappedToOtus(options.Project, '16s')

