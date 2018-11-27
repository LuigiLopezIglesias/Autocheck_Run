import optparse
import PyScripts.GeneralTools    as PyGT
import PyScripts.ExcelGeneration as PyEG
import PyScripts.ReadsAnalysis   as PyRA
import PyScripts.Query           as PyQ
import git
###############################################################################
##  Load arguments
###############################################################################

parser = optparse.OptionParser()

parser.add_option('-p', '--Project', action="store", dest="Project", help="Project name to work", default=False)
parser.add_option('-X', '--FastqPath', action="store", dest="FastqPath", help="Path where was storaged the fastq downloaded", default='/media/yamishakka/Elements/Runs_Good')
parser.add_option('-R', '--ResultPath', action="store", dest="ResultPath", help="Path where result will be stored", default='/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances')
parser.add_option('-G', '--GitPath', action="store", dest="GitPath", help="Path of folder where will be created git repository", default='/home/yamishakka/Escritorio/New_Pipeline/Run4Jenkins')

options, args = parser.parse_args()

print('Project to analyse is \x1b[1;31;10m'+options.Project+'\x1b[0m')

###############################################################################
##  First step: creation of abundance (OTU/Percentaje) files
###############################################################################
PyGT.create_dir(options.ResultPath+'/'+options.Project)
repo = git.Repo(options.GitPath )
for marker in ['16s', 'its']:
  # Change branch to take information about hash
  repo.git.checkout(options.Project+'_'+marker.upper())
  # Download biom files
  PyGT.fileDownloader(options.Project, marker, options.ResultPath, options.GitPath)
  # OTU files creation
  PyEG.mappedToOtus(options.Project, marker, options.ResultPath)
  ## metadata information
  PyRA.informationMerge(options.Project, marker, options.FastqPath, options.ResultPath)
  # PCoA analysis
  PyGT.abundanceAnalysis(options.Project, marker, options.FastqPath, options.ResultPath)
