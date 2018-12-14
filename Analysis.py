import os
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

parser.add_option('-d', '--Date', action="store", dest="Date", help="Project name to work", default=False)
parser.add_option('-X', '--FastqPath', action="store", dest="FastqPath", help="Path where was storaged the fastq downloaded", default=os.getcwd()+'FastqPath')
parser.add_option('-R', '--ResultPath', action="store", dest="ResultPath", help="Path where result will be stored", default=os.getcwd()+'Results')
parser.add_option('-G', '--GitPath', action="store", dest="GitPath", help="Path of folder where will be created git repository", default=os.getcwd()+'Run4Jenkins') #URL del git Run4Jenkins

options, args = parser.parse_args()
print('Project to analyse is \x1b[1;31;10m'+options.Date+'\x1b[0m')

###############################################################################
##  First step: creation of abundance (OTU/Percentaje) files
###############################################################################
PyGT.create_dir(options.ResultPath+'/'+options.Date)
PyGT.create_dir(options.GitPath)
repo = git.Repo.clone_from('git@github.com:BiomeMakers/Run4Jenkins.git', options.GitPath)
#repo = git.Repo(options.GitPath)
#cloned_repo = repo.clone(os.path.join(rw_dir, 'to/this/path'))
#repo.clone_from(git_url, repo_dir)
for marker in ['16s', 'its']:
  # Change branch to take information about hash
  print("Changing branch to \x1b[1;3i6;10m"+options.Date+"_"+marker+"\x1b[0m")
  repo.git.checkout(options.Date+'_'+marker.upper())
  # Download biom files
  print("Downloading mapped biom file to \x1b[1;36;10m"+options.ResultPath+"\x1b[0m")
  PyGT.fileDownloader(options.Date, marker, options.ResultPath, options.GitPath)
  # OTU files creation
  print("Creation of abundance files in \x1b[1;36;10m"+options.ResultPath+"\x1b[0m")
  PyEG.mappedToOtus(options.Date, marker, options.ResultPath)
  ## metadata information
  PyRA.informationMerge(options.Date, marker, options.FastqPath, options.ResultPath)
  # PCoA analysis
  PyGT.abundanceAnalysis(options.Date, marker, options.FastqPath, options.ResultPath)
