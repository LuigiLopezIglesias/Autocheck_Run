#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Several util classes used in the rest of the project
"""

# Includes
import os
import biom
import boto3
import pandas        as pd
import Query         as Q
import ReadsAnalysis as RA
# ---------------------------------------------------------------------------- #
#--# Creation of folder
def create_dir(directory):
  """
  Creates a directory if it do not exists
  :param directory: The name of the directory
  """
  if not os.path.exists(directory):
      os.makedirs(directory)

#--# Creation of files
def FilesMaker(GitPath, FileName, Body):
  """
  Write a file in path
  :param GitPAth : The path of folder
  :param FileName: The name of the file
  :param Body    : Information of file
  """
  file = open(GitPath+'/'+FileName,'w+')
  file.write(Body)

#--# Download files from S3
def fileDownloader(Project, marker, ResultPath, GitPath):
  """
  Knowing the hash of result samples will be uploaded to S3
  :param Project   : project date
  :param marker    : marker asociated to study
  :param ResultPath: Folder where   
  :param GitPath   : Folder where is stoaged the filewith hash result
  """
  gitInfo = pd.read_csv(GitPath+'/output/links.csv', \
                        sep=',', \
                        header=0)
  for index, row in gitInfo.iterrows():
    Marker = row[0]
    S3path = row[1]
    ## Download mapped file
    s3 = boto3.resource('s3')
    print(s3)
    try:
      s3.Bucket(S3path.split("/",3)[2]).download_file(S3path.split("/",3)[3], \
                ResultPath+'/'+Project+'/'+Marker+'_'+Project+'_mapped_reads_tax.biom')
    except botocore.exceptions.ClientError as e:
      if e.response['Error']['Code'] == "404":
        print("The object does not exist.")
      else:
        raise


#--# Output of abundances
def abundanceLoader(Project, marker, ResultPath):
  ## Analisis archivo biom
  table = biom.load_table(ResultPath+'/'+Project+'/'+marker+'_'+Project+'_mapped_reads_tax.biom')
  abunInfo = table.to_dataframe()
  metadataInfo = table.metadata_to_dataframe('observation')
  fullTable = pd.concat([abunInfo, metadataInfo], \
                        axis=1)
  goodCV = fullTable[fullTable['confidence']>= 0.7].copy()
  goodCV.loc[:,'Species'] = goodCV.loc[:, 'taxonomy_6'].str.replace('s:', '', case=False)
  goodCV.taxonomy_6 = goodCV.taxonomy_6.str[2:]
  goodCV = goodCV.drop(['confidence', 'taxonomy_0', 'taxonomy_1', 'taxonomy_2', 'taxonomy_3', 'taxonomy_4', 'taxonomy_5', 'taxonomy_6'], axis=1)
  goodCV = goodCV.rename(columns={'taxonomy_6': 'Species'})
  goodCV = goodCV.sort_values(by=['Species'])
  return(goodCV)

#--# PCoA analysis
def PCoACalculation(File1, ResultPath, Project, marker, Name):
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize = (8,8))
  ax = fig.add_subplot(1,1,1)
  ax.set_xlabel('Principal Component 1', fontsize = 15)
  ax.set_ylabel('Principal Component 2', fontsize = 15)
  ax.set_title('2 component PCA', fontsize = 20)
  targets = list(set(File1['Origin']))
  colors = ['r', 'g']
  for target, color in zip(targets,colors):
      indicesToKeep = File1['Origin'] == target
      ax.scatter(File1.loc[indicesToKeep, 'principal component 1'],
                 File1.loc[indicesToKeep, 'principal component 2'], \
                 c = color, \
                 s = 50)
  ax.legend(targets)
  ax.grid()
  create_dir(ResultPath+'/'+Project+'/'+marker+'_Graphics')
  fig.savefig(ResultPath+'/'+Project+'/'+marker+'_Graphics/'+marker+'_'+Name+'.png')

#--# Analysis of PCoA
def abundanceAnalysis(Project, marker, FastqPath, ResultPath):
  from sklearn.preprocessing import StandardScaler
  from sklearn.decomposition import PCA
  metadata = pd.read_csv(ResultPath+'/'+Project+'/'+marker+'_'+Project+'_Analysis_Metadata.xlsx')
  sampleType = set(metadata['sampleType'])
  for ST in sampleType:
    # Move query here 
    if ST in {'Negative control', 'Soil control', 'Grape control', 'Plant Pine'}:
      print('\x1b[1;31;10m'+ST+'\x1b[0m is not evaluated')
    else:
      print('Query to have abundance of all samples of \x1b[1;32;10m'+ST+'\x1b[0m')
      allSamples = Q.fullAbundancesQuery(ST,marker)
      DBsamplesMatrix = allSamples.pivot(index='genus_specie', \
                                         columns='c_muestra_wineseq', \
                                         values='num_reads').fillna(0)
      ClientsamplesMatrix = abundanceLoader(Project, \
                                            marker, \
                                            ResultPath).groupby(['Species']).sum()
      metadatab = metadata[metadata['sampleType'].str.contains(ST)]['sample']
      ClientsamplesMatrix = ClientsamplesMatrix[ClientsamplesMatrix.columns.intersection(metadatab)]
      allSamples = ClientsamplesMatrix.merge(DBsamplesMatrix, \
                                             left_index=True, \
                                             right_index=True, \
                                             how = 'outer').fillna(0)
      # Standardizing the features
      normAllSamples = StandardScaler().fit_transform(allSamples.T)
      MONormalized = pd.DataFrame(data = normAllSamples, columns = list(allSamples.T))
      pca = PCA(n_components=2)
      principalComponents = pca.fit_transform(allSamples.T)
      principalDf = pd.DataFrame(data = principalComponents, \
                                 columns = ['principal component 1', 'principal component 2'])
      principalDf['Samples'] = allSamples.T.index.tolist()
      principalDf['Origin'] = principalDf['Samples'].isin(metadatab) 
      ## Graphics for sample in run group by sample type against all samples of sample type
      PCoACalculation(principalDf, ResultPath, Project, marker, ST)
      for samp in metadatab:
        principalDf['Origin'] = principalDf['Samples'].str.contains(samp)
        ## Graphic by sample against their sample type
        PCoACalculation(principalDf, ResultPath, Project, marker, samp)
