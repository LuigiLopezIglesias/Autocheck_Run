#!u/sr/bin/env python
# -*- coding: utf-8 -*-

'''
Tools to be called that create Reads analysis files
'''
#··············································································#
## Modules needed
from urllib2 import Request, urlopen, URLError
import json
import math
import os
import gzip
import pandas       as pd
import GeneralTools as GT
import Query        as Q
#··············································································#
def restrequest(rawrequest):
 request = Request(rawrequest)
 try:
  response = urlopen(request)
  json_string = response.read()
  #print(json_string)
  json_obj = json.loads(json_string)
 except URLError, e:
  print 'Got an error code:', e
  sys.exit()
 return json_obj

def samplesIdDetection(Project, marker):
 #Step 1: Find the project ID from project name first, then have biosamples ID
 #assume fewer than 2 projects with name
 ProjectNameRequest = 'https://api.basespace.illumina.com/v1pre3/users/current/projects?Name=%s_%s&Offset=0&Limit=2&SortDir=Asc&access_token=%s' %(Project, marker.upper(),  os.environ['AccessToken'])
 ProjID = restrequest(ProjectNameRequest)['Response']['Items'][0]['Id']
 request = 'https://api.basespace.illumina.com/v2/biosamples?projectid=%s&access_token=%s&limit=1000&SortBy=BioSampleName' %(ProjID, os.environ['AccessToken'])
 json_obj = restrequest(request)
 #print(json_obj)
 nSamples = len(json_obj['Items'])
 samplelist = []
 for sampleindex in range(nSamples):
 # print(json_obj['Items'][sampleindex]['BioSampleName'])
  sampleid = json_obj['Items'][sampleindex]['Id']
  samplelist.append(sampleid)
 samplecsv = ','.join([str(i) for i in samplelist])
 return(samplecsv)
  
def biosamplesDataset(Project, marker):
 #Step 2: Call API to get datasets based on biosample
 samplecsv = samplesIdDetection(Project, marker.upper())
 request = 'https://api.basespace.illumina.com/v2/datasets?inputbiosamples=%s&access_token=%s&limit=1000' %(samplecsv, os.environ['AccessToken'])
 json_obj = restrequest(request)
 totalCount = int(json_obj['Paging']['TotalCount'])
 noffsets = int(math.ceil(float(totalCount)/1000.0))
 IReads = []
 projectRef = []
 sampleName = []
 for index in range(noffsets):
  offset = 1000*index
  request = 'https://api.basespace.illumina.com/v2/datasets?inputbiosamples=%s&access_token=%s&limit=1000&Offset=%s&SortBy=Name' %(samplecsv, os.environ['AccessToken'], offset)
  json_obj = restrequest(request)
  nDatasets = len(json_obj['Items'])
  for fileindex in range(nDatasets):
    SampName = json_obj['Items'][fileindex]['Name'].split("_",1)[0]
    ProjName = Project #json_obj['Items'][fileindex]['Project']['Name']
    InitialReads = json_obj['Items'][fileindex]['Attributes']['common_fastq']['TotalClustersRaw']
    IReads.append(InitialReads)
    sampleName.append(SampName)
    projectRef.append(ProjName)
 fullDataFrame = pd.DataFrame({'initialReads': IReads, 'date': projectRef, 'sample': sampleName})
 projectDataFrame = fullDataFrame[fullDataFrame['date'].str.contains(Project)]
 sortProjectDataFrame = projectDataFrame.sort_values(by=['sample'])
 return(projectDataFrame)

def ReadsAnalysis(Project, marker, ResultPath):
 allSamplesInfo = biosamplesDataset(Project, marker)
 goodCV = pd.DataFrame(GT.abundanceLoader(Project, marker, ResultPath).sum().reset_index())
 goodCV = goodCV.drop(goodCV.index[len(goodCV)-1])
 goodCV.columns=["sample","finalReads"]
 goodCV = goodCV.merge(allSamplesInfo, on='sample', how='inner')
 goodCV['percReadsLost'] = 100-(goodCV['finalReads']*100/goodCV['initialReads'])
 goodCV.finalReads = goodCV.finalReads.astype(int)
 goodCV.percReadsLost = goodCV.percReadsLost.astype(float)
 goodCV['percReadsLost'] = goodCV['percReadsLost'].round(4)
 return(goodCV)

## Calculation of reads number by list of fastq samples in folder
#def readsCount(Project, marker, FastqPath, ResultPath):
#  samples = [i for i in sorted(os.listdir(FastqPath+'/'+Project+'/')) if "_R2_" not in i]
#  if marker == '16s':
#    samples = [i for i in samples if "b" in i]
#  else:  ### else is 'its' marker if is needed other marker add ifelse
#    samples = [i for i in samples if "b" not in i]
#  files = []
#  for fileName in samples:
#    with gzip.open(FastqPath+'/'+Project+'/'+fileName, 'rb') as f:
#      for fastqLines, l in enumerate(f):
#        pass
#    files.append({
#      "sample": fileName.split("_",1)[0],
#      "date": Project,
#      "initialReads": (fastqLines+1)/4
#    })
#    print("File \x1b[1;33;10m{1}\x1b[0m has \x1b[1;31;10m{0}\x1b[0m reads".format((fastqLines+1)/4, fileName.split("_",1)[0]))
#  ### Initial reads
#  allSamplesInfo = pd.DataFrame(files)
#  goodCV = pd.DataFrame(GT.abundanceLoader(Project, marker, ResultPath).sum().reset_index())
#  goodCV = goodCV.drop(goodCV.index[len(goodCV)-1])
#  goodCV.columns=["sample","finalReads"]
#  goodCV = goodCV.merge(allSamplesInfo, on='sample', how='inner')
#  goodCV['percReadsLost'] = 100-(goodCV['finalReads']*100/goodCV['initialReads'])
#  goodCV.finalReads = goodCV.finalReads.astype(int)
#  goodCV.percReadsLost = goodCV.percReadsLost.astype(float)
#  goodCV['percReadsLost'] = goodCV['percReadsLost'].round(4)
#  return goodCV

def samplesMetrics(Project, marker, ResultPath):
  goodCV = pd.DataFrame(GT.abundanceLoader(Project, marker, ResultPath))
  Metrics = []
  for samp in list(goodCV)[:-1]:
    sample = goodCV[[samp,'Species']]
    groupedSpecies = sample.groupby(['Species']).sum()
    cleanGroupedSpecies = groupedSpecies.loc[groupedSpecies[samp] != 0].sort_values(by=samp, ascending=False)
    principalSpPerc = cleanGroupedSpecies.ix[0,0]*100/cleanGroupedSpecies.sum()
    principalSpName = cleanGroupedSpecies.index.tolist()[0]
    principalSpReads = cleanGroupedSpecies.ix[0,0]
    SampSpNumber = len(cleanGroupedSpecies.index)
    Metrics.append({
      "sample": samp,
      "MaxSp": principalSpName,
      "MaxSp%": float(principalSpPerc),
      "MaxSpReads": principalSpReads,
      "SpNumber": SampSpNumber
    })
  goodCV.groupby(['Species']).sum().round(0).to_csv(ResultPath+'/'+Project+'/'+marker.upper()+'_'+Project+'_Abundance.csv', index=True, decimal='.')
  return pd.DataFrame(Metrics)

def samplesDBInfo(Project, marker, ResultPath):
  ReadsInfo = samplesMetrics(Project, marker, ResultPath)
  ReadsInfo['DBnames'] = ReadsInfo['sample'].str.split('b').str.get(0)
  ReadsInfo['DBnames'] = ReadsInfo['DBnames'].str.split('-').str.get(0)
  myString = "','".join(ReadsInfo['DBnames'])
  metadata = Q.DBMetadataQuery(myString)
  metadata.columns = ['DBnames', 'sampleType', 'repeatDNAExtraction', 'repeatPCR16s', 'repeatPCRits', 'stage', 'substage']
  allInfo = metadata.merge(ReadsInfo, on='DBnames', how='inner')
  return allInfo

def informationMerge(Project, marker, ResultPath):
  readsInformation = ReadsAnalysis(Project, marker, ResultPath)
  sampleInformation = samplesDBInfo(Project, marker, ResultPath)
  firstMerge = readsInformation.merge(sampleInformation, on='sample')
  firstMerge['Marker'] = marker
  firstMerge = firstMerge[['sample', 'DBnames', 'date', 'Marker', 'initialReads', 'finalReads', 'percReadsLost', 'sampleType', 'stage', 'substage', 'SpNumber', 'MaxSp', 'MaxSp%', 'MaxSpReads', 'repeatDNAExtraction', 'repeatPCR'+marker]]
  firstMerge.to_csv(ResultPath+'/'+Project+'/'+marker+'_'+Project+'_Analysis_Metadata.xlsx', index=False, decimal=',')
