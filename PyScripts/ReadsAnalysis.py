#!u/sr/bin/env python
# -*- coding: utf-8 -*-

'''
Tools to be called that create Reads analysis files
'''
#··············································································#
## Modules needed
import os
import gzip
import pandas as pd
import GeneralTools as GT
import Query as Q
#··············································································#
def readsCount(Project, marker, FastqPath, ResultPath):
  samples = [i for i in sorted(os.listdir(FastqPath+'/'+Project+'/')) if "_R2_" not in i]
  if marker == '16s':
    samples = [i for i in samples if "b" in i]
  else:  ### else is 'its' marker if is needed other marker add ifelse
    samples = [i for i in samples if "b" not in i]
  files = []
  for fileName in samples:
    with gzip.open(FastqPath+'/'+Project+'/'+fileName, 'rb') as f:
      for fastqLines, l in enumerate(f):
        pass
    files.append({
      "sample": fileName.split("_",1)[0],
      "date": Project,
      "initialReads": (fastqLines+1)/4
    })
    print("File \x1b[1;33;10m{1}\x1b[0m has \x1b[1;31;10m{0}\x1b[0m reads".format((fastqLines+1)/4, fileName.split("_",1)[0]))
  ### Initial reads
  allSamplesInfo = pd.DataFrame(files)
  goodCV = pd.DataFrame(GT.abundanceLoader(Project, marker, ResultPath).sum().reset_index())
  goodCV = goodCV.drop(goodCV.index[len(goodCV)-1])
  goodCV.columns=["sample","finalReads"]
  goodCV = goodCV.merge(allSamplesInfo, on='sample', how='inner')
  goodCV['percReadsLost'] = 100-(goodCV['finalReads']*100/goodCV['initialReads'])
  goodCV.finalReads = goodCV.finalReads.astype(int)
  goodCV.percReadsLost = goodCV.percReadsLost.astype(float)
  goodCV['percReadsLost'] = goodCV['percReadsLost'].round(4)
  #goodCV[['sample', 'initialReads', 'finalReads', 'percReadsLost', 'date']].to_csv(ResultPath+'/'+Project+'/'+marker+'_'+Project+'_Reads_Raw.csv', index=False, decimal=',')
  return goodCV

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
  return pd.DataFrame(Metrics)

def samplesDBInfo(Project, marker, ResultPath):
  ReadsInfo = samplesMetrics(Project, marker, ResultPath)
  ReadsInfo['DBnames'] = ReadsInfo['sample'].str.split('b').str.get(0)
  ReadsInfo['DBnames'] = ReadsInfo['DBnames'].str.split('-').str.get(0)
  myString = "','".join(ReadsInfo['DBnames'])
  metadata = Q.DBMetadataQuery(myString)
  metadata.columns = ['DBnames', 'sampleType', 'repeatDNAExtraction', 'repeatPCR16s', 'repeatPCRits']
  allInfo = metadata.merge(ReadsInfo, on='DBnames', how='inner')
  return allInfo

def informationMerge(Project, marker, FastqPath, ResultPath):
  readsInformation = readsCount(Project, marker, FastqPath, ResultPath)
  sampleInformation = samplesDBInfo(Project, marker, ResultPath)
  firstMerge = readsInformation.merge(sampleInformation, on='sample')
  firstMerge['Marker'] = marker
  firstMerge = firstMerge[['sample', 'DBnames', 'date', 'Marker', 'initialReads', 'finalReads', 'percReadsLost', 'sampleType', 'SpNumber', 'MaxSp', 'MaxSp%', 'MaxSpReads', 'repeatDNAExtraction', 'repeatPCR'+marker]]
  firstMerge.to_csv(ResultPath+'/'+Project+'/'+marker+'_'+Project+'_Analysis_Metadata.csv', index=False, decimal=',')
  print(firstMerge)
