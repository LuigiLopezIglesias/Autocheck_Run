#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Tools to be called in first step of pipeline that create excel files
'''
#··············································································#
## Modules needed
import sys
import os
import pandas as pd
import boto3
import botocore
import biom
import gzip
import GeneralTools as GT
import xlsxwriter
import pandas.io.formats.excel
pandas.io.formats.excel.header_style = None
#··············································································#

def fileDownloader(Project, marker):
  aaa = pd.read_csv('/home/yamishakka/Escritorio/New_Pipeline/RUN_'+Project+'_'+marker+'/output/links.csv', sep=',', header=0)
  for index, row in aaa.iterrows():
    Marker = row[0]
    S3path = row[1]

    ## Download mapped file
    s3 = boto3.resource('s3')    
    try:
      s3.Bucket(S3path.split("/",3)[2]).download_file(S3path.split("/",3)[3], '/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+Marker+'_'+Project+'_mapped_reads_tax.biom')
    except botocore.exceptions.ClientError as e:
      if e.response['Error']['Code'] == "404":
        print("The object does not exist.")
      else:
        raise

def abundanceLoader(Project, marker):
  ## Analisis archivo biom
  table = biom.load_table('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+marker+'_'+Project+'_mapped_reads_tax.biom')
  abunInfo = table.to_dataframe()
  metadataInfo = table.metadata_to_dataframe('observation')
  fullTable = pd.concat([abunInfo, metadataInfo], axis=1)
  goodCV = fullTable[fullTable['confidence']>= 0.7].copy()
  goodCV.loc[:,'Species'] = goodCV.loc[:, 'taxonomy_6'].str.replace('s:', '', case=False)
  goodCV.taxonomy_6 = goodCV.taxonomy_6.str[2:]
  goodCV = goodCV.drop(['confidence', 'taxonomy_0', 'taxonomy_1', 'taxonomy_2', 'taxonomy_3', 'taxonomy_4', 'taxonomy_5', 'taxonomy_6'], axis=1)
  goodCV = goodCV.rename(columns={'taxonomy_6': 'Species'})
  goodCV = goodCV.sort_values(by=['Species'])
  return(goodCV)

def readsCount(Project):
  samples = [i for i in sorted(os.listdir('/media/yamishakka/Elements/Runs_Good/'+Project+'/')) if "_R2_" not in i]
  files = []
  for fileName in samples:
    with gzip.open('/media/yamishakka/Elements/Runs_Good/'+Project+'/'+fileName, 'rb') as f:
      for fastqLines, l in enumerate(f):
        pass
    files.append({
      "sample": fileName.split("_",1)[0],
      "date": Project,
      "initialReads": (fastqLines+1)/4
    })
    print("File \x1b[1;31;10m{1}\x1b[0m has \x1b[1;31;10m{0}\x1b[0m reads".format((fastqLines+1)/4, fileName.split("_",1)[0]))
  ### Initial reads
  allSamplesInfo = pd.DataFrame(files)    
  bac = allSamplesInfo[allSamplesInfo['sample'].str.contains("b", case=True)]
  ### Final reads
  goodCV = abundanceLoader(Project, '16s')
  Total = goodCV.sum() #########continuar desde aqui mergear este con allSamplesInfo by samplename
  print(Total)
  bac[['sample', 'initialReads', 'date']].to_csv('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/16s_'+Project+'_Reads_Raw.csv', index=False)
  fun = allSamplesInfo[~allSamplesInfo['sample'].str.contains("b", case=True)]
  fun[['sample', 'initialReads', 'date']].to_csv('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/its_'+Project+'_Reads_Raw.csv', index=False)
  

def mappedToOtus(Project, marker):
  goodCV = abundanceLoader(Project, marker)
  print('Creation of \x1b[1;31;10m'+Project+'\x1b[0m OTU abundance')
  for sample in list(goodCV.columns.values[:-1]):
    aaa = goodCV[[sample, 'Species']]
    bbb = aaa.loc[aaa[sample] != 0]
    ### Creation quantitative OTU number
    print('Creation of abundance file for \x1b[1;31;10m'+sample+'\x1b[0m')
    GT.create_dir('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+marker+'_OTU_Individual')
    writer = pd.ExcelWriter('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+marker+'_OTU_Individual/'+sample+'_OTU_BAC.xlsx', engine='xlsxwriter')
    bbb.to_excel(writer, sheet_name='Sheet1', startrow=1)
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']
    worksheet.write('B2', sample)
    fmt = writer.book.add_format({"font_name": "Liberation Sans",
                                  "bold": False,
                                  "font_size": 10})
    worksheet.set_column('A:A', 5, fmt)
    worksheet.set_column('B:B', 7, fmt)
    worksheet.set_column('C:C', 25, fmt)
    worksheet.set_row(0, 70)
    worksheet.merge_range('A1:C1', '')
    worksheet.insert_image('A1', '/home/yamishakka/Documentos/Biomemakers/Logo-BM_DL_negro-1024x323-e1521645909584.png', {'x_offset': 15, 'y_offset': 10})
    writer.save()
    ## Creation percent Abundance file
    groupedSpecies = bbb.groupby(['Species']).sum()
    sampIndividual = groupedSpecies*100/groupedSpecies.sum()
    sampInd = sampIndividual.sort_values(by=sample, ascending=False)
    #print(sampInd)
    GT.create_dir('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+marker+'_Ab_Individual')
    writer = pd.ExcelWriter('/home/yamishakka/Escritorio/Biomemakers/00-NP_Abundances/'+Project+'/'+marker+'_Ab_Individual/'+sample+'_Baceria.xlsx', engine='xlsxwriter')
    sampInd.to_excel(writer, sheet_name='Sheet1', startrow=1)
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']
    worksheet.write('B2', sample)
    fmt = writer.book.add_format({"font_name": "Liberation Sans",
                                  "bold": False,
                                  "font_size": 10})
    worksheet.set_column('A:A', 25, fmt)
    worksheet.set_column('B:B', 10, fmt)
    worksheet.set_row(0, 70)
    worksheet.merge_range('A1:B1', '')
    worksheet.insert_image('A1', '/home/yamishakka/Documentos/Biomemakers/Logo-BM_DL_negro-1024x323-e1521645909584.png', {'x_offset': 15, 'y_offset': 10})
    writer.save()

  print('Creation of \x1b[1;31;10m'+Project+'\x1b[0m abundance')
  Total = goodCV.sum()
  print(Total.drop(Total.index[len(Total)-1]))

#  for index ,row in SampleGroup.iterrows():
#    DbName = row[1]
#    SampleType = row[2]
#    RunName = row[3]
#    PoolmixName = row[4]
#    ChainType = row[5]
#    ClientFtpPath = row[6]
#
#    ## Si el usuario tiene OTU_Abundance se hace aqui
#
#    ## Abundancias para cliente (siempre)
#    groupedSpecies = goodCV.groupby(['Species']).sum()
#
#    sampIndividual = groupedSpecies.loc[:,[RunName]].copy()
#    sampIndividual = sampIndividual*100/sampIndividual[RunName].sum()
#    sampInd = sampIndividual[sampIndividual[RunName]>0]
#    sampInd = sampInd.sort_values(by=RunName, ascending=False)
#
