#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Several util classes used in the rest of the project
"""

# Includes
import sqlalchemy as sa
import pandas     as pd
import os
#------------------------------------------------------------------------------#
def DBconexion():
  return sa.create_engine("postgresql://"+os.environ['DB_USER']+":"+os.environ['DB_PASSWORD']+"@"+os.environ['DB_HOST']+"/"+os.environ['DB_NAME'])
#------------------------------------------------------------------------------#
def DBMetadataQuery(samplelist):
  Text = ('select m.c_muestra_wineseq, tm.d_tipo_muestra, mi.repeat_dna_extraction, mi.repeat_pcr_16s, mi.repeat_pcr_its, fs.stage, fss.substage  '
          'from muestra m'
          ' join tipo_muestra tm on tm.c_tipo_muestra = m.c_tipo_muestra'
          ' join muestra_internal mi on m.id = mi.id_muestra'
          ' left join muestra_ferm mf on mf.id_muestra = m.id'
          ' left join ferm_stage fs on fs.id = mf.stage'
          ' left join ferm_substage fss on fss.id = mf.substage'
          ' where m.c_muestra_wineseq in (\''+samplelist+'\')')
  return pd.read_sql(Text, DBconexion())
#------------------------------------------------------------------------------#
def fullAbundancesQuery(sampleType, marker):
  Text= ('select m.c_muestra_wineseq, concat(gs.genus,\' \',gs.species) as genus_specie, a.num_reads'
         '  from muestra m'
         '  join abundances    a  on m.id = a.sample_id'
         '  join chain_types   ct on a.chain_type_id = ct.id'
         '  join genus_species gs on a.gs_id = gs.id'
         '  join tipo_muestra  tm on tm.c_tipo_muestra = m.c_tipo_muestra'
         ' where ct.name = \''+marker.upper()+'\''
         '   and tm.d_tipo_muestra = \''+sampleType+'\''
         '   and m.c_muestra_wineseq in (select c_muestra_wineseq'
         '                               from (select m.c_muestra_wineseq, ct.name, sum(num_reads) '
         '                                      from abundances a'
         '                                      join muestra       m on              m.id = a.sample_id'
         '                                      join chain_types  ct on             ct.id = a.chain_type_id'
         '                                      join tipo_muestra tm on tm.c_tipo_muestra = m.c_tipo_muestra'
         '                                     where tm.d_tipo_muestra = \''+sampleType+'\''
         '                                       and ct.name = \''+marker.upper()+'\' '
         '                                       and m.id_cliente not in (SELECT id '
         '                                                                  FROM cliente' 
         '                                                                 WHERE c_wineseq IN (\'000\', \'999\', \'P00\', \'P01\'))'
         '                                     group by ct.name , m.id, c_muestra_wineseq) as tt'
         '                                     where tt.sum > 10000 '
         '                                     order by c_muestra_wineseq);')
  return pd.read_sql(Text, DBconexion())
