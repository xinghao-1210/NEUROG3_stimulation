#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 11:36:09 2021

@author: user
"""

import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

import gseapy as gp

os.chdir('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction')   #set working path

# read gmt annotation file to annotation table with gene joined by '|'
# take all gene set as background set
gmt_file='inputs/geneLists/gseaDatabases/human/c2.cp.v7.2.symbols.gmt'
geneset_name='CP_all'

with open(f'{gmt_file}','r') as gmt:
    gl=gmt.readlines()
    
gl=[x.split('\t') for x in gl]
geneset=[x[0].split('_')[0]+' '+x[1].split('/')[-1]+' '+'|'.join(x[2:]) for x in gl]
background=set([x if '\n' not in x else x.split('\n')[0] for sublist in gl for x in sublist[2:]])
    
with open(f'{os.path.dirname(gmt_file)}/{geneset_name}.txt','w') as g:
    g.writelines(geneset)

with open(f'inputs/geneLists/gseaDatabases/human/{geneset_name}_bckgrnd.txt','w') as g_background:
    g_background.writelines('\n'.join(background))



