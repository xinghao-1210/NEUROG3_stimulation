#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 11:59:52 2021

@author: user
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns

hioe_path = f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction/graphlet/edgelist'
panc_path = f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/graphlet/edgelist'
 
for tissue_path in [hioe_path, panc_path]:
    for dgcd in Path(tissue_path).rglob('DGCD-129.txt'):
        dgcd_df=pd.read_csv(dgcd, sep='\t', index_col=0)
        if 'comparison' in os.path.dirname(dgcd):
            tissue='Comparison_'
            dgcd_df.columns=['_'.join(os.path.basename(x).split('.')[0].split('_')[0:4]) for x in dgcd_df.columns]
            dgcd_df.index=['_'.join(os.path.basename(x).split('.')[0].split('_')[0:4]) for x in dgcd_df.index]
        elif 'pancreas' in os.path.dirname(dgcd):
            tissue='Pancreas_'
            dgcd_df.columns=['_'.join((tissue+os.path.basename(x)).split('.')[0].split('_')[0:4]) for x in dgcd_df.columns]
            dgcd_df.index=['_'.join((tissue+os.path.basename(x)).split('.')[0].split('_')[0:4]) for x in dgcd_df.index]
        elif 'HIOE' in os.path.dirname(dgcd):
            tissue='HIOE_'
            dgcd_df.columns=['_'.join((tissue+os.path.basename(x)).split('.')[0].split('_')[0:4]) for x in dgcd_df.columns]
            dgcd_df.index=['_'.join((tissue+os.path.basename(x)).split('.')[0].split('_')[0:4]) for x in dgcd_df.index]            
            
        ft=os.path.dirname(dgcd).split('/')[-1]
        
        heatmap=sns.clustermap(dgcd_df, cmap='YlOrRd_r')
        heatmap.fig.suptitle(f'{tissue}{ft}_graphlet_Distance')
        heatmap.savefig(f'{os.path.dirname(dgcd)}/{tissue}{ft}_clustermap.png')

