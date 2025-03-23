#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 16:54:55 2021

@author: user
"""


import os
import math
import random
import glob
import pandas as pd
import numpy as np
import networkx as nx
import dynetx as dnx
import ast
from itertools import islice

import ndlib.models.ModelConfig as mc
import ndlib.models.dynamic as dm
import ndlib.models.epidemics as ep
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend

from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from bokeh.io import show

os.chdir('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction')

filt=['all', 'all_tfs', 'de', 'de_tfs']
graph_dict={}
graph_ts_dict={}
for ft in filt:
    nx_path = f'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}'
    
    Gl={}
    comb_df=pd.DataFrame()
    for f in glob.glob(f'{nx_path}/*all/graphml*None.txt', recursive=True):
        df=pd.read_csv(f, sep='\t', header=None)
        df.columns=['TF', 'Target', 'attr']
        t=int(f'{f.split("/")[-1].split("_")[1]}')
        df['Time']=t
        df_attr=pd.DataFrame([ast.literal_eval(x) for x in df.attr.tolist()])
        #df_attr['weight'] = df_attr.apply(lambda x: -(x['weight']) if x['SignedQuantile']<0 else x['weight'], axis=1)
        df.drop(columns='attr',inplace=True)
        df=df.join(df_attr)
        Gl[t]=df
        comb_df=comb_df.append(df)
    
    
    agg_comb_df=comb_df.groupby(['TF','Target']).agg({'Time':lambda x: sorted(x)}).join(
                comb_df.groupby(['TF','Target']).agg({k:lambda x: tuple(x)[0] for k in comb_df.columns[2:] if k!='Time'}))
    
    
    dynamic_graph = dnx.DynDiGraph()
    for t, df in enumerate([v for k,v in sorted(Gl.items(), key=lambda item: item[0])]):
        edges=tuple(zip(df.TF, df.Target))
        dynamic_graph.add_interactions_from(edges, t=t)
    
    # agg_comb_df.to_csv(nx_path+'/dnx_interaction_df.csv', sep='\t', index=True)
    # dnx.write_interactions(dynamic_graph, path=nx_path+'/dnx_interaction.txt', delimiter='\t', encoding='utf-8')
    
    

    #============================================
    # Combined graph 
    #============================================
    # graph from df
    static_graph=nx.from_pandas_edgelist(agg_comb_df.reset_index(),source='TF', target='Target', edge_attr=True, create_using=nx.DiGraph())
    graph_dict[ft]=static_graph
    
    graph_ts_dict[ft]=sorted([(k,nx.from_pandas_edgelist(g.reset_index(),source='TF', target='Target', edge_attr=True, create_using=nx.DiGraph())) for k,g in Gl.items()], key=lambda x: x[0])



#============================================
# TF-Hormone over combined/temporal graph
#============================================
graph_list_tp=[(x[0], x[1], k) for k,v in graph_ts_dict.items() for x in v ]
graph_list_all=[('static', v, k) for k,v in graph_dict.items()]
graph_list=graph_list_tp + graph_list_all
for graph_ in graph_list:
    graph=graph_[1]
    if  graph_[0]=='static':
        ft=graph_[2] 
        nx_path = f'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}' 
        print(graph_,'\n',nx_path)
    else:
        ft=graph_[2]            
        t=graph_[0]
        nx_path = f'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/{t}hpi all'
        print(graph_,'\n',nx_path)
    
    #============================================
    # Hierarchy graph non-tree 
    #============================================
    if ft in ['all_tfs', 'de_tfs']:
        pos=nx.nx_agraph.graphviz_layout(graph, prog='dot')
        plt.figure(figsize=(20,20))
        nx.draw_networkx(graph, pos=pos, arrowsize=10,  width=0.25, node_size=30,font_size=10, 
                          edge_color=[graph[u][v]['stroke'] for u,v in graph.edges], 
                          with_labels=True)
        plt.title(f'{ft}_hierarchy')
        plt.savefig(f'{nx_path}/{ft}_hierarchy.png',bbox_inches='tight', dpi=300)
        
        #============================================
        # Tree structure TF
        #============================================  
        try:
            os.mkdir(f'{nx_path}/TF hierarchy')
        except FileExistsError:
            pass
        # Tree structure form idicated node
        root_list=['NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB','MAFA','RFX6','ISL1','ARX','LMX1A','INSM1','MYT1','SOX4','ETV1','FEV','GATA4','GATA6','NEUROD1','NEUROD2','PAX4','PAX6']
        for tree_root in root_list:
            if tree_root in graph.nodes():
                tree = nx.bfs_tree(graph, tree_root)
            else:
                continue
            all_edge = [x for x in list(nx.edge_bfs(graph.subgraph(tree.nodes), [tree_root]))]
            tree_dict={'tree':{}, 'complete':{}}
            
            predecessors=[k[0] for k in tree.edges()]
            for t in tree_dict:
                
                # # hierarchy tree 
                # pos = hierarchy_pos(tree, width=1) 
                # plt.figure(figsize=(20,20))
                # if t=='tree':
                #     nx.draw_networkx(tree, pos=pos, width=0.25, node_size=30, font_size=10, arrowsize=10,
                #                      edge_color=[graph[u][v]['stroke'] for u,v in tree.edges],
                #                      with_labels=True)
                # else:
                #     nx.draw_networkx(tree, pos=pos, edgelist=all_edge, width=0.25, node_size=30, font_size=10, arrowsize=7.5, 
                #                      edge_color=[graph[u][v]['stroke'] for u,v in all_edge],
                #                      with_labels=True)
                # plt.title(f'{ft}_{tree_root}_{t}')
                # plt.savefig(f'{nx_path}/TF hierarchy/{ft}_{tree_root}_{t}.png',bbox_inches='tight', dpi=300)
                
                # radial tree 
                pos = hierarchy_pos(tree, tree_root, width = 2*math.pi)
                new_pos = {u:(r*math.cos(theta),r*math.sin(theta)) for u, (theta, r) in pos.items()}
                plt.figure(figsize=(20,20))
                if t=='tree':
                    nx.draw(tree, pos=new_pos, width=0.25, node_size=30, font_size=10, arrowsize=10,
                            edge_color=[graph[u][v]['stroke'] for u,v in tree.edges],
                            with_labels=True)
                else:
                    nx.draw(tree, pos=new_pos, edgelist=all_edge, width=0.25, node_size=30, font_size=10, arrowsize=10,
                            edge_color=[graph[u][v]['stroke'] for u,v in all_edge],
                            with_labels=True)
                nx.draw_networkx_nodes(tree, pos=new_pos, nodelist = [tree_root], node_color = 'blue', node_size=30)
                plt.title(f'{ft}_{tree_root}_radial_{t}')
                plt.savefig(f'{nx_path}/TF hierarchy/{ft}_{tree_root}_radial_{t}.png',bbox_inches='tight', dpi=300)
    
                # dict for tree edges
                for p in predecessors:
                    if t=='tree':
                        tree_dict[t][p]=[k[1] for k in tree.edges() if k[0]==p]
                    else:
                        tree_dict[t][p]=[k[1] for k in graph.edges() if k[0]==p]
    
    
    #============================================
    # tf_hormone pairs direct binding
    #============================================
    if ft in ['all', 'de']:
        try:
            os.mkdir(f'{nx_path}/Hormone hierarchy')
        except FileExistsError:
            pass  
        
        hormone_list = ['INS','GCG','SST','GHRL','GIP','TPH1','TAC1','CCK','PYY','SCT','NTS','GAST']
        tf_hor_list=[x for x in graph.edges if x[1] in hormone_list]
        tf_hor_weight_dict={}
        
        if len(tf_hor_list)!=0:
        #============================================
        # Tree structure hormone
        #============================================  
            for tree_root in list({x[1] for x in tf_hor_list}):
                tree = nx.bfs_tree(graph, tree_root, reverse=True)
                all_edge = [(x[0],x[1]) for x in list(nx.edge_bfs(graph.subgraph([x for x in tree.nodes]), tree_root, orientation='reverse'))]
                hormone_tree_dict={'tree':{}, 'complete':{}}
                
                predecessors=[k[0] for k in tree.edges()]
                for t in ['tree', 'complete']:
                    # radial tree 
                    pos = hierarchy_pos(tree, tree_root, width = 2*math.pi)
                    new_pos = {u:(r*math.cos(theta),r*math.sin(theta)) for u, (theta, r) in pos.items()}
                    plt.figure(figsize=(20,20))
                    if t=='tree':
                        nx.draw(tree, pos=new_pos, width=0.25, node_size=30, font_size=10, arrowsize=10, 
                                edge_color=[graph[v][u]['stroke'] for u,v in tree.edges],
                                with_labels=True)
                    else:
                        nx.draw(tree, pos=new_pos, edgelist=all_edge, width=0.25, node_size=30, font_size=10, arrowsize=10, 
                                edge_color=[graph[u][v]['stroke'] for u,v in all_edge],
                                with_labels=True)
                    nx.draw_networkx_nodes(tree, pos=new_pos, nodelist = [tree_root], node_color = 'blue', node_size=30)
                    plt.title(f'{ft}_{tree_root}_radial_{t}')
                    plt.savefig(f'{nx_path}/Hormone hierarchy/{ft}_{tree_root}_radial_{t}.png',bbox_inches='tight', dpi=300)
                    
                # dict for tree edges
                for p in predecessors:
                    if t=='tree':
                        hormone_tree_dict[t][p]=[k[1] for k in tree.edges() if k[0]==p]
                    else:
                        hormone_tree_dict[t][p]=[k[1] for k in graph.edges() if k[0]==p]   
                
            #============================================
            # tf_hormone weighted binding for correlation
            #============================================
                tree1 = nx.bfs_tree(graph, tree_root, reverse=True, depth_limit=1)
                weight_attr = {(k[1],k[0]):v for k,v in nx.get_edge_attributes(graph.subgraph([x for x in tree1.nodes]), 'SignedQuantile').items() if (k[1],k[0]) in tree1.edges}
                tf_hor_weight_dict.update(weight_attr)
                nx.set_edge_attributes(tree1, weight_attr, 'SignedQuantile')
                
            df_weighted=pd.Series(tf_hor_weight_dict).rename_axis(['Hormone', 'TF']).reset_index(name='SignedQuantile')
            df_weighted=df_weighted.pivot(index='Hormone',columns='TF',values='SignedQuantile').fillna(0)
            df_weighted.to_csv(f'{nx_path}/Hormone hierarchy/tf_hormone_cluster_weighted.csv',sep='\t')
        
            # plot weighted TF-hormone contribution
            fig,ax = plt.subplots(figsize=(12,12))
            ax = sns.heatmap(df_weighted.T, square=True, cmap='bwr', center=0, yticklabels=1, xticklabels=1, 
                             linewidths=.05, linecolor='k',cbar=True)
            plt.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_matrix_weighted.png',bbox_inches='tight', dpi=300)
            
            fig = sns.clustermap(df_weighted, metric="euclidean", method='ward', square=True, cmap='bwr', center=0, yticklabels=1, xticklabels=1, 
                                linewidths=.05, linecolor='k',figsize=(12,4))
            fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_cluster_weighted.png',bbox_inches='tight', dpi=300)
    
            # Correlation plot
            corr_df = df_weighted.T.corr(method='spearman')
            fig = sns.clustermap(corr_df, metric="euclidean", method='ward', square=True,cmap="YlGnBu", yticklabels=1, xticklabels=1, figsize=(15,15))
            fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_corr_weighted.png',bbox_inches='tight', dpi=300)    

            corr_df = df_weighted.corr(method='spearman')
            fig = sns.clustermap(corr_df, metric="euclidean", method='ward', square=True,cmap="YlGnBu", yticklabels=1, xticklabels=1, figsize=(15,15))
            fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_corr_weighted_TF.png',bbox_inches='tight', dpi=300)  
            
            #============================================
            # tf_hormone binding for correlation
            #============================================
            hormone_encode={k:v for v,k in enumerate({x[1] for x in tf_hor_list})}
            tf_encode={k:v for v,k in enumerate({x[0] for x in tf_hor_list})}        
            tf_hor_df = pd.DataFrame(tf_hor_list, columns =['TF', 'Hormone'])
            
            for i,d in enumerate([hormone_encode, tf_encode]):
                if i==0:
                    df=tf_hor_df.pivot(columns='TF', index='Hormone', values='Hormone')
                else:
                    df=tf_hor_df.pivot(columns='Hormone', index='TF', values='TF')
                    
                df.fillna(df.shape[0], inplace=True)
                df.replace(d, inplace=True)
                if i==1:
                    ddf=df.replace({k:(1 if k!=df.shape[0] else 0) for k in range(df.shape[0]+1)})
                    fig,ax = plt.subplots(figsize=(12,12))
                    ax = sns.heatmap(ddf, square=True, cmap=['w','r'], yticklabels=1, xticklabels=1, 
                                     linewidths=.05, linecolor='k',cbar=False)
                    plt.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_matrix.png',bbox_inches='tight', dpi=300)
                    
                
                # Correlation plot
                corr_df = df.corr(method='spearman')
                fig = sns.clustermap(corr_df, metric="euclidean", method='ward', square=True,cmap="YlGnBu", yticklabels=1, xticklabels=1,figsize=(15,15))
                fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_corr_{i}.png',bbox_inches='tight', dpi=300)
            
            # tf_hormone dictionary
            tf_hor_dict={'tf':{}, 'hormone':{}}
            for tf in {x[0] for x in tf_hor_list}:
                tf_hor_dict['tf'][tf]=[x[1] for x in tf_hor_list if x[0]==tf]      
            for h in {x[1] for x in tf_hor_list}:
                tf_hor_dict['hormone'][h]=[x[0] for x in tf_hor_list if x[1]==h]    
          
            #============================================
            # tf_hormone shortest path contribution
            #============================================ 
            tf_list=['NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB','MAFA','RFX6','GATA4','GATA6','NEUROD1','NEUROD2','PAX4','PAX6']+[x.strip() for x in open('inputs/geneSets/Gehart2019_eec.txt')]
            tf_hor_pair=[(x,y) for x in tf_list for y in hormone_list if (x in graph) and (y in graph) and nx.has_path(graph, x, y) and (x not in hormone_list)]
            
            # shortest path by abs(SignedQuantile) and final distance by product of path SignedQuantile
            log_static_graph=graph.copy()
            weight_attr={k:-np.log2(abs(v)) for k,v in nx.get_edge_attributes(graph, 'SignedQuantile').items()}
            nx.set_edge_attributes(log_static_graph, weight_attr, '-log2_weight')
        
            shortest_path={k:{'pos':np.inf, 'pos_path':[], 'pos_node_length':np.inf,
                              'neg':-np.inf, 'neg_path':[], 'neg_node_length':np.inf} for k in tf_hor_pair}
            
            for pair in tf_hor_pair:
                short_path=list(islice(nx.shortest_simple_paths(log_static_graph, pair[0], pair[1], weight='-log2_weight'), 500))
                
                product_dict={}
                for l in range(len(short_path)):
                    product_list=[]
                    for i in range(len(short_path[l])-1):
                        product_list.append(log_static_graph[short_path[l][i]][short_path[l][i+1]]['SignedQuantile'])
                    product_dict['->'.join(short_path[l])]=np.prod(product_list)           
                    
                shortest_path[pair]['pos'],shortest_path[pair]['pos_path'] = max(product_dict.values()),max(product_dict, key=product_dict.get)
                shortest_path[pair]['neg'],shortest_path[pair]['neg_path'] = min(product_dict.values()),min(product_dict, key=product_dict.get)
                shortest_path[pair]['pos_node_length'],shortest_path[pair]['neg_node_length'] = len(shortest_path[pair]['pos_path'].split('->')),len(shortest_path[pair]['neg_path'].split('->'))
                
            shortest_path_df=pd.DataFrame.from_dict(shortest_path,orient='index')
            shortest_path_df.reset_index(inplace=True)
            shortest_path_df.columns=['TF','Hormone'] + shortest_path_df.columns[2:].tolist()
            shortest_path_df.to_csv(f'{nx_path}/Hormone hierarchy/tf_hormone_shortest_path_contribution.csv',sep='\t')
            
            shortest_path_mtx_df=shortest_path_df.pivot(columns='Hormone', index='TF', values=['pos','neg'])
            shortest_path_mtx_df=shortest_path_mtx_df.swaplevel(0, 1, axis=1)
            shortest_path_mtx_df.to_csv(f'{nx_path}/Hormone hierarchy/tf_hormone_shortest_path_contribution_matrix.csv',sep='\t')        
            fig = sns.clustermap(shortest_path_mtx_df.T, method='ward', square=True, cmap='bwr', center=0, xticklabels=1, yticklabels=1, robust=True, figsize=(12,6))
            plt.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_shortest_path_contribution_matrix.png',bbox_inches='tight', dpi=300)        
            
            shortest_path_mtx_df=shortest_path_df.pivot(columns='Hormone', index='TF', values=['pos'])
            shortest_path_mtx_df.columns = shortest_path_mtx_df.columns.droplevel()
            shortest_path_mtx_df = shortest_path_mtx_df.rename_axis(None, axis = 0)
            shortest_path_mtx_df.to_csv(f'{nx_path}/Hormone hierarchy/tf_hormone_shortest_path_positive_contribution_matrix.csv',sep='\t') 
            fig = sns.clustermap(shortest_path_mtx_df.T, method='ward', square=True, cmap='YlOrRd', xticklabels=1, yticklabels=1, robust=True, figsize=(12,6))
            plt.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_shortest_path_positive_contribution_matrix.png',bbox_inches='tight', dpi=300)      
                 
            # Correlation plot
            corr_df = shortest_path_mtx_df.corr(method='spearman')
            fig = sns.clustermap(corr_df, metric="euclidean", method='ward', square=True,cmap="YlGnBu", yticklabels=1, xticklabels=1,figsize=(15,15))
            fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_corr_positive_shortest_path.png',bbox_inches='tight', dpi=300)
                
            corr_df = shortest_path_mtx_df.T.corr(method='spearman')
            fig = sns.clustermap(corr_df, metric="euclidean", method='ward', square=True,cmap="YlGnBu", yticklabels=1, xticklabels=1,figsize=(15,15))
            fig.savefig(f'{nx_path}/Hormone hierarchy/tf_hormone_corr_positive_shortest_path_TF.png',bbox_inches='tight', dpi=300)  
            
#============================================
# merged pancreas+hioe tf-hormone 
#============================================
# tf_hormone weighted binding for correlation
for ft in ['all','de']:
    pan_tfhorm=pd.read_csv(f'../pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/tf_hormone_cluster_weighted.csv', sep='\t', index_col=0)
    hio_tfhorm=pd.read_csv(f'../HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/tf_hormone_cluster_weighted.csv', sep='\t', index_col=0)
    
    pan_tfhorm.index = pan_tfhorm.index.map(lambda x: str(x) + '_pancreas')
    hio_tfhorm.index = hio_tfhorm.index.map(lambda x: str(x) + '_HIOE')
    
    tfhorm_df=pan_tfhorm.T.join(hio_tfhorm.T, how='outer')
    tfhorm_df.fillna(0, inplace=True)
    
    fig,ax = plt.subplots(figsize=(18,18))
    fig = sns.heatmap(tfhorm_df, square=True, cmap='bwr', center=0, yticklabels=1, xticklabels=1, 
                     linewidths=.05, linecolor='k',cbar=True)
    plt.savefig(f'../pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_matrix_weighted.png',bbox_inches='tight', dpi=300)
    plt.savefig(f'../HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_matrix_weighted.png',bbox_inches='tight', dpi=300)
    
    fig = sns.clustermap(tfhorm_df.T, metric="euclidean", method='ward', square=True, cmap='bwr', center=0, yticklabels=1, xticklabels=1, robust=True,
                        linewidths=.05, linecolor='k', figsize=(20,7))
    fig.savefig(f'../pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_cluster_weighted.png',bbox_inches='tight', dpi=300)
    fig.savefig(f'../HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_cluster_weighted.png',bbox_inches='tight', dpi=300)

# tf_hormone shortest path contribution
for ft in ['all','de']:
    pan_tfhorm=pd.read_csv(f'../pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/tf_hormone_shortest_path_positive_contribution_matrix.csv', sep='\t', index_col=0)
    hio_tfhorm=pd.read_csv(f'../HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/tf_hormone_shortest_path_positive_contribution_matrix.csv', sep='\t', index_col=0)

    pan_tfhorm.columns = pan_tfhorm.columns.map(lambda x: str(x) + '_pancreas')
    hio_tfhorm.columns = hio_tfhorm.columns.map(lambda x: str(x) + '_HIOE')
    tfhorm_df=pan_tfhorm.join(hio_tfhorm, how='outer')
    tfhorm_df.fillna(0, inplace=True)

    fig = sns.clustermap(tfhorm_df.T, metric="euclidean", method='ward', square=True, cmap='YlOrRd', yticklabels=1, xticklabels=1, robust=True, z_score=0,
                        linewidths=.05, linecolor='k', figsize=(18,7))
    fig.savefig(f'../pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_shortest_path_positive_contribution_matrix_zscore.png',bbox_inches='tight', dpi=300)
    fig.savefig(f'../HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/Hormone hierarchy/comparison_tf_hormone_shortest_path_positive_contribution_matrix_zscore.png',bbox_inches='tight', dpi=300)

   
#============================================
# Diffusion modeling
#============================================
# statistic modeling
# Model selection
model = ep.ProfileThresholdModel(static_graph)
config = mc.Configuration()
config.add_model_parameter('blocked', 0)
config.add_model_parameter('adopter_rate', 0)
#config.add_model_parameter('fraction_infected', 0.1)
infected_nodes = ['RFX6']
config.add_model_initial_configuration('Infected', infected_nodes)

# Setting nodes parameters
threshold = 0.1
profile = 0.1
for i in static_graph.nodes():
    config.add_node_configuration("threshold", i, threshold)
    config.add_node_configuration("profile", i, profile)

model.set_initial_status(config)

# Simulation
iterations = model.iteration_bunch(5)
trends = model.build_trends(iterations)
status_change_node = {}
for i in [1,-1,0]:
    if i:
        status_change_node[i] = [k for ind in iterations[1:] for k,v in ind['status'].items() if v==i] + [k for k in iterations[0]['status'] if  iterations[0]['status'][k]==1 and i==1]
    else:
        status_change_node[i] = [k for k in iterations[0]['status'] if k not in (status_change_node[1]+status_change_node[-1])]

# Visualization
viz = DiffusionTrend(model, trends)
p = viz.plot(width=400, height=400)
show(p)



# dynamic modeling
# Model selection
model = dm.DynProfileThresholdModel(dynamic_graph)
config = mc.Configuration()
config.add_model_parameter('blocked', 0)
config.add_model_parameter('adopter_rate', 0)
#config.add_model_parameter('fraction_infected', 0.1)
infected_nodes = ['RFX6']
config.add_model_initial_configuration('Infected', infected_nodes)

# Setting nodes parameters
threshold = 0
profile = 0
for i in dynamic_graph.nodes():
    config.add_node_configuration("threshold", i, threshold)
    config.add_node_configuration("profile", i, profile)

model.set_initial_status(config)

# Simulate snapshot based execution
iterations = model.execute_snapshots()
status_change_node = {}
for i in [1,-1,0]:
    if i:
        status_change_node[i] = [k for ind in iterations[1:] for k,v in ind['status'].items() if v==i] + [k for k in iterations[0]['status'] if  iterations[0]['status'][k]==1 and i==1]
    else:
        status_change_node[i] = [k for k in iterations[0]['status'] if k not in (status_change_node[1]+status_change_node[-1])]



#============================================
# hierarchy position
#============================================
def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, leaf_vs_root_factor = 0):

    '''
    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    Based on Joel's answer at https://stackoverflow.com/a/29597209/2966723,
    but with some modifications.  

    We include this because it may be useful for plotting transmission trees,
    and there is currently no networkx equivalent (though it may be coming soon).
    
    There are two basic approaches we think of to allocate the horizontal 
    location of a node.  
    
    - Top down: we allocate horizontal space to a node.  Then its ``k`` 
      descendants split up that horizontal space equally.  This tends to result
      in overlapping nodes when some have many descendants.
    - Bottom up: we allocate horizontal space to each leaf node.  A node at a 
      higher level gets the entire space allocated to its descendant leaves.
      Based on this, leaf nodes at higher levels get the same space as leaf
      nodes very deep in the tree.  
      
    We use use both of these approaches simultaneously with ``leaf_vs_root_factor`` 
    determining how much of the horizontal space is based on the bottom up 
    or top down approaches.  ``0`` gives pure bottom up, while 1 gives pure top
    down.   
    
    
    :Arguments: 
    
    **G** the graph (must be a tree)

    **root** the root node of the tree 
    - if the tree is directed and this is not given, the root will be found and used
    - if the tree is directed and this is given, then the positions will be 
      just for the descendants of this node.
    - if the tree is undirected and not given, then a random choice will be used.

    **width** horizontal space allocated for this branch - avoids overlap with other branches

    **vert_gap** gap between levels of hierarchy

    **vert_loc** vertical location of root
    
    **leaf_vs_root_factor**

    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, leftmost, width, leafdx = 0.2, vert_gap = 0.2, vert_loc = 0, 
                    xcenter = 0.5, rootpos = None, 
                    leafpos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if rootpos is None:
            rootpos = {root:(xcenter,vert_loc)}
        else:
            rootpos[root] = (xcenter, vert_loc)
        if leafpos is None:
            leafpos = {}
        children = list(G.neighbors(root))
        leaf_count = 0
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            rootdx = width/len(children)
            nextx = xcenter - width/2 - rootdx/2
            for child in children:
                nextx += rootdx
                rootpos, leafpos, newleaves = _hierarchy_pos(G,child, leftmost+leaf_count*leafdx, 
                                    width=rootdx, leafdx=leafdx,
                                    vert_gap = vert_gap, vert_loc = vert_loc-vert_gap, 
                                    xcenter=nextx, rootpos=rootpos, leafpos=leafpos, parent = root)
                leaf_count += newleaves

            leftmostchild = min((x for x,y in [leafpos[child] for child in children]))
            rightmostchild = max((x for x,y in [leafpos[child] for child in children]))
            leafpos[root] = ((leftmostchild+rightmostchild)/2, vert_loc)
        else:
            leaf_count = 1
            leafpos[root]  = (leftmost, vert_loc)
#        pos[root] = (leftmost + (leaf_count-1)*dx/2., vert_loc)
#        print(leaf_count)
        return rootpos, leafpos, leaf_count

    xcenter = width/2.
    if isinstance(G, nx.DiGraph):
        leafcount = len([node for node in nx.descendants(G, root) if G.out_degree(node)==0])
    elif isinstance(G, nx.Graph):
        leafcount = len([node for node in nx.node_connected_component(G, root) if G.degree(node)==1 and node != root])
    rootpos, leafpos, leaf_count = _hierarchy_pos(G, root, 0, width, 
                                                    leafdx=width*1./max(leafcount,1), 
                                                    vert_gap=vert_gap, 
                                                    vert_loc = vert_loc, 
                                                    xcenter = xcenter)
    pos = {}
    for node in rootpos:
        pos[node] = (leaf_vs_root_factor*leafpos[node][0] + (1-leaf_vs_root_factor)*rootpos[node][0], leafpos[node][1]) 
#    pos = {node:(leaf_vs_root_factor*x1+(1-leaf_vs_root_factor)*x2, y1) for ((x1,y1), (x2,y2)) in (leafpos[node], rootpos[node]) for node in rootpos}
    xmax = max(x for x,y in pos.values())
    for node in pos:
        pos[node]= (pos[node][0]*width/xmax, pos[node][1])
    return pos



#============================================
# Perturbation simulation on network
#============================================
# initial delta on vsd
rna_vsd = pd.read_table(f'inputs/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt',sep='\t', index_col=0)
for i in range(1,4+1):
    rna_vsd[f'0_{i}']=rna_vsd.iloc[:,(i*6-6):(i*6-3)].mean(axis=1)
    rna_vsd[i]=rna_vsd.iloc[:,(i*6-3):(i*6)].mean(axis=1)
rna_vsd100_mean =  rna_vsd.iloc[:,-7::2]
rna_vsd0_mean =  rna_vsd.iloc[:,-8::2]
rna_vsd0_mean.columns=[1,2,3,4]
df_sub=rna_vsd0_mean.subtract(rna_vsd100_mean, axis=1)

# genes for prediction display
hormone_tf=sorted(list({'NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB','MAFA','RFX6','GATA4','GATA6','NEUROD1','NEUROD2','PAX4','PAX6',
                        'INS','GCG','SST','GHRL','GIP','TPH1','TAC1','CCK','PYY','SCT','NTS','GAST'}))
gene_set={'NEUROG3','PDX1','NKX2-2','NKX6-1','FOXA1','FOXA2','MAFB','MAFA','RFX6','GATA4','GATA6','NEUROD1','NEUROD2','PAX4','PAX6',
          'INS','GCG','SST','GHRL','GIP','TPH1','TAC1','CCK','PYY','SCT','NTS','GAST'}|{x.strip() for x in open(f'inputs/geneSets/Gehart2019_eec.txt')}
gene_set=sorted(list(gene_set))

ft='all'
pred_path=f'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx/{ft}/GRN_prediction'
if not os.path.exists(pred_path):
    os.makedirs(pred_path)
# combined dynamic graph
# grid search for threshold iteration
G_s=graph_dict[ft]
static_eval={}
for tfa in [True, False]:
    for threshold in range(1,21):
        for iteration in range(3,21):
            # combined static graph
            sim_df=perturbation_simulation (G_s, exp_df=rna_vsd100_mean, basis_cond=1, target_ini=[target_ini], exp_type='consistent', perturb_initial='min', 
                                            tfa=tfa, prior=prior_df,
                                            threshold=threshold, iteration=iteration)
            df_st_sub=df_sub[df_sub.index.isin(sim_df.index)].iloc[:,-1].sort_index()
            # ROC_AUC evaluation
            y_true=np.array([1 if x >=0 else 0 for x in df_st_sub])
            y_pred=np.array([1 if x >=0 else 0 for x in sim_df.iloc[:,-1]])
            y_pred_score=np.nan_to_num(sim_df.iloc[:,-1])
            conf_df=pd.DataFrame(confusion_matrix(y_true, y_pred), columns=['Predicted down', "Predicted up"], index=['Actual down', 'Actual up'])
            print(threshold, iteration, conf_df)
            accuray=sum(np.diag(conf_df))/conf_df.sum().sum()
            static_eval[(threshold,iteration)]=accuray
    
    # roc_auc plot
    fig,ax=plt.subplots(figsize=(8,8))
    eval_scatter=ax.scatter(x=[k[0] for k in static_eval], y=[v for k,v in static_eval.items()],c=[k[1] for k in static_eval], s=12, marker='o', cmap='Set1')
    ax.set_ylabel('Accuracy')
    ax.set_xlabel('threshold-tanh')
    plt.legend(handles=eval_scatter.legend_elements()[0], labels=[k[1] for k in static_eval],loc='center right', bbox_to_anchor=(1.15, 0.5))
    ax.set_title(f'iteration-t_threshold-tanh accuracy NEUROG3 KO_static')
    plt.savefig(f'{pred_path}/iteration-t_threshold-tanh accuracy NEUROG3 KO_static.png',bbox_inches='tight')
    # modeling by best hyperparameter
    threshold, iteration=max([(k,v) for k,v in static_eval.items()], key=lambda x: x[1])[0] #8,6
    sim_df=perturbation_simulation (G_s, exp_df=rna_vsd100_mean, basis_cond=1, target_ini=[target_ini], exp_type='consistent', perturb_initial='min', 
                                    tfa=tfa, prior=prior_df,
                                    threshold=threshold, iteration=iteration)
    # heatmap comparison
    for p in ['all', 'subset', 'hormone_TF']:
        if p=='all':
            # heatmap compare all genes
            heat_df=sim_df.iloc[:,-1].to_frame().join(df_st_sub,how='outer').dropna(axis=0)
        elif p=='subset':
            # heatmap compare subset genes
            heat_df=sim_df[sim_df.index.isin(gene_set)].iloc[:,-1].to_frame().join(df_st_sub[df_st_sub.index.isin(gene_set)],how='outer')
        elif p=='hormone_TF':
            heat_df=sim_df[sim_df.index.isin(hormone_tf)].iloc[:,-1].to_frame().join(df_st_sub[df_st_sub.index.isin(hormone_tf)],how='outer')
            
        heat_df.columns=[f'Predicted_thre{threshold}_{iteration}', 'vsd_sub_last']
        # quantile normaliation
        rank_mean = heat_df.stack().groupby(heat_df.rank(method='first').stack().astype(int)).mean()
        heat_df_qn = heat_df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        for c in heat_df_qn.columns:
            heat_df_qn[c]=[heat_df_qn.loc[k,c] if heat_df_qn.loc[k,c]*heat_df.loc[k,c]>=0 else -heat_df_qn.loc[k,c] for k in heat_df.index]
        heat_st=sns.clustermap(heat_df_qn, cmap='bwr', center=0, col_cluster=False)
        heat_st.fig.suptitle(f'HIOE heatmap pre_vsed NEUROG3 KO_static {p}')
        heat_st.savefig(f'{pred_path}/heatmap pre_vsed NEUROG3 KO_static {p}.png')


# dynamic graph
# grid search for threshold iteration
dynamic_eval={}
for tfa in [True, False]:
    for threshold in range(1,21):
        for iteration in range(3,21):
            for i,g in enumerate(graph_ts_dict[ft]):
                G_s=g[1]
                time=g[0]
                target=['NEUROG3']
                threshold, iteration=threshold, iteration
                print(time)
                if i==0:
                    sim_ts_df=perturbation_simulation (G_s, exp_df=rna_vsd100_mean, basis_cond=i+1, target_ini=target, exp_type='consistent', perturb_initial='min', 
                                                       tfa=tfa,
                                                       threshold=threshold, iteration=iteration)
                    target_ini_df=sim_ts_df[sim_ts_df.iloc[:,-1]!=0].iloc[:,-1]
                    target_perturb_initial=sim_ts_df[sim_ts_df.index.isin(target)].iloc[:,-1].values
                    
                else:
                    sim_ts_df=perturbation_simulation(G_s, exp_df=sim_ts_df.iloc[:,0:4], basis_cond=i+1, target_ini=target+[x for x in target_ini_df.index if x in G_s.nodes], 
                                                      tfa=tfa,
                                                      perturb_initial=np.append(target_perturb_initial, [target_ini_df[x] for x in target_ini_df.index if x in G_s.nodes]), 
                                                      exp_type=('consistent',target), dyn='ts', threshold=threshold, iteration=iteration)
                    target_ini_df=sim_ts_df[sim_ts_df.iloc[:,-1]!=0].iloc[:,-1]
                # ROC_AUC evaluation
                y_true=np.array([1 if x >=0 else 0 for x in df_st_sub])
                y_pred=np.array([1 if x >=0 else 0 for x in sim_df.iloc[:,-1]])
                y_pred_score=np.nan_to_num(sim_ts_df.iloc[:,-1])
                conf_df=pd.DataFrame(confusion_matrix(y_true, y_pred), columns=['Predicted down', "Predicted up"], index=['Actual down', 'Actual up'])
                print(threshold, iteration, conf_df)
                accuray=sum(np.diag(conf_df))/conf_df.sum().sum()
                dynamic_eval[(threshold,iteration)]=accuray
            
    # roc_auc plot
    fig,ax=plt.subplots(figsize=(8,8))
    eval_scatter=ax.scatter(x=[k[0] for k in dynamic_eval], y=[v for k,v in dynamic_eval.items()],c=[k[1] for k in dynamic_eval], s=12, marker='o', cmap='Set1')
    ax.set_ylabel('Accuracy')
    ax.set_xlabel('threshold-tanh')
    plt.legend(handles=eval_scatter.legend_elements()[0], labels=[k[1] for k in dynamic_eval],loc='center right', bbox_to_anchor=(1.15, 0.5))
    ax.set_title(f'iteration-t_threshold-tanh accuracy NEUROG3 KO_dynamic')
    plt.savefig(f'{pred_path}/iteration-t_threshold-tanh accuracy NEUROG3 KO_dynamic.png',bbox_inches='tight')
    # modeling by best hyperparameter     
    threshold, iteration=max([(k,v) for k,v in dynamic_eval.items()], key=lambda x: x[1])[0] #2,6
    for i,g in enumerate(graph_ts_dict[ft]):
        G_s=g[1]
        time=g[0]
        target=['NEUROG3']
        threshold, iteration=threshold, iteration
        print(time)
        if i==0:
            sim_ts_df=perturbation_simulation (G_s, exp_df=rna_vsd100_mean, basis_cond=i+1, target_ini=target, exp_type='consistent', perturb_initial='min', 
                                               tfa=tfa,
                                               threshold=threshold, iteration=iteration)
            target_ini_df=sim_ts_df[sim_ts_df.iloc[:,-1]!=0].iloc[:,-1]
            target_perturb_initial=sim_ts_df[sim_ts_df.index.isin(target)].iloc[:,-1].values
            
        else:
            sim_ts_df=perturbation_simulation(G_s, exp_df=sim_ts_df.iloc[:,0:4], basis_cond=i+1, target_ini=target+[x for x in target_ini_df.index if x in G_s.nodes], 
                                              tfa=tfa,
                                              perturb_initial=np.append(target_perturb_initial, [target_ini_df[x] for x in target_ini_df.index if x in G_s.nodes]), 
                                              exp_type=('consistent',target), threshold=threshold, dyn='ts', iteration=iteration, count=1)
            target_ini_df=sim_ts_df[sim_ts_df.iloc[:,-1]!=0].iloc[:,-1]
    # heatmap comparison
    for p in ['all', 'subset', 'hormone_TF']:
        if p=='all':
            # heatmap compare all genes
            heat_df=sim_ts_df.iloc[:,-1].to_frame().join(df_ts_sub,how='outer').dropna(axis=0)
        elif p=='subset':
            # heatmap compare subset genes
            heat_df=sim_ts_df[sim_ts_df.index.isin(gene_set)].iloc[:,-1].to_frame().join(df_ts_sub[df_ts_sub.index.isin(gene_set)],how='outer')
        elif p=='hormone_TF':
            heat_df=sim_ts_df[sim_ts_df.index.isin(hormone_tf)].iloc[:,-1].to_frame().join(df_ts_sub[df_ts_sub.index.isin(hormone_tf)],how='outer')
            
        heat_df.columns=[f'Predicted_thre{threshold}_{iteration}', 'vsd_sub_last']
        # quantile normaliation
        rank_mean = heat_df.stack().groupby(heat_df.rank(method='first').stack().astype(int)).mean()
        heat_df_qn = heat_df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        for c in heat_df_qn.columns:
            heat_df_qn[c]=[heat_df_qn.loc[k,c] if heat_df_qn.loc[k,c]*heat_df.loc[k,c]>=0 else -heat_df_qn.loc[k,c] for k in heat_df.index]
        heat_ts=sns.clustermap(heat_df_qn, cmap='bwr', center=0, col_cluster=False)
        heat_ts.fig.suptitle(f'HIOE heatmap pre_vsed NEUROG3 KO_dynamic {p}')
        heat_ts.savefig(f'{pred_path}/heatmap pre_vsed NEUROG3 KO_dynamic {p}.png')
       

        
def perturbation_simulation (G, exp_df, basis_cond, target_ini, perturb_initial='min', threshold=10, exp_type='consistent', dyn='st', 
                             tfa=False, prior='', iteration=1, count=1):
    
    print(f'\niteration {count}')
    G_=G.copy()
    target_ini=sorted(list(set(target_ini)))
    delta_= []
    if count==1:
        # vsd
        perturb_df = exp_df[exp_df.index.isin(list(G_.nodes))]
        if perturb_initial=='min':
            perturb_initial=[min(perturb_df.min())]*len(target_ini)
        elif type(perturb_initial) in [int, float]:
            perturb_initial=[perturb_initial]
        # delta expression change
        for j,k in enumerate(target_ini):
            delta_.append(2**perturb_initial[j] - 2**perturb_df.loc[k, basis_cond])
        delta_=np.array(delta_)
    else:
        perturb_df= exp_df
        #block_l=list(perturb_df[abs(perturb_df[basis_cond])>=threshold].index)
        # delta expression change
        for k in target_ini:
            delta_.append(perturb_df.loc[k, basis_cond])
        delta_=np.array(delta_)

    # target list
    target_=[]
    for k in target_ini:
        target_.extend(list(G_.successors(k)))
    target_=list(set(target_))
    # block feedback on first round
    if count==1 and dyn=='st':
        target_=[x for x in target_ if x not in target_ini]
        
    # adjacency matrix TF x Target
    adj_df_ini = nx.to_pandas_adjacency(G_, weight='SignedQuantile').sort_index()
    adj_df=adj_df_ini.copy()
    # if count!=1:
    #     for k in block_l:
    #         adj_df.T[k]=[0]*adj_df.shape[1]
       
    # perturbed TFs contribution (delta x edge_weight) to Targets
    # use activation function (tanh) to model input/output saturation
    perturb_df = perturb_df.join(adj_df.loc[adj_df.index.isin(target_ini), adj_df.columns.isin(target_)].apply(lambda x: threshold*np.tanh(np.dot(x.values, threshold*np.tanh(delta_)))).to_frame(name=f'_iter_{count}'), 
                                 how='outer').fillna(0)
    perturb_df.loc[perturb_df.index.isin([x for x in target_ini if x not in target_]),f'_iter_{count}'] = [v for k,v in dict(zip(target_ini,delta_)).items() if k not in target_]
    
    # TFA estimation for each iteration
    if tfa:
        # expression matrix 
        exp_df_=perturb_df.loc[:, [condition[0],f'_iter_{count}']]
        # prior matrix
        prior_df=prior
        # least square solution for TFA
        tfa = np.linalg.lstsq(prior_df.values, exp_df_.values, rcond=-1)[0]
        tfa_df=pd.DataFrame(data=tfa, index=prior_df.columns, columns=exp_df_.columns)
        # update expression matrix TF with TFA
        perturb_df.loc[perturb_df.index.isin(tfa_df.index),f'_iter_{count}']=tfa_df[f'_iter_{count}']        
    
    # consistent on target_ini won't get feed back onto themselves
    if count==1:
        if exp_type[0]=='consistent':
            G_.remove_nodes_from(exp_type[1])
        elif exp_type=='consistent':
            G_.remove_nodes_from(target_ini)
 
    # recursive function iterate over
    if iteration!=1:
        target_=[x for x in target_ if x in G_.nodes]
        perturb_df=perturbation_simulation (G=G_, exp_df=perturb_df, basis_cond=perturb_df.columns[-1], target_ini=target_, threshold=threshold, iteration=iteration-1, count=count+1)
    # sum over delta change at each iteration as final result
    else:
        sum_column=[x for x in perturb_df.columns.tolist() if '_iter' in str(x)]
        perturb_df[f'final_delta_iter{count}']=perturb_df[sum_column].sum(axis=1).apply(lambda x: np.log2(x) if x>0 else -np.log2(abs(x)))
    return perturb_df



#============================================
# Dynamic diiffusion iteration
#============================================
iterations, status_dict = recursive_DynProfileThresholdModel(dynamic_graph, infected=['RFX6'], recursive_round=5)

def recursive_DynProfileThresholdModel(G, blocked=0, adopter_rate=0, threshold=0, profile=0.1, infected=0.1, recursive_round=5):
    model = dm.DynProfileThresholdModel(G)
    config = mc.Configuration()
    config.add_model_parameter('blocked', blocked)
    config.add_model_parameter('adopter_rate', adopter_rate)
    if type(infected)==float:
        config.add_model_parameter('fraction_infected', infected)
    elif type(infected)==list:
        infected_nodes = infected
        config.add_model_initial_configuration('Infected', infected_nodes)
    
    # Setting nodes parameters
    threshold = threshold
    profile = profile
    for i in G.nodes():
        config.add_node_configuration("threshold", i, threshold)
        config.add_node_configuration("profile", i, profile)
    
    model.set_initial_status(config)
    
    # Simulate snapshot based execution
    iterations = model.execute_snapshots()
    status_change_node = {}
    for i in [1,-1,0]:
        if i:
            status_change_node[i] = [k for ind in iterations[1:] for k,v in ind['status'].items() if v==i] + [k for k in iterations[0]['status'] if iterations[0]['status'][k]==1 and i==1]
        else:
            status_change_node[i] = [k for k in iterations[0]['status'] if k not in (status_change_node[1]+status_change_node[-1])]
    
    print(f'\n recursive round: {recursive_round}\ninfected: {status_change_node[1]}\nblocked: {status_change_node[-1]}')
    recursive_round-=1
    
    if recursive_round!=0:
        infected=status_change_node[1]
        print(f'\n{infected}')
        recursive_DynProfileThresholdModel(G, blocked=blocked, adopter_rate=adopter_rate, threshold=threshold, profile=profile, infected=infected, recursive_round=recursive_round)
    else:
        return (iterations, status_change_node)



    

    
        
