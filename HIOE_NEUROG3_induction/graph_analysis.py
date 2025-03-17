#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:44:30 2020

@author: user
"""

import os
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import collections
from itertools import product

import pickle
import community
import networkx as nx
from networkx.algorithms import bipartite
from networkx.algorithms import approximation as apxa

import dynetx as dn

import gseapy as gp
import powerlaw

os.chdir('/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction')   #set working path

tfs=[x.strip() for x in open('inputs/geneSets/human_TF_set.txt')][1:]

##################################################################################################################################################################
# Run the functions with different filters
##################################################################################################################################################################

nx_path = 'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx'
try:
    os.mkdir(nx_path)
except FileExistsError:
    pass

for cond_time in [24,48,72,96]:
    for de_only in [True,False]:
        if de_only:
            out_path_general=f'{nx_path}/allde'
        else:
            out_path_general=f'{nx_path}/de'
        for tf_only in [True,False]:
            if tf_only:
                out_path=out_path_general+'_tfs'
            else:
                out_path=out_path_general
            # set classification
            for classification in [None, 'active', 'repress']:

                if classification==None:
                    c='all'
                else:
                    c=classification

                out_path=out_path+f'/{cond_time}hpi {c}/'
                print(out_path)
                try:
                    os.makedirs(out_path)
                except FileExistsError:
                    pass
                
                print(f'\n\n*****{cond_time}hpi {de_only} {tf_only} {c}*****\n\n')
                
                G_tf=network_construct(cond_time=cond_time)
                G_tf=network_filter(G_tf,classification=classification,de_only=de_only,tf_only=tf_only)
                network_graphml_save(G_tf, path=f'{out_path}/graphml_{cond_time}_{classification}.txt')
                recip_assort_dict=reciprocity_assortivity(G_tf,nodes=None,weight='weight')
            #try:
                cliques,node_clique_max,node_clique=clique_counts(G=G_tf,top_n=50)
                triangles,tri_cluster,squre_cluster=cluster_coeff(G_tf,top_n=50,nodes=None,weight='weight',classification=classification)
                try:
                    dijkstra_path,dijkstra_path_mean=weighted_shortest_path(G_tf,weight='weight',source_target=None,classification=classification)
                    centrality_dict=graph_centrality(G=G_tf,top_n=50,function=None,weight='weight',classification=classification,save_plot=True)
                    for p in [False,True]:
                        partition=draw_network(G_tf, p=p, classification=classification, file_prefix=None, save=True)
                except nx.NetworkXError:
                    print(nx.NetworkXError,f'\n\n====={cond_time}hpi {de_only} {tf_only} {c} is not weakly connected=====\n\n')
                    for weak_c in nx.weakly_connected_components(G_tf):
                        print(weak_c)
                
                power_law_dict=plot_degree_dist(G_tf,classification=classification)
                density,w_connected_comp,s_connected_comp,max_s_connected_nodes,C,C0=directed_connectivity_comp(G_tf,top_n=50,min_s_connected_nodes=3,classification=classification)
            # except:
            #     pass


##################################################################################################################################################################
##################################################################################################################################################################
# set global variable de_only and corresponding path for function outputs
##################################################################################################################################################################
# Outputs path based on de_only and tf_only
de_only=False
tf_only=False
out_path=''


def network_construct(cond_time):
    # find DE not included in meta data from both D9 and D12
    core_trn = pd.read_table(f'outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/{cond_time}hpi_Cores/Core_prior_atac_Miraldi_q_ChIP_bias10_maxComb_fdr5_HIOE_NEUROG3_{cond_time}hpiSet_All_sp.tsv',sep='\t')
    gene_list = set(core_trn.TF.tolist() + core_trn.Target.tolist())
    
    # dict of gene vsd
    rna_vsd = pd.read_table(f'inputs/geneExpression/RNAseq_24_DESeq2_VSDcounts.txt',sep='\t', index_col=0, usecols=[0]+list(range(cond_time//24+1,cond_time//24+6+1)))
    rna_vsd['vsd0'] = rna_vsd.iloc[:,0:3].mean(axis=1)
    rna_vsd['vsd100'] = rna_vsd.iloc[:,3:6].mean(axis=1)
    vsd_dict = {k:{'0_vsd':rna_vsd[rna_vsd.index==k].vsd0.values[0], '100_vsd':rna_vsd[rna_vsd.index==k].vsd100.values[0]} if k in rna_vsd.index else 0 for k in gene_list}
    
    # dict of gene lfc 100vs0
    rna_de = pd.read_table(f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/XZ_07-21-28-08-05-2020_RNA-seq_HIOE/Counts_TPM_mat/analysis/results/WT_HIOE_{cond_time}hpi_filtered.txt',sep='\t', index_col=0)
    #rna_de = pd.read_table(f'inputs/geneExpression/results/WT_HIOE_{cond_time}hpi_filtered.txt',sep='\t', index_col=0)
    lfc_dict = {k:rna_de[rna_de.index==k].log2FoldChange.values[0] if k in rna_de.index else 0 for k in gene_list}
    padj_dict = {k:rna_de[rna_de.index==k].padj.values[0] if k in rna_de.index else 1 for k in gene_list}
    
    # construct graph from dataframe
    # analyze graph from networkx
    df_G = core_trn
    df_G.stroke = df_G.stroke.str.slice(start=4, stop=-1).str.split(',').apply(lambda x: tuple([int(i, 10)/255 for i in x]))
    df_G.weight = df_G.combPCorr.apply(abs)
    
    # dataframe to graph
    G_tf=nx.from_pandas_edgelist(df_G,source='TF', target='Target', edge_attr=True, create_using=nx.DiGraph())
    
    # nodes attribute from pandas
    nx.set_node_attributes(G_tf, vsd_dict, f'tpm_{cond_time}hpi')
    nx.set_node_attributes(G_tf, lfc_dict, f'lfc_{cond_time}hpi')
    nx.set_node_attributes(G_tf, padj_dict, f'padj_{cond_time}hpi')
    return G_tf


def network_filter(G,classification,de_only=de_only,tf_only=tf_only):
    '''
    Parameters
    ----------
    G : graph: graph for centrality
    classification : str/list: [test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud]

    Returns
    -------
    H : graph: G with edges filtered with time/classification
    '''
    #Creates a copy of the graph
    H = G.copy()
    #Checks all the edges and removes some based on time,classification and weight
    for i, j, attr in G.edges(data=True):
        if (i,j) in H.edges():
            #edge colors based on classification
            H[i][j]['color'] = attr['stroke']

        # filter edges with classification
        if classification == None:
            pass
        else:
            if (attr['SignedQuantile']>0 and classification!='active') or (attr['SignedQuantile']<0 and classification!='repress'):
                H.remove_edge(i, j)
                H.remove_nodes_from(list(nx.isolates(H)))
                           
    for i in sorted(H.nodes()):
        if de_only and ((abs(H.nodes[i][f'lfc_{cond_time}hpi']) <1) or (H.nodes[i][f'padj_{cond_time}hpi'] >= 0.05)):
            H.remove_node(i)
        if i in sorted(H.nodes()) and tf_only and not(i in tfs):
            H.remove_node(i)
    H.remove_nodes_from(list(nx.isolates(H)))
    return H


def network_graphml_save(G, path):
    nx.write_edgelist(G, path, comments='#', delimiter='\t', data=True, encoding='utf-8')
    return 'Graphml saved!'
    

def clique_counts(G,top_n):
    '''
    # clique
    # A clique is a maximal subset of the vertices in an undirected network such that every member of the set is connected by an edge to every other.
    # cliques with nodes >=3

    Parameters
    ----------
    G : graph: graph for clique
    top_n : int: top_n cliques with most nodes
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    cliques: list: ordered by lenth of clique
    node_clique_max: dict: Key: nodes, Value: max clique length
    node_clique: dict: Key: nodes, Value: associated cliques in list
    '''

    G_ud=nx.DiGraph.to_undirected(G)
    # identify cliques and sort by size only >=3
    cliques=[x for x in sorted(list(nx.enumerate_all_cliques(G_ud)),key=len,reverse=True) if len(x)>=3]
    clique_element=[item for sublist in cliques for item in sublist]
    # maximal clique size for each/selected node
    node_clique_max={k: v for k, v in sorted(nx.node_clique_number(G_ud).items(),key=lambda item: item[1], reverse=True)}
    # list of cliques for each/selected node
    node_clique={k: v for k, v in sorted(nx.cliques_containing_node(G_ud).items(),key=lambda item: len(item[1]), reverse=True)}

    clique_dict_total={k:v for k,v in sorted(collections.Counter(clique_element).items(),key=lambda item: item[1], reverse=True)}
    clique_dict_top={k:v for k,v in sorted(clique_dict_total.items(),key=lambda item: item[1], reverse=True)[:top_n]}
    clique_dict={k:[] for k in set(clique_element)}
    for i in range(3,max([v for v in node_clique_max.values()])+1):
        for j in clique_dict:
            clique_dict[j].append(np.log1p(collections.Counter([item for sublist in [k for k in cliques if len(k)==i] for item in sublist])[j]))
    clique_df=pd.DataFrame.from_dict(clique_dict,orient='index',columns=list(range(3,max([v for v in node_clique_max.values()])+1)))
    if clique_df.shape!=(0,0):
        print(clique_df.describe())
    
        plt.rcParams.update(plt.rcParamsDefault)
        fig = plt.figure(figsize=(10,5))
        fig.suptitle(f'total_cliques_counts {top_n}', fontsize=15)
        plt.bar(clique_dict_top.keys(), clique_dict_top.values(), color='b')
        plt.gca().get_xticklabels()[1].set
        plt.xticks(fontsize=7,rotation=75)
        for j in [x for x in clique_dict_top if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
            plt.gca().get_xticklabels()[list(clique_dict_top).index(j)].set_color("red")
        for j in [x for x in clique_dict_top if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
            plt.gca().get_xticklabels()[list(clique_dict_top).index(j)].set_color("blue")
        plt.show()
        fig.savefig(f'{out_path}/TF_clique_counts_total_top{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
    
        de_up_dict={key:'up' for key in [x for x in clique_dict if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]}
        de_down_dict={key:'down' for key in [x for x in clique_dict if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]}
        de_updown_dict={**de_up_dict,**de_down_dict}
        de_classification=clique_df.index.map(de_updown_dict)
        row_colors=de_classification.map(dict(zip(['up','down',np.nan],['yellowgreen','red','w'])))
        sns.set(rc={'figure.figsize':(10,10)})
        if clique_df.shape[1]==1:
            col_cluster=False
        else:
            col_cluster=True
        clique_df_heatmap=sns.clustermap(clique_df, square=True,cmap="YlGnBu",row_colors=row_colors,
                                         row_cluster=True, col_cluster=col_cluster)
        clique_df_heatmap.ax_heatmap.set_xticklabels(clique_df_heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 8)
        clique_df_heatmap.ax_heatmap.set_yticklabels(clique_df_heatmap.ax_heatmap.get_ymajorticklabels(), fontsize = 8)
        clique_df_heatmap.fig.suptitle('TF_clique counts log1p', fontsize = 20) # title with fontsize 20
        clique_df_heatmap.ax_heatmap.set_xlabel('n_cliques (>=3)')
        clique_df_heatmap.ax_heatmap.set_ylabel('TF')
        clique_df_heatmap.savefig(f'{out_path}/TF_clique_counts_log1p_{cond_time}_{classification}.png', dpi=300)
        plt.rcParams.update(plt.rcParamsDefault)
    return(cliques,node_clique_max,node_clique)


def cluster_coeff(G,top_n,nodes,weight,classification):
    '''
    # cluster coefficient
    # triangle cluster: the probability that two neighbors of node v are connected with each other
    # square cluster: the probability that two neighbors of node v share a common neighbor different from v
    # transitivity:the fraction of all possible triangles present in graph, favored over average_clustering

    Parameters
    ----------
    G : graph: graph for cluster_coeff
    top_n : int: top_n nodes with highest cluster_coeff
    nodes : list: list of nodes for cluster_coeff, None for all
    weight : str: capacity, edge_counts, hpi_edge
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    triangles: dict: Key: nodes, Value: number of triangles
    tri_cluster: dict: Key: nodes, Value: tri cluster_coeff
    squre_cluster: dict: Key: nodes, Value: squre cluster_coeff
    '''

    triangles=nx.triangles(nx.DiGraph.to_undirected(G), nodes=nodes)
    triangles_top={k:v for k,v in sorted(triangles.items(),key=lambda item: item[1], reverse=True)[:top_n]}
    fig = plt.figure(figsize=(10,5))
    fig.suptitle(f'triangle_counts {top_n}', fontsize=15)
    plt.bar(triangles_top.keys(), triangles_top.values(), color='b')
    plt.gca().get_xticklabels()[1].set
    plt.xticks(fontsize=7,rotation=75)
    for j in [x for x in triangles_top if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(triangles_top).index(j)].set_color("red")
    for j in [x for x in triangles_top if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(triangles_top).index(j)].set_color("blue")
    plt.show()
    fig.savefig(f'{out_path}/TF_triangle_counts_top{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

    tri_cluster=nx.clustering(nx.DiGraph.to_undirected(G),weight=weight)
    n=550
    tri_cluster_top={k:np.log1p(v) for k,v in sorted(tri_cluster.items(),key=lambda item: item[1], reverse=True)[:n]}
    fig = plt.figure(figsize=(20,5))
    fig.suptitle(f'cluster_coeff {n}', fontsize=15)
    plt.bar(tri_cluster_top.keys(), tri_cluster_top.values(), color='b',width=1)
    plt.gca().get_xticklabels()[1].set
    xticks=[x for x in tri_cluster_top][::5]
    plt.xticks(xticks,fontsize=7,rotation=75)
    for j in [x for x in xticks if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[xticks.index(j)].set_color("red")
    for j in [x for x in xticks if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[xticks.index(j)].set_color("blue")
    plt.show()
    fig.savefig(f'{out_path}/TF_cluster_coeff_top{n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

    tri_cluster_de={k:np.log1p(v) for k,v in sorted(tri_cluster.items(),key=lambda item: item[1], reverse=True)[400:550]}
    fig = plt.figure(figsize=(20,5))
    fig.suptitle(f'cluster_coeff DE', fontsize=15)
    plt.bar(tri_cluster_de.keys(), tri_cluster_de.values(), color='b',width=1)
    plt.gca().get_xticklabels()[1].set
    xticks=[x for x in tri_cluster_de if G.nodes[x][f'lfc_{cond_time}hpi']!=0]
    plt.xticks(xticks,fontsize=7,rotation=75)
    for j in [x for x in xticks if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[xticks.index(j)].set_color("red")
    for j in [x for x in xticks if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[xticks.index(j)].set_color("blue")
    plt.show()
    fig.savefig(f'{out_path}/TF_cluster_coeff_top{n}_DE_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

    squre_cluster=nx.square_clustering(G, nodes=nodes)
    squre_cluster_top={k:np.log1p(v) for k,v in sorted(squre_cluster.items(),key=lambda item: item[1], reverse=True)[:top_n]}
    fig = plt.figure(figsize=(10,5))
    fig.suptitle(f'squre_cluster_coeff {top_n}', fontsize=15)
    plt.bar(squre_cluster_top.keys(), squre_cluster_top.values(), color='b')
    plt.ylabel('log1p(squre_cluster_coeff)')
    plt.gca().get_xticklabels()[1].set
    plt.xticks(fontsize=7.5,rotation=75)
    for j in [x for x in squre_cluster_top if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(squre_cluster_top).index(j)].set_color("red")
    for j in [x for x in squre_cluster_top if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(squre_cluster_top).index(j)].set_color("blue")
    plt.show()
    fig.savefig(f'{out_path}/TF_squre_cluster_coeff_top{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

    transitivity=nx.transitivity(G)
    avg_cluster=nx.average_clustering(nx.DiGraph.to_undirected(G), nodes=nodes)
    print(f'transitivity:{transitivity}\navg_cluster_coeff:{avg_cluster}')
    return (triangles,tri_cluster,squre_cluster)


def reciprocity_assortivity(G,nodes,weight):
    '''
    # reciprocity
    # the ratio of the number of edges pointing in both directions to the total number of edges in the graph
    # assortativity
    # measures the similarity of connections in the graph with respect to the node degree/attribute(class/scalar)

    Parameters
    ----------
    G : graph: graph for reciprocity and assortivity
    nodes : list: list of nodes for reciprocity/assortivity, None for all
    weight : str: capacity, hpi_edge
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    reciprocity : int/dict: Key: nodes, Value: reciprocity
    degree_assort : int/dict: Key: nodes, Value: assortivity
    '''

    reciprocity=nx.reciprocity(G, nodes=nodes)
    degree_assort=nx.degree_assortativity_coefficient(G, x='out', y='in', weight=weight, nodes=nodes)
    print(f'reciprocity:{reciprocity}\ndegree_assortivity_coeff:{degree_assort}\n')
    recip_assort_dict={'reciprocity':reciprocity,'degree':degree_assort}

    for i in G.nodes(data=True)[list(G.nodes())[0]]:
        if 'tpm' in i:
            pass
        else: 
            attr_assort=nx.attribute_assortativity_coefficient(G, attribute=i, nodes=nodes)
            print(f'Attribute assortativity for {i}: {attr_assort}')
            recip_assort_dict[i]=attr_assort
    pd.DataFrame(recip_assort_dict.items(),columns=['reciprocity_assortativity', 'value']).to_csv(f'{out_path}/reciprocity_assortivity_{cond_time}_{classification}.txt',sep='\t')
    return recip_assort_dict


def weighted_shortest_path(G,weight,source_target,classification):
    '''
    # weighted shortest path

    Parameters
    ----------
    G : graph: graph for reciprocity and assortivity
    weight : str: capacity, hpi_edge
    source_target: tuple: (source,target), None for all
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    dijkstra_path : dict: Key: nodes, Value: dijkstra_path list
    dijkstra_path_mean : int: average dijkstra_path length
    '''

    dijkstra_path={k:v for k,v in dict(nx.all_pairs_dijkstra(G, cutoff=None, weight=weight)).items() if len(v[0])>1}
    dijkstra_key=list(dijkstra_path.keys())
    dijkstra_target=list(set(k1 for d in [v[0] for k,v in dijkstra_path.items()] for k1,v1 in d.items()))
    dijkstra_index={}
    for e in dijkstra_target:
        index_l=[0]*len(dijkstra_key)
        for k,v in dijkstra_path.items():
            for k1,v1 in v[0].items():
                if k1==e:
                    index_l[dijkstra_key.index(k)]+=v1
        dijkstra_index[e]=index_l
    dijkstra_df=pd.DataFrame.from_dict(dijkstra_index,orient='index',columns=dijkstra_key)
    print(dijkstra_df.describe())

    de_up_dict={key:'up' for key in [x for x in set(dijkstra_index.keys())| set(dijkstra_path.keys()) if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]}
    de_down_dict={key:'down' for key in [x for x in set(dijkstra_index.keys())| set(dijkstra_path.keys()) if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]}
    de_updown_dict={**de_up_dict,**de_down_dict}
    row_classification,col_classification=dijkstra_df.index.map(de_updown_dict),dijkstra_df.columns.map(de_updown_dict)
    row_colors=row_classification.map(dict(zip(['up','down',np.nan],['yellowgreen','red','w'])))
    col_colors=col_classification.map(dict(zip(['up','down',np.nan],['yellowgreen','red','w'])))

    sns.set(rc={'figure.figsize':(10,10)})
    dijkstra_df_heatmap=sns.clustermap(dijkstra_df, square=True,cmap="YlGnBu",row_colors=row_colors,col_colors=col_colors,
                                       row_cluster=True, col_cluster=True)
    dijkstra_df_heatmap.ax_heatmap.set_xticklabels(dijkstra_df_heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
    dijkstra_df_heatmap.ax_heatmap.set_yticklabels(dijkstra_df_heatmap.ax_heatmap.get_ymajorticklabels(), fontsize = 8)
    dijkstra_df_heatmap.fig.suptitle('dijkstra path elements (include target)', fontsize = 20) # title with fontsize 20
    dijkstra_df_heatmap.ax_heatmap.set_xlabel('TF')
    dijkstra_df_heatmap.ax_heatmap.set_ylabel('Targets')
    dijkstra_df_heatmap.savefig(f'{out_path}/TF_dijkstra_path_elements_{cond_time}_{classification}.png', dpi=300)

    plt.rcParams.update(plt.rcParamsDefault)
    dijkstra_path_max={k:max([j for i,j in v[0].items()]) for k,v in dijkstra_path.items()}
    dijkstra_path_max_top={k:v for k,v in sorted(dijkstra_path_max.items(), key=lambda item:item[1],reverse=True)}
    fig = plt.figure(figsize=(20,5))
    fig.suptitle(f'dijkstra path max (include target)', fontsize=15)
    plt.bar(dijkstra_path_max_top.keys(), dijkstra_path_max_top.values(), color='b')
    plt.gca().get_xticklabels()[1].set
    plt.xticks(fontsize=5,rotation=75)
    for j in [x for x in dijkstra_path_max_top if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(dijkstra_path_max_top).index(j)].set_color("red")
    for j in [x for x in dijkstra_path_max_top if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
        plt.gca().get_xticklabels()[list(dijkstra_path_max_top).index(j)].set_color("blue")
    plt.show()
    fig.savefig(f'{out_path}/TF_dijkstra_path_max_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

    dijkstra_path_mean=nx.average_shortest_path_length(G, weight=weight)

    if source_target!=None:
        s_t_path=dijkstra_path[source_target[0]][1][source_target[1]]
        print(f'{source_target[0]} -> {source_target[1]}: shortest_path_length:{len(s_t_path)}\n{s_t_path}')
    return (dijkstra_path,dijkstra_path_mean)


def graph_centrality(G,top_n,function,weight,classification,save_plot):
    '''
    Parameters
    ----------
    G : graph: graph for centrality
    top_n : int: top_n genes for comparison
    function : function/list of functions: centrality methods
    weight : str: capacity, hpi_edge
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud
    save_plot : bool: True to save the plots

    Returns
    -------
    centrality_dict : dict: Keys as centrality method
    '''

    # networkx for vertex centrality
    if function==None:
        functions=[nx.degree_centrality,nx.in_degree_centrality,nx.out_degree_centrality,nx.closeness_centrality,
                   nx.harmonic_centrality,nx.information_centrality,nx.betweenness_centrality,nx.edge_betweenness_centrality,
                   nx.current_flow_betweenness_centrality,nx.edge_current_flow_betweenness_centrality,
                   nx.eigenvector_centrality,nx.katz_centrality,nx.pagerank,nx.hits]
    else:
        if type(function)==list:
            functions=function
        else:
            functions=[function]
    centrality_dict={}
    for i in functions:
        
        if i in [nx.degree_centrality,nx.in_degree_centrality,nx.out_degree_centrality,nx.closeness_centrality]:
            centrality=i(G)
        elif i in [nx.betweenness_centrality,nx.edge_betweenness_centrality]:
            centrality=i(G,weight=weight)
        elif i in [nx.information_centrality,nx.current_flow_betweenness_centrality,nx.edge_current_flow_betweenness_centrality]:
            centrality=i(nx.DiGraph.to_undirected(G),weight=weight)
        elif i==nx.harmonic_centrality:
            centrality=i(G,distance=weight)
        elif i==nx.hits:
            centrality=i(G,max_iter=50000)
        else:
            try:
                centrality=i(G,max_iter=50000, weight=weight)
            except:
                if i==nx.katz_centrality:
                    centrality=nx.katz_centrality_numpy(G, weight=weight)
                elif i==nx.pagerank:
                    centrality=nx.pagerank_numpy(G, weight=weight)
            
        plt.rcParams.update(plt.rcParamsDefault)
        if i==nx.hits:
            centrality_hits={}
            centrality_hits['hubs']={k: v for k, v in sorted(centrality[0].items(),key=lambda item: item[1], reverse=True)[:top_n]}
            centrality_hits['authorities']={k: v for k, v in sorted(centrality[1].items(),key=lambda item: item[1], reverse=True)[:top_n]}
            centrality=centrality_hits.copy()

            fig = plt.figure(figsize=(10,5))
            fig.suptitle(f'{str(i).split(" ")[1]}_hubs top_{top_n}', fontsize=15)
            plt.bar(centrality['hubs'].keys(), centrality['hubs'].values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=7.5,rotation=75)
            for j in [x for x in centrality['hubs'] if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['hubs']).index(j)].set_color("red")
            for j in [x for x in centrality['hubs'] if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['hubs']).index(j)].set_color("blue")
            plt.show()
            print(f'{str(i).split(" ")[1]}_hubs of {str(G)}:\n {centrality["hubs"]}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{str(i).split(" ")[1]}_hubs_centrality top_{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)

            fig = plt.figure(figsize=(10,5))
            fig.suptitle(f'{str(i).split(" ")[1]}_authorities top_{top_n}', fontsize=15)
            plt.bar(centrality['authorities'].keys(), centrality['authorities'].values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=7.5,rotation=75)
            for j in [x for x in centrality['authorities'] if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['authorities']).index(j)].set_color("red")
            for j in [x for x in centrality['authorities'] if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['authorities']).index(j)].set_color("blue")
            plt.show()
            print(f'{str(i).split(" ")[1]}_authorities of {str(G)}:\n {centrality["authorities"]}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{str(i).split(" ")[1]}_authorities_centrality top_{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
        elif i in [nx.edge_betweenness_centrality,nx.edge_current_flow_betweenness_centrality]:
            centrality={k: v for k, v in sorted(centrality.items(),key=lambda item: item[1], reverse=True)[:top_n]}    #str(k) for tuple from edge_betweenness_centrality
            fig = plt.figure(figsize=(10,5))
            fig.suptitle(f'{str(i).split(" ")[1]} top_{top_n}', fontsize=15)
            plt.bar([str(x[0])+'-'+str(x[1]) for x in centrality.keys()], centrality.values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=5,rotation=75)
            plt.show()
            print(f'{str(i).split(" ")[1]} of {str(G)}:\n {centrality}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{str(i).split(" ")[1]}_centrality top_{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
        else:
            centrality={k: v for k, v in sorted(centrality.items(),key=lambda item: item[1], reverse=True)[:top_n]}    #str(k) for tuple from edge_betweenness_centrality
            fig = plt.figure(figsize=(10,5))
            fig.suptitle(f'{str(i).split(" ")[1]} top_{top_n}', fontsize=15)
            plt.bar(centrality.keys(), centrality.values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=7.5,rotation=75)
            for j in [x for x in centrality if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality).index(j)].set_color("red")
            for j in [x for x in centrality if G.nodes[x][f'lfc_{cond_time}hpi']<=-1]:
                plt.gca().get_xticklabels()[list(centrality).index(j)].set_color("blue")
            plt.show()
            print(f'{str(i).split(" ")[1]} of {str(G)}:\n {centrality}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{str(i).split(" ")[1]}_centrality top_{top_n}_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
        centrality_dict[f'{str(i).split(" ")[1]}']=centrality
    return centrality_dict


def draw_network(G, p, classification, file_prefix, save):
    '''
    # default to label nodes with HITS

    Parameters
    ----------
    G : graph: graph for centrality
    p: bool: whether to compute partition
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud
    file_prefix : int: prefix for saved file
    save : bool: True to save the plots

    Returns
    -------
    partition : dict: Key as nodes, Values as partition community
    '''
    H=G.copy()
    #crates a list for edges, weights and colors
    selected_edges=H.edges()
    selected_edges_weight=[H[u][v]['stroke.width'] for u,v in selected_edges]
    edge_colors = [H[u][v]['color'] for u,v in selected_edges]
    #creates list of nodes and a list their degrees that will be used later for their sizes
    nodelist = H.nodes()
    node_sizes =[v for k, v in nx.degree(H)]
    # specific lable based on HITS
    label,label_both,label_hub,label_auth = {},{},{},{}
    centrality_dict=graph_centrality(G,top_n=50,function=nx.hits,weight='stroke.width',classification=classification,save_plot=False)
    for i in H.nodes():
        if (i in centrality_dict['hits']['hubs']) and (i in centrality_dict['hits']['authorities']):
            #set the node name as the key and the label as its value
            label_both[i] = i
        else:
            if i in centrality_dict['hits']['hubs']:
                label_hub[i]=i
            elif i in centrality_dict['hits']['authorities']:
                label_auth[i]=i
            else:
                label[i]=i

    # specific node shape based on DE by spliting node group
    up_nodes = [n for n,lfc in nx.get_node_attributes(G_tf,f'lfc_{cond_time}hpi').items() if lfc>=1]
    down_nodes = [n for n,lfc in nx.get_node_attributes(G_tf,f'lfc_{cond_time}hpi').items() if lfc<=-1]
    nan_nodes = [n for n,lfc in nx.get_node_attributes(G_tf,f'lfc_{cond_time}hpi').items() if -1<lfc<1]

    if p==False:
        partition_community={'Not computed'}
        # k controls the distance between the nodes and varies between 0 and 1
        # iterations is the number of times simulated annealing is run
        # default k =0.1 and iterations=50
        positions = nx.spring_layout(H,weight='stroke.width',k =0.75)
        #positions = nx.circular_layout(H)
        #positions = nx.spectral_layout(H)
        node_color='cyan'
    else:
        # greedy method based louvian alg for modularity maximization and community detection
        partition=community.best_partition(nx.DiGraph.to_undirected(H),weight='stroke.width',random_state=42)
        positions=community_layout(H, partition)
        partition_community,enr_df=gsea_partition_enrich(H,partition=partition,save_plot=True)
        node_color=list(partition.values())
        if file_prefix==None:
            file_prefix='partitioned'
        else:
            file_prefix='partitioned_'+file_prefix
        pd.Series(partition,name='genes').to_csv(f'{out_path}/partition_nodes_{cond_time}_{classification}.txt',sep='\t')
        pd.Series(partition_community,name='genes').to_csv(f'{out_path}/partition_comunity_{cond_time}_{classification}.txt',sep='\t')
    plt.figure(figsize=(20,20))
    #draws nodes
    nx.draw_networkx_nodes(H,positions,node_color=node_color,nodelist=nodelist,
                           #the node size will be now based on its degree
                           node_size=node_sizes,alpha=0.8)
    #Styling for labels
    #nx.draw_networkx_labels(H,positions,label,font_size=5,font_color='k')
    nx.draw_networkx_labels(H,positions,label_both,font_size=7,font_color='red')
    nx.draw_networkx_labels(H,positions,label_hub,font_size=7,font_color='magenta')
    nx.draw_networkx_labels(H,positions,label_auth,font_size=7,font_color='navy')
    #draws the edges
    nx.draw_networkx_edges(H, positions, nodelist=up_nodes, edge_list=selected_edges,style='solid', width=[x/5 for x in selected_edges_weight], edge_color=edge_colors,node_color='lightblue')
    nx.draw_networkx_edges(H, positions, nodelist=down_nodes, edge_list=selected_edges,style='solid', width=[x/5 for x in selected_edges_weight], edge_color=edge_colors,node_color='plum', node_shape='d')
    nx.draw_networkx_edges(H, positions, nodelist=nan_nodes, edge_list=selected_edges,style='solid', width=[x/5 for x in selected_edges_weight], edge_color=edge_colors,node_color='darkgrey', node_shape='s')

    if save==True:
        if file_prefix==None:
            plt.title('HITS_spring_layout')
            plt.savefig(f'{out_path}/HITS_spring_layout_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
        else:
            plt.title(f'{file_prefix}_HITS_spring_layout')
            plt.savefig(f'{out_path}/{file_prefix}_HITS_spring_layout_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
    return partition_community


def directed_connectivity_comp(G,top_n,min_s_connected_nodes,classification):
    '''
    # weakly connected components: a subgraph all vertices are connected to each other by some path regardless of direction
    # stronlgy connected components: a subgraph all vertices are reachable to each other by directed path

    Parameters
    ----------
    G : graph: graph for centrality
    top_n : int: strongly connected components with top length
    min_s_connected_node : int: minimum nodes in a strongly connected components
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    density: int: density of graph
    w_connected_comp: int: number of weakly connected components
    s_connected_comp: int: number of strongly connected components
    max_s_connected_nodes: list: list of tuples containing strongly connected components
    C: graph: condensed graph by mapping all strongly connected nodes into one node
    C0: graph: graph of top strongly connected nodes
    '''

    # density of graph
    density=nx.density(G)
    # check weakly_connected / disconnected components
    nx.is_weakly_connected(G)
    w_connected_comp=nx.number_weakly_connected_components(G)
    # check strongly_connected components
    nx.is_strongly_connected(G)
    s_connected_comp=nx.number_strongly_connected_components(G)
    top_s_connected_len=[len(x) for x in sorted(nx.strongly_connected_components(G_tf),key=len,reverse=True)[:top_n]]
    if [x for x in top_s_connected_len if x>=min_s_connected_nodes]==[]:
        return f'Max connection below {min_s_connected_nodes}',None,None,None,None,None
    else:
        min_s_connected_index=top_s_connected_len.index([x for x in top_s_connected_len if x>=min_s_connected_nodes][-1])
        C=nx.condensation(G)
        max_s_connected_nodes=[(x,y) for x,y in sorted(C.nodes(data=True),key=lambda item: len(item[1]['members']), reverse=True)][:min_s_connected_index+1]
        print(f'Weakly_connected_components:{w_connected_comp}',f'\nstronly_connected_components:{s_connected_comp}')
        print(f'max_strongly_connected_nodes:\n{max_s_connected_nodes}',min_s_connected_index,len(max_s_connected_nodes[0][1]['members']))
    
        C0=G.copy()
        for i in G.nodes():
            if i not in max_s_connected_nodes[0][1]['members']:
                C0.remove_node(i)
        draw_network(C0, p=False, classification=classification, file_prefix='strongly_connected',save=True)
        return (density,w_connected_comp,s_connected_comp,max_s_connected_nodes,C,C0)


def plot_degree_dist(G,classification):
    '''
    # power distribution

    Parameters
    ----------
    G : graph: graph for centrality
    classification : str/list: test_a, test_r, test_ud, ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    power_law_dict: dict: Key: degree type, Values: alpha, sigma, xmin, R, p_val
    '''

    plt.rcParams.update(plt.rcParamsDefault)
    power_law_dict={}
    i_index=0
    for i in [G.degree(),G.in_degree,G.out_degree()]:
        deg_all=[v for k,v in i]
        deg_index=[k for k,v in sorted(collections.Counter(deg_all).items(),key=lambda item: item[0])]
        deg_count=[v for k,v in sorted(collections.Counter(deg_all).items(),key=lambda item: item[0])]
        deg_freq=[x/sum(deg_all) for x in deg_count]

        fit = powerlaw.Fit(deg_all)
        fit_R, fit_p = fit.distribution_compare('power_law', 'lognormal_positive')

        if i_index==1:
            label,width='in_degree',0.5
        elif i_index==2:
            label,width='out_degree',5
        elif i_index==0:
            label,width='degree',5

        i_index+=1
        plt.figure(1, figsize=(12,5))
        plt.subplot(1, 2, 1)
        plt.bar(deg_index,deg_freq,color='b',width=width, log=True)
        plt.title(f'{label} distribution')
        plt.xlabel(f'{label}')
        plt.ylabel('log p(k)')
        plt.subplot(1, 2, 2)
        fit.plot_pdf(marker='.',color='b')
        plt.title('power law distribution')
        plt.xlabel(f'{label}')
        plt.ylabel('log p(k)')
        plt.savefig(f'{out_path}/{label}_distribution_{cond_time}_{classification}.png',bbox_inches='tight', dpi=300)
        plt.show()

        power_law_dict[label]={'alpha':fit.power_law.alpha,'sigma':fit.power_law.sigma,'xmin':fit.xmin,'R':fit_R,'p':fit_p}
        print(f'{label} power law fit\nalpha={fit.power_law.alpha},sigma={fit.power_law.sigma},xmin={fit.xmin}')
        print(f'power_law vs. lognormal_positive: R={fit_R}, p={fit_p}\n')
        with open(f'{out_path}/power_law_{cond_time}_{classification}','a') as f:
            f.write(f'{label} power law fit\nalpha={fit.power_law.alpha},sigma={fit.power_law.sigma},xmin={fit.xmin}\n')
            f.write(f'power_law vs. lognormal_positive: R={fit_R}, p={fit_p}\n\n')
    return power_law_dict


def directed_cycle_search(G,list_path,comp,length,classification):
    '''
    Parameters
    ----------
    G : graph: graph for centrality
    list_path : str: name of genretated simple_path list, None to generate new one
    cycle_list : bool: wheteher to perform simple_sample
    comp : str/list/set: nodes required to present in the cycle, could be long processing time
    length : int: length of cycle
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    simple_cycle: list: list of cycles
    '''

    if list_path==None:
        simple_cycle=[x for x in list(nx.simple_cycles(G))]
        with open(f"{out_path}/simple_cycle_{cond_time}_{classification}.pickle", "wb") as f:   #Pickling
            pickle.dump(simple_cycle, f)
    else:
        with open(f"{out_path}/{list_path}.pickle", "rb") as f:   # Unpickling
            simple_cycle= pickle.load(f)
    print(f'all_cycle_len:{len(simple_cycle)}')
    if comp!=None and length!=None:
        selected_cycle=[x for x in simple_cycle if (set(comp).issubset(x)) and (len(x)==length)]
    elif comp!=None and length==None:
        selected_cycle=[x for x in simple_cycle if set(comp).issubset(x)]
    elif comp==None and length!=None:
        selected_cycle=[x for x in simple_cycle if len(x)==length]
    else:
        selected_cycle=simple_cycle
    print(f'Counts of directed_cycles, comp:{comp}, length:{length}\n',f'result size:{len(simple_cycle)}')
    return (simple_cycle,selected_cycle)

################################################# Position for partition plot #################################################
def community_layout(g, partition):
    """
    Compute the layout for a modular graph.

    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        graph to plot

    partition -- dict mapping int node -> int community
        graph partitions

    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions

    """
    pos_communities = _position_communities(g, partition, scale=7)
    pos_nodes = _position_nodes(g, partition, scale=1.2)
    # combine positions
    pos = dict()
    for node in g.nodes():
        pos[node] = pos_communities[node] + pos_nodes[node]
    return pos

def _position_communities(g, partition, **kwargs):
    '''
    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    '''
    between_community_edges = _find_between_community_edges(g, partition)
    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))
    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]
    return pos

def _find_between_community_edges(g, partition):
    edges = dict()
    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]
        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]
    return edges

def _position_nodes(g, partition, **kwargs):
    """
    Positions nodes within communities.
    """
    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]
    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, **kwargs)
        pos.update(pos_subgraph)
    return pos


################################################# GSEApy_enrich #################################################
def gsea_partition_enrich(G,partition,save_plot):
    '''
    # Enrichment
    Parameters
    ----------
    G : graph: graph for centrality
    partition_dict : dict: Key: nodes, Value: community
    save_plot : bool: wheter to save partition enrichment plot

    Returns
    -------
    partition_nodes: dict: Key: community, Values: nodes in community
    enr.res2d: df: enrichment dataframe for each partitioned community
    '''
    # http://amp.pharm.mssm.edu/Enrichr/
    gsea_lib = gp.get_library_name(database='Human')
    # load downloaded gmt file as gene_set
    # gseapy.parser.gsea_gmt_parser(gmt, min_size=3, max_size=50000)
    gene_sets=['GO_Biological_Process_2018','GO_Molecular_Function_2018','KEGG_2019_Human']

    partition_nodes={v:[] for v in partition.values()}
    for k,v in partition.items():
        partition_nodes[v].append(k)

    for k in partition_nodes:
        for gs in gene_sets:
            enr = gp.enrichr(gene_list=partition_nodes[k],
                             description=f'partition_nodes_{k}',
                             gene_sets=gs,
                             outdir=f'{out_path}/Enrichr_{cond_time}_{classification}',
                             cutoff=0.05,
                             figsize=(8,6),top_term=10,format='png')
            # if save_plot==True:
            #     gp.dotplot(enr.res2d, cutoff=0.05, top_term=10,title=f'community_{k}_gsea',
            #                ofname=f'{out_path}/Enrichr_{cond_time}_{classification}/{gs}_community_{k}_gsea.png',dpi=300)
    return (partition_nodes,enr.res2d)


###############################################################################################################################################

def maxflow_mincut(G,source,target,weight,classification):
    '''
    # maximum flow
    # min cuts

    Parameters
    ----------
    G : graph: graph for centrality
    source : str: source node
    target : str: target node
    weight ： str: capacity, hpi_edge
    classification : str/list: test_a, test_r, test_ud,ChIP_a, ChIP_r, ChIP_ud

    Returns
    -------
    power_law_dict: dict: Key: degree type, Values: alpha, sigma, xmin, R, p_val
    '''
    G=network_filter(G,classification)

    G.remove_edges_from(nx.selfloop_edges(G))
    maxflow_mincut={}
    #flow_func=nx.algorithms.flow.edmonds_karp  shortest_augmenting_path  preflow_push  dinitz  boykov_kolmogorov
    flow_value,flow_dict=nx.maximum_flow(G,_s=source,_t=target,capacity=weight,flow_func=nx.algorithms.flow.boykov_kolmogorov)
    maxflow_mincut['max_flow']={'flow_value':flow_value,'flow_path':flow_dict}

    cutset=set()
    cut_value, partition = nx.minimum_cut(G,_s=source,_t=target,capacity=weight,flow_func=nx.algorithms.flow.boykov_kolmogorov)
    reachable, non_reachable = partition
    for u, nbrs in ((n, G[n]) for n in reachable):
        cutset.update((u, v) for v in nbrs if v in non_reachable)
    maxflow_mincut['min_cut']={'cut_value':cut_value,'cutset':cutset}

    #H = nx.algorithms.connectivity.build_auxiliary_node_connectivity(G)
    #R = nx.algorithms.flow.build_residual_network(H, 'capacity')
    edge_independent_path=list(nx.edge_disjoint_paths(G,s=source,t=target,flow_func=nx.algorithms.flow.boykov_kolmogorov))
    node_independent_path=list(nx.node_disjoint_paths(G,s=source,t=target,flow_func=nx.algorithms.flow.boykov_kolmogorov))
    maxflow_mincut['independent_paths']={'edge_independent_path':edge_independent_path,'node_independent_path':node_independent_path}
    print(f'edge_independevnt_path:\n{maxflow_mincut["independent_paths"]["edge_independent_path"]}')
    print(f'node_independent_path:\n{maxflow_mincut["independent_paths"]["node_independent_path"]}')
    print(f'cutset:\n{maxflow_mincut["min_cut"]["cutset"]}')
    return (maxflow_mincut['independent_paths']['edge_independent_path'],maxflow_mincut['independent_paths']['node_independent_path'],maxflow_mincut['min_cut'])
edge_independent_path,node_independent_path,min_cut=maxflow_mincut(G=G_tf,source='ETV1',target='PDX1',weight='capacity', classification=classification)

