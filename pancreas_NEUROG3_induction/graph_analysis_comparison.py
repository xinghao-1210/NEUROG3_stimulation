#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:30:39 2021

@author: user
"""

import os
import glob
import pandas as pd
import numpy as np
import networkx as nx
import pickle
import community
from matplotlib import pyplot as plt
import ast


hioe_path = f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/HIOE_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_6tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx'
panc_path = f'/Users/user/Desktop/big-data_analysis/NGS_XZ/NGS/mLASSO-StARS/modeling/infTRN_lassoStARS/pancreas_NEUROG3_induction/outputs/networks_targ0p05_SS50_bS5/Network0p05_10tfsPerGene/prior_atac_Miraldi_q_ChIP_bias10_maxComb/networkx'
   
filt=['all', 'all_tfs', 'de', 'de_tfs']
graph_dict_hioe={}
graph_dict_panc={}

for ft in filt:
    nx_hioe_path = f'{hioe_path}/{ft}'
    nx_panc_path = f'{panc_path}/{ft}'
        
    Gl_hioe={}
    Gl_panc={}

    for f in glob.glob(f'{nx_hioe_path}/*/graphml*.gpickle', recursive=True):
        with open(f, 'rb') as ff:
            G_df = pickle.load(ff)
        t=int(f'{f.split("/")[-1].split("_")[1]}')
        Gl_hioe[t]=G_df
    Gl_hioe=dict(sorted(Gl_hioe.items()))
    graph_dict_hioe[ft]=Gl_hioe
    
        
    for f in glob.glob(f'{nx_panc_path}/*/graphml*.gpickle', recursive=True):
        with open(f, 'rb') as ff:
            G_df = pickle.load(ff)
        t=int(f'{f.split("/")[-1].split("_")[1]}')
        Gl_panc[t]=G_df
    Gl_panc=dict(sorted(Gl_panc.items()))
    graph_dict_panc[ft]=Gl_panc

time_pair=tuple(zip([24,48,72,96], [8,24,48,72]))

for ft in filt:
    for t in time_pair:
        print(ft, t)
        
        nx_hioe_path = f'{hioe_path}/comparison/{ft}/{t}'
        nx_panc_path = f'{panc_path}/comparison/{ft}/{t}'
        for p in [nx_hioe_path,nx_panc_path]:
            try:
                os.makedirs(p)
            except FileExistsError:
                pass
        
        H=graph_dict_hioe[ft][t[0]]
        P=graph_dict_panc[ft][t[1]]
        # diff node
        H_node_uni=H.copy()
        H_node_uni.remove_nodes_from(n for n in H if n in P)
        P_node_uni=P.copy()
        P_node_uni.remove_nodes_from(n for n in P if n in H)
        # sahred node
        H_node_share=H.copy()
        H_node_share.remove_nodes_from(n for n in H if n not in P)
        P_node_share=P.copy()
        P_node_share.remove_nodes_from(n for n in P if n not in H)
        # shared edge in shared node
        HP_edge_share=nx.intersection(H_node_share, P_node_share)
        # diff edge in shared node
        H_edge_uni=nx.difference(H_node_share, P_node_share)
        P_edge_uni=nx.difference(P_node_share, H_node_share)
        g_dict={'HIOE_node_uni':H_node_uni,'Pancreas_node_uni':P_node_uni,'HIOE_node_share':H_node_share,'Pancreas_node_share':P_node_share,
                'HIOE_Pancreas_edge_share':HP_edge_share,'HIOE_edge_uni':H_edge_uni,'Pancreas_edge_uni':P_edge_uni}
        
        for k,G in g_dict.items():
            summary_df=pd.DataFrame(index=g_dict.keys(), columns=['count','file']).fillna(0)
            
            if 'HIOE' in k:
                cond_time=t[0]
            else:
                cond_time=t[1]
            
            for out_path in [nx_hioe_path, nx_panc_path]:
                network_graphml_save(G,f'{out_path}/{k}.txt')
                if 'edge' not in k:
                    draw_network(G, True, ft, None, True)
                    draw_network(G, False, ft, None, True)
                
                # set of node and edge
                node_edge_set=sorted(list({n for n in G}))
                with open(f'{out_path}/{k}_list.txt','w') as f:
                    f.write('\n'.join(node_edge_set))
             
            summary_df.loc[k,'count']=len(node_edge_set)
            summary_df.loc[k,'file']=f'{k}_list.txt'
            summary_df.to_csv(f'Summary_node_edge_list_HIOE_Pancreas_{t}.txt')

# ====================================================================================
# ====================================================================================
def network_graphml_save(G, path):
    nx.write_weighted_edgelist(G, path=path.replace('.txt','_weighted_edgelist.txt'))
    nx.write_edgelist(G, path, comments='#', delimiter='\t', data=True, encoding='utf-8')
    nx.write_gpickle(G, path=path.replace('.txt','.gpickle'))
    return 'Graphml saved!'        
        


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
            fig.suptitle(f'{k}_{str(i).split(" ")[1]}_hubs top_{top_n}', fontsize=15)
            plt.bar(centrality['hubs'].keys(), centrality['hubs'].values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=7.5,rotation=75)
            for j in [x for x in centrality['hubs'] if (G.nodes[x][f'lfc_{cond_time}hpi']>=1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['hubs']).index(j)].set_color("red")
            for j in [x for x in centrality['hubs'] if (G.nodes[x][f'lfc_{cond_time}hpi']<=-1) and (G.nodes[x][f'padj_{cond_time}hpi'] < 0.05)]:
                plt.gca().get_xticklabels()[list(centrality['hubs']).index(j)].set_color("blue")
            plt.show()
            print(f'{k}_{str(i).split(" ")[1]}_hubs of {str(G)}:\n {centrality["hubs"]}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{k}_{str(i).split(" ")[1]}_hubs_centrality top_{top_n}_{classification}.png',bbox_inches='tight', dpi=300)

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
            print(f'{k}_{str(i).split(" ")[1]}_authorities of {str(G)}:\n {centrality["authorities"]}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{k}_{str(i).split(" ")[1]}_authorities_centrality top_{top_n}_{classification}.png',bbox_inches='tight', dpi=300)
        elif i in [nx.edge_betweenness_centrality,nx.edge_current_flow_betweenness_centrality]:
            centrality={k: v for k, v in sorted(centrality.items(),key=lambda item: item[1], reverse=True)[:top_n]}    #str(k) for tuple from edge_betweenness_centrality
            fig = plt.figure(figsize=(10,5))
            fig.suptitle(f'{str(i).split(" ")[1]} top_{top_n}', fontsize=15)
            plt.bar([str(x[0])+'-'+str(x[1]) for x in centrality.keys()], centrality.values(), color='b')
            plt.gca().get_xticklabels()[1].set
            plt.xticks(fontsize=5,rotation=75)
            plt.show()
            print(f'{k}_{str(i).split(" ")[1]} of {str(G)}:\n {centrality}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{k}_{str(i).split(" ")[1]}_centrality top_{top_n}_{classification}.png',bbox_inches='tight', dpi=300)
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
            print(f'{k}_{str(i).split(" ")[1]} of {str(G)}:\n {centrality}\n')
            if save_plot==True:
                fig.savefig(f'{out_path}/{k}_{str(i).split(" ")[1]}_centrality top_{top_n}_{classification}.png',bbox_inches='tight', dpi=300)
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
                
    # specific node shape based on DE by spliting node group
    up_nodes = [n for n,lfc in nx.get_node_attributes(G,f'lfc_{cond_time}hpi').items() if lfc>=1]
    down_nodes = [n for n,lfc in nx.get_node_attributes(G,f'lfc_{cond_time}hpi').items() if lfc<=-1]
    nan_nodes = [n for n,lfc in nx.get_node_attributes(G,f'lfc_{cond_time}hpi').items() if -1<lfc<1]
    # specific lable dict for font color based on DE
    up_nodes_dcit = {n:n for n in up_nodes}
    down_nodes_dcit = {n:n for n in down_nodes}
    nan_nodes_dcit = {n:n for n in nan_nodes}
    
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
        try:
            partition_community,enr_df=gsea_partition_enrich(H,partition=partition,save_plot=True)
            partition_group.Series(partition_community,name='partition_group').sort_index().sort_values(kind = 'mergesort')
            partition_group.to_csv(f'{out_path}/{k}_partition_comunity_{classification}.txt',sep='\t')
        except:
            print('Connection aborted no GSEA analysis')
            partition_community=None
        node_color=list(partition.values())
        if file_prefix==None:
            file_prefix='partitioned'
        else:
            file_prefix='partitioned_'+file_prefix
        partition_group=pd.Series(partition,name='partition_group').sort_index().sort_values(kind = 'mergesort')
        partition_group.to_csv(f'{out_path}/{k}_partition_nodes_{classification}.txt',sep='\t')

    plt.figure(figsize=(20,20))
    #draws nodes
    nx.draw_networkx_nodes(H,positions,node_color=node_color,nodelist=nodelist,
                           #the node size will be now based on its degree
                           node_size=node_sizes,alpha=0.8)
    #Styling for labels
    #nx.draw_networkx_labels(H,positions,label,font_size=5,font_color='k')
    nx.draw_networkx_labels(H,positions,up_nodes_dcit,font_size=7,font_color='m')
    nx.draw_networkx_labels(H,positions,down_nodes_dcit,font_size=7,font_color='green')
    nx.draw_networkx_labels(H,positions,nan_nodes_dcit,font_size=7,font_color='k')
    #draws the edges
    nx.draw_networkx_edges(H, positions, nodelist=up_nodes, edgelist=selected_edges,style='solid', width=[x/25 for x in selected_edges_weight], edge_color=edge_colors)
    nx.draw_networkx_edges(H, positions, nodelist=down_nodes, edgelist=selected_edges,style='solid', width=[x/25 for x in selected_edges_weight], edge_color=edge_colors)
    nx.draw_networkx_edges(H, positions, nodelist=nan_nodes, edgelist=selected_edges,style='solid', width=[x/25 for x in selected_edges_weight], edge_color=edge_colors)

    if save==True:
        if file_prefix==None:
            plt.title('HITS_spring_layout')
            plt.savefig(f'{out_path}/{k}_HITS_spring_layout_{classification}.png',bbox_inches='tight', dpi=300)
        else:
            plt.title(f'{file_prefix}_HITS_spring_layout')
            plt.savefig(f'{out_path}/{file_prefix}_{k}_HITS_spring_layout_{classification}.png',bbox_inches='tight', dpi=300)
    return partition_community     
        

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
        
        
        
        
    