# # -*- coding: utf-8 -*-
# """
# Created on Fri Nov 29 16:03:12 2019
#
# @author: samic
# """
#
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from networkx.algorithms import community
import community
# import networkx as nx
# from networkx.algorithms import community

# block = pd.read_csv('block1.csv')
# block = pd.read_csv('block2.csv')
# block = pd.read_csv('block3.csv')
block = pd.read_csv('block4.csv')

block.drop(block.columns[[1]], axis=1, inplace=True)

toGraph = block.stack().reset_index()

del block

toGraph.columns = ['var1', 'var2','value']

toGraph_filtered = toGraph.loc[ (toGraph['value'] > .8) & (toGraph['var1'] != toGraph['var2']) ]
G = nx.from_pandas_edgelist(toGraph_filtered, 'var1', 'var2')

# g = nx.draw(G, node_size=10)
# parts = community.best_partition(G)
# values = [parts.get(node) for node in G.nodes()]

#this works???
communities_generator = community.girvan_newman(G)
top_level_communities = next(communities_generator)
next_level_communities = next(communities_generator)
# print(sorted(map(sorted, next_level_communities)))
print(next_level_communities)
