# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:03:12 2019

@author: samic
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#block = pd.read_csv('block1.csv')
#block = pd.read_csv('block2.csv')
#block = pd.read_csv('block3.csv')
block = pd.read_csv('block4.csv')

block.drop(block.columns[[1]], axis=1, inplace=True)

toGraph = block.stack().reset_index()

del block

toGraph.columns = ['var1', 'var2','value']

toGraph_filtered = toGraph.loc[ (toGraph['value'] > .8) & (toGraph['var1'] != toGraph['var2']) ]
G = nx.from_pandas_edgelist(toGraph_filtered, 'var1', 'var2')

#g = nx.draw(G, node_size=10)

# centrality
dictNodes = nx.eigenvector_centrality(G, max_iter=1000)
deg = nx.degree(G)