# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 17:40:04 2019

@author: Ethan
"""


import networkx as nx
#import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt


blocks = ['Salt Data//salt-block4.csv']

for x in blocks:
    print(x)
    block = pd.read_csv(x)
    block = block.drop(['Unnamed: 0'], axis =1)
    
    for i in block.columns:
        block.loc[block[i] < .3, i] = 0
        block.loc[block[i] >= .3, i] = 1
    
    print("done with for loop")
    block = block.reset_index(drop = True)
    block.columns = range(block.shape[1])
    G = nx.from_pandas_adjacency(block)
    print(nx.info(G))
    nx.write_gpickle(G, str(x)[:-4] + '.p') 


#G.remove_nodes_from(list(nx.isolates(G)))
#g = nx.draw(G, node_size=1, width = .01)
 