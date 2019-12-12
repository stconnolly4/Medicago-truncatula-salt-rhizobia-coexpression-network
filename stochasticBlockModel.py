# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:28:05 2019

@author: samic
"""

import networkx as nx
from scipy.stats import kstest
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from networkx.generators.community import stochastic_block_model as sbm

path = "C:\\Users\\samic\\Documents\\Coexpression Network\\Pickled Networks\\"
graph = nx.read_gpickle(path + "rhizobia-block1.p")
deg_dist = nx.degree_histogram(graph)
ecdf = ECDF(deg_dist)

cluster_sizes = [100, 200]
probs = [[.1], [.2]]
results = []
#for p in probs:
#print(p)
model = sbm(cluster_sizes, probs)
k = kstest(deg_dist, model)
results += [(probs, k)]
 




#import pandas as pd
#df = pd.read_csv("C:\\Users\\samic\\Documents\\Coexpression Network\\Coexpression - CSVs\\no-rhizobia-block4.csv", nrows=3)

# rhizobia block 1: 9996
# rhizobia block 2: 9491
# rhizobia block 3: 8978
# rhizobia block 4: 2479

# non-rhizobia block 1: 9615
# non-rhizobia block 2: 7337
# non-rhizobia block 3: 7020
# non-rhizobia block 4: 6802


# salt block 1: 
# salt block 2: 
# salt block 3: 
# salt block 4: 


# non-salt block 1: 
# non-salt block 2: 
# non-salt block 3: 
# non-salt block 4: 