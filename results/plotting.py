import pandas as pd
from functools import partial
import timeit
import numpy as np
import matplotlib.pyplot as plt
import utree
import networkx as nx
from math import log
from tqdm import tqdm
#import tree


SETUP = '''
import utree
import networkx as nx

G = nx.erdos_renyi_graph({}, 0.5)

'''


CODE = 'utree.unweighted_lst(G, 1)'

ts=[]
xs = []
es=[]
for n in tqdm(range(100, 1600, 100)):
    xs.append(n)
    sn = str(n)
    es.append(n*log(n,2))
    times = timeit.timeit(setup = SETUP.format(sn),
                          stmt = CODE,
                          number = 10)
    ts.append(times)

vs = [x*x*log(x,2) for x in xs]

coeffs = np.polyfit(vs, ts, 1)
fit = np.poly1d(coeffs)


fig, ax = plt.subplots()
plt.scatter(xs, ts)
plt.plot(xs, fit(vs), color='m')
plt.xlabel("Number of vertices")
plt.ylabel("Run time (s)")
plt.title("Erdos Renyi Random Graph, p = 0.5")
plt.show()
