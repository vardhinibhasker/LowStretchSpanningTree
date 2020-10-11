import networkx as nx
import utree
from tqdm import tqdm
import matplotlib.pyplot as plt
import random
import tree

def bfs_tree(G, source):
    T = nx.Graph()
    T.add_node(source)
    es = nx.bfs_edges(G, source)
    for u,v in es:
        T.add_edge(u,v,weight=G[u][v]["weight"])
    return T

def ave_stretch(m, edges, T):
    edge_num = G.size()
    sum = 0
    for u, v, d in edges:
        weight = d.get('weight', 1)
        sum += nx.shortest_path_length(T, source=u, target=v)*weight
    return (sum/edge_num)

def weighted_ave_stretch(G, T):
    edge_num = G.size()
    sum = 0
    for u, v, d in G.edges(data=True):
        weight = d.get('weight', 1)
        sum += nx.shortest_path_length(T, source=u, target=v)*weight
    return (sum/edge_num)

def weighted_ave_stretch_lst(G, T):
    edge_num = G.size()
    sum = 0
    missed=0
    for u, v, d in G.edges(data=True):
        weight = d.get('weight', 1)
        try:
            sum += nx.shortest_path_length(T, source=u, target=v)/weight
        except nx.NetworkXNoPath:
            missed+=1
            continue
    return (sum/(edge_num - missed))

rand_weights = [x for x in range(1,6)] #[1,2,3,4,5,6,7,8,9,10]


x=[]
LST_stretches=[]
BFS_stretches=[]
MST_stretches=[]
#os = []
for n in tqdm(range(40, 1000, 100)):
    x.append(n)
    connected=False
    while not(connected):
        G=nx.random_regular_graph(4, n)
        if nx.is_connected(G): connected=True
    H = nx.Graph()
    H.add_nodes_from(G)
    for u, v in G.edges:
        w = random.choice(rand_weights)
        G[u][v]["weight"] = w
        H.add_edge(u,v,weight=1/w)

    MST = nx.maximum_spanning_tree(G)
    source = 0
    BFS = bfs_tree(G, source)
    LST = tree.weighted_lst(H, source)

    for u,v in BFS.edges:
        BFS[u][v]["weight"] = 1/BFS[u][v]["weight"]
    for u,v in MST.edges:
        MST[u][v]["weight"] = 1/MST[u][v]["weight"]


    LST_stretches.append(weighted_ave_stretch_lst(H, LST))
    MST_stretches.append(weighted_ave_stretch(G, MST))
    BFS_stretches.append(weighted_ave_stretch(G, BFS))


plt.plot(x, LST_stretches, '-*', color='m', label='Elkin et al')
plt.plot(x, MST_stretches, '-*', color='c', label='Kruksal')
plt.plot(x, BFS_stretches, '-*', color='g', label='BFS')
plt.ylabel('Average Stretch')
#plt.ylim(1.9, 2.05)
plt.legend()
plt.xlabel('Number of vertices')
plt.title("Random 4-Regular Graph for edges of weight 1 to 5")
#plt.title('Plot showing the average stretch against number of vertices in a complete graph.')
plt.show()
