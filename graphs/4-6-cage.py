# 4 REG GIRTH 6
import networkx as nx

d={1 : {2, 3, 4, 5}, 2 : {6, 7, 8}, 3 : {9, 10, 11}, 4 : {12, 13, 14}, 5 : {15, 16, 17}, 6 : {18, 19, 20}, 7 : {21, 22, 23}, 8 : {24, 25, 26}, 9 : {18, 21, 24}, 10 : {19, 22, 25}, 11 : {20, 23, 26}, 12 : {18, 22, 26}, 13 : {19, 23, 24}, 14 : {20, 21, 25}, 15 : {18, 23, 25}, 16 : {19, 21, 26}, 17 : {20, 22, 24}}
G=nx.Graph(d)
