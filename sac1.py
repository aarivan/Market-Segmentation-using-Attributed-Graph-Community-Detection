import sys
import csv
import random
import math
import pandas as pa
import numpy
from igraph import *
from scipy import spatial

#Check for valid arguments
if (len(sys.argv) != 2):
    print('Invalid Arguments Error: Enter - python sac1.py <alpha>')
    sys.exit(1)


def getEdgeList(edges):
	edgelist_file = open('./data/fb_caltech_small_edgelist.txt')
	edge_list = edgelist_file.read().split("\n")

	for edge in edge_list:
		v = edge.split(' ')
		if v[0] != '' and v[1] != '':
			edges.append((int(v[0]),int(v[1])))
	return edges


def getAttrList(attrList):
	attrList = pa.read_csv('./data/fb_caltech_small_attrlist.csv')
	return attrList


def graphFormation(graph, num_vertices, edges, attributes_names):
	graph.add_vertices(num_vertices)
	graph.add_edges(edges)
	graph.es["weight"] = [1 for x in range(len(edges))]
	for attribute in attributes_names :
		graph.vs[attribute] = list(attributes_ds[attribute])
	return graph


def compute_cossimilarity() :
	global num_vertices, cossim
	cossim = [[0 for x in range(num_vertices)] for x in range(num_vertices)] 
	for i in range(0, num_vertices) :
		vert_i = graph.vs.select(i)[0].attributes().values()
		for j in range(i, num_vertices) :
			vert_j = graph.vs.select(j)[0].attributes().values()
			distance = spatial.distance.cosine(vert_i, vert_j) + 1.0
			# cossim[j][i] = 1.0 / (distance)
			# cossim[i][j] = cossim[i][j]
			cossim[i][j] = 1.0 / (distance)
			cossim[j][i] = cossim[i][j]


def qneuman(x, comm, nod, num_edges):
	return x - sum(graph.degree(comm)) * graph.degree(nod) / (2 * num_edges)

def qattr(g_attr, cossim, comm, nod):
	for item in comm:
		g_attr = g_attr + cossim[item][nod]
	return g_attr / len(comm) / len(comm)

def compute_modularity_gain(nod, comm):
	x = 0
	deg = 0
	num_edges = len(graph.es)
	comm = list(set(comm))

	for item in comm:
		if graph.are_connected(nod, item):
			ind = graph.get_eid(nod, item)
			x += graph.es["weight"][ind]

	g_neuman = qneuman(x, comm, nod, num_edges)
	# g_neuman = x - sum(graph.degree(comm)) * graph.degree(nod) / (2 * num_edges)
	g_neuman = g_neuman / (2.0 * num_edges)
	g_attr = 0.0
	# g_attr = qattr(g_attr, cossim, comm, nod)
	for item in comm:
		g_attr = g_attr + cossim[item][nod]
	g_attr = g_attr / len(comm) / len(comm)
	# print g_attr
	return alpha * g_neuman + (1 - alpha) * g_attr


def find_community(community, vertex):
    for item in community:
        if vertex in item:
            return item
    return []


def form_community(graph, community):
    count = 0
    for item in range(num_vertices):
        gains = []
        comm_i = find_community(community, item)

        max_gain = -1
        max_comm = []
        
        for com in community:
            gain = compute_modularity_gain(item, com)
            if gain > 0:
            	if gain > max_gain:
            	    max_gain = gain
            	    max_comm = com
                
        if set(comm_i) != set(max_comm):
        	if max_gain > 0:
	            comm_i.remove(item)
    	        max_comm.append(item)

    	        count += 1
            	if len(comm_i) == 0:
                	community.remove([])
    return count


def phaseone(graph, cossim, community):
	#n = 1
	cossim = compute_cossimilarity()
	cnt = form_community(graph, community)
	i = 0
	while cnt > 0 and i < 15:
		print "Iteration" ,  i
		i+=1
		cnt = form_community(graph, community)


def phasetwo(graph, cossim, mapped_communities, mapped_vertices):
	global num_vertices

	nv = 0
	for comm in mapped_communities :
		for v in comm :
			mapped_vertices[v] = nv
		nv += 1

	graph.contract_vertices(mapped_vertices, combine_attrs = "mean")
	graph.simplify(multiple = True, loops = True)

	num_vertices = nv
	mapped_communities = [[x] for x in range(num_vertices)]
	graph.es["weight"] = [0 for x in range(len(graph.es))]
	
	for edge in edges :
		left_comm = mapped_vertices[edge[0]]
		right_comm = mapped_vertices[edge[1]]
		
		if left_comm != right_comm:
			id = graph.get_eid(left_comm, right_comm)
			graph.es["weight"][id] += 1
	
	cossim = compute_cossimilarity()
	phaseone(graph, cossim, mapped_communities)


alpha = float(sys.argv[1])

attributes_ds = []
attributes_ds = getAttrList(attributes_ds)

# store the number of vertices(rows) and attributes(cols)
num_vertices = attributes_ds.shape[0]
num_attributes = attributes_ds.shape[1]

edges = []
edges = getEdgeList(edges)

graph = Graph()

attributes_names = list(attributes_ds.columns.values)

graph = graphFormation(graph, num_vertices, edges, attributes_names)
cossim = [[0 for x in range(num_vertices)] for x in range(num_vertices)] 
# cossim = compute_cossimilarity()

communities = [[x] for x in range(num_vertices)]
mapped_vertices = [0 for x in range(num_vertices)]

phaseone(graph, cossim, communities)
phasetwo(graph, cossim, communities, mapped_vertices)

op_file = open("./communities.txt", "w")
for item in communities:
    for i in range(len(item)):
        if i != 0:
            op_file.write(",")
        op_file.write(str(item[i]))
    op_file.write("\n")
op_file.close()