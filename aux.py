#Auxiliary Functions
from networkx import nx
import sys
from collections import defaultdict

def _get_com_wise_nodes(dictionary):
	#m = max(dictionary.values())
	louvain_p = defaultdict(set)
	for l in dictionary.keys():
	    louvain_p[dictionary[l]].add(l)
	return louvain_p

def printsomeinformation(node, com_node, best_node, incr,node_l,node_c):
	pass; return;
	print("Printing Information: {0} {1} {2} {3} {4} {5}".format(node, com_node, best_node, incr,node_l,node_c))

def build_network(layer, node_l, node_c,weightednode_l,weightednode_c, top, bot, couple, edge_l, edge_c,weighted):
    G = nx.Graph()
    if (weighted==1):
        for n in node_l:
            for n2 in node_l[n]:
                G.add_edge(n, n2,weight = weightednode_l[n][n2])
        for n in node_c:
            for n2 in node_c[n]:
                G.add_edge(n, n2,weight = weightednode_c[n][n2])
    
    else:
        for n in node_l:
            for n2 in node_l[n]:
                G.add_edge(n, n2)
        for n in node_c:
            for n2 in node_c[n]:
                G.add_edge(n, n2)
    
    return G

def read_raw_network(filename,weighted):
    print(filename)
    fp=open(filename,'r')
    line=fp.readline()
    line=line.rstrip()
    n_layer=int(float(line))
    layer={}
    node_l={}
    l_ID=1
    edge_l={}
    edge_c={}
    weightednode_l ={}
    weightednode_c ={}
    # f_el = open(filename+'_edges_list_commod'+str(g), 'w')
    for i in range(0,n_layer):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        layer[l_ID]=set()
        #print line
        for n in line:
            layer[l_ID].add(int(float(n)))
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        #print n_edge
        edge_l[l_ID]=n_edge
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(float(line[0]))
            n2=int(float(line[1]) )
            if n1 not in node_l:
                node_l[n1]=set()
            node_l[n1].add(n2)    
            if n2 not in node_l:
                node_l[n2]=set()
            node_l[n2].add(n1)
            
            if(weighted==1):
                weightednode_l[n1]  = weightednode_l.get(n1,dict())
                weightednode_l[n2]  = weightednode_l.get(n2,dict())
                weightednode_l[n1][n2] = int(float(line[2]))
                weightednode_l[n2][n2] = int(float(line[2]))

            
        l_ID+=1
        
    line=fp.readline()
    line=line.rstrip()
    n_couple=int(float(line))
    #print n_couple
    node_c={}      
    top={}
    bot={}
    c_ID=1
    couple={}

    for i in range(0,n_couple):
        line=fp.readline()
        #print line
        line=line.rstrip()
        line=line.split()
        top[c_ID]=int(float(line[0]))
        bot[c_ID]=int(float(line[1]))
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        #print n_edge
        edge_c[c_ID]=n_edge
        count_edge = 0
        for j in range(0,n_edge):
            line=fp.readline()
            line=line.rstrip()
            line=line.split()
            n1=int(float(line[0]))
            n2=int(float(line[1]))
            if n1 not in node_c:
                node_c[n1]=set()
            node_c[n1].add(n2)
            if n2 not in node_c:
                node_c[n2]=set()
            node_c[n2].add(n1)  
            count_edge += 1

            if(weighted==1):
                weightednode_c[n1]  = weightednode_c.get(n1,dict())
                weightednode_c[n2]  = weightednode_c.get(n2,dict())
                weightednode_c[n1][n2] = int(float(line[2]))
                weightednode_c[n2][n2] = int(float(line[2]))

            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')
        edge_c[c_ID] = count_edge
        c_ID=c_ID+1
    commu={}
    print(filename,filename.find('yelp'))
    if(filename.find('yelp')!=-1 ):
        mu=0
        ml_network =build_network(layer, node_l, node_c,weightednode_l,weightednode_c, top, bot, couple, edge_l, edge_c,weighted)
        return ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu,commu

    line=fp.readline()
    line=line.rstrip()
    #print line
    n_comm=int(float(line))
    commu={}
    com_ID=1
    for i in range(0,n_comm):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        commu[com_ID]=set()
        for n in line:
            commu[com_ID].add(int(float(n)))
        com_ID+=1      
    mu=0

    ml_network =build_network(layer, node_l, node_c,weightednode_l,weightednode_c, top, bot, couple, edge_l, edge_c,weighted)

    return ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu,commu
    
