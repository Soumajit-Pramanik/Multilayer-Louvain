#-------------------------------------------------------------------------------------
#
#Generating networks by adding weights to the edges
#
#-------------------------------------------------------------------------------------
import os, sys, time
import subprocess
import networkx as nx
import random
import pickle
import random
import math
from sklearn.metrics import *
from collections import defaultdict
import matplotlib.pyplot as plt
from multiprocessing import Pool
from copy import deepcopy

desired_dm = 1.0
N=100
num_layers =2
mu = 0.05
#filepath = "./netsForDtDmDb/_networks/netsByGenerateNetsv3/alpha0.7/"

#file  ="network_0.7_0.7_0.05_0.7_0.0_3"
#file  ="network_0.3_0.7_0.05_0.7_0.0_10"

def applyweighttoDM(weightednode_c,desired_dm,initial_dm,weighttoadd_dm,node_c):
	loop = weighttoadd_dm
	while(loop>0):
		n1 = random.choice(weightednode_c.keys())
		n2 = random.choice(weightednode_c[n1].keys())
		weightednode_c[n1][n2] +=1
		weightednode_c[n2][n1] +=1
		loop-=1
	return weightednode_c

def sumupweights(weightednode_l):
	sum = 0
	for k in weightednode_l.keys():
		for k1 in weightednode_l[k].keys():
			sum+=weightednode_l[k][k1]
	#print("sum of weightednote_l: {0}".format(sum));

def build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c):
	G = nx.Graph()
	for n in node_l:
		for n2 in node_l[n]:
			G.add_edge(n, n2)
	for n in node_c:
		for n2 in node_c[n]:
			G.add_edge(n, n2)
	return G

def read_raw_network(filename):
    ##print(filename)
    fp=open(filename,'r')
    line=fp.readline()
    line=line.rstrip()
    n_layer=int(float(line))
    layer={}
    node_l={}
    l_ID=1
    edge_l={}
    edge_c={}
    # f_el = open(filename+'_edges_list_commod'+str(g), 'w')
    for i in range(0,n_layer):
        line=fp.readline()
        line=line.rstrip()
        line=line.split()
        layer[l_ID]=set()
        ###print line
        for n in line:
            layer[l_ID].add(int(float(n)))
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ###print n_edge
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
            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')

        l_ID+=1
        
    line=fp.readline()
    line=line.rstrip()
    n_couple=int(float(line))
    ###print n_couple
    node_c={}      
    top={}
    bot={}
    c_ID=1
    couple={}

    for i in range(0,n_couple):
        line=fp.readline()
        ###print line
        line=line.rstrip()
        line=line.split()
        top[c_ID]=int(float(line[0]))
        bot[c_ID]=int(float(line[1]))
        
        couple[c_ID]=layer[top[c_ID]].union(layer[bot[c_ID]])
        
        line=fp.readline()
        line=int(float(line.rstrip()))
        n_edge=line
        ###print n_edge
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
            # f_el.write(str(n1-1)+' '+str(n2-1)+'\n')
        edge_c[c_ID] = count_edge
        c_ID=c_ID+1

    line=fp.readline()
    line=line.rstrip()
    ###print line
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
    mu=0.05

    ml_network =build_network(layer, node_l, node_c, top, bot, couple, edge_l, edge_c)

    return ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu

def writenet(dt,db, dm ,weightednode_c, weightednode_l, ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12,filepath):
	filename = filepath + str(float(dm)) +"_"+ str(float(dt)) +"_" + str(float(db))
	fp = open(filename,'w')

	fp.write(str(num_layers)+"\n")

	#write all nodes in each layer
	for l in layer:
		towrite =""
		for n in layer[l]:
			towrite= towrite+ str(n)+" "
		fp.write(towrite+'\n')

		fp.write(str(int(2*E[l]))+'\n')
		#print(2*E[l])
		count=0
		for n1 in layer[l]:
			for n2 in node_l[n1]:
				count+=1
				fp.write(str(n1) + ' ' + str(n2) +' ' + str(weightednode_l[n1][n2]) + '\n')
		#print("count=" ,count)
	fp.write("1\n")
	fp.write("1 2\n")
	fp.write(str(int(2*E12))+'\n')
	for n1 in node_c:
		for n2 in node_c[n1]:
			fp.write(str(n1) + ' ' + str(n2) +' ' + str(weightednode_c[n1][n2]) + '\n')

	fp.write(str(len(commu))+'\n')
	for c in commu:
		towrite =""
		for n in commu[c]:
			towrite= towrite+ str(n)+" "
		fp.write(towrite+' \n')

	fp.close()

#def create(dt,db,ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12):

def getinfo(graph,layer,commu,node_l, node_c):
	#Construct a dict = {node1: layer, node2: layer......}
    nodelayer = {}
    nodecomm ={}
    for l in layer:
        for node in layer[l]:
            nodelayer[node] = l

    #calculate Intra_inter ----------------------------------------------------------
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)

    #--------------------------------------------------------------------------------

    #calculate |E1| , |E2| , |E12|---------------------------------------------------
    E={}
    E12=0
    for l in layer:
        E[l]=0
        for n in layer[l]:
            for nei in node_l.get(n,set()):
                E[l]+= 1
        E[l] = E[l]/2

    for n in node_c:
        for nei in node_c[n]:
            E12+=1
    E12 = E12/2
    #--------------------------------------------------------------------------------

    for c in commu:
    	for n in commu[c]:
    		nodecomm[n] = c

    return nodecomm , nodelayer, intra_inter,E, E12

def applyremoval(d,layer,node_l,weightednode_l,commu,nodelayer,nodecomm,E,l):
	toremove = E[l] - (N*d/2.0)
	removeoutofcommu = int(round(mu*toremove))
	removewithincommu  =int(round(toremove - removeoutofcommu))

	#print("toremove: {0} , removewithincommu: {1} ,removeoutofcommu: {2}".format(toremove,removewithincommu,removeoutofcommu))
	#print("E[1]: ",E[l])

	loop = removewithincommu
	removed = 0
	while(loop>0):
		n1 = random.choice(node_l.keys())
		if n1 not in layer[l]:
			continue
		templist = list(node_l[n1])
		if len(templist)<=0:
			continue
		n2 = random.choice(templist)
		if n2 not in layer[l]:
			continue
		if(nodecomm[n1]!=nodecomm[n2] or n1==n2 or weightednode_l[n1][n2]<=0):
			continue
		weightednode_l[n1][n2]-=1
		weightednode_l[n2][n1]-=1
		loop-=1
		removed+=1


	loop = removeoutofcommu
	while(loop>0):
		n1 = random.choice(node_l.keys())
		if n1 not in layer[l]:
			continue
		templist = list(node_l[n1])
		if len(templist)<=0:
			continue
		n2 = random.choice(templist)
		if n2 not in layer[l]:
			continue
		if(nodecomm[n1]==nodecomm[n2] or weightednode_l[n1][n2]<=0):
			continue
		weightednode_l[n1][n2]-=1
		weightednode_l[n2][n1]-=1
		loop-=1
		removed +=1

	#print("weight removed: {0}".format(2*removed))
	#print("for layer",l," density achieved: ",2.0*(E[l]+removed*1.0)/N," desired: ",d)
	return E,weightednode_l

def applyadditional(d,layer,node_l,weightednode_l,commu,nodelayer,nodecomm,E,l):
	edgestoadd = ((N*d)/2.0) - E[l]
	edgesoutofcommu  =int(round(mu*edgestoadd))
	edgesincommu = int(round(edgestoadd - edgesoutofcommu))

	#print("edgestoadd: {0} , edgesincommu: {1} ,edgesoutofcommu: {2} to achieve d: {3}".format(edgestoadd,edgesincommu,edgesoutofcommu,d))
	##print("E[1]: ",E[l])

	added = 0
	loop = edgesincommu
	while(loop>0):
		n1 = random.choice(node_l.keys())
		if n1 not in layer[l]:
			continue
		templist = list(node_l[n1])
		if len(templist)<=0:
			continue
		n2 = random.choice(templist)
		if n2 not in layer[l]:
			continue
		if(nodecomm[n1]!=nodecomm[n2] or n1==n2):
			continue
		weightednode_l[n1][n2]+=1
		weightednode_l[n2][n1]+=1
		loop-=1
		added+=1

	loop = edgesoutofcommu
	while(loop>0):
		n1 = random.choice(node_l.keys())
		if n1 not in layer[l]:
			continue
		templist = list(node_l[n1])
		if len(templist)<=0:
			continue
		n2 = random.choice(templist)
		if n2 not in layer[l]:
			continue
		if(nodecomm[n1]==nodecomm[n2] or n1==n2):
			continue
		weightednode_l[n1][n2]+=1
		weightednode_l[n2][n1]+=1
		loop-=1
		added +=1

	#E[l] += (edgesoutofcommu+edgesincommu)
	#print("weight added: {0}".format(2*added))
	#print("for layer",l," density achieved: ",2.0*(E[l]+added*1.0)/N," desired: ",d)
	return E,weightednode_l

def checkforNegweights(weightednode_l):
	for k in weightednode_l.keys():
		for k1 in weightednode_l[k].keys():
			if(weightednode_l[k][k1]<0):
				return 1 
	return 0;

def getSeries(filename,newfilepath):
	ml_network, layer, node_l, node_c, top, bot, couple, edge_l, edge_c, mu, commu =read_raw_network(filename)
	nodecomm,nodelayer, intra_inter,E, E12 = getinfo(ml_network, layer,commu, node_l, node_c)

	#initial_dm= int(round(E12*1.0/N))
	initial_dm= (E12*1.0/N)
	initial_dt = int(round(2.0*E[1]/N))
	initial_db = int(round(2.0*E[2]/N))
	##print("initial_db: {0} , initial Dt: {1} , initial Dm: {2}".format(initial_db,initial_dt,initial_dm))

	totalweightofdmreq = desired_dm*N
	weighttoadd_dm = totalweightofdmreq - E12

	##print("Weight to add in dm: {0} , so as to get Dm: {1}".format(weighttoadd_dm,desired_dm))

	#-------Initialise weighted node_l and node_c---------------------------
	weightednode_l = dict()
	for n1 in node_l:
		weightednode_l[n1]=weightednode_l.get(n1,dict())
		for n2 in node_l[n1]:
			weightednode_l[n1][n2] =1
			weightednode_l[n2] = weightednode_l.get(n2,dict())
			weightednode_l[n2][n1] =1

			#weightednode_l[n1].add(tempdict[n2])
	weightednode_c = dict()
	for n1 in node_c:
		weightednode_c[n1]=weightednode_c.get(n1,dict())
		for n2 in node_c[n1]:
			weightednode_c[n1][n2] =1
			weightednode_c[n2] = weightednode_c.get(n2,dict())
			weightednode_c[n2][n1] =1
	#-----------------------------------------------------------------------
	
	weightednode_c = applyweighttoDM(weightednode_c,desired_dm,initial_dm,weighttoadd_dm,node_c)
	initial_dm = desired_dm

	required_fractions= [0.1,0.3,0.5,0.7,0.9,1.0,3.0,5.0,7.0,9.0,10.0]
	required_dts = [v*initial_dm for  v in reversed(required_fractions)]

	print(len(set(required_dts)))

	node_lcopy = deepcopy(node_l)
	Ecopy= deepcopy(E)
	weightednode_lcopy = deepcopy(weightednode_l)
	weightednode_ccopy = deepcopy(weightednode_c)

	sumupweights(weightednode_l)

	for dt in required_dts:
		#dt=5.0
		node_l = deepcopy(node_lcopy)
		E = deepcopy(Ecopy)
		weightednode_l = deepcopy(weightednode_lcopy)
		weightednode_c = deepcopy(weightednode_ccopy)

		#print("Before applying dt")
		sumupweights(weightednode_l)

		if(dt < initial_dt):
			E,weightednode_l = applyremoval(dt,layer,node_l,weightednode_l ,commu,nodelayer,nodecomm,E,1)
		else:
			E,weightednode_l = applyadditional(dt,layer,node_l,weightednode_l,commu,nodelayer,nodecomm,E,1)
		
		node_lcopydt = deepcopy(node_l)
		Ecopydt= deepcopy(E)
		
		#print("After applying dt")
		sumupweights(weightednode_l)

		weightednode_lcopydt = deepcopy(weightednode_l)
		weightednode_ccopydt = deepcopy(weightednode_c)

		for db in required_dts:
			#db=10.0
			#print("for bd: ",db," dt:",dt)
			node_l = deepcopy(node_lcopydt)
			E = deepcopy(Ecopydt)
			weightednode_l = deepcopy(weightednode_lcopydt)
			weightednode_c = deepcopy(weightednode_ccopydt)
			
			##print("Before applying db")
			#sumupweights(weightednode_l)
			
			if(db < initial_db):
				E,weightednode_l = applyremoval(db,layer,node_l,weightednode_l,commu,nodelayer,nodecomm,E,2)
			else:
				E,weightednode_l = applyadditional(db,layer,node_l,weightednode_l,commu,nodelayer,nodecomm,E,2)
			
			##print("After applying db")
			#sumupweights(weightednode_l)
			
			'''
			check = checkforNegweights(weightednode_l)
			if(check==1):
				#print(dt,db)
				sys.exit()
			'''
			writenet(dt,db, initial_dm ,weightednode_c, weightednode_l, ml_network, layer, node_l, node_c,commu,nodecomm,nodelayer, intra_inter,E, E12,newfilepath)
			##print("******************************************************************")
		
		
			
		

#getSeries("./netsForDtDmDb/_networks/baseNetworks/" + file)

