#Generates synthetic networks
import os, sys, time
import subprocess
from lfr_multilayer_v3 import create_layers

new_dir = sys.argv[1]
dir_lfr = sys.argv[2]

if not os.path.exists(dir_lfr):
	os.makedirs(dir_lfr)
if not os.path.exists(new_dir):
	os.makedirs(new_dir)
if not os.path.exists(new_dir + "_networks"):
	os.makedirs(new_dir + "_networks")


#n=1000
#k=11
#maxk=30
n=100
k=6
maxk=10
nb_layers=2

# config1
list_alpha = [0.1, 0.2]
#list_alpha = [0.1, 0.2]
list_mu = [0.05]
#list_p = [1.0,0.95,0.9,0.85,0.8,0.75,0.7]
list_p = [0.6,0.7]
list_p1 = [0.3]
list_p2 = [0.0]

mu=0.05

intra_params={}
intra_params['n']=n
intra_params['k']=k
intra_params['maxk']=maxk
intra_params['mu']=0.05

create_layers(dir_lfr, nb_layers, intra_params)

#nb_com_layer , nb_link_layer = benchmark_lfr(0, dir_lfr, n, k, maxk, mu)
#nb_com_layer , nb_link_layer = benchmark_lfr(1, dir_lfr, n, k, maxk, mu)

# config2
# list_alpha = [0.6]
# list_mu = [0.05]
# list_p = [0.6]
# list_p1 = [0.1, 0.3, 0.6, 0.8]
# list_p2 = [0.3]

# config3
# list_alpha = [0.6]
# list_mu = [0.05]
# list_p = [0.6]
# list_p1 = [0.7]
# list_p2 = [0.1, 0.3, 0.6, 0.8]

 
# config4
# list_alpha = [0.6]
# list_mu = [0.05]
# list_p = [0.2,0.4,0.6,0.8]
# list_p1 = [0.7]
#list_p2 = [0.3]

# config5
# list_alpha = [0.6]
# list_mu = [0.05, 0.2, 0.4, 0.55]
# list_p = [0.6]
# list_p1 = [0.7]
# list_p2 = [0.3]


i=0
nb_iteration = len(list_alpha) * len(list_mu) * len(list_p) * len(list_p1) * len(list_p2)

for alpha in list_alpha:
	new_dir_alpha = new_dir + "/alpha-" + str(alpha)
	if not os.path.exists(new_dir_alpha):
		os.makedirs(new_dir_alpha)

	for p in list_p:
		new_dir_p = new_dir_alpha + "/p-" + str(p)
		if not os.path.exists(new_dir_p):
			os.makedirs(new_dir_p)
			
		for mu in list_mu:
			new_dir_mu = new_dir_p + "/mu-" + str(mu)
			if not os.path.exists(new_dir_mu):
				os.makedirs(new_dir_mu)	
		
			for p1 in list_p1:
				new_dir_p1 = new_dir_mu + "/p1-" + str(p1)
				if not os.path.exists(new_dir_p1):
					os.makedirs(new_dir_p1)
				
				for p2 in list_p2:
					new_dir_p2 = new_dir_p1 + "/p2-" + str(p2)
					if not os.path.exists(new_dir_p2):
						os.makedirs(new_dir_p2)

					i += 1
					print "%i/%i" % (i, nb_iteration)
					
					cmd= "cp -r %s/layer0 %s" % (dir_lfr,new_dir_p2)
					print cmd
					os.popen(cmd)
					cmd= "cp -r %s/layer1 %s" % (dir_lfr,new_dir_p2)
					print cmd
					os.popen(cmd)
					cmd = "python2.7 lfr_multilayer_v3.py %s %s -p %f -a %f -p1 %f -p2 %f" % (new_dir_p2, new_dir_p2, p, alpha, p1, p2)
					print cmd
					#sys.exit(1)
					os.popen(cmd)

					if os.path.isfile(new_dir_p2 + "/new_format"):
						cmd = "cp " + new_dir_p2 + "/new_format " + new_dir + "_networks/network_" + str(alpha)  + "_" + str(p) + "_" + str(mu) + "_" + str(p1) + "_" + str(p2)
		 				os.popen(cmd)
						print "File copied ----------"
					time.sleep(0.15)
