#Generates synthetic networks
import os, sys, time
import subprocess
from lfr_multilayer_v3 import create_layers

nb_layers=2
new_dir = sys.argv[1]
dir_lfr = sys.argv[2]

if not os.path.exists(dir_lfr):
	os.makedirs(dir_lfr)
if not os.path.exists(new_dir):
	os.makedirs(new_dir)
if not os.path.exists(new_dir + "_networks"):
	os.makedirs(new_dir + "_networks")


n=100
k=6
maxk=10
mu=0.05

# config1
list_alpha = [0.1, 0.2]
list_mu = [mu]
list_p = [0.6,0.7]
list_p1 = [0.3]
list_p2 = [0.0]



intra_params={}
intra_params['n']=n
intra_params['k']=k
intra_params['maxk']=maxk
intra_params['mu']=mu

create_layers(dir_lfr, nb_layers, intra_params)



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
