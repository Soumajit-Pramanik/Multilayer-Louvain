import sys, os, subprocess, random, math, copy
from random import randint, random, uniform, sample

separator = '\t'

fgh = open("test.data", 'a')

# Depth First Search of communities starting from the first layer
def recursive_cross_com_gathering(cross_communities_l1_l2, index_com1):
	if index_com1 not in cross_communities_l1_l2:
		return [index_com1]

	list_coms = [index_com1]

	for index_com2 in cross_communities_l1_l2[index_com1]:		
		list_coms.extend(recursive_cross_com_gathering(cross_communities_l1_l2, index_com2))

	return list_coms

def load_layers(dir_lfr):
	layers = {}

	layers[0] = {}
	layers[1] = {}

	nb_nodes, nb_links, nb_coms = load_com_nodes(0, dir_lfr, layers)
	layers[0]["from_id_node"] = 0
	layers[0]["from_id_com"] = 0
	layers[0]['nb_coms'] = nb_coms
	layers[0]['nb_links'] = nb_links
	layers[0]['nb_nodes'] = nb_nodes
	nb_nodes, nb_links, nb_coms = load_com_nodes(1, dir_lfr, layers)
	layers[1]["from_id_node"] = layers[0]['nb_nodes']
	layers[1]["from_id_com"] = layers[0]['nb_coms']
	layers[1]['nb_coms'] = nb_coms
	layers[1]['nb_links'] = nb_links
	layers[1]['nb_nodes'] = nb_nodes

	return layers

def create_layers(outdir, nb_layers, dict_intra_params):
	layers = {}

	# create every layer based on LFR benchmark 
	for i in range(0, nb_layers):
		nb_com_layer , nb_link_layer = benchmark_lfr(i, outdir, dict_intra_params['n'], dict_intra_params['k'], dict_intra_params['maxk'], dict_intra_params['mu'])
		layers.setdefault(i, {})
		layers[i]['nb_coms'] = nb_com_layer
		layers[i]['nb_links'] = nb_link_layer
		layers[i]['nb_nodes'] = dict_intra_params['n']
		print "Layer %i created (with %i intra layer coms)" % (i, nb_com_layer)

	index_node = 0
	index_com = 0
	for i in range(0, nb_layers):
		rewrite_nodes_and_coms(i, outdir, index_node, index_com)
		
		layers[i]["from_id_node"] = index_node
		layers[i]["from_id_com"] = index_com

		index_node += layers[i]['nb_nodes']
		index_com += layers[i]['nb_coms']

	return layers


def rewrite_nodes_and_coms(num, directory, from_node, from_com):
	dir_layer = directory + "/layer" + str(num)
	filename_nodes = dir_layer + "/network.dat"
	filename_coms = dir_layer + "/community.dat"

	f_nodes = open(filename_nodes + "bis", 'w')
	with open(filename_nodes, 'r') as file:
		for line in file:
			line = line.replace("\n", "").replace(" ", "").split(separator)

			new_id1 = str(int(line[0]) + from_node)
			new_id2 = str(int(line[1]) + from_node)

			f_nodes.write(new_id1 + separator + new_id2 + '\n')
	f_nodes.close()

	os.popen("rm " + filename_nodes)
	os.popen("mv " + filename_nodes + "bis " + filename_nodes)

	f_coms = open(filename_coms + "bis", 'w')
	with open(filename_coms, 'r') as file:
		for line in file:
			line = line.replace("\n", "").replace(" ", "").split(separator)

			new_id1 = str(int(line[0]) + from_node)
			new_id2 = str(int(line[1]) + from_com)

			f_coms.write(new_id1 + separator + new_id2 + '\n')
	f_coms.close()

	os.popen("rm " + filename_coms)
	os.popen("mv " + filename_coms + "bis " + filename_coms)

def load_com_nodes(num, directory, layers):
	dir_layer = directory + "/layer" + str(num)
	filename_coms = dir_layer + "/community.dat"

	layers[num]['com_nodes'] = {}
	layers[num]['node_com'] = {}
	layers[num]['nodes'] = []

	nb_link = 0
	with open(filename_coms, 'r') as file:
		for line in file:
			line = line.replace("\n", "").replace(" ", "").split(separator)
			nb_link += 1

			index_com = int(line[1]) - 1
			layers[num]['com_nodes'].setdefault(index_com, [])
			layers[num]['com_nodes'][index_com].append(line[0])
			layers[num]['nodes'].append(line[0])
			layers[num]['node_com'][line[0]] = index_com

	return len(layers[num]['nodes']), nb_link, len(layers[num]['com_nodes'])

def benchmark_lfr(num, directory, n, k, maxk, mu):
	dir_layer = directory + "/layer" + str(num)
	if not os.path.exists(dir_layer):
		os.makedirs(dir_layer)

	cmd = "./benchmark -N %i -k %f -maxk %i -mu %f" % (n, k, maxk, mu)
	print cmd
	os.popen(cmd)
	cmd = "mv *.dat %s" % (dir_layer)
	os.popen(cmd)
	
	cmd = "cat %s/community.dat | cut -f 2 | sort | uniq | wc -l" % (dir_layer)
	process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	errcode = process.returncode
	nb_com_layer = int(out)

	cmd = "wc -l %s/network.dat | cut -d\" \" -f1" % (dir_layer)
	process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	errcode = process.returncode
	nb_link_layer = int(out)

	return nb_com_layer, nb_link_layer

def create_cross_layer_communities(nb_layers, layers, alpha):	
	cross_communities_l1_l2 = {}
	intra_coms = {}
	nb_intra_coms = 0

	if alpha == 0:
		for i in range(0, nb_layers):
			intra_coms[i] = range(layers[i]["from_id_com"], layers[i]["from_id_com"]+layers[i]['nb_coms'])
			nb_intra_coms += len(intra_coms[i])			

		return [], intra_coms, nb_intra_coms

	for i in range(0, nb_layers-1):
		nb_min = min(layers[i]['nb_coms'], layers[i+1]['nb_coms'])
		nb_cross_coms = int(math.ceil(nb_min * alpha))

		list_id_coms_layer1 = range(layers[i]["from_id_com"], layers[i]["from_id_com"]+layers[i]['nb_coms'])
		list_id_coms_layer2 = range(layers[i+1]["from_id_com"], layers[i+1]["from_id_com"]+layers[i+1]['nb_coms'])
		
		j = 0
		while j < nb_cross_coms:
			j+= 1

			i_l1 = randint(0, len(list_id_coms_layer1)-1)
			i_l2 = randint(0, len(list_id_coms_layer2)-1)

			index_coms_l1 = list_id_coms_layer1[i_l1]
			index_coms_l2 = list_id_coms_layer2[i_l2]

			del list_id_coms_layer1[i_l1]
			del list_id_coms_layer2[i_l2]

			cross_communities_l1_l2.setdefault(index_coms_l1, [])
			cross_communities_l1_l2[index_coms_l1].append(index_coms_l2)

		intra_coms[i] = list_id_coms_layer1
		intra_coms[i+1] = list_id_coms_layer2
		nb_intra_coms += len(list_id_coms_layer1) + len(list_id_coms_layer2)

	cross_communities = []
	for index_com1 in cross_communities_l1_l2:
		list_coms = recursive_cross_com_gathering(cross_communities_l1_l2, index_com1)
		cross_communities.append(list_coms)

	return cross_communities, intra_coms, nb_intra_coms

def apply_p1_parameter(nb_layers, layers, couplings_links, cross_communities, p1):
	if len(cross_communities) == 0:
		return 0

	nb_link_inside_cross_coms = 0

	for i in range(0, nb_layers-1):
		couplings_links.setdefault(i, {})
		couplings_links[i].setdefault(i+1, [])

		list_id_coms_layer1 = range(layers[i]["from_id_com"], layers[i]["from_id_com"]+layers[i]['nb_coms'])
		list_id_coms_layer2 = range(layers[i+1]["from_id_com"], layers[i+1]["from_id_com"]+layers[i+1]['nb_coms'])
		list_cross_l1_l2 = []
		for id_com_l1, id_com_l2 in cross_communities:
			if id_com_l1 in list_id_coms_layer1 and id_com_l2 in list_id_coms_layer2:
				list_cross_l1_l2.append((id_com_l1, id_com_l2))

		for id_com_l1, id_com_l2 in list_cross_l1_l2:
			list_node_com_l1 = layers[i]['com_nodes'][id_com_l1]
			list_node_com_l2 = layers[i+1]['com_nodes'][id_com_l2]

			nb_nodes_coupling_l1 = int(round(len(list_node_com_l1) * p1))
			nb_nodes_coupling_l2 = int(round(len(list_node_com_l2) * p1))

			sample_node_l1 = sample(list_node_com_l1, nb_nodes_coupling_l1)
			sample_node_l2 = sample(list_node_com_l2, nb_nodes_coupling_l2)

			s_l1 = copy.copy(sample_node_l1)
			s_l2 = copy.copy(sample_node_l2)

			while len(s_l1) > 0 or len(s_l2) > 0:
				n_min = min(len(s_l1), len(s_l2))
				for j in range(n_min):
					i_l1 = randint(0, len(s_l1)-1)
					i_l2 = randint(0, len(s_l2)-1)
					
					node_l1 = s_l1[i_l1]
					node_l2 = s_l2[i_l2]

					couplings_links[i][i+1].append((node_l1, node_l2))
					nb_link_inside_cross_coms += 1	

					s_l1.remove(node_l1)
					s_l2.remove(node_l2)

				if len(s_l1) > 0:
					s_l2 = sample(sample_node_l2, min(len(s_l1), len(sample_node_l2)))
				elif len(s_l2) > 0:
					s_l1 = sample(sample_node_l1, min(len(s_l2), len(sample_node_l1)))

	return nb_link_inside_cross_coms

def apply_p2_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, p2):
	if len(intra_coms) == 0:
		return 0

	nb_link_outside_cross_coms = 0

	for i in range(0, nb_layers-1):
		couplings_links.setdefault(i, {})
		couplings_links[i].setdefault(i+1, [])

		list_id_coms_layer1 = range(layers[i]["from_id_com"], layers[i]["from_id_com"]+layers[i]['nb_coms'])
		list_id_coms_layer2 = range(layers[i+1]["from_id_com"], layers[i+1]["from_id_com"]+layers[i+1]['nb_coms'])
		list_cross_l1_l2 = []
		for id_com_l1, id_com_l2 in cross_communities:
			if id_com_l1 in list_id_coms_layer1 and id_com_l2 in list_id_coms_layer2:
				list_cross_l1_l2.append((id_com_l1, id_com_l2))

		list_node_l1 = []
		list_node_l2 = []

		for id_com_l1 in intra_coms[i]:
			 list_node_l1.extend(layers[i]['com_nodes'][id_com_l1])
		for id_com_l2 in intra_coms[i+1]:
			 list_node_l2.extend(layers[i+1]['com_nodes'][id_com_l2])
		for id_com_cross_l1, id_com_cross_l2 in list_cross_l1_l2:
			list_node_l1.extend(layers[i]['com_nodes'][id_com_cross_l1])
			list_node_l2.extend(layers[i+1]['com_nodes'][id_com_cross_l2])

		for id_com_l1 in intra_coms[i]:
			list_node_com_l1 = layers[i]['com_nodes'][id_com_l1]
			nb_nodes_coupling_l1 = int(round(len(list_node_com_l1) * p2))
			sample_node_l1 = sample(list_node_com_l1, nb_nodes_coupling_l1)

			for j in range(len(sample_node_l1)):
				n_l2 = randint(0, len(list_node_l2)-1)
				couplings_links[i][i+1].append((sample_node_l1[j], list_node_l2[n_l2]))
				nb_link_outside_cross_coms += 1

		for id_com_l2 in intra_coms[i+1]:
			list_node_com_l2 = layers[i+1]['com_nodes'][id_com_l2]
			nb_nodes_coupling_l2 = int(round(len(list_node_com_l2) * p2))
			sample_node_l2 = sample(list_node_com_l2, nb_nodes_coupling_l2)

			for j in range(len(sample_node_l2)):
				while True:
					n_l1 = randint(0, len(list_node_l1)-1)
					if (list_node_l1[n_l1], sample_node_l2[j]) not in couplings_links[i][i+1]:
						break

				couplings_links[i][i+1].append((list_node_l1[n_l1], sample_node_l2[j]))
				nb_link_outside_cross_coms += 1

	return nb_link_outside_cross_coms

def apply_p_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, nb_link_inside_cross_coms, nb_link_outside_cross_coms, p):
	if len(cross_communities) == 0:
		return nb_link_inside_cross_coms, nb_link_outside_cross_coms

	list_com_node = {}

	for i in range(0, nb_layers-1):
		list_id_coms_layer1 = range(layers[i]["from_id_com"], layers[i]["from_id_com"]+layers[i]['nb_coms'])
		list_id_coms_layer2 = range(layers[i+1]["from_id_com"], layers[i+1]["from_id_com"]+layers[i+1]['nb_coms'])
		
		for id_com_l1, id_com_l2 in cross_communities:
			if id_com_l1 in list_id_coms_layer1 and id_com_l2 in list_id_coms_layer2:
				list_com_node[(id_com_l1, id_com_l2)] = {}
				list_com_node[(id_com_l1, id_com_l2)][i] = set()
				list_com_node[(id_com_l1, id_com_l2)][i+1] = set()
				
				for id_node_l1, id_node_l2 in couplings_links[i][i+1]:
					if id_node_l1 in layers[i]['com_nodes'][id_com_l1] and id_node_l2 in layers[i+1]['com_nodes'][id_com_l2]:
						list_com_node[(id_com_l1, id_com_l2)][i].add(id_node_l1)
						list_com_node[(id_com_l1, id_com_l2)][i+1].add(id_node_l2)

	p_temp = (1.0 * nb_link_inside_cross_coms) / (nb_link_inside_cross_coms + nb_link_outside_cross_coms)
	print "actual p value : %f" % p_temp
	fgh.write("actual p value : %f\n" % p_temp)
	nb_new_coupling = 0

	if p > p_temp:
		total_possible_nb_coupling = 0
		for id_com_l1, id_com_l2 in list_com_node:
			for i in range(0, len(list_com_node[(id_com_l1, id_com_l2)])-1):
				total_possible_nb_coupling += len(list_com_node[(id_com_l1, id_com_l2)][i]) * len(list_com_node[(id_com_l1, id_com_l2)][i+1])
		total_possible_nb_coupling -= nb_link_inside_cross_coms
	
		list_id_cross_coms = list_com_node.keys()
		nb_new_coupling = 0
		while p > p_temp and nb_new_coupling < total_possible_nb_coupling:
			i_com = randint(0, len(list_id_cross_coms)-1)
			i_layer = randint(0, nb_layers-2)
		
			i_node_l1 = randint(0, len(list_com_node[list_id_cross_coms[i_com]][i_layer])-1)
			i_node_l2 = randint(0, len(list_com_node[list_id_cross_coms[i_com]][i_layer+1])-1)

			node_l1 = list(list_com_node[list_id_cross_coms[i_com]][i_layer])[i_node_l1]
			node_l2 = list(list_com_node[list_id_cross_coms[i_com]][i_layer+1])[i_node_l2]

			if (node_l1, node_l2) not in couplings_links[i_layer][i_layer+1]:
				couplings_links[i_layer][i_layer+1].append((node_l1, node_l2))
				nb_new_coupling += 1
				p_temp = (1.0 * nb_link_inside_cross_coms + nb_new_coupling) / (nb_link_inside_cross_coms + nb_link_outside_cross_coms)

		nb_link_inside_cross_coms += nb_new_coupling
		print "Nb of coupling links created (inside cross-layer communities) : %i" % nb_new_coupling
	
	elif p < p_temp:
		list_id_cross_coms = list_com_node.keys()

		total_possible_nb_coupling = 0
		for i in range(0, len(list_id_cross_coms)):
			for j in range(0, len(list_id_cross_coms)):
				if i == j:
					continue

				com_l1_1, com_l2_1 = list_id_cross_coms[i]
				com_l1_2, com_l2_2 = list_id_cross_coms[j]
				
				total_possible_nb_coupling += len(list_com_node[(com_l1_1, com_l2_1)][0]) * len(list_com_node[(com_l1_2, com_l2_2)][1])

		nb_new_coupling = 0
		while p < p_temp and nb_new_coupling < total_possible_nb_coupling:
			# for initialisation
			com_l1_1, com_l2_2 = list_id_cross_coms[0]
			com_l1_2 = com_l2_1 = 0

			while (com_l1_1, com_l2_2) in list_id_cross_coms:
				i_com1 = randint(0, len(list_id_cross_coms)-1)
				i_com2 = randint(0, len(list_id_cross_coms)-1)

				if i_com1 != i_com2:					
					com_l1_1, com_l2_1 = list_id_cross_coms[i_com1]
					com_l1_2, com_l2_2 = list_id_cross_coms[i_com2]

			i_layer = randint(0, nb_layers-2)
		
			i_node_l1 = randint(0, len(list_com_node[(com_l1_1, com_l2_1)][i_layer])-1)
			i_node_l2 = randint(0, len(list_com_node[(com_l1_2, com_l2_2)][i_layer+1])-1)

			node_l1 = list(list_com_node[(com_l1_1, com_l2_1)][i_layer])[i_node_l1]
			node_l2 = list(list_com_node[(com_l1_2, com_l2_2)][i_layer+1])[i_node_l2]

			if (node_l1, node_l2) not in couplings_links[i_layer][i_layer+1]:
				couplings_links[i_layer][i_layer+1].append((node_l1, node_l2))
				nb_new_coupling += 1
				p_temp = (1.0 * nb_link_inside_cross_coms) / (nb_link_inside_cross_coms + nb_link_outside_cross_coms + nb_new_coupling)

		nb_link_outside_cross_coms += nb_new_coupling
		print "Nb of coupling links created (outside cross-layer communities) : %i" % nb_new_coupling

	print "p value reached : %f" % p_temp
	fgh.write("p value reached : %f\n\n" % p_temp)

	return nb_link_inside_cross_coms, nb_link_outside_cross_coms

def write_soum_format(dir_lfr, outdir, nb_layers, layers, nb_couplings, couplings_links, nb_com, cross_communities, intra_coms):
	f = open(outdir + "/new_format", 'w')

	#nb of layers
	f.write(str(nb_layers) + '\n')

	for i in range(0, nb_layers):
		# nodes of layer i
		f.write(" ".join(layers[i]['nodes']) + '\n')

		filelayer = dir_lfr + "/layer" + str(i) + "/network.dat"
		# edges of layer i
		edges = []
		with open(filelayer, 'r') as file:
			for line in file:
				line = line.replace("\n", "").split(separator)

				edges.append((line[0], line[1]))


		f.write(str(len(edges)) + '\n')
		for e1, e2 in edges:
			f.write(e1 + ' ' + e2 + '\n')

	#nb couplings
	f.write(str(len(couplings_links)) + '\n')
	
	# coupling edges
	for l1 in couplings_links:
		for l2 in couplings_links[l1]:
			f.write(str(l1+1) + ' ' + str(l2+1) + '\n')

			f.write(str(len(couplings_links[l1][l2])) + '\n')
			for e1, e2 in couplings_links[l1][l2]:
				f.write(e1 + ' ' + e2 + '\n')

	# nb of community 
	f.write(str(nb_com) + '\n')
	
	for list_coms in cross_communities:
		for index_com in list_coms:
			for i in range(0, nb_layers):
				if index_com in layers[i]['com_nodes']:
					f.write(" ".join(layers[i]['com_nodes'][index_com]))
					break
			f.write(' ')
		f.write('\n')

	for i in range(0, nb_layers):
		for index_com in intra_coms[i]:
			f.write(" ".join(layers[i]['com_nodes'][index_com]) + '\n')

	f.close()

def usage():
	print "Usage : ./lfr_multilayer <outdir> <nb_layers> <INTRA> <INTER>"
	print "--------------------------------------"
	print "outdir : choose a directory where will be stored the multilayer network"
	print "nb_layers : number of layers of the multilayer network"
	print "---------------INTRA---------------"
	print "Intralayer parameters (over all layers) :"
	print "-n : number of nodes"
	print "-k : average degree"
	print "-maxk : maximum degree"
	print "-mu : mixing parameter"
	print "---------------INTER---------------"
	print "-p : set p to the control fraction of interlayer edges (coupling links) within community"
	print "-a : set alpha parameter to control the fraction of multilayer communities"
	print "-p1 : controls the fraction of nodes connected within coupling layer, from cross-layer communities"
	print "-p2 : controls the fraction of nodes connected within coupling layer, from intra-layer communities"
	print "-d : set density of interlayers edges"

#Main function
def main(argv):
	if len(argv) != 10:
		print "here",len(argv)
		usage()
		sys.exit(1)

	dict_intra_params = {}
	dict_inter_params = {}
	outdir = ""
	nb_layers = 2

	skip = False
	for i in range(0, len(argv)):
		if skip:
			skip = False
			continue

		if i == 0:
			outdir = argv[0]
		elif i == 1:
			dir_lfr = argv[1]
		else:
			skip = True
			'''
			if argv[i] == "-n":
				dict_intra_params['n'] = float(argv[i+1])
			elif argv[i] == "-k":
				dict_intra_params['k'] = float(argv[i+1])
			elif argv[i] == "-maxk":
				dict_intra_params['maxk'] = float(argv[i+1])
			elif argv[i] == "-mu":
				dict_intra_params['mu'] = float(argv[i+1])
			'''
			if argv[i] == "-p":
				dict_inter_params['p'] = float(argv[i+1])
			elif argv[i] == "-a":
				dict_inter_params['a'] = float(argv[i+1])
			elif argv[i] == "-p1":
				dict_inter_params['p1'] = float(argv[i+1])
			elif argv[i] == "-p2":
				dict_inter_params['p2'] = float(argv[i+1])
			elif argv[i] == "-d":
				dict_inter_params['d'] = float(argv[i+1])
			elif argv[i] == "-h":
				usage()
				sys.exit(1)
			else:
				print "%s parameter is unknown" % argv[i]
				sys.exit(1)

	layers = load_layers(dir_lfr)
	#layers = create_layers(outdir, nb_layers, dict_intra_params)

	cross_communities, intra_coms, nb_intra_coms = create_cross_layer_communities(nb_layers, layers, dict_inter_params['a'])	
	nb_com = len(cross_communities) + nb_intra_coms
	
	print "Setting network communities (alpha = %f)" % (dict_inter_params['a'])
	print "%i communities created - %i single layer - %i cross-layer" % (nb_com, nb_intra_coms, len(cross_communities))

	couplings_links =  {}
	nb_link_inside_cross_coms = apply_p1_parameter(nb_layers, layers, couplings_links, cross_communities, dict_inter_params['p1'])
	nb_link_outside_cross_coms =  apply_p2_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, dict_inter_params['p2'])
	nb_link_inside_cross_coms, nb_link_outside_cross_coms = apply_p_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, nb_link_inside_cross_coms, nb_link_outside_cross_coms, dict_inter_params['p'])

	nb_couplings = nb_link_inside_cross_coms + nb_link_outside_cross_coms 

	f = open(outdir + "/couplings", 'w')
	for i in range(0, nb_layers-1):
		for (index_node1, index_node2) in couplings_links[i][i+1]:
			f.write(index_node1 + ' ' + index_node2 + '\n')
	f.close()

	write_soum_format(dir_lfr, outdir, nb_layers, layers, nb_couplings, couplings_links, nb_com, cross_communities, intra_coms)
	
	print "Statistics :"
	for i in range(0, nb_layers):
		print "Layer %i" % (i)
		print "  - Nb of links %i" % layers[i]["nb_links"]

	for i in range(0, nb_layers-1):
		print "Cross-layer %i-%i" % (i, i+1)
		print "  - Nb of coupling links : %i " % (len(couplings_links[i][i+1]))
		print "  - Nb inside coms : %i " % (nb_link_inside_cross_coms)
		print "  - Nb outside coms : %i " % (nb_link_outside_cross_coms)

if __name__  == "__main__":
	main(sys.argv[1:])
	fgh.close()
