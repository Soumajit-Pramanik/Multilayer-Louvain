#Class definition
class Status :
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}

    layer={}
    node_l={}
    node_c={}       
    top={}
    bot={}
    edge_l={}
    edge_c={}
    couple={}
    mu = 0

    in_layer_in_comm ={}
    in_layer_out_comm ={}
    out_layer_in_comm ={}
    out_layer_out_comm ={}

    def __init__(self) :
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])

        self.layer=dict([])
        self.node_l=dict([])
        self.node_c=dict([])
        self.top=dict([])
        self.bot=dict([])
        self.edge_l=dict([])
        self.edge_c=dict([])
        self.couple=dict([])
        self.mu = 0

        self.in_layer_in_comm   = dict()
        self.in_layer_out_comm  = dict()
        self.out_layer_in_comm  = dict()
        self.out_layer_out_comm = dict()

    def __str__(self) :
        return ("node2com : " + str(self.node2com) + " degrees : "
            + str(self.degrees) + " internals : " + str(self.internals)
            + " total_weight : " + str(self.total_weight)) 

    def copy(self) :
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight
        new_status.layer=self.layer.copy()
        new_status.node_l=self.node_l.copy()
        new_status.node_c=self.node_c.copy()
        new_status.top=self.top.copy()
        new_status.bot=self.bot.copy()
        new_status.edge_l=self.edge_l.copy()
        new_status.edge_c=self.edge_c.copy()
        new_status.couple=self.couple.copy()
        new_status.mu = self.mu
        new_status.in_layer_in_comm = self.in_layer_in_comm.copy()
        new_status.in_layer_out_comm = self.in_layer_out_comm.copy()
        new_status.out_layer_in_comm = self.out_layer_in_comm.copy()
        new_status.out_layer_out_comm = self.out_layer_out_comm.copy()
        return new_status

    def updatelists(self, graph):
        self.in_layer_in_comm   = dict()
        self.in_layer_out_comm  = dict()
        self.out_layer_in_comm  = dict()
        self.out_layer_out_comm = dict()

        node2layer = dict()
        for l in self.layer:     #nodeset is the set of nodes in each layer
            nodeset = self.layer[l]
            for node in nodeset:
                node2layer[node] = l
        #print("node2layer: ",node2layer)
        #update in_layer_in_comm
        datakey = 'weight'                  #this is the key used to get edge weight from weightdict

        for node in graph.nodes():
            node_neighbours = graph[node]  #node_neighbours will be a dict {2: {'weight':3}, 3: {'weight':4}}
            self.in_layer_in_comm[node] = 0
            self.in_layer_out_comm[node] = 0
            self.out_layer_in_comm[node] = 0
            self.out_layer_out_comm[node] = 0
            
            for dest,edge_data in node_neighbours.items():
                edge_weight = edge_data.get("weight", 1)
                if(self.node2com[dest] == self.node2com[node]): #both nodes are in the same community
                    if(node2layer[node] == node2layer[dest]): #both nodes in same layer
                        self.in_layer_in_comm[node] += edge_weight
                        if(node == dest): self.in_layer_in_comm[node] += edge_weight     #Add self loop twice
                    else:
                        self.out_layer_in_comm[node] += edge_weight

                else:                                 #both nodes are in different community
                    if(node2layer[node] == node2layer[dest]): #both nodes in same layer
                        self.in_layer_out_comm[node] += edge_weight
                    else:
                        self.out_layer_out_comm[node] += edge_weight
        #print("In layer in comm: ",self.in_layer_in_comm)


    def init(self, graph, part = None) :

        """Initialize the status of a graph with every node in one community"""
        count = 1
        #count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight = 'weight')

        if part == None :
            for node in graph.nodes() :
                self.node2com[node] = count
                deg = float(graph.degree(node, weight = 'weight'))
                if deg < 0 :
                    raise ValueError("Bad graph type, use positive weights")
                self.degrees[count] = deg
                self.gdegrees[node] = deg
                self.loops[node] = float(graph.get_edge_data(node, node,
                                                 {"weight":0}).get("weight", 1))
                self.internals[count] = self.loops[node]
                count = count + 1
        else :
            for node in graph.nodes() :
                com = part[node]
                self.node2com[node] = com
                deg = float(graph.degree(node, weight = 'weight'))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg
                inc = 0.
                for neighbor, datas in graph[node].items() :
                    weight = datas.get("weight", 1)
                    if weight < 0 :
                        raise ValueError("Bad graph type, use positive weights")
                    if part[neighbor] == com :
                        if neighbor == node :
                            inc += float(weight)
                        else :
                            inc += float(weight) / 2.
                self.internals[com] = self.internals.get(com, 0) + inc

        self.updatelists(graph)
