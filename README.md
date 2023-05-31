# Multilayer Louvain
This repository contains the code for our proposed "Multilayer Louvain" algorithm which is able to detect communities in multilayer networks. Simultaneously, we also provide the code to generate synthetic multilayer networks with ground truth communities. This kind of networks can be used for evaluating the performance of the community detection algorithms on multilayer networks.

## Setup -
  + We require python 2.7 and a handful of other supporting libraries. To install the other dependencies, use: 
  ```pip install -r requirements.txt```

## Multilayer Louvain algorithm 
  + How to run? -
    + Parameters - Update the following two paths in Multilayer-Louvain.py as required,
      + Pathtosave - directory path to save result file 
      + networkpath - directory path from where to read network files
    + Command - 
      ```python Multilayer-Louvain.py``` 

## Synthetic Network Generation
  + How to run? -
    + Parameters - Update the following variables in generate_networks.py,
      + list_alpha - Set of 'alpha' values for which networks need to be generated 
      + list_mu - Set of 'mu' values for which networks need to be generated
      + list_p - Set of 'p' values for which networks need to be generated
      + list_p1 - Set of 'p1' values for which networks need to be generated 
      + list_p2 - Set of 'p2' values for which networks need to be generated 
      + n - Number of nodes in each layer
      + k - Average degree of the individual network layers
      + maxk - Maximum degree of the individual network layers
    + Command -
      ```python generate_networks.py <folder_name1> <folder_name2>```
    + This generates 2-layer multilayer networks in the folder named 'folder_name1_Networks' with all possible combinations of the specified parameter values. 'folder_name2' would contain the LFR networks.
    + The format of the generated network files is the following -
      + number_of_layers
      
        layer1 vertices
        
        number_of_layer1_edges
        
        layer1_edge1
        
        layer1_edge2
        
        ...
        
        layer2 vertices
        
        number_of_layer2_edges
        
        layer2_edge1
        
        layer2_edge2
        
        ...
        
        number_of_couplings
        
        coupling1_top_layer
        
        coupling1_bot_layer
        
        number_of_coupling1_edges
        
        coupling1_edge1
        
        coupling1_edge2
        
        ...
        
        coupling2_top_layer
        
        coupling2_bot_layer
        
        number_of_coupling2_edges
        
        coupling2_edge1
        
        coupling2_edge2
        
        ...
        
        number_of_communities
        
        community1_vertices
        
        community2_vertices
        
        ...
        

