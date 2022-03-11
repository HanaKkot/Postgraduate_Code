import networkx as nx
import matplotlib.pyplot as plt
from sympy import true
import numpy as np
import random

# Dataset path
pwd = "C:/Users/HanaKkot\Desktop/karate_dataset.csv"
pwd2 = "C:/Users/HanaKkot\Desktop/G1_dataset.csv"
# Open source data file
karate_dataset = open(pwd, "rb")
G1_dataset = open(pwd2, "rb")
#Gx = nx.read_edgelist(karate_dataset, delimiter=',', nodetype=int)
Gx = nx.Graph()
Gx = nx.read_adjlist(G1_dataset, delimiter=',', nodetype=int)
G1_dataset.close()
print('Node number of G1 is {}'.format(len(Gx.nodes())))
#karate_dataset.close()
#print(Gx.nodes())
SIR_S = []
SIR_I = []
SIR_R = []
# Plot structure figure
def draw_struct_fig(G):
    plt.figure(figsize=(10,7))
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=100)
    nx.draw_networkx_edges(G, pos = pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size = 5,font_color='white', font_family='sans-serif')
    plt.axis('off')
    plt.title('G1 network structure layout')
    plt.show()

# Plot adjacent matrix

def draw_adj_matrix(G):
    plt.imshow(nx.to_numpy_matrix(G))
    cbar = plt.colorbar()
    cbar.set_ticks([0,1])
    cbar.ax.set_yticklabels(['0','1'],)

    plt.title('Adjacency Matrix')
    cbar.set_label('link', rotation=270)
    plt.xlabel('nodes')
    plt.ylabel('nodes')
    plt.show()
    return None

def network_analysis(G):
    average_shortest_path_length = average_shortest_path_length(G)
    print("average_shortest_path_length is :", average_shortest_path_length)
    
draw_struct_fig(Gx)

draw_adj_matrix(Gx)

print("The population of karate network is {}".format(len(Gx.nodes())))

 
# Return the neighbour of the nodes in the social network.


# Define the SIR dynamic considering opinion dynamic(Polarization)

def SIR_dynamic_model_Opinion_Dynamic(G, I_init, beta, time, recoveredTime):
    
    # Using 'J-A' theory to update opinion value
    Initial_Infected_Nodes =I_init
    # Create timelist from time 0 to formal paramater "time" 
    timeList = []
    S_data = []
    I_data = []
    R_data = []
    node_num = len(G)

    # Return infected node in the network
    def getInfected(graph):
        return [x for x,y in graph.nodes(data=True) if y['infected'] == True and y['recovered'] == False]

    # Return neighbour node of the givenNode
    def getNeighbors(graph, givenNode):
        return [x for x in graph.neighbors(givenNode)]
    
    # Return recovered node in the network
    def getRecovered(graph):
        return [x for x,y in graph.nodes(data=True) if y['recovered'] == True]
    
    nodes = G.nodes()
    
    #setting initial conditions
    for node in nodes:
        if node in Initial_Infected_Nodes:
            G.nodes[node]['infected'] = True # If nodes are in the list, adding attribute to those nodes.
            print("Node {} is infected node".format(node))
        else:
            G.nodes[node]['infected'] = False
        
        G.nodes[node]['recovered'] = False # All people in the initial status are not recovered. 
        G.nodes[node]['Opinion'] = 1.0 - 2.0*random.random() #random value in [-1,1]
        G.nodes[node]['recoveredTime'] = 0 # Recovery Timing
    
    # Simulation
    for t in range(time):
        timeList.append(t) # Adding time in each iteration
        infected_list = getInfected(G)

        # Infected rules
        for i in infected_list:
            neighbors = getNeighbors(G, i) # Find neighbour
            for j in neighbors:
                randomValue = np.random.rand() # Generate a random number from [0,1)
                if randomValue <= beta:
                    G.nodes[j]['infected'] = True # Random infected
                    
                    G.nodes[j]["Opinion"]
        
        # Recovery rules
        for k in infected_list:
            if G.nodes[k]['recoveredTime'] > recoveredTime:
                G.nodes[k]['recovered'] = True
            G.nodes[k]['recoveredTime'] += 1
        
        # Collect simulation data 
        infected_list = getInfected(G)
        recovered_list = getRecovered(G)
        
        # Converted to the percentage by divided total nodes number
        # S+I+R = total nodes number
        I_data.append(len(infected_list)/node_num)
        S_data.append((node_num-len(infected_list)-len(recovered_list))/node_num)
        R_data.append(len(recovered_list)/node_num)
    
    # Plot figure
    plt.title('SIR opinion dynamic model')
    plt.plot(timeList, S_data, label = "Susceptible")
    plt.plot(timeList, I_data, label = "Infected")
    plt.plot(timeList, R_data, label = "Recovered")
    plt.legend()
    plt.show()

    return S_data, I_data, R_data

# Parameter setting
beta = 0.001
simulation_time = 400 
recoveryTime = 50 
Initial_infected_index = [5] 

SIR_S, SIR_I, SIR_R = SIR_dynamic_model_Opinion_Dynamic(Gx, Initial_infected_index, beta, simulation_time, recoveryTime)

#print(SIR_S)
#print(SIR_I)
#print(SIR_R)

