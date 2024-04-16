import networkx as nx

def are_graphs_isomorphic(graph1_edges, graph2_edges):
    G1 = nx.Graph(graph1_edges)
    G2 = nx.Graph(graph2_edges)
    return nx.is_isomorphic(G1, G2)


def read_graph_from_file(filename):
    graph = []
    sum = 0
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) == 3:
                try:
                    u, v, w = map(int, numbers)
                    graph.append((u, v, w))
                    sum+=w
                except ValueError:
                    print(f"Ignoring line '{line.strip()}' as it doesn't contain three integers.")
    return graph,sum

# Read graph from file
filename1 = "randomGraphLarge.txt_out.txt"  
graph_edges1 ,sum1= read_graph_from_file(filename1)
# Create graph from edge list
ourGraph = nx.Graph()
ourGraph.add_weighted_edges_from(graph_edges1)



filename2 = "randomGraphLarge.txt"
graph_edges2,sum2 = read_graph_from_file(filename2)

refGraph = nx.Graph()
refGraph.add_weighted_edges_from(graph_edges2)

# Compute the minimum spanning tree
mst = nx.minimum_spanning_tree(refGraph)

sum_mst = 0

# Iterate over the edges of the minimum spanning tree
for u, v, data in mst.edges(data=True):
    # Add the weight of the current edge to the sum
    sum_mst += data['weight']



if(nx.is_connected(ourGraph) and sum1 == sum_mst):
     print("Our graph is a minimum spanning tree!")
else:
     print("ERROR: Our graph is NOT a minimum spanning tree!")



