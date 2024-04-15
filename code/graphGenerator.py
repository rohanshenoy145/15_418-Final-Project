import networkx as nx
import random



def graphGenerator(numVertices, numEdges):
    if numEdges < numVertices - 1:
        raise ValueError("Number of edges must be at least num_vertices - 1 to ensure connectivity")
    
    if(numEdges > (numVertices)*(numVertices-1)/2):
        raise ValueError("Number of edges must be at most n(n-1)/2 to prevent redundancy")

    r = random.randint(0,numVertices -1 )
    inTree = [False for _ in range(numVertices)]
    next = [ -1 for _ in range(numVertices)]
    for i in range (numVertices):
        inTree[i] = False
        next[r] = -1
        inTree[r] = True
    
    for i in range(numVertices):
        u = i
        while not inTree[u]:
            next[u] = randomSuccessor(u,numVertices)
            u = next[u]
        u = i
        while not inTree[u]:
            inTree[u] = True
            u = next[u]


    return next

def randomSuccessor(u,numVertices):
     
    while (True):
        randomSuccessor = random.randint(0,numVertices - 1)
        if(randomSuccessor != u):
            return randomSuccessor


numVertices = 100000
numEdges = 300000
print(graphGenerator(numVertices,numEdges))


    


# numNodes = 100000
# numEdges = 300000
# random_weighted_graph = graphGenerator(numNodes, numEdges)
# file_path = "randomGraph.txt"
# with open(file_path,'w') as file:
#     for u, v, data in random_weighted_graph.edges(data=True):
#         weight = data['weight']
#         file.write(f"{u} {v} {weight}\n")


                
    


    


    




