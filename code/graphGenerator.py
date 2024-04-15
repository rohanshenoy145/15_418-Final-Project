import networkx as nx
import random



def graphGenerator(numVertices, numEdges):
    if numEdges < numVertices - 1:
        raise ValueError("Number of edges must be at least num_vertices - 1 to ensure connectivity")
    
    if(numEdges > (numVertices)*(numVertices-1)/2):
        raise ValueError("Number of edges must be at most n(n-1)/2 to prevent redundancy")

    masterList = []
    numEdgesAdd = 0
    edgeExists = set()
    S = {i for i in range(numVertices)}
    nodes = {i for i in range(numVertices)}
    T = set()
    current_node = random.sample(list(S), 1).pop()
    S.remove(current_node)
    T.add(current_node)


    while len(S)>0:
        neighbor_node = random.sample(list(nodes), 1).pop()
        #print(neighbor_node)
        if neighbor_node not in T:
            #edge = (current_node, neighbor_node,random.randint(0,100))
            weight = random.randint(0,100)
            masterList.append((current_node, neighbor_node,weight))
            edgeExists.add((current_node,neighbor_node))
            numEdgesAdd+=1
            S.remove(neighbor_node)
            T.add(neighbor_node)
        current_node = neighbor_node


    while(numEdgesAdd < numEdges):
        firstVert = random.randint(0,numVertices-1)
        secondVert = random.randint(0,numVertices -1)
        if firstVert!=secondVert and (firstVert,secondVert)not in edgeExists  and (secondVert,firstVert) not in edgeExists:
            numEdgesAdd+=1
            weight = random.randint(0,100)
            masterList.append((firstVert,secondVert,weight))
            edgeExists.add((firstVert,secondVert))

    return masterList

#     r = random.randint(0,numVertices -1 )
#     inTree = [False for _ in range(numVertices)]
#     next = [ -1 for _ in range(numVertices)]
#     for i in range (numVertices):
#         inTree[i] = False
#         next[r] = -1
#         inTree[r] = True
    
#     for i in range(numVertices):
#         u = i
#         while not inTree[u]:
#             next[u] = randomSuccessor(u,numVertices)
#             u = next[u]
#         u = i
#         while not inTree[u]:
#             inTree[u] = True
#             u = next[u]


#     return next

# def randomSuccessor(u,numVertices):
     
#     while (True):
#         randomSuccessor = random.randint(0,numVertices - 1)
#         if(randomSuccessor != u):
#             return randomSuccessor


numVertices = 100000
numEdges = 200000
x= graphGenerator(numVertices,numEdges)
file_path = "randomGraphLarge.txt"
with open(file_path,'w') as file:
    for data in x:
       
        file.write(f"{data[0]} {data[1]} {data[2]}\n")





# numNodes = 100000
# numEdges = 300000
# random_weighted_graph = graphGenerator(numNodes, numEdges)
# 


                
    


    


    




