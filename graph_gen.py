"""
Code to generate random networks. 
The concepts are taken from 'How bad is naive multicasting?' (Doar, Leslie) and
'How to model an internetwork' (Zegura et al.)

The edge probability uses Doar-Leslie's probability
   P(u,v) = ke/n * Bexp(-d/a*L)
   k, B, a are parameters (default values given by paper)
   e is average nodal degree. n is number of nodes, L is maximum distance
   between any two nodes, d is the distance between u and v
(The values of parameters k, B, a used in this script are from Zegura et al.)

The nodes are placed uniformly within a grid and for each pair the above
probability is calculated to determine if there is an edge between the nodes.
Euclidean distances are used as distance between two nodes. The edges are
assumed to be bidirectional.

Note, the script allows you to specify a nodal degree, but the graphs are
randomly generated, so any particular run probably won't have the exact nodal
degree specified.

author: Neal Charbonneau
"""

from sys import argv, exit
from math import sqrt, fabs, exp
from random import randint, random


class Node:
    def __init__(self, x, y, id):
        self._x = x
        self._y = y
        self._id = id
    def distance(self, other):
        x = fabs(self._x - other._x)
        y = fabs(self._y - other._y)
        return sqrt(x**2 + y**2)



def placeNodes(dimension, numNodes):
    """Uniformly distributes numNodes over a dimension x dimension sized grid
    """

    nodes = []
    n = 0
    while n < numNodes:
        nodes.append(Node(randint(0, dimension-1), randint(0, dimension-1), n))
        n += 1

    return nodes


def createGraph(nodes, k, alpha, beta, aveDegree, maxLen):
    """Use edge probability to create adjacency list of graph

    The adjacency list is a dictionary keyed on the node IDs.
    Each element is a list of the neighbor nodes and the distances
    to the node stored as tuples.
    e.g. {3, [(5, 66.3), (1, 33.1), ...], ...}
    """

    adjList = {}
    for i in range(0, len(nodes)):
        adjList[i] = []
    numEdges = 0
    nodes.sort(lambda x,y: x._id - y._id)

    #bidirectional links, so only check pair (u,v) for all v > u
    for i in nodes:
        for j in nodes[i._id+1:]:

            dist = i.distance(j)
            prob = ((k*aveDegree *1.0 / len(nodes)) * beta *
                    exp(-dist *1.0/(alpha*maxLen)))
            if random() < prob:
                adjList[i._id].append( (j._id, dist) )
                adjList[j._id].append( (i._id, dist) )
                numEdges += 1
                
    return adjList, numEdges

def isConnected(adjList):
    """Do a BFS of the graph to make sure it is connected
    """

    numNodes = len(adjList.keys())
    visited = {}
    for i in range(0, numNodes):
        visited[i] = False
    queue = [0]
    visited[0] = True

    while len(queue) > 0:
        next = queue[0]
        queue = queue[1:]
        for n in adjList[next]:
            # tuples are (nodeID, dist), we just want nodeID, so n[0]
            if not visited[n[0]]:
                queue.append(n[0])
                visited[n[0]] = True

    for i in range(0, numNodes):
        if not visited[i]:
            return False
    return True

def printMatrix(adjList, out):
    """Write the graph as an adjacency matrix to file out
    """
    
    numNodes = len(adjList.keys())
    for i in range(0, numNodes):
        neighbors = adjList[i]
        neighbors.sort(lambda x, y: x[0] - y[0])
        for j in range(0, numNodes):
            if len(neighbors) > 0 and j == neighbors[0][0]:
                out.write("%d " % neighbors[0][1])
                neighbors = neighbors[1:]
            else:
                out.write("-1 ")
        out.write("\n")
    out.close()


if __name__ == '__main__':

    if len(argv) != 4:
        print "Usage: graph_gen.py <number of nodes> <average nodal degree>\
        <output file>"
        print "Example: graph_gen.py 75 3.5"
        exit(1)


    gridDimen = 100
    numNodes = int(argv[1])
    nodalDegree = float(argv[2])
    maxLen = gridDimen * sqrt(2)

    if numNodes > gridDimen**2:
        print "Max of %d nodes" % (gridDimen**2)
        exit(1)

 
    #parameters defined by Zegura et al. in 'How to model an internetwork'
    k = 27
    alpha = 0.3
    beta = 0.1

    count = 1
    while True:
        nodes = placeNodes(gridDimen, numNodes)
        adjList, numEdges = createGraph(nodes, k, alpha, beta,
                                        nodalDegree, maxLen)
        if isConnected(adjList):
            break
        print "Graph (%d) not connected...trying next random graph" % (count)
        count+=1
    
    printMatrix(adjList, open(argv[3], "w"))
    print "Graph generation complete. Written to '%s'" % (argv[3])
    print "Graph contains %d nodes, %d edges. Average nodal degree: %.2f"\
          % (numNodes, numEdges, 2.0*numEdges/numNodes)

