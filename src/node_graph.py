#!/usr/bin/env python

# Python program to print connected  
# components in an undirected graph 

# This code is contributed by Abhishek Valsan
# Updated by Saad Ansari for a directed graph application

import pdb

class Graph: 
      
    # init function to declare class variables 
    def __init__(self,V): 
        self.V = V 
        self.next = [[] for i in range(V)] 
        self.prev = [[] for i in range(V)] 
  
    def VisitNext(self, temp, v, visited): 
  
        # Visited this node
        visited[v] = True

        # Store the vertex to list 
        try:
            if temp[-1] != v:
                temp.append(v) 
        except:
            temp.append(v) 
  
        # Repeat for all vertices adjacent 
        # to this vertex v 
        for i in self.next[v]: 
            if visited[i] == False: 
                  
                # Update the list 
                temp= self.VisitNext(temp, i, visited) 

        return temp
  
    def VisitPrev(self, temp, v, visited): 
  
        # Visited this node
        visited[v] = True

        # Store the vertex to list 
        try:
            if temp[0] != v:
                temp.insert(0, v) 
        except:
            temp.insert(0, v) 
  
        # Repeat for all vertices adjacent 
        # to this vertex v 
        for i in self.prev[v]: 
            if visited[i] == False: 
                  
                # Update the list 
                temp= self.VisitPrev(temp, i, visited) 

        return temp

    # method to add an directed edge 
    def addEdge(self, v, w): 
        self.next[v].append(w) 
        self.prev[w].append(v) 
  
    # Method to retrieve connected components 
    # in a directed graph 
    def connectedComponents(self): 
        visited = [] 
        cc = [] 
        for i in range(self.V): 
            visited.append(False) 
        for v in range(self.V): 
            if visited[v] == False: 

                temp = [] 
                temp= self.VisitNext( temp, v, visited) 
                cc.append( self.VisitPrev( temp, v, visited) )

        return cc 
  
# Driver Code 
if __name__=="__main__": 
      
    # Create a graph given in the above diagram 
    # 5 vertices numbered from 0 to 4 
    g = Graph(7); 
    g.addEdge(1, 0) 
    g.addEdge(3, 4) 
    g.addEdge(0, 6) 
    g.addEdge(5, 1) 
    cc = g.connectedComponents() 
    print("Following are connected components") 
    print(cc) 
  
