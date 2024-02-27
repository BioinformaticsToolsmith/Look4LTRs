/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * DirectedGraph.h
 *
 *  Created on: Oct 21, 2022
 *      Author: Anthony B. Garza.
 *
 * Purpose: Directed graph that contains nodes with edges pointing to other nodes. All inner workings
 *          should be hidden from the user with methods to retrieve the needed values
 *
 * Academic use: Affero General Public License version 1.
 *
 * Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 *
 * Copyright (C) 2022 by the authors.
 */
#pragma once

#include "Node.h"

#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <stack>
#include <string>
#include <fstream>

template <typename V>
class DirectedGraph
{
private:
    /**
     * Variables
     */

    // If a node is allowed to have a node to itself
    bool allowCircular;

    // Helps with identifying subgraphs
    int subGraphIndex = 0;

    // In the form of: a pointer to Element -> a pointer of the corresponding NodeDg
    std::unordered_map<V *, std::shared_ptr<NodeDg<V>>> nodeTable;

    // In the form of: a pointer to Element -> a vector of pointers to nodes that this Elements points to
    std::unordered_map<V *, std::shared_ptr<std::vector<std::shared_ptr<NodeDg<V>>>>> forwardTable;

    // Table contains all backward connections to the node given by the value
    // Allows for quick deletion of a node and all backward connections to it
    // In the form of: a pointer to Element -> a vector of pointers to Nodes pointing back to the corresponding NodeDg of this Element.
    std::unordered_map<V *, std::shared_ptr<std::vector<std::shared_ptr<NodeDg<V>>>>> backwardTable;

    // Table that contains an index of each graph component and which nodes belong to it
    std::unordered_map<V *, int> subGraphIndexTable;

    // Table that contains a vector of each node belonging to each subgraph
    std::unordered_map<int, std::shared_ptr<std::unordered_set<V *>>> subGraphTable;

    /**
     * Methods
     */
    // Helper method; removes nodes pointed by a vector of values
    void removeValPtrs(std::vector<V *> &valPtrVec);

    // Helper method for asserting on a value
    void valCheck(V &aNode);

    V &getValue(std::shared_ptr<NodeDg<V>> aNodePtr);

public:
    // For backend of separating subgraphs. Returns a set pointer containing all of the connected components
    // If the break value is found, stop searching for the rest of the connected components and return.
    std::unordered_set<V *> *retrieveConnectedBreak(V &src, V &breakValue);

    void checkGraph();

public:
    DirectedGraph(bool allowCircular_ = false);
    ~DirectedGraph();

    /**
     * Getters
     */

    // Get number of nodes in the graph
    int getNodeCount() const;

    // Get a vector of pointers to all values in the graph
    std::vector<V *> getValueVec();

    // Get the number of subgraphs the graph can be broken into
    int getSubGraphCount() const;

    // Get the subgraph that the given node is a part of
    DirectedGraph<V> getSubGraph(V &aNode);


    /**
     * Methods
     */
    // adds a node holding a reference to the given value
    void addNode(V &value);

    // adds an edge between the source and destination node with the given weight
    void addEdge(V &src, V &dest, double weight = 1.);

    void addGraph(DirectedGraph<V> &aGraph);

    // Updates the weight of an edge from src to dest
    void updateWeight(V &src, V &dest, double weight);

    // Removes the node from the graph as well as all outgoing and incoming connection
    void removeNode(V &src);

    // Removes the edge from src to dest
    void removeEdge(V &src, V &dest);

    // Remove all of the connections of a node
    void clearConnections(V &src);

    // If the given graph shares nodes that have connections, copy the connections
    void copyConnections(DirectedGraph<V> &aGraph);

    // Retrieves the weight of the edge from src to destination
    double retrieveWeight(V &src, V &dest);

    // Retrieve pointers to values that are connected to the given node; either forward or backward
    std::vector<V *> retrieveConnectedValues(V &aNode, bool isForward = true);

    // Retrieve pointers to all values within the same subgraph
    std::vector<V *> retrieveSubGraphValueVec(V &aNode);

    // Is nodeOne part of the same subgraph as nodeTwo?
    int isSameSubGraph(V &nodeOne, V &nodeTwo);

    // Does the given node exist in the graph?
    bool isExist(V &aNode);

    // Is the src connected to destination (Asks for one way)
    bool isConnected(V &src, V &dest);

    // Retrieve all nodes that have no outgoing or incoming edges
    std::vector<V *> retrieveAlone();

    // Removes all nodes that have no outgoing or incoming edges
    void removeAlone();

    // Return all nodes that have no outgoing or incoming edges and then remove them
    std::vector<V *> popAlone();

    // Does the node have no outgoing or incoming edges?
    bool isAlone(V &aNode);

    // Remove all nodes that have no incoming edges
    void removeUnseen();

    // Return all nodes that have no incoming edges and then remove them
    std::vector<V *> popUnseen();

    // Retrieve all nodes that have no incoming edges
    std::vector<V *> retrieveUnseen();

    // Does the node have no incoming eges?
    bool isUnseen(V &aNode);

    // Do both nodes connect to each other?
    bool isBidirectional(V &nodeOne, V &nodeTwo);

    // Do both nodes have only one edge, both connecting to each other? (They can have other incoming edges)
    // A perfect couple that only see each other
    // O <--> O
    bool isCouple(V &nodeOne, V &nodeTwo);
    // Is the node in a couple relationship with some unknown node?
    bool isCouple(V &aNode);

    // Are the weights of the edges of the given node above the threshold
    // bool isAllAboveThreshold(V& aNode, double thresh);

    // Is the weight of the edge from src to dest above the threshold
    // bool isAboveThreshold(V& src, V& dest, double thresh);

    // Can the source node reach the destination by traveling?  If direction doesnt matter, pass false to isDirected
    bool isReachable(V &src, V &dest, bool isDirected = true);

    // Retrieve a connected component (all of the nodes that have a connection) by going forward (only outgoing connections) or backwards (incoming connections)
    std::unordered_set<V *> retrieveConnectedComponents(V &src, bool isForward = true);

    int retrieveSubGraphIndex(V &aVal) const;


    // Breaks the graph into smaller graphs
    std::vector<DirectedGraph<V>> breakGraph();

    // Writes to os; gives access to private variables
    void write(std::ostream &os) const;

    // A debuggin method
    void printSubGraphTable();

    void print(std::string filePath);

};


template <typename V>
std::ostream &operator<<(std::ostream &os, const DirectedGraph<V> &graph);

#include "DirectedGraph.cpp"