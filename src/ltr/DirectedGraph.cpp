/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * DirectedGraph.cpp
 *
 *  Created on: Oct 21, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer: Hani Z. Girgis.
 *
 * Purpose: NodeDg in a graph; contains a value and a vector of edges; templatized to hold whatever
 *          value is needed
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

#include "DirectedGraph.h"

template <typename V>
DirectedGraph<V>::DirectedGraph(bool allowCircular_)
{
    allowCircular = allowCircular_;
}

// [ok]
template <typename V>
DirectedGraph<V>::~DirectedGraph()
{

    // for (auto i : backwardTable)
    // {
    //     delete i.second;
    // }

    // for (auto const [_, setPtr] : subGraphTable)
    // {
    //     delete setPtr;
    // }
}

// [ok]
template <typename V>
void DirectedGraph<V>::addNode(V &value)
{
    assert(nodeTable.count(&value) == 0);
    assert(forwardTable.count(&value) == 0);
    assert(backwardTable.count(&value) == 0);
    assert(subGraphIndexTable.count(&value) == 0);

    nodeTable[&value] = std::make_shared<NodeDg<V>>(value);
    forwardTable[&value] = std::make_shared<std::vector<std::shared_ptr<NodeDg<V>>>>();
    backwardTable[&value] = std::make_shared<std::vector<std::shared_ptr<NodeDg<V>>>>();
    subGraphIndexTable[&value] = subGraphIndex;

    subGraphTable[subGraphIndex] = std::make_shared<std::unordered_set<V *>>();
    subGraphTable[subGraphIndex]->insert(&value);

    subGraphIndex++;
}

// [ok]
template <typename V>
void DirectedGraph<V>::addEdge(V &src, V &dst, double weight)
{
    // Preconditions
    valCheck(src);
    valCheck(dst);
    assert(weight >= 0.0);

    if (!allowCircular)
    {
        assert(&src != &dst);
    }

    auto srcPtr = nodeTable[&src];
    auto dstPtr = nodeTable[&dst];
    srcPtr->addEdge(dstPtr, weight);
    forwardTable[&src]->push_back(dstPtr);
    backwardTable[&dst]->push_back(srcPtr);

    int srcIndex = subGraphIndexTable[&src];
    int destIndex = subGraphIndexTable[&dst];

    if (srcIndex != destIndex) // They are not part of the same subgraph
    {
        for (auto &ptr : *subGraphTable.at(destIndex))
        {
            subGraphIndexTable.at(ptr) = srcIndex;
            subGraphTable.at(srcIndex)->insert(ptr);
        }
        subGraphTable.erase(subGraphTable.find(destIndex));
    }
}

/**
 * Purpose: This method is usually called while constructing a subgraph
 */
// [OK]
template <typename V>
void DirectedGraph<V>::copyConnections(DirectedGraph<V> &aGraph)
{
    // Ask if node exists in aGraph;
    for (auto const [key, ptr] : nodeTable)
    {
        if (aGraph.isExist(*key))
        {
            // Check if connected nodes of aGraph exist in this graph; if so, add edge
            for (auto val : aGraph.retrieveConnectedValues(*key))
            {
                if (nodeTable.count(val) == 1)
                {
                    addEdge(*key, *val, aGraph.retrieveWeight(*key, *val));
                }
            }
        }
    }
}

// [OK]
/**
 * Merge a graph with this graph
 */
template <typename V>
void DirectedGraph<V>::addGraph(DirectedGraph<V> &aGraph)
{
    auto valueVec = aGraph.getValueVec();
    if (valueVec.size() > 0)
    {
        for (auto &valuePtr : valueVec)
        {
            if (!isExist(*valuePtr))
            {
                addNode(*valuePtr);
            }
        }

        copyConnections(aGraph);
    }
}

// [OK]
template <typename V>
void DirectedGraph<V>::updateWeight(V &src, V &dst, double weight)
{
    valCheck(src);
    valCheck(dst);
    assert(weight >= 0); // 0 means a vertial connection, i.e., the same segment is found in the forward and backward scores.

    auto srcPtr = nodeTable.at(&src);
    auto destPtr = nodeTable.at(&dst);
    assert(srcPtr->isConnected(destPtr));

    srcPtr->updateWeight(destPtr, weight);
}

// [OK]
template <typename V>
void DirectedGraph<V>::removeNode(V &aValue)
{
    valCheck(aValue);

    // Remove all connections
    clearConnections(aValue);

    // Make sure that this node makes its own one-node component
    assert(subGraphTable.at(subGraphIndexTable.at(&aValue))->size() == 1);

    nodeTable.erase(nodeTable.find(&aValue));
    forwardTable.erase(forwardTable.find(&aValue));
    backwardTable.erase(backwardTable.find(&aValue));
    subGraphTable.erase(subGraphIndexTable.at(&aValue));
    subGraphIndexTable.erase(&aValue);

    checkGraph();
}

// [OK]
template <typename V>
void DirectedGraph<V>::removeEdge(V &src, V &dst)
{
    valCheck(src);
    valCheck(dst);
    assert(isConnected(src, dst));

    nodeTable.at(&src)->removeEdge(nodeTable.at(&dst));

    // The fan-out connections of the src
    auto forwardVecPtr = forwardTable.at(&src);
    // The fan-in connections of the dst
    auto backVecPtr = backwardTable.at(&dst);

    forwardVecPtr->erase(std::find(forwardVecPtr->begin(), forwardVecPtr->end(), nodeTable.at(&dst)));
    backVecPtr->erase(std::find(backVecPtr->begin(), backVecPtr->end(), nodeTable.at(&src)));

    // connectedSet may (i) include the src and the dst or (ii) the src -- not the dst -- and its complete subgraph
    auto connectedSet = retrieveConnectedBreak(src, dst);
    if (connectedSet->count(&dst) == 0)
    {
        auto srcSet = subGraphTable[subGraphIndexTable[&src]];
        // Calculate set difference
        for (auto valuePtr : *connectedSet)
        {
            subGraphIndexTable[valuePtr] = subGraphIndex;
            srcSet->erase(valuePtr);
        }
        subGraphTable[subGraphIndex] = std::make_shared<std::unordered_set<V *>>(*connectedSet);

        subGraphIndex++;
    }
}

// [OK]
template <typename V>
void DirectedGraph<V>::clearConnections(V &src)
{
    valCheck(src);
    // Remove all fan-out connections
    for (auto &[fanOutPtr, _] : nodeTable.at(&src)->getEdgeTable())
    {
        removeEdge(src, fanOutPtr->getValue());
    }

    // Remove all fan-in connections
    auto fanInVec = backwardTable.at(&src);
    while (!fanInVec->empty())
    {
        removeEdge(fanInVec->at(0)->getValue(), src);
    }
}

// [OK]
template <typename V>
double DirectedGraph<V>::retrieveWeight(V &src, V &dst)
{
    valCheck(src);
    valCheck(dst);
    assert(nodeTable.at(&src)->isConnected(nodeTable.at(&dst)));

    return nodeTable.at(&src)->retrieveWeight(nodeTable.at(&dst));
}

// [OK]
template <typename V>
std::vector<V *> DirectedGraph<V>::retrieveConnectedValues(V &aValue, bool isForward)
{
    valCheck(aValue);

    std::vector<V *> r;

    auto vec = isForward == true ? forwardTable[&aValue] : backwardTable[&aValue];
    for (auto nodePtr : *vec)
    {
        r.push_back(&nodePtr->getValue());
    }

    return r;
}

// [OK]
/**
 * Get all singleton values (i.e., no fan-out and no fan-in edges)
 */
template <typename V>
std::vector<V *> DirectedGraph<V>::retrieveAlone()
{
    std::vector<V *> r;
    for (auto [valuePtr, nodePtr] : nodeTable)
    {
        if (isAlone(*valuePtr))
        {
            r.push_back(valuePtr);
        }
    }

    return r;
}

// [OK]
template <typename V>
void DirectedGraph<V>::removeAlone()
{
    removeValPtrs(retrieveAlone());
}

// [OK]
template <typename V>
std::vector<V *> DirectedGraph<V>::popAlone()
{
    std::vector<V *> r = retrieveAlone();
    removeValPtrs(r);
    return r;
}

// [OK]
template <typename V>
bool DirectedGraph<V>::isAlone(V &aValue)
{
    valCheck(aValue);

    return (nodeTable.at(&aValue)->getEdgeCount() == 0 && backwardTable.at(&aValue)->size() == 0) ? true : false;
}

// [OK]
template <typename V>
bool DirectedGraph<V>::isBidirectional(V &nodeOne, V &nodeTwo)
{
    valCheck(nodeOne);
    valCheck(nodeTwo);

    return isConnected(nodeOne, nodeTwo) && isConnected(nodeTwo, nodeOne) ? true : false;
}

// [OK]
template <typename V>
std::vector<V *> DirectedGraph<V>::retrieveSubGraphValueVec(V &aValue)
{
    valCheck(aValue);

    auto aSet = subGraphTable.at(subGraphIndexTable.at(&aValue));
    std::vector<V *> r{aSet->begin(), aSet->end()};

    return r;
}

// [ok]
template <typename V>
int DirectedGraph<V>::isSameSubGraph(V &valueOne, V &valueTwo)
{
    valCheck(valueOne);
    valCheck(valueTwo);

    return subGraphIndexTable.at(&valueOne) == subGraphIndexTable.at(&valueTwo) ? true : false;
}

// [ok]
template <typename V>
bool DirectedGraph<V>::isExist(V &aNode)
{
    return nodeTable.count(&aNode) == 1 ? true : false;
}

// [OK]
template <typename V>
bool DirectedGraph<V>::isConnected(V &src, V &dst)
{
    valCheck(src);
    valCheck(dst);

    return nodeTable.at(&src)->isConnected(nodeTable.at(&dst));
}

// [OK]
/**
 *
 * If src and breakValue are connected, it will stop prematurely.
 * If src and breakValue are unconnected, it will return a whole subgraph including the src.
 *
 */
template <typename V>
std::unordered_set<V *> *DirectedGraph<V>::retrieveConnectedBreak(V &src, V &breakValue)
{
    valCheck(src);
    valCheck(breakValue);

    auto r = new std::unordered_set<V *>;
    if (&src == &breakValue || retrieveSubGraphValueVec(src).size() == 1)
    {
        r->insert(&src);
    }
    else
    {
        std::stack<V *> unvisited;
        unvisited.push(&src);

        do
        {
            V *top = unvisited.top();
            r->insert(top);
            unvisited.pop();
            // Visit forward connectoins (fan-out nodes)

            bool breakOuter = false;

            for (auto valPtr : retrieveConnectedValues(*top, true))
            {
                if (valPtr == &breakValue)
                {
                    r->insert(valPtr);
                    breakOuter = true;
                    break;
                }
                // If it is not marked visited in the set
                if (r->count(valPtr) == 0)
                {
                    // Push it to the stack
                    unvisited.push(valPtr);
                }
            }

            if (breakOuter)
            {
                break;
            }

            // Visit backwaord connections (fan-in nodes)
            for (auto valPtr : retrieveConnectedValues(*top, false))
            {
                if (valPtr == &breakValue)
                {
                    r->insert(valPtr);
                    breakOuter = true;
                    break;
                }

                if (r->count(valPtr) == 0)
                {
                    unvisited.push(valPtr);
                }
            }

            if (breakOuter)
            {
                break;
            }
        } while (!unvisited.empty());
    }

    return r;
}

// [OK]
template <typename V>
bool DirectedGraph<V>::isReachable(V &src, V &dst, bool isDirected)
{
    // Using DFS algorithm

    valCheck(src);
    valCheck(dst);

    bool r = false;

    if (&src == &dst)
    {
        r = true;
    }
    else if (isSameSubGraph(src, dst))
    {
        // If the nodes aren't part of the same subgraph, then it is impossible for them to be reachable
        // But being part of the same subgraph does not imply that they are reachable.
        if (!isDirected)
        {
            r = true;
        }
        else
        {
            std::unordered_set<V *> visited;
            std::stack<V *> unvisited;
            unvisited.push(&src);

            do
            {
                // visiting top
                V *top = unvisited.top();
                visited.insert(top);
                unvisited.pop();

                for (auto fanOutPtr : retrieveConnectedValues(*top))
                {
                    if (fanOutPtr == &dst)
                    {
                        r = true;
                        break;
                    }

                    if (visited.count(fanOutPtr) == 0)
                    {
                        unvisited.push(fanOutPtr);
                    }
                }

                if (r == true)
                {
                    break;
                }

            } while (!unvisited.empty());
        }
    }

    return r;
}

// [OK]
/**
 * Find all connected components in this graph
 */
template <typename V>
std::vector<DirectedGraph<V>> DirectedGraph<V>::breakGraph()
{
    std::vector<DirectedGraph<V>> r;
    for (auto const &[_, nodeSet] : subGraphTable)
    {
        DirectedGraph<V> g;
        for (auto aValuePtr : *nodeSet)
        {
            g.addNode(*aValuePtr);
        }
        g.copyConnections(*this);
        r.push_back(g);
    }

    return r;
}

// [OK]
template <typename V>
void DirectedGraph<V>::removeValPtrs(std::vector<V *> &valPtrVec)
{
    for (auto &key : valPtrVec)
    {
        removeNode(*key);
    }
}

// [OK]
template <typename V>
void DirectedGraph<V>::valCheck(V &aValue)
{
    assert(nodeTable.count(&aValue) == 1);
    assert(forwardTable.count(&aValue) == 1);
    assert(backwardTable.count(&aValue) == 1);
    assert(subGraphIndexTable.count(&aValue) == 1);
    assert(subGraphTable.count(subGraphIndexTable.at(&aValue)) == 1);
    assert(subGraphTable.at(subGraphIndexTable.at(&aValue))->count(&aValue) == 1);
}

// [OK]
template <typename V>
std::vector<V *> DirectedGraph<V>::getValueVec()
{
    std::vector<V *> r;
    for (auto const &[valuePtr, _] : nodeTable)
    {
        r.push_back(valuePtr);
    }

    return r;
}

// [OK]
// Use call after each test case
// Put in private afterwards
template <typename V>
void DirectedGraph<V>::checkGraph()
{

    // 1: Check that number of nodes in the graph are the same for every table.
    int size = nodeTable.size();
    assert(backwardTable.size() == size);
    assert(forwardTable.size() == size);
    assert(subGraphIndexTable.size() == size);
    int sum = 0;
    for (auto const &[_, setPtr] : subGraphTable)
    {
        sum += setPtr->size();
    }
    assert(sum == size);

    for (auto const &[valPtr, _] : nodeTable)
    {
        assert(backwardTable.count(valPtr) == 1);
        assert(forwardTable.count(valPtr) == 1);
        assert(subGraphIndexTable.count(valPtr) == 1);
    }

    for (auto const &[_, setPtr] : subGraphTable)
    {
        for (auto const valuePtr : *setPtr)
        {
            assert(nodeTable.count(valuePtr) == 1);
        }
    }
}

// [OK]
template <typename V>
int DirectedGraph<V>::getNodeCount() const
{
    return backwardTable.size();
}

// [OK]
template <typename V>
int DirectedGraph<V>::getSubGraphCount() const
{
    return subGraphTable.size();
}

// [OK]
/**
 * Get the connected component including the input node
 */
template <typename V>
DirectedGraph<V> DirectedGraph<V>::getSubGraph(V &aNode)
{
    valCheck(aNode);

    DirectedGraph<V> r;

    // Getting all elements pointers from the subgraph, adding to result graph
    for (auto ptr : *subGraphTable.at(subGraphIndexTable.at(&aNode)))
    {
        r.addNode(*ptr);
    }

    r.copyConnections(*this);

    return r;
}

template <typename V>
int DirectedGraph<V>::retrieveSubGraphIndex(V &aVal) const
{
    valCheck(aVal);
    return subGraphIndexTable.at(&aVal);
}

// [OK]
/**
 * Print adjacency lists for each node, i.e., an edge and its weight
 */
template <typename V>
void DirectedGraph<V>::write(std::ostream &os) const
{
    for (auto [key, nodePtr] : nodeTable)
    {
        os << "> " << nodePtr->getValue() << " --> ";
        for (auto [destPtr, weight] : nodePtr->getEdgeTable())
        {
            os << "{" << destPtr->getValue() << " W: " << weight << "}, ";
        }
        os << std::endl;
    }
}

// [OK]
/**
 * index -> {node_1, node_2, .... , node_n}
 */
template <typename V>
void DirectedGraph<V>::printSubGraphTable()
{
    for (auto const [anIndex, aSet] : subGraphTable)
    {
        std::cout << anIndex << ": ";
        for (auto const aVal : *aSet)
        {
            std::cout << *aVal << " ";
        }
        std::cout << std::endl;
    }
}

// [OK]
template <typename V>
std::ostream &operator<<(std::ostream &os, const DirectedGraph<V> &graph)
{
    graph.write(os);
    return os;
}

template <typename V>
void DirectedGraph<V>::print(std::string filePath)
{
    std::cout << filePath << std::endl;
    std::ofstream out(filePath);
    out << *this;
    out.close();
}