/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * NodeDg.cpp
 *
 *  Created on: Oct 21, 2022
 *      Author: Anthony B. Garza.
 *      Author: Hani Z. Girgis, PhD
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
template <typename V>
NodeDg<V>::NodeDg(V &aValue) : value(aValue)
{
}

template <typename V>
void NodeDg<V>::addEdge(std::shared_ptr<NodeDg<V>> aNode, double w)
{
    assert(edgeTable.count(aNode) == 0);
    edgeTable[aNode] = w;
}

template <typename V>
void NodeDg<V>::removeEdge(std::shared_ptr<NodeDg<V>> aNode)
{
    assert(edgeTable.count(aNode) == 1);
    auto itr = edgeTable.find(aNode);
    edgeTable.erase(itr);
}

template <typename V>
void NodeDg<V>::updateWeight(std::shared_ptr<NodeDg<V>> aNode, double w) {
    edgeTable.at(aNode) = w;
}

template <typename V>
double NodeDg<V>::retrieveWeight(std::shared_ptr<NodeDg<V>> aNode) const {
    assert(edgeTable.count(aNode) == 1);

    return edgeTable.at(aNode);
}

template <typename V>
bool NodeDg<V>::isConnected(std::shared_ptr<NodeDg<V>> aNode)
{
    return edgeTable.count(aNode) == 0 ? false : true;
}

template <typename V>
V& NodeDg<V>::getValue() const
{
    return value;
}


template <typename V>
std::unordered_map<std::shared_ptr<NodeDg<V>>, double> NodeDg<V>::getEdgeTable() const
{
    return edgeTable;
}

template <typename V>
int NodeDg<V>::getEdgeCount() const {
    return edgeTable.size();
}