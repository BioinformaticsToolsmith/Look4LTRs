/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * NodeDg.h
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
#pragma once

#include <iostream>
#include <assert.h>
#include <unordered_map>
#include <memory>

template <typename V>
class NodeDg
{
public:
    /**
     * Constructor
     */
    // store a value in node to access.
    NodeDg(V &aValue);

    /**
     * Getters
     */
    V& getValue() const;

    std::unordered_map<std::shared_ptr<NodeDg<V>>, double> getEdgeTable() const;
    int getEdgeCount() const;

    /**
     * Methods
     */
    void addEdge(std::shared_ptr<NodeDg<V>> aNode, double w);
    void removeEdge(std::shared_ptr<NodeDg<V>> aNode);
    void updateWeight(std::shared_ptr<NodeDg<V>> aNode, double w);
    double retrieveWeight(std::shared_ptr<NodeDg<V>> aNode) const;

    bool isConnected(std::shared_ptr<NodeDg<V>> aNode);

private:
    /**
     * Variables
     */
    V &value;

    std::unordered_map<std::shared_ptr<NodeDg<V>>, double> edgeTable;
};

#include "Node.cpp"