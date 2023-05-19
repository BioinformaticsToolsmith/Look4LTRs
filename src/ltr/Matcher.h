/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Matcher.h
 *
 *  Created on: Oct 19, 2022
 *      Author: Anthony B. Garza.
 *    Reviewer: 
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

#include "RT.h"
#include "RTSolo.h"
#include "Element.h"
#include "Stretch.h"
#include "DirectedGraph.h"
#include "Node.h"
#include "LtrParameters.h"

#include "CaseMatcher.h"
#include "CaseSingle.h"
#include "CaseRecent.h"
#include "CaseSolo.h"
#include "CaseRecentComplex.h"

#include "../IdentityCalculator.h"
#include "../red/Red.h"
#include "../utility/LocAlign.h"

#include <vector>
#include <set>
#include <iterator>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <tuple>
#include <string>
#include <stack>
#include <unordered_set>

class Matcher
{
private:


    /**
     * Variables
     */
    // graph nodes are references to the Elements given in the constructor
    DirectedGraph<Element> graph;

    // A vector of all connected sub components
    std::vector<DirectedGraph<Element>> subGraphVec;

    // A map of all connected sub components: key (RT) -> value (subgraph)
    std::unordered_map<RT*, DirectedGraph<Element>*> rtSubGraphMap;

    // Passed in during constructor; forward and backward elements
    std::vector<Element> &fElementVec;
    std::vector<Element> &bElementVec;

    Red &red;
    IdentityCalculator<int32_t> &ic;
    IdentityCalculator<int32_t> &icRecent;

    // The sequence of a chromosome
    const std::string *seq;

    // The result vector, which includes complete RTs or solo RTs with matched elements
    std::vector<RT *> rtVec;

    // A result vector that contains the starts and ends of complex regions
    std::vector<RT *> complexVec;
    
    // Overlapping regions in the forward and the backword elements have been matched
    bool isVertical;
    
    // The chromosome graph has been broken to connected subcomponents
    bool isBroken;

    /**
     * Methods
     */

    // Check validity of passed in arguments
    void checkValid();

    // Check that elements of a ordered vector of elements dont merge
    void checkNoMerge(std::vector<Element> &eVec);

    // Add elements to graph
    void addElementVec(std::vector<Element> &elementVec);

    // Filter out MITEs using local alignment
    void filterMites(std::vector<Element> &elementVec);

    // Match forward and backward elements with each other
    void matchElements(DirectedGraph<Element> &g);

    // Not currently used
    void matchElementsForward(DirectedGraph<Element> &g);

    // Not corrently used
    void matchElementsBackward(DirectedGraph<Element> &g);

    // Helper method; adds an edge from src to dest if no edge exists; otherwise, update the weight of the existing edge
    // with the size of the current matching stretch and the previously stored stretch size retrieved from sizeTable
    void addEdge(Element &src, Element &dest, Stretch &matchStretch, std::unordered_map<Element *, int> &sizeTable, DirectedGraph<Element> &g);

    // Remove weak connections
    void filterMatches();

    // Not currently used; not very useful
    void mergeCrosses();

    // Finds forward elements that overlap with backward elements, i.e., elements at the same exact positions but one
    // is pointing forward and the other backwards; build edge between them with weight 0
    void matchOverlaps();


    std::vector<Element*> retrieveVertCons(Element &father);

    // Merges elements together if they overlap the same element.
    void mergeVerticalOverlaps(std::vector<Element*> &eleVec);

    void validateVerticals();

    // Simple method that calls break on the graph and assigns the vector of graphs to subGraphVec
    void breakGraph(); 

    // Creating our candidate graphs into further graphs! Graph-ception!
    // void breakCandidates();

    // Call multiple analysis cases
    void runCases();

    // Helper method ???????????
    bool checkRank(RT *r1, RT *r2);

    bool isSameRT(RT *r1, RT *r2);

public:
    /**
     * Constructor
     * 
     */
    Matcher(std::vector<Element> &_fElementVec, std::vector<Element> &_bElementVec, Red &_red, 
            IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent, const std::string *_seq);

    // Destructor; destroys the RTs in the RtVec
    ~Matcher();

    /**
     * Getters
     */
    // Retrieves a sorted (based on the start of the candidates) vector of pointers to the elements belonging to a graph
    std::tuple<std::vector<Element *>, std::vector<Element *>> retrieveForwardBackward(DirectedGraph<Element> &g);
    std::tuple<std::vector<Element *>, std::vector<Element *>> retrieveForwardBackward();
    std::vector<Element *> retrieveElementVec(bool isForward, DirectedGraph<Element> &g);

    DirectedGraph<Element> *getGraph() const;
    std::vector<DirectedGraph<Element>> getSubGraphVec() const;
    DirectedGraph<Element> *getSubGraph(RT *rt) const;

    // Get the results
    std::vector<RT *> *getRtVec();
    std::vector<RT *> *getComplexVec();
};

// Not used at this time
// bool compareElementGraph(DirectedGraph<Element> &a, DirectedGraph<Element> &b);

// 
bool compareValuePtrs(Element *a, Element *b);
