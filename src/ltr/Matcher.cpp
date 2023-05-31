/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 *
 * Matcher.cpp
 *
 *  Created on: Oct 19, 2022
 *      Author: Anthony B. Garza.
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
#include "Matcher.h"

Matcher::Matcher(std::vector<Element> &_fElementVec, std::vector<Element> &_bElementVec, Red &_red, IdentityCalculator<int32_t> &_ic, IdentityCalculator<int32_t> &_icRecent,  const std::string *_seq)
    : fElementVec(_fElementVec), bElementVec(_bElementVec), red(_red), ic(_ic), icRecent(_icRecent), seq(_seq){

    checkValid();

    // Vertical connections have not been added yet    
    isVertical = false;

    // Adding element vectors to directed graph; every node contains a pointer to the element
    addElementVec(fElementVec);
    addElementVec(bElementVec);

    //filterMites(fElementVec);
    //filterMites(bElementVec);

    // Matching Elements
    matchElements(graph);

    // Filter matches based on low weights
    // filterMatches();

    // Matching overlaps with weight of 0 for edge
    // [OK]
    matchOverlaps();

    // Merge overlapping elements
    auto fEleVec = retrieveElementVec(true, graph);
    mergeVerticalOverlaps(fEleVec);
    auto bEleVec = retrieveElementVec(false, graph);
    mergeVerticalOverlaps(bEleVec);

    // Match elements again after merging vertical connections
    matchElements(graph);

    // Validate that each node only has one vertical connection
    validateVerticals();

    //graph.print("/home/transposons/Data/Look4LTRs/Graphs/" + *filepath);

    // Break the graph into subgraphs with only matches
    breakGraph();

    // Creating candidate graphs (breaking up possible hyperextensions)
    runCases();

    rtVec.shrink_to_fit();
    complexVec.shrink_to_fit();
    // std::sort(rtVec.begin(), rtVec.end(), [](RT *r1, RT *r2)
    //         { return r1->getStart() < r2->getStart(); });

}

Matcher::~Matcher() {
    // destroying RTs assigned on heap that weren't destroyed later as well as any elements assigned on heap.
    // for (auto r : rtVec) {
    //     if (r != nullptr) {
    //         delete r;
    //     }
    // }
}

void Matcher::checkValid() {
    checkNoMerge(fElementVec);
    checkNoMerge(bElementVec);
}

void Matcher::checkNoMerge(std::vector<Element> &eVec) {
    for (int i = 1; i < eVec.size(); i++) {
        assert(eVec[i].getStart() >= eVec[i - 1].getEnd());
    }
}


// [OK]
void Matcher::addElementVec(std::vector<Element> &elementVec)
{
    for (int i = 0; i < elementVec.size(); i++)
    {
        graph.addNode(elementVec[i]);
    }
}


/**
 * Filters out MITEs from the graph. These MITEs are transposable elements with inverted repeats.
 * Local alignment will be used to detect them.
*/
void Matcher::filterMites(std::vector<Element> &elementVec) {
    for (auto &ele : elementVec) {
        auto start = ele.getStart();
        auto end = ele.getEnd();

        // If the element can't hold at least two inverted repeats, then it is not a MITE
        if (ele.getSize() < LtrParameters::MIN_MITE_TIR_SIZE * 2 || ele.getSize() > LtrParameters::MAX_MITE_SIZE) {
            continue;
        }

        std::string leftTirSeq;
        std::string rightTirSeq;
        if (ele.getSize() > LtrParameters::TIR_SEARCH_RANGE * 2) {
            leftTirSeq = seq->substr(start, LtrParameters::TIR_SEARCH_RANGE);
            rightTirSeq = seq->substr(end - LtrParameters::TIR_SEARCH_RANGE, LtrParameters::TIR_SEARCH_RANGE);
        }
        else {
            int half = ele.getSize()/2;
            leftTirSeq = seq->substr(start, half);
            rightTirSeq = seq->substr(end - half, half);
        }

        // Reverse complement the right TIR sequence. This is done to align the left TIR sequence with the right TIR sequence.
        std::string rightTirSeqRev = LtrUtility::reverseComplement(rightTirSeq);

 
        // Perform local alignment
        // Arguments are as follow: sequence 1, start, end, sequence 2, start, end, match, mismatch, gap open, gap continue
        // leftTirSeq and rightTirSeqRev need to be converted to a const char * for LocAlign
        int length = LtrUtility::getAlignmentLength(leftTirSeq, rightTirSeqRev, 2, -3, -5, -2);
        
        // If the alignment length is greater than LtrParameters::MIN_MITE_TIR_ALIGN, then this is a MITE
        if (length >= LtrParameters::MIN_MITE_TIR_SIZE) {
            std::cout << "MITE found at " << start << " " << end << " " << length << std::endl;
            // Remove the node from the graph
            graph.removeNode(ele);
        }

    }
}



// [OK]
void Matcher::matchElements(DirectedGraph<Element> &g)
{
    auto [fEleVec, bEleVec] = retrieveForwardBackward(g);

    int bIndex = 0;
    // Looping through all forward elements
    for (auto fElement : fEleVec)
    {
        // Looping through until you find the begining of the potential matching region in the backward vector
        for (int i = bIndex; i < bEleVec.size(); i++)
        {
            // Find the first backward element after the start of the forward element, loop should start here
            if (bEleVec.at(i)->getStart() >= fElement->getStart())
            {
                bIndex = i; // Should it be i+1?
                break;
            }
        }
            
        // Contains the size of the current forward element connection to the backward element
        // Key is the backward element being pointed to, value is the size
        // Element ptr backward (or forward) -> total size of stretches overlapping with forward (or backward)
        std::unordered_map<Element *, int> sizeTable;

        // Looping through all hypothetical stretches that all forward stretches are pointing to
        for (auto fStretchPtr : fElement->getStretchVec())
        {
            // Where this element should be matching (a possible match), i.e., the matching region
            Stretch mStretch = fStretchPtr->buildMatch();
            // Looping through backward elements to find potential matches
            for (int i = bIndex; i < bEleVec.size(); i++)
            {
                // if the mStretch is before the backward element, then we should break
                if (mStretch.getEnd() <= bEleVec[i]->getStart())
                {
                    break;
                }

                // Otherwise, add edge if there is overlap
                if (bEleVec[i]->calcOverlap(mStretch) > 0 && 
                        (!graph.isConnected(*fElement, *bEleVec.at(i)) || 
                        !LtrUtility::isEqual(graph.retrieveWeight(*fElement, *bEleVec.at(i)), 0.0)))
                {
                    addEdge(*fElement, *bEleVec.at(i), mStretch, sizeTable, g);
                }
            }
        }

        // Creating edges from backward matches to forward
        std::vector<Element*> bElementPtrVec;
        for (auto const &[bElementPtr, _] : sizeTable)
        {
            bElementPtrVec.push_back(bElementPtr);
        }
        for (auto bElementPtr : bElementPtrVec) 
        {
            sizeTable[fElement] = 0;
            for (auto bStretchPtr : bElementPtr->getStretchVec())
            {
                Stretch bMatchStretch = bStretchPtr->buildMatch();

                if (fElement->calcOverlap(bMatchStretch) > 0 && 
                    (!graph.isConnected(*bElementPtr, *fElement) || 
                    !LtrUtility::isEqual(graph.retrieveWeight(*bElementPtr, *fElement), 0.0)))
                {
                    addEdge(*bElementPtr, *fElement, bMatchStretch, sizeTable, g);
                }
            }
        }
    }
}

// Not currently used
// Match elements forward
void Matcher::matchElementsForward(DirectedGraph<Element> &g) {

    auto [fEleVec, bEleVec] = retrieveForwardBackward(g);

    int bIndex = 0;
    // Looping through all forward elements
    for (auto &fElement : fEleVec)
    {
        // If index (in the bElementVec) is assigned then we can start looking at the backward elements starting at bIndex
        bool isIndexAssigned = false;

        // Contains the size of the current forward element connection to the backward element
        // Key is the backward element being pointed to, value is the size
        // Element ptr backward (or forward) -> total size of stretches overlapping with forward (or backward)
        std::unordered_map<Element *, int> sizeTable;

        // Looping through all hypothetical stretches that all forward stretches are pointing to
        for (auto fStretchPtr : fElement->getStretchVec())
        {
            // Where this element should be matching (a possible match), i.e., the matching region
            Stretch mStretch = fStretchPtr->buildMatch();

            // Looping through until you find the begining of the potential matching region
            if (!isIndexAssigned)
            {
                for (int i = bIndex; i < bEleVec.size(); i++)
                {
                    // Find the first backward element after the start of the forward element, loop should start here
                    if (bEleVec[i]->getStart() >= fElement->getStart())
                    {
                        bIndex = i; // Should it be i+1?
                        isIndexAssigned = true;
                        break;
                    }
                }
            }

            // Looping through backward elements to find potential matches
            for (int i = bIndex; i < bEleVec.size(); i++)
            {
                // if the mStretch is before the backward element, then we should break
                if (mStretch.getEnd() <= bEleVec[i]->getStart())
                {
                    break;
                }

                // Otherwise, add edge if there is overlap
                if (bElementVec[i].calcOverlap(mStretch) > 0)
                {
                    addEdge(*fElement, *bEleVec[i], mStretch, sizeTable, g);
                }
            }
        }
    }
}

// Not currently used
// Match elements backward
void Matcher::matchElementsBackward(DirectedGraph<Element> &g) {

    auto [fEleVec, bEleVec] = retrieveForwardBackward(g);

    int fIndex = fEleVec.size() - 1;
    // Looping through all backward elements
    for (auto bElement = bElementVec.rbegin(); bElement != bElementVec.rend(); ++bElement)
    {
        // If index (in the bElementVec) is assigned then we can start looking at the backward elements starting at bIndex
        bool isIndexAssigned = false;

        // Contains the size of the current forward element connection to the backward element
        // Key is the backward element being pointed to, value is the size
        // Element ptr backward (or forward) -> total size of stretches overlapping with forward (or backward)
        std::unordered_map<Element *, int> sizeTable;

        // Looping through all hypothetical stretches that all forward stretches are pointing to
        for (auto bStretchPtr : bElement->getStretchVec())
        {
            // Where this element should be matching (a possible match), i.e., the matching region
            Stretch mStretch = bStretchPtr->buildMatch();

            // Looping through until you find the begining of the potential matching region
            if (!isIndexAssigned)
            {
                for (int i = fIndex; i < fEleVec.size(); i--)
                {
                    if (fEleVec[i]->getStart() >= bElement->getEnd())
                    {
                        fIndex = i - 1;
                        isIndexAssigned = true;
                        break;
                    }
                }
            }

            // Looping through backward elements to find potential matches
            for (int i = fIndex; i < fEleVec.size(); i--)
            {
                // if the mStretch is before the backward element, then we should break
                if (mStretch.getStart() >= fEleVec[i]->getEnd())
                {
                    break;
                }

                // Otherwise, add edge if there is overlap
                if (fEleVec[i]->calcOverlap(mStretch) > 0)
                {
                    addEdge(*bElement, *fEleVec[i], mStretch, sizeTable, g);
                }
            }
        }
    }
}

// [OK]
/**
 * sizeTable: backward element ptr -> size of overlap with a forward element
 */
void Matcher::addEdge(Element &src, Element &dest, Stretch &matchStretch, std::unordered_map<Element *, int> &sizeTable, DirectedGraph<Element> &g)
{
    sizeTable[&dest] += dest.calcOverlap(matchStretch);
    double w = sizeTable[&dest] / static_cast<double>(src.getSize());
    if (g.isConnected(src, dest))
    {
        g.updateWeight(src, dest, w);
    }
    else
    {
        g.addEdge(src, dest, w);
    }
}

/**
 * This method is not working
 * Assumption: Vertical connections have not been added 
 */
// void Matcher::filterMatches() {
    
//     assert(!isVertical);

//     for (auto ele : graph.getValueVec()) {
//         // fan-out connectoins
//         auto fConVec = graph.retrieveConnectedValues(*ele);

//         for (auto fEle : fConVec) {
//             double fWeight = graph.retrieveWeight(*ele, *fEle);
//             if (fWeight < LtrParameters::MIN_WEIGHT) {
//                 graph.removeEdge(*ele, *fEle);
//             }
//         }

//     }
// }

void Matcher::filterMatches() {
    for (auto ele : graph.getValueVec()) {
        auto fConVec = graph.retrieveConnectedValues(*ele);

        for (auto fEle : fConVec) {
            double fWeight = graph.retrieveWeight(*ele, *fEle);
            if (fWeight != 0.0 && fWeight < LtrParameters::MIN_WEIGHT) {
                if (!graph.isConnected(*fEle, *ele)) {
                    graph.removeEdge(*ele, *fEle);
                } 
                else if (graph.retrieveWeight(*fEle, *ele) < LtrParameters::MIN_WEIGHT) {
                    graph.removeEdge(*ele, *fEle);
                    if (graph.isConnected(*fEle, *ele)) {
                        graph.removeEdge(*fEle, *ele);
                    }
                }

            }
        }

    }
}

/**
 * Retrieves forward or backward elements (dependent on the bool given)
 * and sorts them
 */
std::vector<Element *> Matcher::retrieveElementVec(bool isForward, DirectedGraph<Element> &g) {
    std::vector<Element *> valueVec = g.getValueVec();
    std::vector<Element *> r;
    for (auto &valuePtr : valueVec)
    {
        if (valuePtr->getIsForward() == isForward)
        {
            r.push_back(valuePtr);
        }
    }
    std::sort(r.begin(), r.end(), compareValuePtrs);
    r.shrink_to_fit();
    return r;
}

/**
 * Returns forward and backward elements (sorted) as a tuple
*/
std::tuple<std::vector<Element *>, std::vector<Element *>> Matcher::retrieveForwardBackward(DirectedGraph<Element> &g)
{
    return {retrieveElementVec(true, g), retrieveElementVec(false, g)};
}

std::tuple<std::vector<Element *>, std::vector<Element *>> Matcher::retrieveForwardBackward()
{
    return {retrieveElementVec(true, this->graph), retrieveElementVec(false, this->graph)};
}

/**
 * Will merge crosses in the graph; i.e., fragmented elements that belong to the same element, going forward
*/
void Matcher::mergeCrosses() {

    // Start with the forward elements
    // Moving in windows of two and extending when candidate found
    auto [fEleVec, bEleVec] = retrieveForwardBackward(graph);

    int i = 0; // Forward Element iterator
    int m = 0; // Backward Element iterator
    while (i < fEleVec.size() - 1) {
        auto firstConVec = graph.retrieveConnectedValues(*fEleVec[i]);
        auto secondConVec = graph.retrieveConnectedValues(*fEleVec[i + 1]);

        // Moving backward element iterator to first backward element of the possible connected backward elements;
        for (; m < bEleVec.size() - 1; m++) {
            if (bEleVec.at(m)->getStart() > fEleVec.at(i)->getEnd()) {
                m = m > 0? m - 1 : 0;
                break;
            }
        }
        for (; m > 0; m--) {
            if (bEleVec.at(m)->getStart() <= fEleVec.at(i)->getStart()) {
                break;
            }
        }

        // Both forward elements in the window need to be connected to at least one element (looking forward)
        if (firstConVec.size() > 0 && secondConVec.size() > 0) {
            std::sort(firstConVec.begin(), firstConVec.end(), compareValuePtrs);
            std::sort(secondConVec.begin(), secondConVec.end(), compareValuePtrs);

            Element *fEle = fEleVec[i];
            Element *bEle = firstConVec[0];


            // Finding where the bEle starts
            assert (bEleVec[m]->getStart() <= bEle->getStart());
            int n = m;
            for (; n < bEleVec.size(); n++) {
                if (bEleVec[n]->getStart() == bEle->getStart()) {
                    break;
                }
            }


            /**
             * Assuming that f is our window of forward element and b is the first connected element of each forward element in f
             * And assuming that t is the total length of the f window and y is the total length of the b window, apply these constraints:
             * Ft (the final forward element) must be before B1 (the first backward element) && b must be ordered by ascending position
             * Additionally, b must be ordered the same way as the bEleVec, such that no elements are skipped
             */ 
            const int b1 = firstConVec[0]->getStart(); // Won't change
            int ft = fEleVec[i + 1]->getEnd();
            int bj = firstConVec[0]->getEnd();
            int by = secondConVec[0]->getStart(); 

            std::vector<Element*> fMergeVec;
            std::vector<Element*> bMergeVec;


            int k = i + 2;
            // Satisfying the rule Ft < B1 && that Bj < By && that B is the same order as bEleVec
            // Extend the window until rule is not satisifed
            while (ft < b1 && (bj < by || bj == secondConVec[0]->getEnd()) && bEleVec[n + 1]->getStart() == by) {
                if (k < fEleVec.size()) {
                    fMergeVec.push_back(fEleVec[k - 1]);
                    bMergeVec.push_back(secondConVec[0]);

                    n = bj == secondConVec[0]->getEnd() ? n : n + 1;

                    firstConVec = secondConVec;
                    secondConVec = graph.retrieveConnectedValues(*fEleVec[k]);
                    if (secondConVec.size() > 0) {
                        std::sort(secondConVec.begin(), secondConVec.end(), compareValuePtrs);

                        ft = fEleVec[k]->getEnd();
                        bj = firstConVec[0]->getEnd();
                        by = secondConVec[0]->getStart();
                        k++;
                    }
                    else {
                        k++;
                        break;
                    }

                }
                else {
                    break;
                }
            }

            // If the condition was satisfied at least once
            if (fMergeVec.size() > 0) {
                // Merging the elements and removing old elements from graph
                for (int j = 0; j < fMergeVec.size(); j++) {
                    fEle->merge(*fMergeVec[j]);
                    graph.removeNode(*fMergeVec[j]);
                }
                // Merging in the backward elements and removing old elements from graph
                for (int j = 0; j < bMergeVec.size(); j++) {
                    bEle->merge(*bMergeVec[j]);
                    graph.removeNode(*bMergeVec[j]);
                }


                // Optimize this to a single method that removes all edges fanning in and out
                graph.removeNode(*fEle);
                graph.addNode(*fEle);

                graph.removeNode(*bEle);
                graph.addNode(*bEle);
            }

            i = k - 1;
        }
        else {
            i++;
        }
    }

    matchElements(graph);
}


/**
 * Match vertical nodes, i.e., forward elements that have overlaps with backword elements
 * Happen in recently nested repeats or sequential.
 */
void Matcher::matchOverlaps()
{
    isVertical = true;
    
    auto [forwardVec, backwardVec] = retrieveForwardBackward(graph);

    int i = 0;
    int j = 0;
    while (i < forwardVec.size() && j < backwardVec.size())
    {
        Element &f = *forwardVec[i];
        Element &b = *backwardVec[j];

        if (graph.isConnected(f, b) || graph.isConnected(b, f)) {
            i++;
            continue;
        }

        // If the forward and backward element have overlap of at least one nucleotide
        if (f.calcOverlap(b) > 0) {
            // Weights of zero mean overlap
            graph.addEdge(f, b, 0.0);
            graph.addEdge(b, f, 0.0);
        }

        // If the forward element is before the backward, the next forward has potential for overlap
        if (f.getEnd() < b.getEnd()) {
            i++;
        }
        // If the backward element is before the forward, the next backward has potential for overlap
        else if (f.getEnd() > b.getEnd()) {
            j++;
        }
        // Otherwise, their ends are equal to each other and both should be advanced.
        else {
            i++;
            j++;
        }
    }
}

/**
 * Find connected nodes with zero weights
 * [OK]
 */
std::vector<Element*> Matcher::retrieveVertCons(Element &father) {
    std::vector<Element *> r;
    for (auto cPtr : graph.retrieveConnectedValues(father)) {
        if (LtrUtility::isEqual(graph.retrieveWeight(father, *cPtr), 0.0)) {
            r.push_back(cPtr);
        }
    }
    r.shrink_to_fit();
    return r;
}

/**
 * Merge elements overlapping with the same forward (or backward) element. 
 * [OK]
 */
void Matcher::mergeVerticalOverlaps(std::vector<Element*> &eleVec) {
    assert (isVertical == true);

    for (auto &elePtr : eleVec) {
        Element &father = *elePtr;
        std::vector<Element*> childVec = retrieveVertCons(father); // Vertical connections of the father
        std::unordered_set<Element*> zeroSet;   // Will contain the vertical connections of the merged element

        if (childVec.size() > 1) {
            // Get all of the vertical connections of the to-be-merged children
            for (auto childPtr : childVec) {
                Element &child = *childPtr;
                std::vector<Element*> zeroVec = retrieveVertCons(child);
                zeroSet.insert(zeroVec.begin(), zeroVec.end());
            }

            auto &newChild = *childVec[0];

            // Merge the elements in childVec
            for (int i = 1; i < childVec.size(); i++) {
                auto &mergedChild = *childVec[i];
                newChild.merge(mergedChild);
                graph.removeNode(mergedChild);
            }

            graph.clearConnections(newChild);

            // Re-adding vertical connections
            for (auto ptr : zeroSet) {
                graph.addEdge(newChild, *ptr, 0.0);
                graph.addEdge(*ptr, newChild, 0.0);
            }
        }
    }
}

void Matcher::validateVerticals() {
    for (auto fPtr : graph.getValueVec()) {
        auto vertVec = retrieveVertCons(*fPtr);
        assert(vertVec.size() <= 1);
        if (!vertVec.empty()) {
            assert (LtrUtility::isEqual(graph.retrieveWeight(*vertVec.front(), *fPtr), 0.0));
        }

    }
}



void Matcher::breakGraph() {
    subGraphVec = graph.breakGraph();
    isBroken = true;
}


/**
 * Given that the graph has already been broken into subcomponents, break these subcomponents into further subcomponents based on a certain criteria
 * 1: If the weight between two nodes is high in one direction but low in the other, that is a sign of possible hyperextension. If the node with the lower
 *      weight can be broken into a smaller element such that the weight is within acceptable limits
 * 2: If this applies, the overlap connected should be revisited:
 *       O  O              O  OO                O  O
 *        \ |\      --->    \ |\        --->     \ |\
 *         \| \              \| \                 \| \
 *          O  O              O  O                 O  O
 *      The second element in the forward (first in the backward) breaks into two pieces based on weight. The second part of this piece is not part of the graph
 *      anymore and should thus be removed. 
*/
// void Matcher::breakCandidates() {
//     assert(isBroken == true);
    
//     std::vector<DirectedGraph<Element>> r;
//     for (auto &subGraph : subGraphVec) {
//         auto [forwardVec, backwardVec] = retrieveForwardBackward(subGraph);

//     }



//     subGraphVec.insert(subGraphVec.end(), r.begin(), r.end());

// }


void Matcher::runCases()
{
    // For every subgraph
    for (auto &subGraph : subGraphVec)
    {
        // std::cout << *filePath << " Matcher Loop" << std::endl;
        // std::vector<CaseMatcher *> caseVec{new CaseSingle{ic, red, seq}, new CaseRecent{ic, red, seq}, 
        //                                     new CaseSequential{ic, red, seq}, new CaseSolo{ic, red, seq}, 
        //                                     new CaseAllSolo{ic, red, seq}, new CaseDegenRecent{ic, red, seq}};
        // std::vector<CaseMatcher *> caseVec{new CaseSingle{ic, red, seq}, new CaseRecent{ic, red, seq}, 
        //                                     new CaseSequential{ic, red, seq}, new CaseSolo{ic, red, seq}};

        std::vector<CaseMatcher *> caseVec{new CaseSolo{ic, icRecent, red, seq}, new CaseSingle{ic, icRecent, red, seq}, new CaseRecent{ic, icRecent,red, seq}};
        CaseRecentComplex * complex = new CaseRecentComplex{ic, icRecent, red, seq};
                                            
        auto [forward, backward] = retrieveForwardBackward(subGraph);
        assert(!forward.empty() || !backward.empty());

        int subIndex = !forward.empty()? graph.retrieveSubGraphIndex(*forward.at(0)) : graph.retrieveSubGraphIndex(*backward.at(0));
        std::vector<RT *> r;

        // Getting all RTs from case analysis tests
        for (auto &caseTest : caseVec) {
            caseTest->apply(subGraph, forward, backward, subIndex);
            auto result = caseTest->getRTVec();
            r.insert(r.end(), result.begin(), result.end());
        }

        complex->apply(subGraph, forward, backward, subIndex);
        auto complexResult = complex->getRTVec();
        if (!complexResult.empty()) {
            complexVec.push_back(complexResult.front());
        }


        // LtrUtility::sortRTs(r);
        LtrUtility::sortRTs(complexVec);
        // // Ranking RTs
        LtrUtility::rankRTs(r);

        for (auto rt : r) {
            assert(rt != nullptr);
            rtVec.push_back(rt);
            rtSubGraphMap[rt] = &subGraph;
        }

        for (auto p : caseVec) {
            delete p;
        }
        delete complex;

    }
}

bool Matcher::checkRank(RT *r1, RT *r2) {
    int cr1 = r1->getCaseRank();
    int cr2 = r2->getCaseRank();

    return cr1 >= cr2;
}

bool Matcher::isSameRT(RT *r1, RT *r2) {
    // bool r = false;
    // // MAKE LTRPARAMETER HERE
    // if (LtrUtility::isGreaterEqual(r1->calcOverlap(r2), 0.8) {
    //     r = true;
    // }
    // return r;
    return true;
}



DirectedGraph<Element> *Matcher::getGraph() const
{
    return &graph;
}

std::vector<DirectedGraph<Element>> Matcher::getSubGraphVec() const
{
    return subGraphVec;
}

std::vector<RT *> *Matcher::getRtVec()
{
    return &rtVec;
}

std::vector<RT *> *Matcher::getComplexVec() {
    return &complexVec;
}

DirectedGraph<Element> *Matcher::getSubGraph(RT *rt) const
{
    return rtSubGraphMap[rt];
}

// Most likely is unused
// bool compareElementGraph(DirectedGraph<Element> &a, DirectedGraph<Element> &b)
// {
//     auto aFVec = a.getValueVec();
//     auto bFVec = b.getValueVec();
//     std::sort(aFVec.begin(), aFVec.end(), compareValuePtrs);
//     std::sort(bFVec.begin(), bFVec.end(), compareValuePtrs);

//     return (aFVec[0]->getStart() < bFVec[0]->getStart());
// }

// [OK]
// Most likely is unused
bool compareValuePtrs(Element *a, Element *b)
{
    return (a->getStart() < b->getStart());
}