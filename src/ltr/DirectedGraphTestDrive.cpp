#include "DirectedGraph.h"
#include "Node.h"

#include <vector>
#include <iostream>
#include <unordered_map>

DirectedGraph<int> makeGraph(std::vector<int> &v)
{
    DirectedGraph<int> g;
    for (auto &x : v)
    {
        g.addNode(x);
    }
    return g;
}

int main()
{
    std::vector<int> vec{0, 1, 2, 3, 4, 5, 6, 7, 8};

    // // #### TESTING SUBGRAPH ####
    // std::cout << "Testing subgraphs!" << std::endl;

    // DirectedGraph<int> graph_1 = makeGraph(vec);

    // graph_1.addEdge(vec[0], vec[1], 0.5);
    // graph_1.addEdge(vec[0], vec[2], 0.4);
    // graph_1.addEdge(vec[3], vec[0], 0.8);
    // graph_1.addEdge(vec[7], vec[2], 0.2);

    // graph_1.printSubGraphTable();
    // std::cout << std::endl << std::endl << std::endl << "---------------------------" << std::endl;

    // // #### TESTING DFS ALGORITHM ####
    // std::cout << "Testing DFS Algorithm!" << std::endl;

    // DirectedGraph<int> graph_2 = makeGraph(vec);

    // graph_2.addEdge(vec[0], vec[1]);
    // graph_2.addEdge(vec[1], vec[2]);
    // graph_2.addEdge(vec[2], vec[3]);
    // graph_2.addEdge(vec[3], vec[4]);
    // graph_2.addEdge(vec[4], vec[5]);
    // graph_2.addEdge(vec[5], vec[6]);
    // graph_2.addEdge(vec[6], vec[7]);
    // graph_2.addEdge(vec[7], vec[0]);

    // std::cout << graph_2 << std::endl;
    // std::cout << "Is 0 able to reach 7?  Answer: " << graph_2.isReachable(vec[0], vec[7]) << std::endl;
    // std::cout << "Is 0 able to reach 0?  Answer: " << graph_2.isReachable(vec[0], vec[0]) << std::endl;

    // std::cout << std::endl << std::endl << std::endl << "---------------------------" << std::endl;

    // DirectedGraph<int> graph_3 = makeGraph(vec);

    // graph_3.addEdge(vec[0], vec[1]);
    // graph_3.addEdge(vec[1], vec[2]);
    // graph_3.addEdge(vec[2], vec[3]);
    // graph_3.addEdge(vec[3], vec[4]);
    // graph_3.addEdge(vec[4], vec[5]);
    // graph_3.addEdge(vec[5], vec[6]);
    // // Dropped connection from 6 to 7
    // graph_3.addEdge(vec[7], vec[0]);

    // std::cout << graph_3 << std::endl;
    // std::cout << "Is 0 able to reach 7?  Answer: " << graph_3.isReachable(vec[0], vec[7]) << std::endl;
    // std::cout << "Is 0 able to reach 7 if we ignore direction? Answer: " << graph_3.isReachable(vec[0], vec[7], false) << std::endl;
    // std::cout << std::endl << std::endl << std::endl << "---------------------------" << std::endl;

    // DirectedGraph<int> graph_4 = makeGraph(vec);

    // graph_4.addEdge(vec[0], vec[1]);
    // graph_4.addEdge(vec[1], vec[2]);
    // graph_4.addEdge(vec[2], vec[3]);
    // // Adding connection from 2 to 6
    // graph_4.addEdge(vec[2], vec[6]);
    // graph_4.addEdge(vec[3], vec[4]);
    // graph_4.addEdge(vec[4], vec[5]);
    // graph_4.addEdge(vec[5], vec[6]);
    // graph_4.addEdge(vec[6], vec[7]);
    // graph_4.addEdge(vec[7], vec[0]);

    // std::cout << graph_4 << std::endl;
    // std::cout << "Is 0 able to reach 7?  Answer: " << graph_4.isReachable(vec[0], vec[7]) << std::endl;
    // std::cout << std::endl << std::endl << std::endl << "---------------------------" << std::endl;

    // ############## TESTING REMOVAL OF EDGE ################ //

    // Case One; just one node; should return a set with one value
    // DirectedGraph<int> graph_5;
    // graph_5.addNode(vec[0]);
    // std::cout << graph_5 << "Subgraph: ";
    // for (auto i : graph_5.retrieveSubGraphValueVec(vec[0])) {
    //     std::cout << *i << " ";
    // }
    // std::cout << std::endl << std::endl << std::endl << "---------------------------" << std::endl;
    // graph_5.checkGraph();

    // Case Two; two nodes A and B
    // DirectedGraph<int> graph_6;
    // graph_6.addNode(vec[0]);
    // graph_6.addNode(vec[1]);
    // graph_6.addEdge(vec[1], vec[0], 0.9);
    // std::cout << graph_6 << "This graph should all be connected!" << std::endl;
    // std::cout << "Subgraph: ";
    // for (auto i : graph_6.retrieveSubGraphValueVec(vec[0]))
    // {
    //     std::cout << *i << " ";
    // }
    // std::cout << std::endl
    //           << std::endl;
    // graph_6.checkGraph();

    // std::cout << "If we remove the edge between 1 and 0..." << std::endl;
    // graph_6.removeEdge(vec[1], vec[0]);
    // std::cout << "This graph should not be connected!" << std::endl;
    // graph_6.printSubGraphTable();
    // graph_6.checkGraph();

    // std::cout << graph_6 << "Subgraph of 0: ";
    // for (auto i : graph_6.retrieveSubGraphValueVec(vec[0]))
    // {
    //     std::cout << *i << " ";
    // }
    // std::cout << std::endl
    //           << "Subgraph of 1: ";
    // for (auto i : graph_6.retrieveSubGraphValueVec(vec[1]))
    // {
    //     std::cout << *i << " ";
    // }
    // std::cout << std::endl
    //           << std::endl
    //           << std::endl
    //           << "---------------------------" << std::endl;

    // Case Three; complex graph
    DirectedGraph<int> graph_7;
    graph_7.addNode(vec[1]);
    graph_7.addNode(vec[2]);
    graph_7.addNode(vec[3]);
    graph_7.addNode(vec[4]);
    graph_7.addNode(vec[5]);
    graph_7.addNode(vec[6]);
    graph_7.addNode(vec[7]);
    graph_7.addNode(vec[8]);

    graph_7.addEdge(vec[1], vec[7]);
    graph_7.addEdge(vec[6], vec[7]);
    graph_7.addEdge(vec[6], vec[8]);
    graph_7.addEdge(vec[5], vec[6]);
    graph_7.addEdge(vec[4], vec[6]);
    graph_7.addEdge(vec[4], vec[5]);
    graph_7.addEdge(vec[5], vec[4]);
    graph_7.addEdge(vec[3], vec[4]);
    graph_7.addEdge(vec[3], vec[1]);
    graph_7.addEdge(vec[3], vec[2]);

    graph_7.checkGraph();

    std::cout << graph_7 << "This graph should all be connected!" << std::endl;
    graph_7.printSubGraphTable();
    std::cout << std::endl
              << std::endl;

    std::cout << "If we remove the edge between 3 and 4..." << std::endl;
    graph_7.removeEdge(vec[3], vec[4]);
    graph_7.checkGraph();
    graph_7.printSubGraphTable();
    std::cout << std::endl
              << std::endl;

    std::cout << "If we remove the edge between 1 and 7..." << std::endl;
    graph_7.removeEdge(vec[1], vec[7]);
    graph_7.checkGraph();
    graph_7.printSubGraphTable();

    std::cout << std::endl
              << std::endl
              << std::endl
              << "---------------------------" << std::endl;
}