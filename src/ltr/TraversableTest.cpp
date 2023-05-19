#include "DirectedGraph.h"
#include <vector>
#include "LtrUtility.h"
#include <unordered_set>

std::vector<int*> getDiagonals(int *e, DirectedGraph<int> & graph) {
    auto conVec = graph.retrieveConnectedValues(*e);
    std::vector<int *> r;
    for (auto ele : conVec) {
        if (!LtrUtility::isEqual(graph.retrieveWeight(*e, *ele), 0.0)) {
            r.push_back(ele);
        }
    }

    std::sort(r.begin(), r.end());
    return r;
}

int * CaseMatcher::getVertical(int *e, DirectedGraph<int> &graph) {
    auto eleVec = graph.retrieveConnectedValues(*e);
    int *r = nullptr;
    for (auto ele : eleVec) {
        if (LtrUtility::isEqual(graph.retrieveWeight(*e, *ele), 0.0)) {
            r = ele;
            break;
        } 
    }
    return r;
}

std::vector<int> getTraversable(int *e, DirectedGraph<int> &graph) {
    std::unordered_set<int*> visistedSet;
}


int main() {
    DirectedGraph<int> dg;

    int node1 = 1;
    int node2 = 2;
    int node3 = 3;
    int node4 = 4;
    int node5 = 5;
    int node6 = 6;
    int node7 = 7;
    int node8 = 8;

    dg.addNode(node1);
    dg.addNode(node2);
    dg.addNode(node3);
    dg.addNode(node4);
    dg.addNode(node5);
    dg.addNode(node6);
    dg.addNode(node7);
    dg.addNode(node8);

    dg.addEdge(node1, node2, 0.5);
    dg.addEdge(node2, node1, 0.5);

    dg.addEdge(node1, node4, 0.5);
    dg.addEdge(node4, node1, 0.5);

    dg.addEdge(node3, node4, 0.5);
    dg.addEdge(node4, node3, 0.5);

    dg.addEdge(node5, node4, 0.5);
    dg.addEdge(node4, node5, 0.5);

    dg.addEdge(node6, node7, 0.5);
    dg.addEdge(node7, node5, 0.5);

    node.addEdge(node6, node8, 0.5);
    node.addEdge(node8, node6, 0.5);

    node.addEdge(node2, node3, 0.0);
    node.addEdge(node3, node2, 0.0);

    node.addEdge(node4, node6, 0.0);
    node.addEdge(node6, node4, 0.0);



}