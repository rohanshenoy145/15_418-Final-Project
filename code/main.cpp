#include <iostream>
#include <unordered_map>
#include <vector>
#include <set>
#include <random>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>


struct Edge {
    int u;
    int v;
    int w;
    int l;
};

std::unordered_map<int, std::vector<Edge>> make_graph(std::vector<Edge>&input_edges) {
    std::unordered_map<int, std::vector<Edge>> G;
    for (size_t i = 0; i < input_edges.size(); i ++) {
        Edge e = input_edges[i];
        auto it = G.find(e.u);
        if (it == G.end()){
            G[e.u] = {e};
        }
        else{
            G[e.u].push_back(e);
        }
    }
    return G;
}

std::unordered_map<int, int> flip_coins (std::unordered_map<int, std::vector<Edge>> &G){
    std::unordered_map<int, int> flips;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    for (const auto &item: G ){
        int u = item.first;
        flips[u] = dis(gen);
    }
    return flips;
}

std::vector<int> MST(std::unordered_map<int, std::vector<Edge>> &G, std::vector<int> &T){

    std::unordered_map<int, int> flips = flip_coins(G);
    while (G.size() > 1) {
        // For each vertex, find its minimum cost edge
        std::unordered_map<int, int> P;
        for (const auto &item: G ){
            int src = item.first;
            std::vector<Edge> neighbors = item.second;
            int min_weight = INT_MAX;
            int min_dest = -1;
            int min_label = INT_MAX;
            for (size_t i = 0; i < neighbors.size(); i++){
                Edge cur_edge = neighbors[i];
                if (cur_edge.w < min_weight){
                    min_dest = cur_edge.v;
                    min_weight = cur_edge.w;
                    min_label = cur_edge.l;
                }
                else if (cur_edge.w == min_weight && cur_edge.l < min_label){
                    min_dest = cur_edge.v;
                    min_weight = cur_edge.w;
                    min_label = cur_edge.l;
                }
            }
            Edge e;
            e.u = src;
            e.w = min_weight;
            e.v = min_dest;
            e.l = min_label;
            // Add the label to the MST vector and contract
            if ( !flips[e.u] && flips[e.v] ) {
                T.push_back(e.l);
                P[e.u] = e.v;
            }
        }
        for (const auto &item: G ){
            auto it = P.find(item.first);
            if (it == P.end()){
                P[item.first] = item.first;
            }
        }

        std::unordered_map<int, std::vector<Edge>> new_G;
        for (const auto &item: G) {
            std::vector<Edge> neighbors = item.second;
            for (size_t i = 0; i < neighbors.size(); i ++ ) {
                int u = neighbors[i].u;
                int v = neighbors[i].v;
                if (P[u] != P[v]) {
                    auto it = new_G.find(P[u]);
                    Edge new_edge;
                    new_edge.u = P[u];
                    new_edge.v = P[v];
                    new_edge.w = neighbors[i].w;
                    new_edge.l = neighbors[i].l;
                    if (it == new_G.end()){
                        new_G[P[u]] = {new_edge};
                    }
                    else{
                        new_G[P[u]].push_back(new_edge);
                    }
                }
            }
        }
        G = new_G;
    }
    return T;
}

int main(int argc, char *argv[]) {
    char* inputFilePath = nullptr;
    int opt;
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        switch (opt) {
            case 'f':
                // -f option is used to specify the input file path
                inputFilePath = optarg;
                break;
            default:
                // Print usage information if an invalid option is provided
                std::cerr << "Usage: " << argv[0] << " -f <input_file_path>" << std::endl;
                return 1;
        }
    }
    FILE * fd = fopen(inputFilePath, "rt");
    std::vector<Edge> input_edges;
    if (fd){
        int i = 0;
        char linebuf[50];
        while (fgets (linebuf, 50, fd)) {
            int u;
            int v;
            int w;
            sscanf(linebuf, "%d %d %d", &u, &v, &w);
            Edge e;
            e.u = u;
            e.v = v;
            e.w = w;
            e.l = i;
            input_edges.push_back(e);
            i ++;
        }
    }
    fclose(fd);

    std::unordered_map<int, std::vector<Edge>> G = make_graph(input_edges);

    std::vector<int> T;
    std::vector<int> res = MST(G, T);

    std::vector<Edge> final_edges;
    for (size_t i = 0; i < T.size(); i ++) {
        final_edges.push_back(input_edges[T[i]]);
    }

    // std::ofstream outputFile();
    // for (size_t i = 0; i < final_edges.size(); i ++ ) {
        
    // }


    
    return 0;
}
