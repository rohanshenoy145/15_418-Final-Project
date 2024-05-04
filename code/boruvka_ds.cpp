#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <set>
#include <random>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <DisjointSets.h>
#include <chrono>

using namespace std;


struct Edge {
    size_t u;
    size_t v;
    int w;
};

struct Graph {
    vector<Edge> edges;
    vector<size_t> nodes;
};

Graph make_graph(std::vector<Edge>&input_edges) {
    Graph G;
    vector<Edge> edge_list;
    size_t max_node = 0;
    for (size_t i = 0; i < input_edges.size(); i ++) {
        Edge e = input_edges[i];
        edge_list.push_back(e);
        Edge back;
        back.u = e.v;
        back.v = e.u;
        back.w = e.w;
        edge_list.push_back(back);
        max_node = max(max_node, max(e.u,e.v));
    }

    vector<size_t> node_list;
    for (size_t i=0; i < max_node + 1; i ++ ) {
        node_list.push_back(i);
    }
    G.edges = edge_list;
    G.nodes = node_list;
    return G;
}

std::unordered_map<int, int> flip_coins (std::unordered_map<int,std::vector<Edge>> &G){
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

vector<Edge> MST(Graph &G){
    vector<Edge> mst_edges;

    size_t init_size = G.nodes.size();

    ds::DisjointSets union_find(init_size);
    double timeFindShortestEdges = 0.0;
    double  addMSt = 0.0;
    double mapNewEdges = 0.0;
    double mapNewNodes = 0.0;

    
    while (G.nodes.size() > 1) {
        auto start = std::chrono::high_resolution_clock::now();
        // Find shortest edges for each node.
        vector<pair<int, int>> shortest_edges(init_size, { 0, INT_MAX});
        for (size_t i = 0; i < G.edges.size(); i ++) {
            Edge cur = G.edges[i];
            pair<int, int> shortest = shortest_edges[cur.u];
            
            if (cur.w < shortest.second) {
                shortest_edges[cur.u] = {i, cur.w };
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration<double>(end - start).count();
        timeFindShortestEdges+=duration;

        auto start2 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < G.nodes.size(); i ++ ) { 
            size_t u = G.nodes[i];
            Edge shortest_from_u = G.edges[shortest_edges[u].first];

            size_t v = shortest_from_u.v;
            Edge shortest_from_v = G.edges[shortest_edges[v].first];

            // Ensure not a duplicate edge.
            if (u != shortest_from_v.v || (u == shortest_from_v.v && u < v)){
                mst_edges.push_back(shortest_from_u);
                union_find.unite(u, v);
            }
        }
        auto end2 = std::chrono::high_resolution_clock::now();
        double duration2 = std::chrono::duration<double>(end2 - start2).count();
        addMSt+=duration2;

        auto start3 = std::chrono::high_resolution_clock::now();
        vector<Edge> new_edges;
        for (size_t i = 0; i < G.edges.size(); i ++ ){
            Edge cur = G.edges[i];
            // Cross edges only
            if (!union_find.same(cur.u, cur.v)) {
                // Map endpoints to their new representative nodes
                cur.u = union_find.find(cur.u);
                cur.v = union_find.find(cur.v);
                new_edges.push_back(cur);
            }
        }
        auto end3 = std::chrono::high_resolution_clock::now();
        double duration3 = std::chrono::duration<double>(end3 - start3).count();
        mapNewEdges+=duration3;

        auto start4 = std::chrono::high_resolution_clock::now();
        // Only keep nodes who are the representative vectors.
        vector<size_t> new_nodes;
        for (size_t i = 0; i < G.nodes.size(); i ++){
            size_t cur_node = G.nodes[i];
            if (union_find.find(cur_node) == cur_node) {
                new_nodes.push_back(cur_node);
            }
        }
        auto end4 = std::chrono::high_resolution_clock::now();
        double duration4 = std::chrono::duration<double>(end4 - start4).count();
        mapNewNodes+=duration4;


        G.nodes = new_nodes;
        G.edges = new_edges;
    }
    std::cout << "Total time taken to find minEdges seq: " << timeFindShortestEdges << " seconds" << std::endl;
    std::cout << "Total time taken to add MST edges seq: " << addMSt << " seconds" << std::endl;
    std::cout << "Total time taken to map new edges seq: " << mapNewEdges << " seconds" << std::endl;
    std::cout << "Total time taken to map new nodes seq: " << mapNewNodes << " seconds" << std::endl;




    

    return mst_edges;
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
        char linebuf[50];
        while (fgets (linebuf, 50, fd)) {
            size_t u;
            size_t v;
            size_t w;
            sscanf(linebuf, "%lu %lu %lu", &u, &v, &w);
            Edge e;
            e.u = u;
            e.v = v;
            e.w = w;
            input_edges.push_back(e);
        }
        fclose(fd);
    }


    Graph G = make_graph(input_edges);
    size_t init_N = G.nodes.size();
    auto start = std::chrono::high_resolution_clock::now();
    vector<Edge> res = MST(G);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;

    std::string out_name = inputFilePath;
    std::string output_file_path = out_name + "_out.txt";
    std::ofstream outputFile(output_file_path);

    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return 1;
    }

    for (size_t i = 0; i < res.size(); i ++) {
        Edge e = res[i];
        outputFile << e.u << " " << e.v << " " << e.w << std::endl;
    }
     
    cout << "N = " << init_N << " M = " << res.size() << endl;
    
    outputFile.close();



    return 0;
}

