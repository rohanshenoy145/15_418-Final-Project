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
#include <ParallelDisjointSets.h>
#include <omp.h>
#include <chrono>
#include <mutex>
#include <algorithm>
#include <parallel/algorithm>
#include <parallel/numeric>



using namespace std;
int number_of_threads;



struct Edge {
    size_t u;
    size_t v;
    int w;

   // Edge() : u(0), v(0), w(0) {}
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

// struct shortEdge {
//     int index;
//     int weight;
//     std::mutex mtx; // mutex for locking access to the struct's data

//     // Constructor to initialize the struct members
//     shortEdge() : index(0), weight(INT_MAX) {}
// };

bool compareBySource(const Edge& a, const Edge& b) {
    return a.u < b.u;
}


vector<Edge> MST(Graph &G){
    double timeFindShortestEdges = 0.0;
    double  addMSt = 0.0;
    double mapNewEdges = 0.0;
    double mapNewNodes = 0.0;

    size_t init_size = G.nodes.size();

    vector<Edge> mst_edges(init_size - 1);
    // size_t rounded_size = 4096 * ((init_size + 4095) / 4096);
    ds::DisjointSets union_find(init_size);
    // size_t mstIndexOffset = 0;
    // pair<int, int>*local_shortest = new pair<int,int>[number_of_threads*rounded_size];
    vector< int> select_edges(G.edges.size());
    vector< int> prefix_sum2(G.edges.size());
    vector< int> prefix_sum3(G.nodes.size());
    vector<int> selectNewEdges(G.edges.size());


    while (G.nodes.size() > 1) {

        // Find shortest edges for each node.

        // vector<shortEdge> shortest_edges(init_size);
        auto start = std::chrono::high_resolution_clock::now();

        __gnu_parallel::sort(G.edges.begin(), G.edges.end(), compareBySource);
        vector<int> shifts(G.edges.size());
        shifts[0] = 0;
        omp_set_num_threads(number_of_threads);
        #pragma omp parallel for
        for (size_t i = 1; i < G.edges.size(); i ++){
            shifts[i] = shifts[i-1] != shifts[i];
        }
        vector<int> pSum(G.edges.size());
        __gnu_parallel::partial_sum(shifts.begin(), shifts.end(), pSum.begin()); 
        vector<int> offsets(G.nodes.size()+1);
        offsets[0] = 0;
        omp_set_num_threads(number_of_threads);
        #pragma omp parallel for 
        for (size_t i = 1; i < G.edges.size(); i++){
            if (shifts[i]){
                offsets[pSum[i]] = i;
            }
        }
        offsets[G.nodes.size()] = G.edges.size();
        vector<pair<int, int>> shortest_edges(init_size, { 0, INT_MAX});
        omp_set_num_threads(number_of_threads);
        #pragma omp parallel for
        for (size_t i = 0; i < G.nodes.size(); i ++ ){
            for (int j = offsets[i]; j < offsets[i+1]; j ++){
                select_edges[j] = 0;
                Edge cur = G.edges[j];
                pair<int, int> shortest = shortest_edges[cur.u];
                if (cur.w < shortest.second) {
                    shortest_edges[cur.u] = {j, cur.w };
                }
            }
        }

       
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration<double>(end - start).count();
        timeFindShortestEdges+=duration;

        //vector< int> select_edges(G.edges.size());
        auto start2 = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<Edge>> mst_private_edge(number_of_threads);
        std::vector<int> privateCount(number_of_threads, 0); // Initialize all elements to 0
        
        #pragma omp parallel for 
        for (size_t i = 0; i < G.nodes.size(); i ++ ) {
            
            size_t u = G.nodes[i];
            Edge shortest_from_u = G.edges[shortest_edges[u].first];
            size_t v = shortest_from_u.v;
            Edge shortest_from_v = G.edges[shortest_edges[v].first];
            // Ensure not a duplicate edge.
            if (u != shortest_from_v.v || (u == shortest_from_v.v && u < v)){
                mst_private_edge[omp_get_thread_num()].push_back(shortest_from_u);// select edge 
                privateCount[omp_get_thread_num()]++; 
                union_find.unite(u, v);  
                
            }
        }
    
        
        vector<int> prefixSum(number_of_threads);
        __gnu_parallel::partial_sum(privateCount.begin(), privateCount.end(), prefixSum.begin()); //prefix sum
        vector <Edge> mst_edges_add(prefixSum[number_of_threads - 1]);

         #pragma omp parallel
         {
            int threadId = omp_get_thread_num();
            int startIdx;
            if(threadId!=0)
            {
                startIdx = prefixSum[threadId - 1];
            }
            else{
                startIdx = 0;
            }
            for(int i = startIdx; i < startIdx + privateCount[threadId]; i++)
            {
                mst_edges_add[i] = mst_private_edge[threadId][i-startIdx];
            }


         }
         mst_edges.insert(mst_edges.end(), mst_edges_add.begin(), mst_edges_add.end());

        auto end2 = std::chrono::high_resolution_clock::now();
        double duration2 = std::chrono::duration<double>(end2 - start2).count();
        addMSt+=duration2;
        auto start3 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for
        for(size_t i = 0; i < G.edges.size(); i++ )
        {
            Edge cur = G.edges[i];
            selectNewEdges[i] = !union_find.same(cur.u, cur.v); //Cross edges only
        }

        
       
   
         size_t offset = G.edges.size();
         //end_it = std::next(selectNewEdges.begin(), offset);

        __gnu_parallel::partial_sum(selectNewEdges.begin(), selectNewEdges.begin() + offset, prefix_sum2.begin()); //prefix sum
        vector<Edge> new_edges(prefix_sum2[G.edges.size() - 1]);
            
        #pragma omp parallel for
        for (size_t i = 0; i < G.edges.size(); i ++ ){
            if(selectNewEdges[i])
            {

                Edge cur = G.edges[i];
                cur.u = union_find.find(cur.u);
                cur.v = union_find.find(cur.v);
                new_edges[prefix_sum2[i] - 1] = cur;
            }

        }

        auto end3 = std::chrono::high_resolution_clock::now();
        double duration3 = std::chrono::duration<double>(end3 - start3).count();
        mapNewEdges+=duration3;
            
    
        
        auto start4 = std::chrono::high_resolution_clock::now();
        
        //Only keep nodes who are the representative nodes.
        vector<int> selectNewNodes(G.nodes.size());
        #pragma omp parallel for
        for (size_t i = 0; i < G.nodes.size(); i ++){
            size_t cur_node = G.nodes[i];
            selectNewNodes[i] = (union_find.find(cur_node) == cur_node);
        }
        
        offset = G.nodes.size();
         //end_it = std::next(selectNewNodes.begin(), offset);
        __gnu_parallel::partial_sum(selectNewNodes.begin(), selectNewNodes.begin() + offset, prefix_sum3.begin()); //prefix sum
        vector<size_t> new_nodes(prefix_sum3[G.nodes.size() - 1]);

        #pragma omp parallel for 
        for (size_t i = 0; i < G.nodes.size(); i ++){
            
            if(selectNewNodes[i])
            {
                new_nodes[prefix_sum3[i] - 1 ] = G.nodes[i];
            }
                
        }
        auto end4 = std::chrono::high_resolution_clock::now();
        double duration4 = std::chrono::duration<double>(end4 - start4).count();
        mapNewNodes+=duration4;


        G.nodes = new_nodes;
        G.edges = new_edges;
    
    
    }

    std::cout << "Total time taken to find minEdges parallel: " << timeFindShortestEdges << " seconds" << std::endl;
    std::cout << "Total time taken to add MST edges parallel: " << addMSt << " seconds" << std::endl;
    std::cout << "Total time taken to map new edges paralllel: " << mapNewEdges << " seconds" << std::endl;
    std::cout << "Total time taken to map new nodes parallel: " << mapNewNodes << " seconds" << std::endl;
    
    return mst_edges;
}

int main(int argc, char *argv[]) {
    char* inputFilePath = nullptr;
    int opt;
    while ((opt = getopt(argc, argv, "f:n:")) != -1) {
        switch (opt) {
            case 'f':
                // -f option is used to specify the input file path
                inputFilePath = optarg;
                break;
            case 'n':
                number_of_threads = atoi(optarg);
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
