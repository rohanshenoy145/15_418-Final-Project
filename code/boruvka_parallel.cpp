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
#include <LockFreeDisjointSets.h>
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


vector<Edge> MST(Graph &G){

    size_t init_size = G.nodes.size();
    size_t rounded_size = 4096 * ((init_size + 4095) / 4096);
    vector<Edge> mst_edges;
    
    ds::DisjointSets union_find(init_size);
    while (G.nodes.size() > 1) {

        // Find shortest edges for each node.

        vector<pair<int, int>> shortest_edges(init_size, { 0, INT_MAX});
        pair<int, int>*local_shortest = new pair<int,int>[number_of_threads*rounded_size];
        

        omp_set_num_threads(number_of_threads);
        #pragma omp parallel
        {

            const int nthreads = omp_get_num_threads();
            const int ithread = omp_get_thread_num();

            

        //    #pragma omp single 
        //     {
        //         local_shortest = new pair<int,int>[nthreads*rounded_size];
        //         for (size_t i = 0; i < nthreads*rounded_size; i ++) local_shortest[i] = make_pair(0, INT_MAX);
        //     }
            #pragma omp for schedule(static, rounded_size)
            for (size_t i = 0; i < nthreads*rounded_size; i ++) {
                local_shortest[i] = make_pair(0, INT_MAX);
            }
            #pragma omp barrier

            #pragma omp for
            for (size_t i = 0; i < G.edges.size(); i++) {
                Edge cur = G.edges[i];
                if (cur.w < local_shortest[(ithread*rounded_size) + cur.u].second) {
                    local_shortest[(ithread*rounded_size) + cur.u] = make_pair(i, cur.w);
                }
            }
            
            #pragma omp for
            for (size_t i = 0; i < init_size; i ++ ) {
                for (int t = 0 ; t < nthreads; t++ ) {
                    if (local_shortest[(rounded_size*t) + i].second < shortest_edges[i].second) {
                        shortest_edges[i] = local_shortest[(rounded_size*t) + i];
                    }
                }
            }

        }
        ////////////////////////////////////////////////////////////////////////

        omp_set_num_threads(number_of_threads);
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
    ///////////////////////////////////////////////////////////////////////////////////////////////////

        
    // Only keep edges who are the representative edges.
        std::vector<std::vector<Edge>> mst_private_select_edge(number_of_threads);
        vector<int>edgeCount(number_of_threads);
        
        #pragma omp parallel for 
        for (size_t i = 0; i < G.edges.size(); i ++ ){
            Edge cur = G.edges[i];
            if(!union_find.same(cur.u, cur.v))
            {
                cur.u = union_find.find(cur.u);
                cur.v = union_find.find(cur.v);
                mst_private_select_edge[omp_get_thread_num()].push_back(cur);
                edgeCount[omp_get_thread_num()]++;
            }
        }
        vector< int> prefix_sum2(number_of_threads);
        __gnu_parallel::partial_sum(edgeCount.begin(), edgeCount.end(), prefix_sum2.begin()); //prefix sum
        vector<Edge> new_edges(prefix_sum2[number_of_threads - 1]);
            

        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            int startIdx;
            if(threadId!=0)
            {
                startIdx = prefix_sum2[threadId - 1];
            }
            else{
                startIdx = 0;
            }

            for (int i = startIdx; i < startIdx + edgeCount[threadId]; i++){
                
                new_edges[i] = mst_private_select_edge[threadId][i - startIdx];

            }
        }    
    
        
    /////////////////////////////////////////////////////////////////////////////    
        // Only keep nodes who are the representative nodes.

        std::vector<std::vector<size_t>> mst_private_select_nodes(number_of_threads);
        vector<int>nodeCount(number_of_threads);

        #pragma omp parallel for
        for (size_t i = 0; i < G.nodes.size(); i ++){
            size_t cur_node = G.nodes[i];
            if(union_find.find(cur_node) == cur_node)
            {
                mst_private_select_nodes[omp_get_thread_num()].push_back(cur_node);
                nodeCount[omp_get_thread_num()]++;
            }
            
        }
        
        vector< int> prefix_sum3(number_of_threads);
        __gnu_parallel::partial_sum(nodeCount.begin(), nodeCount.end(), prefix_sum3.begin()); //prefix sum
        vector<size_t> new_nodes(prefix_sum3[number_of_threads - 1]);

        #pragma omp parallel 
        {
            int threadId = omp_get_thread_num();
            int startIdx;
            if(threadId!=0)
            {
                startIdx = prefix_sum3[threadId - 1];
            }
            else{
                startIdx = 0;
            }

            for (int i = startIdx; i < startIdx + nodeCount[threadId]; i++){
                
                new_nodes[i] = mst_private_select_nodes[threadId][i - startIdx];

            }
        }

        G.nodes = new_nodes;
        G.edges = new_edges;
    
    }
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
