//
//  main.cpp
//  ClusteringCoefficient
//
//  Created by Abdullah Giray YAGLIKCI on 9/8/15.
//  Copyright (c) 2015 Abdullah Giray YAGLIKCI. All rights reserved.
//

#include <iostream>
#include "main.h"

extern "C" {
    extern void mcsim_skip_instrs_begin();
    extern void mcsim_skip_instrs_end();
    extern void mcsim_spinning_begin();
    extern void mcsim_spinning_end();
    int32_t log_2(uint64_t);
}

using namespace std;
pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t wrback_mutex = PTHREAD_MUTEX_INITIALIZER;



Graph my_graph(edgefile_graph_size(TIMESTAMP));

int main(int argc, char* argv[]) {
    //    int ts = TIMESTAMP;
    create_graph(TIMESTAMP);
    common_neighbors_master_thread(graph_size.second);
    clustering_coefficient_master_thread(graph_size.first);
    /*
     ts++;
     create_graph(ts);
     clustering_coefficient_master_thread(graph_size.first);
     
     ts++;
     create_graph(ts);
     clustering_coefficient_master_thread(graph_size.first);
     //*/
     cout <<"Done!"<< endl;
     
    return 0;
}


//============================================================================
// Functions About Common Neighbours PreProcessing

void * common_neighbors_slave_thread(void *arg){
    
    cn_thread_args * args = (cn_thread_args *) arg;
    vector<pair<vertex_desc, vertex_desc> > vi = *(args->vertex_pairs_ptr);
    
    vector< pair<vertex_desc, vertex_desc> > vertex_pairs = *(args->vertex_pairs_ptr) ;
    pthread_mutex_lock(&cout_mutex);
    cout << args->thread_id << ": started with the following pairs";
    
    while (!vertex_pairs.empty()) {
        vertex_desc v1 = (vertex_pairs.back().first), v2 = (vertex_pairs.back().second);
        vertex_pairs.pop_back();
        cout << " (" << v1 << "," << v2 << ")";
    
        adjac_iter ai , ai_end;
        int num_of_common_neighbors = 0;
        tie(ai,ai_end) = adjacent_vertices(v1,my_graph);
        for (  ; ai != ai_end ; ++ai )
            if (boost::edge(v2, *ai, my_graph).second)
                num_of_common_neighbors ++;
    
        if (num_of_common_neighbors > CN_THRESHOLD){
            pair<int, vertex_desc> v1rel, v2rel;
            v1rel.first=num_of_common_neighbors;
            v2rel.first=num_of_common_neighbors;
            v1rel.second = v2; v2rel.second = v1;
    
            pthread_mutex_lock(&wrback_mutex);
            my_graph[v1].relatedNodes.push_back(v1rel);
            my_graph[v2].relatedNodes.push_back(v2rel);
            push_heap(my_graph[v1].relatedNodes.begin(),my_graph[v1].relatedNodes.end());
            push_heap(my_graph[v2].relatedNodes.begin(),my_graph[v2].relatedNodes.end());
            pthread_mutex_unlock(&wrback_mutex);
        }
    }
    cout << endl;
    pthread_mutex_unlock(&cout_mutex);
    return NULL;
}


void common_neighbors_master_thread(int number_of_edges){
    pthread_t threads[NUM_OF_PTHREADS];
    cn_thread_args *th_args = (cn_thread_args*) malloc((NUM_OF_PTHREADS) * sizeof(cn_thread_args));
    int num_of_pairs = (graph_size.first * (graph_size.first - 1) / 2 ) - graph_size.second; //Number of non-neighbor pairs
    std::pair<vertex_iter , vertex_iter > vp = vertices(my_graph);
    vertex_iter source = vp.first;
    vertex_iter target = source + 1;
    vector<pair<vertex_desc, vertex_desc> > vertex_pairs [NUM_OF_PTHREADS];
//    adjac_iter tai, tai_end, sai, sai_end;
//    tie(sai,sai_end) = adjacent_vertices(*source,my_graph);
//    tie(tai,tai_end) = adjacent_vertices(*target,my_graph);
    
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        int chunk_size = ceil(((double) num_of_pairs ) / (NUM_OF_PTHREADS-i));
        num_of_pairs -= chunk_size;
        
        while (chunk_size > 0 && source < vp.second )  {
//            if (sai!=sai_end){
                while (target != vp.second && chunk_size > 0 ) {
                    if (/*tai != tai_end && */!boost::edge(*source, *target, my_graph).second) {
                        pair<vertex_desc,vertex_desc> vertex_pair;
                        vertex_pair.first = *source;
                        vertex_pair.second = *target;
                        vertex_pairs[i].push_back(vertex_pair);
                        chunk_size --;
                    }
                    target ++;
//                    tie(tai,tai_end) = adjacent_vertices(*target,my_graph);
                }
 //           }
            if (chunk_size > 0){
                source ++;
                target = source + 1;
//                tie(sai,sai_end) = adjacent_vertices(*source,my_graph);
//                tie(tai,tai_end) = adjacent_vertices(*target,my_graph);
            }
        }
        th_args[i].thread_id = i+1;
        th_args[i].vertex_pairs_ptr = &(vertex_pairs[i]);

    }
    
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_create(&threads[i], NULL, common_neighbors_slave_thread, (void*) &th_args[i]);
    }
    
    for(int j = 0; j < NUM_OF_PTHREADS; j++){
        pthread_join(threads[j], NULL);
    }
    free(th_args);
    
    // Log results
    pair<vertex_iter, vertex_iter> vi = vertices(my_graph);
    for (vertex_iter i = vi.first; i != vi.second; i++) {
        cout << *i << ":";
        vector<pair<int, vertex_desc> > relateds = my_graph[*i].relatedNodes;
        for (vector<pair<int, vertex_desc> >::iterator j = relateds.begin(); j != relateds.end(); j++) {
            cout << " " << (*j).first << "," << (*j).second ;
        }
        cout << endl;
    }
    // */

}

//============================================================================
// Functions About Clustering Coefficient Process
void * clustering_coefficient_slave_thread(void *arg){
    cc_thread_args * args = (cc_thread_args *) arg;
    int id = args->thread_id;
    vector<vertex_desc> vertices = *(args->vertices_pointer);
    
    while (!vertices.empty()) {
        vertex_desc vd = vertices.back();
        vertices.pop_back();
        out_edge_iter oei, oei_end;
        adjac_iter ai, ai_end;
        hint_access(vd,id);
        int num_of_neighbors = 0, num_of_triangles = 0;
        
        for ( tie(ai,ai_end) = adjacent_vertices(vd,my_graph) ; ai != ai_end ; ++ai ) {
            hint_access(*ai,id);
            num_of_neighbors++;
            adjac_iter ai2, ai2_end;
            
            for ( tie(ai2,ai2_end) = adjacent_vertices(*ai,my_graph) ; ai2 != ai2_end ; ++ai2 ) {
                if (boost::edge(vd, *ai2, my_graph).second) // 2nd order neighbor has an edge to current vertex
                    num_of_triangles++;
            }
        }
        my_graph[vd].clustering_coefficient = (num_of_neighbors > 1)? ((double)num_of_triangles) / ((double)(num_of_neighbors * (num_of_neighbors - 1))) : 0.0;
    }
    return NULL;
}


void clustering_coefficient_master_thread(int number_of_vertices){
    std::pair<vertex_iter, vertex_iter> vp = vertices(my_graph);
    vertex_iter vi = vp.first;
    
    pthread_t threads[NUM_OF_PTHREADS]; //array to hold thread information
    cc_thread_args *th_args = (cc_thread_args*) malloc((NUM_OF_PTHREADS) * sizeof(cc_thread_args));
    
    vector<vertex_desc> verticesArr [NUM_OF_PTHREADS];
    
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        int chunk_size = ceil(((double) number_of_vertices ) / (NUM_OF_PTHREADS-i));
        number_of_vertices -= chunk_size;
        
        while (chunk_size > 0) {
            verticesArr[i].push_back(*vi);
            my_graph[*vi].assigned = true;
            chunk_size --;
            
            /* Maximum common neighbor
             if (chunk_size > 0 && !my_graph[*vi].relatedNodes.empty()){
             vertex_iter related_iter = my_graph[*vi].relatedNodes.front().second;
             th_args[i].vertices.push_back(related_iter);
             my_graph[*related_iter].assigned = true;
             chunk_size --;
             }
             // */
            
            adjac_iter ai, ai_end;
            tie(ai,ai_end) = adjacent_vertices(*vi,my_graph);
            while (chunk_size > 0 && ai != ai_end) {
                if (!my_graph[*ai].assigned){
                    verticesArr[i].push_back(*ai);
                    my_graph[*ai].assigned = true;
                    chunk_size --;
                }
                ai ++;
            }
            while (my_graph[*vi].assigned && vi <= vp.second) {
                vi++;
            }
        }
        th_args[i].thread_id = i+1;
        th_args[i].vertices_pointer =  &(verticesArr[i]);
    }
    
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_create(&threads[i], NULL, clustering_coefficient_slave_thread, (void*) &th_args[i]);
    }
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_join(threads[i], NULL);
    }
    free(th_args);
    
    // Log results
    pair<vertex_iter, vertex_iter> viter = vertices(my_graph);
    for (vertex_iter i = viter.first; i != viter.second; i++) {
        cout << *i << ":" << my_graph[*i].clustering_coefficient<<endl;
    }
    // */
}



//============================================================================
// File Reading and Graph Creating Functions
int edgefile_graph_size(short timestamp){
    string file_prefix = FILE_PREFIX;
    string extension = FILE_EXTENSION;
    ifstream file(DIRECTORY+file_prefix + to_string(timestamp) + extension);
    //pthread_mutex_lock(&cout_mutex);
    //cout << "File: " << DIRECTORY <<  file_prefix <<  to_string(timestamp) <<  extension << endl;
    //pthread_mutex_unlock(&cout_mutex);
    short max_id = 0;
    if (file){
        for (string line; getline(file, line); ) {
            int comma_pos = (int) line.find(',');
            int temp_id = stoi(line.substr(0,comma_pos));
            if ( temp_id > max_id)
                max_id = temp_id;
            temp_id = stoi(line.substr(comma_pos+1,line.length()));
            if ( temp_id > max_id )
                max_id = temp_id;
        }
        max_id ++;
        //pthread_mutex_lock(&cout_mutex);
        //cout << max_id << " vertices are there in node pairs file." << endl;
        //pthread_mutex_unlock(&cout_mutex);
    }
    else {
        //pthread_mutex_lock(&cout_mutex);
        //cout << "File doesn't exist." << endl;
        //pthread_mutex_unlock(&cout_mutex);
    }
    graph_size.first = max_id;
    return max_id;
}

int create_graph(short timestamp){
    string file_prefix = FILE_PREFIX;
    string extension = FILE_EXTENSION;
    ifstream file(DIRECTORY+file_prefix + to_string(timestamp) + extension);
    int num_of_edges = 0;
    if (file){
        for (string line; getline(file, line); ) {
            int comma_pos = (int) line.find(',');
            int from_id = stoi(line.substr(0,comma_pos));
            int to_id = stoi(line.substr(comma_pos+1,line.length()-(comma_pos+1)));
            edge e;
            add_edge(from_id, to_id , e , my_graph);
            num_of_edges++;
        }
        
        std::pair<vertex_iter, vertex_iter> vp = vertices(my_graph);
        int i = 0;
        for (vertex_iter it = vp.first ; it < vp.second ; it++){
            my_graph[*it].index = i;
            i++;
        }
    }
    else{
        pthread_mutex_lock(&cout_mutex);
        cout << "File doesn't exist" << endl;
        pthread_mutex_unlock(&cout_mutex);
    }
    
    //pthread_mutex_lock(&cout_mutex);
    //cout << num_of_edges << " edges are there in node pairs file." << endl;
    //pthread_mutex_unlock(&cout_mutex);
    graph_size.second = num_of_edges;
    
    return num_of_edges;
}


void hint_access(vertex_desc vd , int thread_id = -1) {
    while ( !my_graph[vd].relatedNodes.empty() ){
        _mm_prefetch((char *) &my_graph.out_edge_list( my_graph[vd].relatedNodes.front().second) , 0);
        pop_heap(my_graph[vd].relatedNodes.begin(),my_graph[vd].relatedNodes.end());
        my_graph[vd].relatedNodes.pop_back();
    }
}