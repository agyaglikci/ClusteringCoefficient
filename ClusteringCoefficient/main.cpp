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
pthread_mutex_t args_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t edgelist_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t my_graph_mutex = PTHREAD_MUTEX_INITIALIZER;


Graph my_graph(edgefile_graph_size(TIMESTAMP));

int main(int argc, char* argv[]) {
    int ts = TIMESTAMP;
    create_graph(TIMESTAMP);
    common_neighbors_master_thread(graph_size.second);
    clustering_coefficient_master_thread(graph_size.first);

    ts++;
    create_graph(ts);
    clustering_coefficient_master_thread(graph_size.first);
    
    ts++;
    create_graph(ts);
    clustering_coefficient_master_thread(graph_size.first);
    
    //cout <<"Done!"<< endl;
    return 0;
}


//============================================================================
// Functions About Common Neighbours PreProcessing

void * common_neighbors_slave_thread(void *arg){
    cn_thread_args * args = (cn_thread_args *) arg;
    edge_iter ei = args->ei;
    //int id = args->thread_id;
    for (int chunk_size = args->chunk_size ; chunk_size > 0 ; chunk_size --) {
        vertex_desc v1 = source(*ei, my_graph);
        vertex_desc v2 = target(*ei, my_graph);
        adjac_iter ai , ai_end;
        int num_of_common_neighbors = 0;
        tie(ai,ai_end) = adjacent_vertices(v1,my_graph);
        for (  ; ai != ai_end ; ++ai ){
            if (boost::edge(*ai, v2, my_graph).second){
                num_of_common_neighbors ++;
            }
        }
        pthread_mutex_lock(&wrback_mutex);
        my_graph[*ei].common_neighbors = num_of_common_neighbors;
        
        int v1index = my_graph[v1].related_node_top_index;
        my_graph[v1].relatedNodes[v1index].pointer = &my_graph.out_edge_list(v2);
        my_graph[v1].relatedNodes[v1index].weight = num_of_common_neighbors;
        my_graph[v1].relatedNodes[v1index].v = v2;
        my_graph[v1].related_node_top_index = v1index + 1;
        
        int v2index = my_graph[v2].related_node_top_index;
        my_graph[v2].relatedNodes[v2index].pointer = &my_graph.out_edge_list(v1);
        my_graph[v2].relatedNodes[v2index].weight = num_of_common_neighbors;
        my_graph[v2].relatedNodes[v2index].v = v1;
        my_graph[v2].related_node_top_index = v2index + 1;
        pthread_mutex_unlock(&wrback_mutex);

        
        ei++;
    }
    
    return NULL;
}


void common_neighbors_master_thread(int number_of_edges){
    pthread_t threads[NUM_OF_PTHREADS];
    cn_thread_args *th_args = (cn_thread_args*) malloc((NUM_OF_PTHREADS) * sizeof(cn_thread_args));
    edge_iter ei, ei_end;
    tie(ei, ei_end) = edges(my_graph);
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        int chunk_size = ceil(((double) number_of_edges ) / (NUM_OF_PTHREADS-i));
        th_args[i].ei = ei;
        th_args[i].chunk_size = chunk_size;
        th_args[i].thread_id = i+1;
        number_of_edges -= chunk_size;
        for (; ei != ei_end && chunk_size > 0; chunk_size--) {
            ei++;
        }
    }
    
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_create(&threads[i], NULL, common_neighbors_slave_thread, (void*) &th_args[i]);
    }
    
    for(int j = 0; j < NUM_OF_PTHREADS; j++){
        pthread_join(threads[j], NULL);
    }
    free(th_args);
    //cout << "Common neighbors done." << endl;
}

//============================================================================
// Functions About Clustering Coefficient Process
void * clustering_coefficient_slave_thread(void *arg){
    cc_thread_args * args = (cc_thread_args *) arg;
    int id = args->thread_id;
    int vector_size = args->vector_size;
    vertex_desc * vertices = (args->vertices);
    
    for (int i = 0; i < vector_size ; i++) {
        vertex_desc vd = vertices[i];
        out_edge_iter oei, oei_end;
        adjac_iter ai, ai_end;
        hint_access(vd,id);
        int num_of_neighbors = 0;
        int num_of_triangles = 0;
        //pthread_mutex_lock(&cout_mutex);
        //cout << "Thread " << id << " Accessed place: " << hex << &my_graph.out_edge_list(vd) << endl;
        //pthread_mutex_unlock(&cout_mutex);
        for ( tie(ai,ai_end) = adjacent_vertices(vd,my_graph) ; ai != ai_end ; ++ai ) {
            hint_access(*ai,id);
            num_of_neighbors++;
            adjac_iter ai2, ai2_end;
            //pthread_mutex_lock(&cout_mutex);
            //cout << "Thread " << id << " Accessed place: " << hex << &my_graph.out_edge_list(*ai) << endl;
            //pthread_mutex_unlock(&cout_mutex);
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
    
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        int chunk_size = ceil(((double) number_of_vertices ) / (NUM_OF_PTHREADS-i));
        number_of_vertices -= chunk_size;
        th_args[i].vector_size = chunk_size;
        int k = 0;
        while (chunk_size > 0) {
            th_args[i].vertices[k] = *vi; //vertices.push_back(*vi);
            k++;
            my_graph[*vi].assigned = true;
            chunk_size --;
            int num_of_relateds = my_graph[*vi].related_node_top_index;
            for (int j = 0 ; chunk_size > 0 && j < num_of_relateds; j++) {
                vertex_desc vd = my_graph[*vi].relatedNodes[j].v;
                if (my_graph[vd].assigned == false){
                    th_args[i].vertices[k] = vd; //vertices.push_back(vd);
                    my_graph[vd].assigned = true;
                    chunk_size --;
                    k++;
                }
            }
            while (my_graph[*vi].assigned && vi <= vp.second) {
                vi++;
            }
        }
        th_args[i].thread_id = i+1;
    }
    
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_create(&threads[i], NULL, clustering_coefficient_slave_thread, (void*) &th_args[i]);
    }
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_join(threads[i], NULL);
    }
    free(th_args);
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
    for (int i = 0 ; i < my_graph[vd].related_node_top_index ; i++){
        _mm_prefetch((char *) my_graph[vd].relatedNodes[i].pointer, 0);
    }
}