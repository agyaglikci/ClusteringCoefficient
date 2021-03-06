//
//  main.cpp
//  ClusteringCoefficient
//
//  Created by Abdullah Giray YAGLIKCI on 9/8/15.
//  Copyright (c) 2015 Abdullah Giray YAGLIKCI. All rights reserved.
//

#include <iostream>
#include "main.h"

using namespace std;
pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t cntr_mutex = PTHREAD_MUTEX_INITIALIZER;

Graph my_graph;

int main(int argc, char* argv[]) {
    my_graph = Graph(edgefile_graph_size(START_TIMESTAMP,DATA_SIZE));
    
    for (int ts = START_TIMESTAMP; ts <= STOP_TIMESTAMP; ts++) {
        create_graph(START_TIMESTAMP,DATA_SIZE);
        //if (ts == START_TIMESTAMP)
            //common_neighbors_master_thread(graph_size.second,CN_THRESHOLD);
        clustering_coefficient_master_thread(graph_size.first - numof_noedge_vertices, ts);
    }
    
    //cout << "NUM_OF_THREADS:" << NUM_OF_PTHREADS << endl;
    //cout << "DATASET_SIZE:" << DATA_SIZE << endl;
    //cout << "START_TIME_STAMP:" << START_TIMESTAMP << endl;
    //cout << "STOP_TIME_STAMP:" << STOP_TIMESTAMP << endl;
    log_resource_stats("CNMASTER ",cn_master_stats);
    for(int i = 0 ; i < NUM_OF_PTHREADS ; i++){
        ostringstream title(""); title << "CNSLAVE." << i;
        log_resource_stats(title.str() ,cn_slave_stats[i]);
    }
    for (int i = START_TIMESTAMP ; i <= STOP_TIMESTAMP; i++) {
        ostringstream title(""); title << "CCMASTER" << i;
        log_resource_stats(title.str(),cc_master_stats[i]);
        for(int j = 0 ; j < NUM_OF_PTHREADS ; j++){
            ostringstream title(""); title << "CCSLAVE" << i << "." << j;
            log_resource_stats(title.str() ,cc_slave_stats[i][j]);
        }
    }
    
    return 0;
}

void log_resource_stats( string title, pair<rusage,rusage> stats){
    
    long UT_SEC = (stats.second.ru_utime.tv_sec - stats.first.ru_utime.tv_sec);
    long UT_USEC = (stats.second.ru_utime.tv_usec - stats.first.ru_utime.tv_usec);
    cout << title << ":";
    if (UT_USEC >= 0)
        cout << UT_SEC <<"."<< UT_USEC;
    else
        cout << UT_SEC - 1 <<"."<< 1000000 + UT_USEC; // */
    
    long ST_SEC = (stats.second.ru_stime.tv_sec - stats.first.ru_stime.tv_sec);
    long ST_USEC = (stats.second.ru_stime.tv_usec - stats.first.ru_stime.tv_usec);
    cout << " , ";
    if (ST_USEC >= 0)
        cout << ST_SEC <<"."<< ST_USEC;
    else
        cout << ST_SEC - 1 <<"."<< 1000000 + ST_USEC; // */

    cout << endl;
}

//============================================================================
// Functions About Common Neighbours PreProcessing

void * common_neighbors_slave_thread(void *arg){
    cn_thread_args * args = (cn_thread_args *) arg;
    int thread_id = args->thread_id;
    getrusage(RUSAGE_SELF,&(cn_slave_stats[thread_id].first));
    
    int cn_threshold = args->cn_threshold;
    vertex_iter source = args->source_vertex;
    vertex_iter target = source + 1;
    bool lock_is_mine = false;
    pair<vertex_iter , vertex_iter > vp = vertices(my_graph);
    
    //pthread_mutex_lock(&cout_mutex);
    //cout << args->thread_id << ": started with the following pairs";
    
    while ( source < vp.second )  {
        pthread_mutex_lock(&(my_graph[*source].wrback_mutex));
        if (!my_graph[*source].paired){
            my_graph[*source].paired = true;
            lock_is_mine = true;
        }
        else
            lock_is_mine = false;
        pthread_mutex_unlock(&(my_graph[*source].wrback_mutex));
        
        if (!my_graph[*source].noedge && lock_is_mine){
            while (target < vp.second ) {
                if (!my_graph[*target].noedge && !boost::edge(*source, *target, my_graph).second) {
                    
                    adjac_iter ai , ai_end;
                    int num_of_common_neighbors = 0;
                    tie(ai,ai_end) = adjacent_vertices(*source,my_graph);
                    for (  ; ai != ai_end ; ++ai )
                        if (boost::edge(*target, *ai, my_graph).second)
                            num_of_common_neighbors ++;
                    
                    if (num_of_common_neighbors > cn_threshold){
                        pair<int, vertex_desc> v1rel, v2rel;
                        v1rel.first=num_of_common_neighbors;
                        v2rel.first=num_of_common_neighbors;
                        v1rel.second = *target; v2rel.second = *source;
                        
                        pthread_mutex_lock(&(my_graph[*source].wrback_mutex));
                        my_graph[*source].relatedNodes.push_back(v1rel);
                        push_heap(my_graph[*source].relatedNodes.begin(),my_graph[*source].relatedNodes.end());
                        pthread_mutex_unlock(&(my_graph[*source].wrback_mutex));
                        
                        pthread_mutex_lock(&(my_graph[*target].wrback_mutex));
                        my_graph[*target].relatedNodes.push_back(v2rel);
                        push_heap(my_graph[*target].relatedNodes.begin(),my_graph[*target].relatedNodes.end());
                        pthread_mutex_unlock(&(my_graph[*target].wrback_mutex));
                        
                        pthread_mutex_lock(&cntr_mutex);
                        numof_link_predictions[num_of_common_neighbors] ++;
                        if (num_of_common_neighbors > max_link_predictions)
                            max_link_predictions = num_of_common_neighbors;
                        pthread_mutex_unlock(&cntr_mutex);
                    }
                }
                target ++;
            }
        }
        
        source ++;
        target = source;
        target ++;
    }

    
   
    getrusage(RUSAGE_SELF,&(cn_slave_stats[thread_id].second));
    //cout << endl;
    //pthread_mutex_unlock(&cout_mutex);
    return NULL;
}


void common_neighbors_master_thread(int number_of_edges, int cn_threshold){
    getrusage(RUSAGE_SELF,&(cn_master_stats.first));
    
    pthread_t threads[NUM_OF_PTHREADS];
    cn_thread_args *th_args = (cn_thread_args*) malloc((NUM_OF_PTHREADS) * sizeof(cn_thread_args));

    pair<vertex_iter , vertex_iter > vp = vertices(my_graph);
    vertex_iter source = vp.first;
    
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        th_args[i].thread_id = i;
        th_args[i].cn_threshold = cn_threshold;
        th_args[i].source_vertex = source;
        pthread_create(&threads[i], NULL, common_neighbors_slave_thread, (void*) &th_args[i]);
        source ++;
    }
    
    getrusage(RUSAGE_SELF,&(cn_master_stats.second));
    
    for(int j = 0; j < NUM_OF_PTHREADS; j++){
        pthread_join(threads[j], NULL);
    }
    free(th_args);
}

//============================================================================
// Functions About Clustering Coefficient Process
void * clustering_coefficient_slave_thread(void *arg){
    cc_thread_args * args = (cc_thread_args *) arg;
    int thread_id = args->thread_id;
    int iter = args->iter;
    getrusage(RUSAGE_SELF,&(cc_slave_stats[iter][thread_id].first));
    vector<vertex_desc> vertices = *(args->vertices_pointer);
    
    while (!vertices.empty()) {
        vertex_desc vd = vertices.back();
        vertices.pop_back();
        out_edge_iter oei, oei_end;
        adjac_iter ai, ai_end;
        hint_access(vd);
        int num_of_neighbors = 0, num_of_triangles = 0;
        
        for ( tie(ai,ai_end) = adjacent_vertices(vd,my_graph) ; ai != ai_end ; ++ai ) {
            hint_access(*ai);
            num_of_neighbors++;
            adjac_iter ai2, ai2_end;
            
            for ( tie(ai2,ai2_end) = adjacent_vertices(*ai,my_graph) ; ai2 != ai2_end ; ++ai2 ) {
                if (boost::edge(vd, *ai2, my_graph).second) // 2nd order neighbor has an edge to current vertex
                    num_of_triangles++;
            }
        }
        my_graph[vd].clustering_coefficient = (num_of_neighbors > 1)? ((double)num_of_triangles) / ((double)(num_of_neighbors * (num_of_neighbors - 1))) : 0.0;
    }
    
    getrusage(RUSAGE_SELF,&(cc_slave_stats[iter][thread_id].second));
    return NULL;
}

void clustering_coefficient_master_thread(int number_of_vertices, int iter){
    getrusage(RUSAGE_SELF,&(cc_master_stats[iter].first));
    std::pair<vertex_iter, vertex_iter> vp = vertices(my_graph);
    vertex_iter vi = vp.first;
    
    pthread_t threads[NUM_OF_PTHREADS]; //array to hold thread information
    cc_thread_args *th_args = (cc_thread_args*) malloc((NUM_OF_PTHREADS) * sizeof(cc_thread_args));
    
    vector<vertex_desc> verticesArr [NUM_OF_PTHREADS];
    
    for ( int i = 0 ;  i < NUM_OF_PTHREADS; i++){
        int chunk_size = ceil(((double) number_of_vertices ) / (NUM_OF_PTHREADS-i));
        number_of_vertices -= chunk_size;
        //cout << "Thread " << i+1 << " is being started by " << chunk_size << " vertices." << endl;
        while (chunk_size > 0 && vi < vp.second) {
            if (!my_graph[*vi].noedge){
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
                    if (!my_graph[*ai].assigned && !my_graph[*ai].noedge){
                        verticesArr[i].push_back(*ai);
                        my_graph[*ai].assigned = true;
                        chunk_size --;
                    }
                    ai ++;
                }
            }
            
            while ((my_graph[*vi].assigned || my_graph[*vi].noedge) && vi < vp.second) {
                vi++;
            }
        }
        th_args[i].thread_id = i;
        th_args[i].vertices_pointer =  &(verticesArr[i]);
        th_args[i].iter = iter;
    }

    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_create(&threads[i], NULL, clustering_coefficient_slave_thread, (void*) &th_args[i]);
    }
    getrusage(RUSAGE_SELF,&(cc_master_stats[iter].second));
    for(int i = 0; i < NUM_OF_PTHREADS; i++){
        pthread_join(threads[i], NULL);
    }
    free(th_args);
    
    /* Log results
    pair<vertex_iter, vertex_iter> viter = vertices(my_graph);
    for (vertex_iter i = viter.first; i != viter.second; i++) {
        if (my_graph[*i].clustering_coefficient > 0.1)
            cout << *i << ":" << my_graph[*i].clustering_coefficient<<endl;
    }
    // */
}



//============================================================================
// File Reading and Graph Creating Functions
int edgefile_graph_size(int timestamp, int data_size){
    
    ostringstream full_path("");
    full_path << DIRECTORY << "enron-" << timestamp << "-" << data_size << ".txt";
    ifstream file(full_path.str());
    //pthread_mutex_lock(&cout_mutex);
    //cout << "File: " << full_path.str() << endl;
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

int create_graph(int timestamp, int data_size){
    string file_prefix = FILE_PREFIX;
    ostringstream full_path("");
    full_path << DIRECTORY << "enron-" << timestamp << "-" << data_size << ".txt";
    ifstream file(full_path.str());
    
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
        
        numof_noedge_vertices = 0;
        for (vertex_iter it = vp.first ; it < vp.second ; it++){
            my_graph[*it].assigned = false;
            //my_graph[*it].paired.store(false,std::memory_order_relaxed);
            my_graph[*it].paired = false;
            adjac_iter ai, ai_end;
            tie(ai,ai_end) = adjacent_vertices(*it,my_graph);
            if (ai == ai_end){
                my_graph[*it].noedge = true;
                numof_noedge_vertices ++;
            }
            else
                my_graph[*it].noedge = false;
            
        }
    }
    else{
        pthread_mutex_lock(&cout_mutex);
        cout << "File doesn't exist" << endl;
        pthread_mutex_unlock(&cout_mutex);
    }
    
    graph_size.second = num_of_edges;
    
    return num_of_edges;
}


void hint_access(vertex_desc vd) {
    while ( !(my_graph[vd].relatedNodes.empty()) ){
        _mm_prefetch((char *) &my_graph.out_edge_list( my_graph[vd].relatedNodes.front().second) , _MM_HINT_T1);
        pop_heap(my_graph[vd].relatedNodes.begin(),my_graph[vd].relatedNodes.end());
        my_graph[vd].relatedNodes.pop_back();
    }
}