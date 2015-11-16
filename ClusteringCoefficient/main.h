#ifndef __ClusteringCoefficient__ClusteringCoefficient__
#define __ClusteringCoefficient__ClusteringCoefficient__

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <pthread.h>
#include <stdint.h>
#include <getopt.h>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <list>
#include <xmmintrin.h>


#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>

#define DIRECTORY "/Users/agyaglikci/Workspace/ClusteringCoefficient/ClusteringCoefficient/"
#define FILE_PREFIX "enron-"
#define FILE_EXTENSION ".txt"
#define NUM_OF_VERTICES 22430
#define TIMESTAMP 1
#define NUM_OF_PTHREADS 4
#define VERTEX_CHUNK_SIZE 5608//1402 // ceil(NUM_OF_VERTICES/NUM_OF_PTHREADS)
#define MAX_NEIGHBORS 608
//22430 vertices are there in node pairs file.
//608 edges are there in node pairs file.
struct vertex;
struct edge;

//typedef boost::adjacency_matrix<boost::undirectedS, vertex, edge> Graph;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, vertex, edge > Graph;
std::pair<int,int> graph_size;

template <typename G>
struct graph_traits {
    typedef typename G::vertex_descriptor   vertex_descriptor;
    typedef typename G::vertex_iterator     vertex_iterator;
    typedef typename G::out_edge_iterator   out_edge_iterator;
    
    typedef typename G::edge_descriptor     edge_descriptor;
    typedef typename G::edge_iterator       edge_iterator;
    typedef typename G::adjacency_iterator  adjacency_iterator;
};

typedef graph_traits<Graph>::vertex_descriptor  vertex_desc;
typedef graph_traits<Graph>::vertex_iterator    vertex_iter;
typedef graph_traits<Graph>::adjacency_iterator adjac_iter;

typedef graph_traits<Graph>::out_edge_iterator  out_edge_iter;
typedef graph_traits<Graph>::edge_descriptor    edge_desc;
typedef graph_traits<Graph>::edge_iterator      edge_iter;

struct __attribute__ ((packed)) hint {
    int weight = 0;
    vertex_desc v;
    void * pointer;
} ;

struct vertex   {
public:
    bool assigned = false;
    int index = 0;
    int related_node_top_index = 0;
    double clustering_coefficient;
    hint relatedNodes [MAX_NEIGHBORS];
};

struct edge {
    int  common_neighbors;
};

int  edgefile_graph_size(short timestamp);
int  create_graph(short timestamp);

void * common_neighbors_slave_thread(void *arg);
void common_neighbors_master_thread(int number_of_edges);

void * clustering_coefficient_slave_thread(void *arg);
void clustering_coefficient_master_thread(int number_of_vertices);

void hint_access(vertex_desc , int);

//============================================================================
// THREAD ARGUMENTS
class cn_thread_args {
public:
    edge_iter ei;
    int chunk_size;
    int thread_id;
};

class cc_thread_args {
public:
    vertex_desc vertices [VERTEX_CHUNK_SIZE];
    int vector_size;
    int thread_id;
};

#endif /* defined(__ClusteringCoefficient__ClusteringCoefficient__) */