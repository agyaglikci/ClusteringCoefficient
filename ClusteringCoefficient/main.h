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
#include "gtest/gtest.h"


#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>

#define DIRECTORY "/Users/agyaglikci/Workspace/ClusteringCoefficient/ClusteringCoefficient/"
#define FILE_PREFIX "enron-"

#define DATA_SIZE atoi(argv[1])
#define START_TIMESTAMP atoi(argv[2])
#define STOP_TIMESTAMP atoi(argv[3])
#define NUM_OF_PTHREADS 4
#define CN_THRESHOLD atoi(argv[4])

struct vertex;
struct edge;

//typedef boost::adjacency_matrix<boost::undirectedS, vertex, edge> Graph;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, vertex, edge > Graph;
std::pair<int,int> graph_size;
int numof_noedge_vertices;
int max_link_predictions;
int numof_link_predictions [57];

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
    bool noedge = false;
    double clustering_coefficient = 0 ;
    std::vector< std::pair<int,vertex_desc> > relatedNodes;
};

struct edge {
    int  common_neighbors;
};

int  edgefile_graph_size(int timestamp, int data_size);
int  create_graph(int timestamp, int data_size);

void * common_neighbors_slave_thread(void *arg);
void common_neighbors_master_thread(int number_of_edge, int cn_threshold);

void * clustering_coefficient_slave_thread(void *arg);
void clustering_coefficient_master_thread(int number_of_vertices, int iter);

void hint_access(vertex_desc);

//============================================================================
// THREAD ARGUMENTS

class cn_thread_args {
public:
    int thread_id;
    int cn_threshold;
    std::vector<std::pair<vertex_desc, vertex_desc> > * vertex_pairs_ptr ;
};



class cc_thread_args {
public:
    std::vector<vertex_desc> * vertices_pointer;
    int thread_id;
    int iter;
};

//
// STATISTICS STRUCTURES
std::pair<rusage,rusage> cn_master_stats;
std::pair<rusage,rusage> cn_slave_stats [NUM_OF_PTHREADS];
std::pair<rusage,rusage> cc_master_stats [3];
std::pair<rusage,rusage> cc_slave_stats [3][NUM_OF_PTHREADS];

// STATISTICS FUNCTIONS
void log_resource_stats( std::string title, std::pair<rusage,rusage> stats);

/*
 struct rusage {
    struct timeval ru_utime; // user time used
    struct timeval ru_stime; // system time used
    long ru_maxrss;          // max resident set size
    long ru_ixrss;           // integral shared text memory size
    long ru_idrss;           // integral unshared data size
    long ru_isrss;           // integral unshared stack size
    long ru_minflt;          // page reclaims
    long ru_majflt;          // page faults
    long ru_nswap;           // swaps
    long ru_inblock;         // block input operations
    long ru_oublock;         // block output operations
    long ru_msgsnd;          // messages sent
    long ru_msgrcv;          // messages received
    long ru_nsignals;        // signals received
    long ru_nvcsw;           // voluntary context switches
    long ru_nivcsw;          // involuntary context switches
 };
*/

#endif /* defined(__ClusteringCoefficient__ClusteringCoefficient__) */