#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/fuzzy_pattern_matching.hpp>
#include <havoqgt/graph.hpp>
//#include <havoqgt/label_propagation_pattern_matching.hpp>
#include <havoqgt/label_propagation_pattern_matching_bsp.hpp> 
#include <havoqgt/pattern_util.hpp>
#include <havoqgt/token_passing_pattern_matching.hpp>
#include <havoqgt/vertex_data_db.hpp>
#include <havoqgt/vertex_data_db_degree.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
using namespace havoqgt;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;

template<typename T>
  using SegmentAllocator = bip::allocator<T, segment_manager_t>;

int main(int argc, char** argv) {
  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
  havoqgt::havoqgt_init(&argc, &argv);
  {

  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);  

  std::string graph_input = argv[1];
//  std::string backup_filename; 

  // for fuzzy pattern matching
  std::string vertex_data_input_filename = argv[2];
  //std::string pattern_input_filename = argv[3];
  std::string pattern_dir = argv[3];
  std::string vertex_rank_output_filename = argv[4];
  std::string backup_filename = argv[5];
  std::string result_dir = argv[6];
  bool use_degree_as_vertex_data = std::stoull(argv[7]); // 1 - yes, 0 - no   

  // TODO: parse commandline

  MPI_Barrier(MPI_COMM_WORLD);

  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  //std::cout << "Graph input source : " << graph_input << std::endl; // Test

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  //segment_manager_t* segment_manager = ddb.get_segment_manager();
  //  bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  //graph_type *graph = segment_manager->
  //  find<graph_type>("graph_obj").first;
  //assert(graph != nullptr);
  
  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr); 

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }

  //graph->print_graph_statistics(); // causes MPI error
  MPI_Barrier(MPI_COMM_WORLD);  

  // fuzzy pattern matching
  {

  // types used by the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  //typedef typename graph_type::edge_iterator eitr_type;

  // user defined types
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData; // for string hash
  //typedef uint8_t VertexData; // for log binning  
  typedef uint64_t VertexRankType;

  //typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetaData;
  //typedef graph_type::vertex_data<VertexRankType, SegmentAllocator<VertexRankType> > VertexRank;
  //typedef graph_type::vertex_data<bool, SegmentAllocator<bool> > VertexActive;
  //typedef graph_type::vertex_data<uint64_t, SegmentAllocator<uint64_t> > VertexIteration;

  typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetaData;
  typedef graph_type::vertex_data<VertexRankType, std::allocator<VertexRankType> > VertexRank;
  typedef graph_type::vertex_data<bool, std::allocator<bool> > VertexActive;
  typedef graph_type::vertex_data<uint64_t, std::allocator<uint64_t> > VertexIteration;

  typedef vertex_state<uint64_t> VertexState;
  typedef std::unordered_map<Vertex, VertexState> VertexStateMap;

  if(mpi_rank == 0) {
    std::cout << "Distributed fuzzy pattern matching." << std::endl;
  }

  double time_start = MPI_Wtime();
  double time_end = MPI_Wtime();

  // per-rank containers
  VertexStateMap vertex_state_map; 

  // vertex containers
//  VertexMetaData vertex_metadata(*graph, alloc_inst);
//  VertexRank vertex_rank(*graph, alloc_inst);
//  VertexActive vertex_active(*graph, alloc_inst);
//  VertexIteration vertex_iteration(*graph, alloc_inst);
// TODO: need a new alloc_inst to use bip/mmap
 
  VertexMetaData vertex_metadata(*graph); 
  VertexRank vertex_rank(*graph);
  VertexActive vertex_active(*graph);
  VertexIteration vertex_iteration(*graph);

  // application parameters
 
  // token passing types
  size_t token_passing_algo = 0; // TODO: use enum if this stays

  MPI_Barrier(MPI_COMM_WORLD);
 
  ///////////////////////////////////////////////////////////////////////////// 

  // build the distributed vertex data db
  // each rank reads 10K lines at a time
  time_start = MPI_Wtime();

  if (use_degree_as_vertex_data) {
    vertex_data_db_degree<graph_type, VertexMetaData, Vertex, VertexData>
      (graph, vertex_metadata);   
  } else { 
    vertex_data_db<graph_type, VertexMetaData, Vertex, VertexData>
      (graph, vertex_metadata, vertex_data_input_filename, 10000);
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | Vertex Data DB : " << time_end - time_start << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////////// 

  // test print
  //if(mpi_rank == 0) {
/*  int set_mpi_rank = 0;
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    //if (mpi_rank == set_mpi_rank)
      std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_metadata[vertex] << std::endl; 
  }
  
  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    //if (mpi_rank == set_mpi_rank)
      std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_metadata[vertex] << std::endl; 
  } 
  //}

  if (use_degree_as_vertex_data) {
    return 0;
  }*/
  // test print
  
  // result
  std::string pattern_set_result_filename = result_dir + "/result_pattern_set"; 
  std::ofstream pattern_set_result_file(pattern_set_result_filename, std::ofstream::out);
    
  // setup pattern set 
  // a pattern set is collection of directories containing pattern files  
   
  // loop over pattern set
  for(size_t ps = 0; ps < 1; ps++) { // TODO: for now only reading from pattern_dir/0 
  // beginning of an elemnet in the pattern set
   
  // setup pattern to search
  if(mpi_rank == 0) { 
    std::cout << "Setting up pattern to search ... " << std::endl;
  }
   
  // setup pattern - sub-graph
  std::string pattern_input_filename = pattern_dir + "/" + std::to_string(ps) + "/pattern";

  typedef ::graph<Vertex, Edge, VertexData> PatternGraph; // TODO: fix graph class name conflict
  PatternGraph pattern_graph(pattern_input_filename + "_edge",
    pattern_input_filename + "_vertex", 
    pattern_input_filename + "_vertex_data", 
    pattern_input_filename + "_stat",
    false, false);

  // test print
  if(mpi_rank == 0) {
  std::cout << "Searching pattern : " << std::endl;
  for (auto v = 0; v < pattern_graph.vertex_count; v++) {
    std::cout << v << " : off-set " << pattern_graph.vertices[v] << " vertex_data " 
    << pattern_graph.vertex_data[v] << " vertex_degree " 
    << pattern_graph.vertex_degree[v] << std::endl;
    std::cout << " neighbours : "; 
    for (auto e = pattern_graph.vertices[v]; e < pattern_graph.vertices[v + 1]; e++) {
    auto v_nbr = pattern_graph.edges[e];
      std::cout << v_nbr << ", " ;
    }
    std::cout << std::endl; 
  }
  std::cout << "diameter : " << pattern_graph.diameter << std::endl; 
  }
  // test print

  // TODO: remove from here    
  // setup pattern - token passing
  pattern_util<VertexData, Vertex> ptrn_util_two(pattern_input_filename, true);
  auto pattern = std::get<0>(ptrn_util_two.input_patterns[0]);
  auto pattern_indices = std::get<1>(ptrn_util_two.input_patterns[0]);

  // initialize containers - per-pattern
  vertex_state_map.clear(); // important
  vertex_rank.reset(0);
  vertex_active.reset(true);
  vertex_iteration.reset(0); // TODO: -1 ?

  // std::cout << "Vertex state map size (initially): " << vertex_state_map.size() << std::endl; // Test 
  
  bool global_initstep = true;  
  bool global_not_finished = false;

  uint64_t global_itr_count = 0;

  // result
  std::string itr_result_filename = result_dir + "/" + std::to_string(ps) + "/result_itr";
  std::ofstream itr_result_file(itr_result_filename, std::ofstream::out);

  std::string step_result_filename = result_dir + "/" + std::to_string(ps) + "/result_step";
  std::ofstream step_result_file(step_result_filename, std::ofstream::out); 

  std::string superstep_result_filename = result_dir + "/" + std::to_string(ps) + "/result_superstep";
  std::ofstream superstep_result_file(superstep_result_filename, std::ofstream::out);

  std::string active_vertices_count_result_filename = result_dir + "/" + 
    std::to_string(ps) + "/all_ranks_active_vertices_count/active_vertices_" + std::to_string(mpi_rank); 
  std::ofstream active_vertices_count_result_file(active_vertices_count_result_filename, std::ofstream::out);

  std::string active_vertices_result_filename = result_dir + "/" +
    std::to_string(ps) + "/all_ranks_active_vertices/active_vertices_" + std::to_string(mpi_rank);
  std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);

  MPI_Barrier(MPI_COMM_WORLD);

  double pattern_time_start = MPI_Wtime();

  // run application

  do {

  global_not_finished = false;

  double itr_time_start = MPI_Wtime();

  /////////////////////////////////////////////////////////////////////////////
   
  double label_propagation_time_start = MPI_Wtime();

  // clone pattern matchng
  //fuzzy_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank);

  // label propagation pattern matching 
//  label_propagation_pattern_matching<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    VertexRank, VertexActive, VertexIteration, VertexStateMap, PatternGraph>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph);

  // label propagation pattern matching bsp 
  label_propagation_pattern_matching_bsp<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
    VertexRank, VertexActive, VertexIteration, VertexStateMap, PatternGraph>
    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
    vertex_iteration, vertex_state_map, pattern_graph, global_initstep, global_not_finished, 
    global_itr_count, superstep_result_file, active_vertices_count_result_file);

  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
  double label_propagation_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | Label Propagation : " 
      << label_propagation_time_end - label_propagation_time_start << std::endl;
  }

  if(mpi_rank == 0) {
    step_result_file << global_itr_count << ", LP, "
      << (label_propagation_time_end - label_propagation_time_start) << "\n";
  }

  if (global_initstep) {
    global_initstep = false;
  }

  // verify global termination condition
  //std::cout << "Fuzzy Pattern Matching | Global Not Finished status (local) : " << global_not_finished << std::endl; // Test
  // global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work
  global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching | Global Finished status : "; 
    if (global_not_finished) { 
      std::cout << "Continue" << std::endl;
    } else {
      std::cout << "Stop" << std::endl;
    } 
  }

  // std::cout << "Vertex state map size (finally): " << vertex_state_map.size() << std::endl; // Test 

  // test print
  uint64_t vertex_active_count = 0;
  uint64_t vertex_inactive_count = 0;
  uint64_t vertex_iteration_valid_count = 0;  
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_iteration[vertex] >= 2*pattern_graph.vertex_data.size() + 1) { 
      //std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_iteration[vertex] << std::endl;
      vertex_iteration_valid_count++;
    }   
    if (vertex_rank[vertex] > 0) {
      std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl;
    }
    if (!vertex_active[vertex]) {
      //std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_active[vertex] << std::endl;
      vertex_inactive_count++;
    } else {
      vertex_active_count++;
    }
  }
  
  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_iteration[vertex] >= 2*pattern_graph.vertex_data.size() + 1) {
      //std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_iteration[vertex] << std::endl;
      vertex_iteration_valid_count++;
    }
    if (vertex_rank[vertex] > 0) {
      std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl; 
    }  
    if (!vertex_active[vertex]) {   
      //std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_active[vertex] << std::endl;
      vertex_inactive_count++;
    } else {
      vertex_active_count++;
    }
  }
  //  std::cout << mpi_rank << " # active vertices " << vertex_active_count << std::endl;
  //  std::cout << mpi_rank << " # inactive vertices " << vertex_inactive_count << std::endl; 
  //  std::cout << mpi_rank << " # vertices reached max-iterations " << vertex_iteration_valid_count << std::endl;

  // test print
  
  
  // test print
  //for (auto& v : vertex_state_map) {
  //  auto v_locator = graph->label_to_locator(v.first); 
  //  std::cout << v.first << " " << v.second.vertex_pattern_index << " " << vertex_metadata[v_locator] << std::endl; 
  //}
  // test print
 
  /////////////////////////////////////////////////////////////////////////////
 
  // toekn passing
  double token_passing_time_start = MPI_Wtime();   

  if ((token_passing_algo == 0) /*&& global_not_finished*/) { // do token passing ? 
  global_not_finished = false;  

  typedef std::unordered_map<Vertex, bool> TokenSourceMap; 
  TokenSourceMap token_source_map;

  //std::vector<bool> pattern_found(ptrn_util_two.input_patterns.size(), false); 
  // TODO: bool does not work with mpi_all_reduce_inplace   
  std::vector<uint8_t> pattern_found(ptrn_util_two.input_patterns.size(), 0); 

  // loop over the patterns and run token passing
  for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) 
  {
 
  auto pattern_tp = std::get<0>(ptrn_util_two.input_patterns[pl]);
  auto pattern_indices_tp = std::get<1>(ptrn_util_two.input_patterns[pl]);
  auto pattern_cycle_length_tp = std::get<2>(ptrn_util_two.input_patterns[pl]); // uint
  auto pattern_valid_cycle_tp = std::get<3>(ptrn_util_two.input_patterns[pl]); // boolean
  
  if(mpi_rank == 0) {
    std::cout << "Token Passing [" << pl << "] | Searching pattern: ";
    pattern_util<VertexData>::output_pattern(pattern_tp);
  }
    
  time_start = MPI_Wtime();

  token_passing_pattern_matching(graph, vertex_metadata, pattern_tp, 
    pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map, 
    token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp, pattern_found[pl]);
 
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?    
  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | Token Passing [" << pl << "] : " << time_end - time_start << std::endl;
  }

  if(mpi_rank == 0) {
    superstep_result_file << global_itr_count << ", TP, "
      << pl << ", "
      << time_end - time_start << "\n";      
  } 

  } // for - loop ove the patterns

  // pattren found
  havoqgt::mpi::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
  if(mpi_rank == 0) {   
    for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
      std::string s = pattern_found[pl] == 1 ? "true" : "false";
      std::cout << "Token Passing [" << pl << "] | Found pattern : " << s << std::endl;
    }
  }

  // remove the invalid vertices from the vertex_state_map

  // TODO: In the case, a vertex is on multiple cycles/chains (not as the token source)
  // only invalidate it as a token source, but do not remove it from the vertex_state_map 

  //std::cout << "token_source_map size " << token_source_map.size() << std::endl; // Test
  uint64_t remove_count = 0; 
  for (auto& s : token_source_map) {
    if (!s.second) {
      if (vertex_state_map.erase(s.first) < 1) { // s.first is the vertex
        // also remove s.first from the token_source_map
        std::cerr << "Error: failed to remove an element from the map."
          << std::endl;
      } else {
        //std::cout << s.first << " was removed from token_source_map" << std::endl; // Test
        remove_count++;
        vertex_active[graph->label_to_locator(s.first)] = false; 
        if (!global_not_finished) {
          global_not_finished = true;
        } 
      }
    } else {
      //std::cout << "token_source_map " << s.first << " is " << s.second << std::endl; // Test
    }
  } // for - token_source_map

  if(mpi_rank == 0) {
    std::cout << "Token Passing | Removed " << remove_count << " vertices."<< std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

  // result
  // Important : This may slow down things -only for presenting results
  uint64_t active_vertices_count = 0;
  for (auto& v : vertex_state_map) {
    auto v_locator = graph->label_to_locator(v.first);
    if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
      active_vertices_count++;
    } else if (!v_locator.is_delegate()) {
      active_vertices_count++;
    }
  }

  active_vertices_count_result_file << global_itr_count << ", TP, "
    << "0, "  
    << active_vertices_count << "\n";    
 
  } else {
    std::cout << "Fuzzy Pattern Matching | Skipping Token Passing." << std::endl;
  } // do token passing ? 

  // verify global termination condition  
  //std::cout << "Fuzzy Pattern Matching | Global Not Finished status (local) : " << global_not_finished << std::endl; // Test
  //global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work  
  global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching | Global Finished status : ";
    if (global_not_finished) {
      std::cout << "Continue" << std::endl;
    } else {
      std::cout << "Stop" << std::endl;
    }
  }

  double token_passing_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | Token Passing : " 
      << (token_passing_time_end - token_passing_time_start) << std::endl;
  }

  if(mpi_rank == 0) {
    step_result_file << global_itr_count << ", TP, "  
      << (token_passing_time_end - token_passing_time_start) << "\n"; 
  }

  /////////////////////////////////////////////////////////////////////////////
  double itr_time_end = MPI_Wtime(); 
  if(mpi_rank == 0) { 
    std::cout << "Fuzzy Pattern Matching Time | Iteration [" 
      << global_itr_count << "] " << itr_time_end - itr_time_start << std::endl;  
  }
  
  if(mpi_rank == 0) {
    // iteration number, time
    itr_result_file << global_itr_count << ", "
      << (itr_time_end - itr_time_start) << "\n";
  }

  global_itr_count++;

  // verify global termination condition
  //if (global_itr_count >= 2) { // Test
  //  global_not_finished = false;
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  } while (global_not_finished); // application loop

  MPI_Barrier(MPI_COMM_WORLD);  
  double pattern_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | Pattern " << ps << " : " << pattern_time_end - pattern_time_start << std::endl;
  }
   
  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching | # Iterations : " << global_itr_count << std::endl;
  }

  // result
  if(mpi_rank == 0) {  
    // pattern set element ID, number of ranks, total number of iterations (lp-tp), time, 
    // #edges in the pattern, #vertices in the pattern,  #token passing paths 
    pattern_set_result_file << ps << ", " 
      << mpi_size << ", "
      << global_itr_count << ", " 
      << (pattern_time_end - pattern_time_start) << ", "
      << pattern_graph.edge_count << ", " 
      << pattern_graph.vertex_count << ", "
      << ptrn_util_two.input_patterns.size() << "\n"; 
  }

  // Important : This may slow things down -only for presenting results

  for (auto& v : vertex_state_map) {
    auto v_locator = graph->label_to_locator(v.first);
    if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
      active_vertices_result_file << mpi_rank << ", "
        << v.first << ", " 
        << v.second.vertex_pattern_index << "\n";
    } else if (!v_locator.is_delegate()) {
      active_vertices_result_file << mpi_rank << ", "
        << v.first << ", "
        << v.second.vertex_pattern_index << "\n";
    }
  }

  // cleanup memeory
  //vertex_rank.clear(); // TODO: add clear() method to vertex_data.cpp   
  //vertex_active.clear(); // TODO: add clear() method to vertex_data.cpp
  //vertex_state_map.clear();

  // close files
  itr_result_file.close();
  step_result_file.close();
  superstep_result_file.close();
  active_vertices_count_result_file.close(); 
  active_vertices_result_file.close();

  // end of an elemnet in the pattern set
  } // for - loop over pattern set

  pattern_set_result_file.close(); // close file

  } // fuzzy pattern matching  

  } // havoqgt_init
  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;  
} // main
