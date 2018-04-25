#ifndef HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_NONUNIQUE_HPP_INCLUDED
#define HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_NONUNIQUE_HPP_INCLUDED

#include <bitset>

#include <boost/dynamic_bitset.hpp>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

# define output_result

namespace havoqgt { namespace mpi {

//static constexpr size_t max_bit_vector_size = 16;

template<typename IntegralType>
class vertex_state {
public:
  vertex_state() :
  //is_active(false),
  vertex_pattern_index(0) {}

  //bool is_active;
  size_t vertex_pattern_index; // TODO: change type
  std::unordered_map<size_t, IntegralType> pattern_vertex_itr_count_map;
  // TODO: not itr_count anymore, more like true / false 
};

/*template<typename VertexData, typename DynamicBitSet, typename IntegralType>
class vertex_state_generic_A {
public:
  vertex_state_generic_A(boost::dynamic_bitset<>::size_type bit_set_size) : 
  vertex_pattern_index(0) {
    template_vertices.resize(bit_set_size);
    template_neighbors.resize(bit_set_size);  
  }

  DynamicBitSet template_vertices; 
  DynamicBitSet template_neighbors;
  std::unordered_map<VertexData, IntegralType> 
    template_neighbor_metadata_count_map;
  size_t vertex_pattern_index; // TODO: dummy 
};*/

template<typename VertexData, typename IntegralType, typename BitSet>
class vertex_state_generic {
public:
  vertex_state_generic() :
  vertex_pattern_index(0) {}

  BitSet template_vertices;
  BitSet template_neighbors;
  std::unordered_map<VertexData, IntegralType>
    template_neighbor_metadata_count_map;
  size_t vertex_pattern_index; // TODO: dummy 
};

template<typename Visitor>
class lppm_queue {

public:
  lppm_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    data.pop_back();
  }
 
  Visitor const& top() {
    return data.back();
  } 
  
  size_t size() const {
    return data.size();;
  }

  bool empty() const {
    return data.empty();
  }

  void clear() {
    data.clear();
  }

protected:
  std::vector<Visitor> data;

};

// label propagation pattern matching visitor class
template<typename Graph, typename Vertex, typename VertexData, typename BitSet, typename DynamicBitSet>
class lppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  lppm_visitor() : 
    msg_type(0) {}

  lppm_visitor(vertex_locator _vertex, uint8_t _msg_type = 0) : 
    vertex(_vertex), 
    msg_type(_msg_type) {}

  lppm_visitor(vertex_locator _vertex, size_t _parent_pattern_index,
    uint8_t _msg_type) :
    vertex(_vertex),
    parent_pattern_index(_parent_pattern_index),
    msg_type(_msg_type) {}
 
/*  template<typename ByteVector>
  lppm_visitor(vertex_locator _vertex, ByteVector _parent_template_vertices, 
    uint8_t _msg_type) :
    vertex(_vertex),
    msg_type(_msg_type) {

    std::copy(_parent_template_vertices.begin(), 
      _parent_template_vertices.end(), parent_template_vertices.begin());

    //for (size_t i = 0; i < pattern_graph.vertex_count; i++) { // Test // dummy
    //  if (parent_template_vertices[i] == (uint8_t)1) {
    //    parent_pattern_index = i;
    //    break;
    //  }   
    //}
 
  }*/

  //template<typename BitSet>
  lppm_visitor(vertex_locator _vertex, BitSet _parent_template_vertices, 
    uint8_t _msg_type) : 
    vertex(_vertex), 
    parent_template_vertices_bitset(_parent_template_vertices),  
    msg_type(_msg_type) {}  

  ~lppm_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep 
    auto g = std::get<10>(alg_data); 
    // std::get<11>(alg_data) - template_vertices
 
    bool match_found = false;
    bool valid_parent_found = false;

//--    size_t vertex_pattern_index = 0;

//    DynamicBitSet vertex_template_vertices(pattern_graph.vertex_count);
//    DynamicBitSet parent_template_vertices_bitset(pattern_graph.vertex_count);
    
//    for (size_t i = 0; i < pattern_graph.vertex_count; i++) { // TODO: temporary fix
//      if (parent_template_vertices[i] == (uint8_t)1) {
//        parent_template_vertices_bitset.set(i);
//      }
//    } // for   

//--    BitSet vertex_template_vertices; // initialized with zeros

    int mpi_rank = havoqgt_env()->world_comm().rank();    
 
      // a delegate but not the controller
      if (vertex.is_delegate() && g->master(vertex) != mpi_rank && msg_type == 1) { 
        // the vertex_state is only maintained on the controller

        // first LP superstep of the first iteration 
        if (std::get<8>(alg_data) == 0  && std::get<9>(alg_data)) {
   
        BitSet vertex_template_vertices; // initialized with zeros
        // TODO: avoid figuring out vertex_pattern_index everytime

        // does vertex_data match any entry in the query pattern
        for (size_t vertex_pattern_index = 0;
          vertex_pattern_index < pattern_graph.vertex_data.size();
          vertex_pattern_index++) {
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
            assert(vertex_pattern_index < 16); // Test 
            vertex_template_vertices.set(vertex_pattern_index);
            match_found = true;   

            // verify if heard from a valid neighbor (parent)
            if (msg_type == 1 && match_found) {
              // unique metadata
              //match_found = false;
//              for (auto e = pattern_graph.vertices[vertex_pattern_index];
//                e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
//                if (pattern_graph.edges[e] == parent_pattern_index) {
//                  //match_found = true;
//                  valid_parent_found = true; 
//                  break; 
//                }
//              } // for

              // repeated metadata
              if (parent_template_vertices_bitset.none()) {
                return false;
              } 

              for (size_t i = 0; i < parent_template_vertices_bitset.size(); i++) {
                if (parent_template_vertices_bitset.test(i)) {
                  for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
                    e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {       
                    if (pattern_graph.edges[e] == i) {
                      valid_parent_found = true;    
                      break;
                    }   
                  } // for
                  if (valid_parent_found) {
                    break;  
                  } 
                } // if  
              } // for       
              // repeated metadata 

//              if (!valid_parent_found) { // TODO: This is wrong! Moved it outside the for loop.
//                return false; 
//              }    

            } // if            
 
          } // if

          if (valid_parent_found) {
            break;
          }
        } // for 

        if (match_found && !valid_parent_found) {
          return false;   
        }  

        // initial case - return true to handle delegates // TODO: try this one too
        //if (std::get<4>(alg_data)[vertex] && msg_type == 0 && !match_found) {
        //  return true;  
        //}

        if (!match_found) {
          std::get<4>(alg_data)[vertex] = false; 
          //return true; // send to the controller
          return false;
        } else {
          std::get<11>(alg_data)[vertex] = vertex_template_vertices.to_ulong();  
          return true; // send to the controller
        }

        } else { // if // redundent code is to avoid memory copy

          if (parent_template_vertices_bitset.none()) {
              return false;
          }

          BitSet vertex_template_vertices(std::get<11>(alg_data)[vertex]);

          if (vertex_template_vertices.none()) {
            //std::get<4>(alg_data)[vertex] = false; // TODO: ?
            return false;
          } else {
            match_found = true; // TODO: 
          }

          // verify if heard from a valid neighbor (parent)
          if (msg_type == 1 && match_found) {

            for (size_t vertex_pattern_index = 0; // TODO: change the name vertex_pattern_index to template_vertex
              vertex_pattern_index < vertex_template_vertices.size();
              vertex_pattern_index++) {
              if (vertex_template_vertices.test(vertex_pattern_index)) {
                match_found = true;
    
                // repeated metadata
                for (size_t i = 0; i < parent_template_vertices_bitset.size(); i++) {
                  if (parent_template_vertices_bitset.test(i)) {
                    for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
                      e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {       
                      if (pattern_graph.edges[e] == i) {
                        valid_parent_found = true;    
                        break;
                      }   
                    } // for
                    if (valid_parent_found) {
                      break;  
                    } 
                  } // if  
                } // for       
                // repeated metadata 
                 
              } // if
              if (valid_parent_found) {
                break;
              }
            } // for 

          } // if
           
          if (match_found && !valid_parent_found) {
            return false;
          }
 
          if (!match_found) {
            //std::get<4>(alg_data)[vertex] = false; // TODO: ?
            //return true; // send to the controller
            return false;
           } else {
            // TODO: std::get<11>(alg_data)[vertex].reset() // if it is a bit vector 
            //std::get<11>(alg_data)[vertex] = 0; // TODO: how could you avoid this store operation? // does not work
            //std::get<11>(alg_data)[vertex] must be updated by the value stored on the controller
            //Probably need a new visitor 
            return true; // send to the controller
           }
        } // else  
 
      } 

      // for local vertex and controller only
 
      // first LP superstep of the first iteration 
      if (std::get<8>(alg_data) == 0  && std::get<9>(alg_data)) {
        BitSet vertex_template_vertices; // initialized with zeros

        for (size_t vertex_pattern_index = 0; 
          vertex_pattern_index < pattern_graph.vertex_data.size(); 
          vertex_pattern_index++) { 
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
             assert(vertex_pattern_index < 16); // Test  
             vertex_template_vertices.set(vertex_pattern_index); 
             match_found = true;
             //break; // Important 
          }
        }

        if (!match_found) { // TODO: controller return true?
          std::get<4>(alg_data)[vertex] = false;
          if (vertex.is_delegate() && g->master(vertex) == mpi_rank) { // controller
            //return true; // to invalidate the delegates
            return false;
          }           
          return false; 
        }

        if (msg_type == 1) {
          if (parent_template_vertices_bitset.none()) {
            return false;
          }

          if (vertex_template_vertices.none()) {
               std::cerr << "Error: no bit is set." << std::endl;
               return false;
          }

          std::get<11>(alg_data)[vertex] = vertex_template_vertices.to_ulong(); 

          verify_and_update_vertex_state(alg_data, vertex_template_vertices);
          return false;        
        }
      } // if

      //if (msg_type == 1 && parent_template_vertices_bitset.none()) {
      //  return false;
      //}
     
      // if the vertex is not in the global map after the first LP 
      // superstep of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) { 
        auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false;
        } else {
//--          vertex_pattern_index = find_vertex->second.vertex_pattern_index;
          //vertex_template_vertices.set(vertex_pattern_index);
//--          vertex_template_vertices = find_vertex->second.template_vertices;         
          if (msg_type == 1) {
            if (parent_template_vertices_bitset.none()) {
              return false;
            }
  
            if (find_vertex->second.template_vertices.none()) {
              std::cerr << "Error: no bit is set." << std::endl;
              return false;
            }

            //verify_and_update_vertex_state(alg_data, find_vertex->second.template_vertices);
            BitSet vertex_template_vertices(std::get<11>(alg_data)[vertex]);
            if (vertex_template_vertices.none()) {
              std::cerr << "Error: no bit is set." << std::endl;
              return false;
            }             
         
            verify_and_update_vertex_state(alg_data, vertex_template_vertices);
            return false;  
          } 
        }   
      } // if
 
//      if (msg_type == 1) {
        //verify_and_update_vertex_state(alg_data, vertex_pattern_index);

        //for (size_t p = 0; p < patter_graph.vertex_count; p++) { // Test
        //  if (parent_template_vertices[p] == (uint8_t)1) {
        //    std::cout << (uint32_t)parent_template_vertices[p] << " "; //<< std::endl; 
        //  }   
        //} 
        //std::cout << std::endl;  

//        verify_and_update_vertex_state(alg_data, vertex_template_vertices);        
//      }
   
      return false;     
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {    
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank(); 

    // Important : skip this verification for the delegates as the 
    // vertex_state is only maintained on the contrller
    if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
      // if the vertex is not in the global map after the first LP superstep 
      // of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) {
        auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false; 
        }
      } 
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    // std::get<6>(alg_data) // vertex state map
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep
    // std::get<10>(alg_data) - g
    // std::get<11>(alg_data) - template_vertices 

    // does vertex_data match an entry in the query pattern
    bool match_found = false;

    // TODO: do you want to compute this every time or store in the memory? Overhead is not noticable though.
//    std::vector<size_t> vertex_pattern_indices(0); // a vertex label could be a match for multiple pattern labels

//    std::array<uint8_t, max_bit_vector_size> vertex_template_vertices_array;
//    vertex_template_vertices_array.fill((uint8_t)0);

//    DynamicBitSet vertex_template_vertices(pattern_graph.vertex_count);

//--    BitSet vertex_template_vertices; // initialized with zeros
//    size_t vertex_pattern_index_tmp = 0; // Test  

    // figure out the template vertices (of 'vertex') in the first LP superstep 
    // of the first iteration and store 
    if (std::get<8>(alg_data) == 0  && std::get<9>(alg_data)) {

    BitSet vertex_template_vertices;

    for (size_t vertex_pattern_index = 0; 
      vertex_pattern_index < pattern_graph.vertex_data.size(); 
      vertex_pattern_index++) { 
      if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
//        vertex_pattern_indices.push_back(vertex_pattern_index);
        // TODO: compare with the entry in pattern_indices to detect loop or 
        // use token passing
//        vertex_template_vertices.set(vertex_pattern_index);
        assert(vertex_pattern_index < 16); // Test    
        vertex_template_vertices.set(vertex_pattern_index);
        match_found = true;
//        vertex_pattern_index_tmp = vertex_pattern_index; // Test 
        //break; 
      }       
    } // for

    // TODO: if all the bits in vertex_template_vertices are zero then the 
    // vertex is invalid; may not need match_found or vertex_active
    if (!match_found) { 
      std::get<4>(alg_data)[vertex] = false;
      return false;
      //return true; // TODO: ask Roger?
    } else {
      std::get<11>(alg_data)[vertex] = vertex_template_vertices.to_ulong();

      //std::cout << vertex_data << " " //<< vertex_pattern_index_tmp
      //    << " " << vertex_template_vertices << " " 
      //    << std::get<11>(alg_data)[vertex] << std::endl; // Test
  
      if (msg_type == 0 && match_found) {
        // send to all the neighbors
        for(eitr_type eitr = g.edges_begin(vertex);
          eitr != g.edges_end(vertex); ++eitr) {
          vertex_locator neighbor = eitr.target();
          lppm_visitor new_visitor(neighbor, vertex_template_vertices, 1);
          vis_queue->queue_visitor(new_visitor); 
        } // for 
        return true;
      } else if (msg_type == 1 && match_found) {
        // must go all the way to the controller 
        //return true; // false?
        return false; // TODO: ask Roger?
      } else {
        return false;
      }      
    }

    } else { // if // redundent code is to avoid memory copy // TODO: vertex_template_vertices pointer ? 
     
    BitSet vertex_template_vertices(std::get<11>(alg_data)[vertex]);   

    if (vertex_template_vertices.none()) {
      //std::get<4>(alg_data)[vertex] = false; // TODO: ?
      return false; 
    } else {
      match_found = true; // TODO 
    }  

    if (msg_type == 0 && match_found) {
      // send to all the neighbors
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();

        // TODO: only handling undirected grpahs

        //for (auto vertex_pattern_index : vertex_pattern_indices) {
          // do this for all the pattern indices for this vertex
          //lppm_visitor new_visitor(neighbor, vertex_pattern_index, 1);
          //vis_queue->queue_visitor(new_visitor);
        //} // for

//        for (size_t i = 0; i < vertex_template_vertices.size(); i++) { // TODO: temporary fix
//          if (vertex_template_vertices.test(i)) {
            //lppm_visitor new_visitor(neighbor, i, 1);
            //vis_queue->queue_visitor(new_visitor);
//            vertex_template_vertices_array[i] = (uint8_t)1;
//          } else {vertex_template_vertices_array[i] = (uint8_t)0;}      
//        }
        //std::cout << "\n"; 
        
        //for (size_t i = 0; i < vertex_template_vertices_array.size(); i++) { // Test
        //  std::cout << (uint32_t)vertex_template_vertices_array[i] << " ";
        //}
        //std::cout << "\n";       
        
        //std::cout << vertex_data << " " //<< vertex_pattern_index_tmp
        //  << " " << vertex_template_vertices << " "
        //  << std::get<11>(alg_data)[vertex] << std::endl; // Test
 
        //lppm_visitor new_visitor(neighbor, vertex_template_vertices_array, 1);
        lppm_visitor new_visitor(neighbor, vertex_template_vertices, 1); 
        vis_queue->queue_visitor(new_visitor);

      } // for
      return true;
    } else if (msg_type == 1 && match_found) {        
      // must go all the way to the controller 
      //return true; // false?
      return false; // TODO: ask Roger?
    } else {
      return false; 
    } 

    } //else

    return true;
  }

  friend inline bool operator>(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  /*template<typename AlgData>
  uint64_t verify_and_update_vertex_state(AlgData& alg_data, 
    BitSet& vertex_template_vertices) const {
    if (vertex_template_vertices.none()) {
      std::cerr << "Error: no bit is set." << std::endl;       
      return 0;       
    }

    auto g = std::get<10>(alg_data);
    
    //for (auto i = 0; i < vertex_template_vertices.size(); i++) { // Test
    //  if (vertex_template_vertices.test(i)) {
    //     //verify_and_update_vertex_state(alg_data, i); // works         
    //  }          
    //}

    verify_and_update_vertex_state_impl(alg_data, vertex_template_vertices);
      
    return 0;  
  }*/

  template<typename AlgData>
  uint64_t verify_and_update_vertex_state(AlgData& alg_data,
    BitSet& vertex_template_vertices) const {

//    typedef vertex_state_generic<VertexData, DynamicBitSet, uint64_t> VertexState;
    typedef vertex_state_generic<VertexData, uint8_t, BitSet> VertexState;     
    
    // std::get<6>(alg_data) - vertex_state_map 
    auto& pattern_graph = std::get<7>(alg_data); 
    auto g = std::get<10>(alg_data);   
    // std::get<11>(alg_data) - template_vertices 
 
    bool match_found = false; 
    
    // verify if vertex heard from a valid neighbor (parent)
    bool valid_parent_found = false;
    
//    DynamicBitSet parent_template_vertices_bitset(pattern_graph.vertex_count);

//    for (size_t i = 0; i < pattern_graph.vertex_count; i++) { // TODO: temporary fix
//      if (parent_template_vertices[i] == (uint8_t)1) {
//        parent_template_vertices_bitset.set(i);
//      }
      //std::cout << (uint32_t)parent_template_vertices[i] << " "; // Test
//    } // for     
    //std::cout << "\n"; // Test
    
        for (size_t vertex_pattern_index = 0; // TODO: change the name vertex_pattern_index to template_vertex
          vertex_pattern_index < vertex_template_vertices.size();
          vertex_pattern_index++) {
          if (vertex_template_vertices.test(vertex_pattern_index)) {
            match_found = true;  

	    //TODO: represent the pattern_graph using a bitset
	    // t = parent_template_vertices_bitset & pattern_graph.vertices[vertex_pattern_index].edges_bitset
	    // if (t.any()) valid_parent_found = true;  

            // verify if heard from a valid neighbor (parent)
            if (match_found) {
              // repeated metadata
              for (size_t i = 0; i < parent_template_vertices_bitset.size(); i++) {
                if (parent_template_vertices_bitset.test(i)) {

                  for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
                    e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {       
                    if (pattern_graph.edges[e] == i) {
                      valid_parent_found = true;   
                      //std::cout << "Valid parent found." << std::endl; // Test   
                      break;
                    }   
                  } // for

                  if (valid_parent_found) {
                    break;  
                  } 
                } // if  
              } // for       
              // repeated metadata 

//              if (!valid_parent_found) { // TODO: This is wrong! Moved it outside the for loop.
//                return 0; //false; 
//              } 

            } // if            
 
          } // if

          if (valid_parent_found) {
            //std::cout << "Valid parent found." << std::endl; // Test
            break;
          }
        } // for 
    
        if (match_found && !valid_parent_found) {
          //std::cout << "Ignorign invalid parent." << std::endl; // Test
          return 0; //false;
        }

    if (!match_found) {
      return 0;
    }

    //std::cout << "Valid parent found." << std::endl; // Test

    // vertex heard from a valid neighbor (parent) 
    // create an entry for this vertex in the vertex_state_map or
    // update it, if exists already 

    auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
    if (find_vertex == std::get<6>(alg_data).end()) {
      auto insert_status = std::get<6>(alg_data)
        .insert({g->locator_to_label(vertex), 
        VertexState()}); 
      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return 0;
      }     	
      find_vertex = insert_status.first;
      //find_vertex->second.vertex_pattern_index = vertex_pattern_index; // ID of the vertex in the pattern_graph 
      find_vertex->second.template_vertices = vertex_template_vertices;       
    }

    if (std::get<6>(alg_data).size() < 1) {
      return 0;
    }

    // set template_vertices    
    //for (size_t i = 0; i < find_vertex->second.template_vertices.size(); i++) { // Test
    //  std::cout << (uint32_t)find_vertex->second.template_vertices[i] << " ";
    //}
    //std::cout << "\n"; 

    // set template_neighbors 
    //for (size_t i = 0; i < pattern_graph.vertex_count; i++) {
    //  if (parent_template_vertices[i] == (uint8_t)1) {
    //    if (!find_vertex->second.template_neighbors.test(i)) {
    //      find_vertex->second.template_neighbors.set(i);    
    //    } 
    //  } 
    //}

/*    for (size_t i = 0; i < find_vertex->second.template_neighbors.size(); i++) {
      if (!find_vertex->second.template_neighbors.test(i)) { 
        if (parent_template_vertices_bitset.test(i)) {
          find_vertex->second.template_neighbors.set(i);
        }
      } 
    }*/
     
    find_vertex->second.template_neighbors|=parent_template_vertices_bitset; // bitwise or on the bitset
    //std::cout << find_vertex->second.template_neighbors << std::endl; // Test 

    //for (size_t i = 0; i < find_vertex->second.template_neighbors.size(); i++) { // Test
    //  std::cout << (uint32_t)find_vertex->second.template_neighbors[i] << " ";
    //}
    //std::cout << "\n";                   
    
    // TODO: populate template_neighbor_metadata_count_map. This should contain what is expected and/or the actual count? 
    // It is possible to validate it in the post traversal step. You could collect here. You could provide what is expected 
    // as an input / or precompute on a per rank basis.
    // If heard from a valid parent, even if it thinks it has multiple IDs, it still have only one label, just add the label to 
    // the map and increase the count by one. Do verification in the post processing step. E.g., I need thee gov and two net etc.
    // may use 1-byte label ID for each label in the template and map it back to the original label. 

    return 1;  
  }  

/*  template<typename AlgData>
  uint64_t verify_and_update_vertex_state(AlgData& alg_data, 
    size_t vertex_pattern_index) const {

    typedef vertex_state<uint8_t> VertexState; // TODO: use Vertex type

    //auto& pattern = std::get<1>(alg_data);  
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data); 
    auto g = std::get<10>(alg_data);
 
    bool match_found = false;

    // verify if parent_pattern_index is valid   
    for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
      e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
      if (pattern_graph.edges[e] == parent_pattern_index) {
        match_found = true;
        break; 
      }  
    }

    if (!match_found) {
      return 0;
    }

    // vertex heard from a valid neighbor 
    // create an entry for this vertex in the vertex_state_map or 
    // update, if exists already 
    auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
    if (find_vertex == std::get<6>(alg_data).end()) {
      auto insert_status = std::get<6>(alg_data).insert({g->locator_to_label(vertex), VertexState()});
      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return 0;
      }     	
      find_vertex = insert_status.first;
      find_vertex->second.vertex_pattern_index = vertex_pattern_index; // ID of the vertex in the pattern_graph 
    }

    if (std::get<6>(alg_data).size() < 1) {
      return 0;
    }

    // figure out what pattern indices are expected and add them to pattern_vertex_itr_count_map
    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      for (auto e = pattern_graph.vertices[vertex_pattern_index];
        e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
              
        auto pattern_index = pattern_graph.edges[e];   
        
          auto find_pattern_vertex =  find_vertex->second.pattern_vertex_itr_count_map.find(pattern_index);
          if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
            auto insert_status = find_vertex->second.pattern_vertex_itr_count_map.insert({pattern_index, 0});
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl;
              return 0;
            }
          } 

      } // for
      
    } // if

    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      return 0;
    }

    // set status of parent_pattern_index to 1
    auto find_pattern_vertex = find_vertex->second.pattern_vertex_itr_count_map.find(parent_pattern_index);  
    if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
      std::cerr << "Error: did not find the expected item in the map." << std::endl;
      return 0;
    }      
   
    // update status of the pattern vertex 
    if (find_pattern_vertex->second < 1) {
      find_pattern_vertex->second = 1;
    }  
 
    return 1; 
  }*/ 

  vertex_locator vertex;
  size_t parent_pattern_index; // TODO: pass type as template argument
  //std::array<uint8_t, max_bit_vector_size> parent_template_vertices; // TODO: remove this 
  //std::bitset<max_bit_vector_size> parent_template_vertices_bitset; 
  BitSet parent_template_vertices_bitset;
  uint8_t msg_type; // 0 - init, 1 - alive
};

template <typename TGraph, typename AlgData, typename VertexStateMap, 
  typename PatternGraph, typename VertexActive, typename VertexIteration, typename BitSet, typename TemplateVertex>
void verify_and_update_vertex_state(TGraph* g, AlgData& alg_data, 
  VertexStateMap& vertex_state_map, PatternGraph& pattern_graph, 
  VertexActive& vertex_active, 
  VertexIteration& vertex_iteration, uint64_t superstep, bool initstep, bool& global_not_finished, TemplateVertex& template_vertices) {

  typedef typename TGraph::vertex_iterator vertex_iterator;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_rank = havoqgt_env()->world_comm().rank();

  // Important : invalidate vertices that have valid labels but were not added 
  // to the vertex_state_map (becuase they did not hear from a valid parent)
  if (superstep == 0 && initstep) { // Important
    for (vertex_iterator vitr = g->vertices_begin(); 
      vitr != g->vertices_end(); ++vitr) {  
      vertex_locator vertex = *vitr; 
      if (vertex_active[vertex]) {
        auto find_vertex = vertex_state_map.find(g->locator_to_label(vertex));
        if (find_vertex == vertex_state_map.end()) { 
          vertex_active[vertex] = false;     
        } 
      }  
    }

    for(vertex_iterator vitr = g->delegate_vertices_begin();
      vitr != g->delegate_vertices_end(); ++vitr) {
      vertex_locator vertex = *vitr;
      if (vertex.is_delegate() && (g->master(vertex) == mpi_rank)) {
        auto find_vertex = vertex_state_map.find(g->locator_to_label(vertex));
        if (find_vertex == vertex_state_map.end()) {
          vertex_active[vertex] = false;
        }
      }
      // skip the delegates, reduction on vertex_active will take care of them 
    }
  }

  //auto vertex_temp = vertex_state_map.begin()->first;
  //std::vector<decltype(vertex_temp)> vertex_remove_from_map_list(0); // hack
  std::vector<uint64_t> vertex_remove_from_map_list; // TODO: use Vertex type

//  for (auto& v : vertex_state_map) { // TODO: use C++11 approach to remove item from map
//    auto v_locator = g->label_to_locator(v.first);
    
//    for (auto& p : v.second.pattern_vertex_itr_count_map) {
//      if (p.second < 1) {
//        vertex_remove_from_map_list.push_back(v.first);
//        vertex_active[v_locator] = false; 
//        break;
//      } else {
//        p.second = 0; // reset for next iteration   
//      }   
//    }   
//  } // for

  for (auto& v : vertex_state_map) {
    auto v_locator = g->label_to_locator(v.first);
    // v.second.template_vertices // bitset 
    // v.second.template_neighbors // bitset
    
    //TODO: Read from template_vertices[v_locator] instead of v.second.template_vertices.
    // At the moment they are synced. 
    
    // Test
    // reddit patterns_labels_A_3
    //if (v.second.template_vertices.test(0) && v.second.template_vertices.test(1)) {
    //  std::cout << v.first << " " << v.second.template_vertices << std::endl;  
    //}  
    // Test 
     
    for (size_t vertex_pattern_index = 0; // TODO: change the name vertex_pattern_index to template_vertex
      vertex_pattern_index < v.second.template_vertices.size();
      vertex_pattern_index++) {
      if (v.second.template_vertices.test(vertex_pattern_index)) {

        // if did not receive from all the required neighbours, set the vertex_pattern_index bit in v.second.template_vertices to 0
        // if all the bits in v.second.template_vertices are 0, v is added to the vertex_remove_from_map_list
         
        // there is a catch, now we need to maintain template_vertices on the delegates as well, 
        // and in addition to vertex_active.all_min_reduce(), we need to sync all the set bits on the delegates.
        // to_ullong(), all_max_reduce, vertex_data<bitset> but when you do reduction pass to_ullong()? or 
        BitSet tmp_pattern_vertex_edges;  // TODO: this should be done once per MPI rank  
        for (auto e = pattern_graph.vertices[vertex_pattern_index];
          e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
          assert(pattern_graph.edges[e] < 16); // Test
          tmp_pattern_vertex_edges.set(pattern_graph.edges[e]);  
        }

        assert(!tmp_pattern_vertex_edges.none()); // Important 

        //BitSet tb (std::string("1111001010111100")); // Test 
        //auto tmp_bitset = tmp_pattern_vertex_edges & tb;//v.second.template_neighbors;
        auto tmp_bitset = tmp_pattern_vertex_edges & v.second.template_neighbors;

        if (tmp_pattern_vertex_edges == tmp_bitset && !tmp_bitset.none()) {
          // heard from all required neighbours 
          //std::cout << tmp_pattern_vertex_edges << " same as " << tmp_bitset << " " << tb << std::endl;//<< v.second.template_neighbors << std::endl; // Test
          //std::cout << tmp_pattern_vertex_edges << " same as " << tmp_bitset << " " << v.second.template_neighbors << std::endl; // Test
        } else {
          // did not hear from all required neighbours
          v.second.template_vertices.reset(vertex_pattern_index);  
          //std::cout << tmp_pattern_vertex_edges << " not same as " << tmp_bitset << " " << tb << std::endl;//<< v.second.template_neighbors << std::endl; // Test
          //std::cout << tmp_pattern_vertex_edges << " not same as " << tmp_bitset << " " << v.second.template_neighbors << std::endl; // Test          
          //std::cout << v.second.template_vertices << " - " << v.first << std::endl;
        } 

      } // if
    } // for

    if (v.second.template_vertices.none()) {
      vertex_remove_from_map_list.push_back(v.first);
      vertex_active[v_locator] = false;
      //template_vertices[v_locator] = 0;  
      //std::cout << "removing " << v.second.template_vertices << " - " << v.first << std::endl;
    } else {
      template_vertices[v_locator] = v.second.template_vertices.to_ulong();
      v.second.template_neighbors.reset(); // reset for next iteration
      assert(v.second.template_neighbors.none());
    } 
 
  }

  if (vertex_remove_from_map_list.size() > 0) {
    global_not_finished = true;
  }

  for (auto v : vertex_remove_from_map_list) {
    if (vertex_state_map.erase(v) < 1) {
      std::cerr << "Error: failed to remove an element from the map." 
        << std::endl;  
    }    
  }

  // TODO: this is a temporary patch, forcing all the delegates to have no identity
  for(vertex_iterator vitr = g->delegate_vertices_begin();
      vitr != g->delegate_vertices_end(); ++vitr) {
      vertex_locator vertex = *vitr;
      if (vertex.is_delegate() && (g->master(vertex) == mpi_rank)) {
        continue;  // skip the controller 
      } 
      else { 
        auto find_vertex = vertex_state_map.find(g->locator_to_label(vertex));
        if (find_vertex == vertex_state_map.end()) {
          template_vertices[vertex] = 0;
        }
      }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  vertex_active.all_min_reduce();
  //template_vertices.all_bor_reduce();  
  template_vertices.all_max_reduce(); // ensure all the delegates have the value as the controller 
  MPI_Barrier(MPI_COMM_WORLD); 
}   

template <typename TGraph, typename VertexMetaData, typename VertexData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename VertexActive, 
  typename VertexIteration, typename VertexStateMapGeneric, typename PatternGraph, typename BitSet, typename TemplateVertex>
void label_propagation_pattern_matching_bsp(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, VertexRank& vertex_rank,
  VertexActive& vertex_active, VertexIteration& vertex_iteration, VertexStateMapGeneric& vertex_state_map_generic, 
  PatternGraph& pattern_graph, bool initstep, bool& global_not_finished, size_t global_itr_count, 
  std::ofstream& superstep_result_file, std::ofstream& active_vertices_count_result_file, TemplateVertex& template_vertices) {

  int mpi_rank = havoqgt_env()->world_comm().rank();
  uint64_t superstep_var = 0;
  uint64_t& superstep_ref = superstep_var;

  typedef uint64_t Vertex;
  //typedef std::bitset<max_bit_vector_size> BitSet;
  typedef boost::dynamic_bitset<> DynamicBitSet; 

  //typedef vertex_state_generic<VertexData, DynamicBitSet, uint64_t> VertexStateGeneric;
//  typedef vertex_state_generic<VertexData, uint64_t, BitSet> VertexStateGeneric;
//  typedef std::unordered_map<Vertex, VertexStateGeneric> VertexStateMapGeneric;
//  VertexStateMapGeneric vertex_state_map_generic;  

  typedef lppm_visitor<TGraph, Vertex, VertexData, BitSet, DynamicBitSet> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank,
    vertex_active, vertex_iteration, vertex_state_map_generic, pattern_graph, superstep_var, initstep, g, template_vertices);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // beiginning of BSP execution
  // TODO: change for loop to use a local termination detection at the end of a seperstep
  //bool not_finished = false;
 
  //for (uint64_t superstep = 0; superstep < 1; superstep++) {
  for (uint64_t superstep = 0; superstep < pattern_graph.diameter; superstep++) {
    superstep_ref = superstep;
    if (mpi_rank == 0) { 
      //std::cout << "Superstep #" << superstep << std::endl;
      std::cout << "Label Propagation | Superstep #" << superstep;
    }

    //MPI_Barrier(MPI_COMM_WORLD); 
    double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); 
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {    
      //std::cout << "Superstep #" << superstep <<  " Synchronizing ... " << std::endl;
      std::cout <<  " | Synchronizing ...";
    }

    //vertex_active.all_min_reduce(); // do not need this here
    ///MPI_Barrier(MPI_COMM_WORLD);
    //std::cout << "MPI Rank : " << mpi_rank <<  " - vertex_state_map_generic.size() : " << vertex_state_map_generic.size() << std::endl; // Test

    // TODO: TBA new verify_and_update_vertex_state_generic( ... )
 
    verify_and_update_vertex_state<TGraph, decltype(alg_data), VertexStateMapGeneric,
      PatternGraph, VertexActive, VertexIteration, BitSet, TemplateVertex>(g, alg_data, vertex_state_map_generic, pattern_graph, 
      vertex_active, vertex_iteration, superstep, initstep, global_not_finished, template_vertices);
    //MPI_Barrier(MPI_COMM_WORLD);
    //std::cout << "MPI Rank : " << mpi_rank <<  " - vertex_state_map_generic.size() : " << vertex_state_map_generic.size() << std::endl; // Test
   
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      //std::cout << "Superstep #" << superstep <<  " Time " << time_end - time_start << std::endl;
      std::cout << " | Time : " << time_end - time_start << std::endl;
    }

#ifdef output_result
    // result
    if (mpi_rank == 0) { 
      superstep_result_file << global_itr_count << ", LP, "
        << superstep << ", "
        << time_end - time_start << "\n"; 
    }

    // Important : This may slow things down -only for presenting results
    uint64_t active_vertices_count = 0; 
    for (auto& v : vertex_state_map_generic) {
      auto v_locator = g->label_to_locator(v.first);
      if (v_locator.is_delegate() && (g->master(v_locator) == mpi_rank)) {
        active_vertices_count++;  
      } else if (!v_locator.is_delegate()) {
        active_vertices_count++;
      }
    }
 
    active_vertices_count_result_file << global_itr_count << ", LP, "
      << superstep << ", " 
      << active_vertices_count << "\n";
#endif

    // TODO: global reduction on global_not_finished before next iteration

  } // for 
  // end of BSP execution  
}  

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_NONUNIQUE_HPP_INCLUDED 
