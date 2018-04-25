#ifndef HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_PATH_CHECKING_HPP_INCLUDED
#define HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_PATH_CHECKING_HPP_INCLUDED

#include <deque>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class tppm_queue {

public:
  tppm_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    //data.pop_back();
    data.pop_front();
  }

  Visitor const& top() {
    //return data.back();
    return data.front();
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
  //std::vector<Visitor> data;
  std::deque<Visitor> data;
};

// token passing pattern matching visitor class
template<typename Graph>
class tppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  tppm_visitor() : 
  itr_count(0), 
  do_pass_token(false), 
  is_init_step(true), 
  source_index_pattern_indices(0), 
  parent_pattern_index(0) {}

  tppm_visitor(vertex_locator _vertex) :  
    vertex(_vertex), 
    itr_count(0),
    do_pass_token(false), 
    is_init_step(true), 
    source_index_pattern_indices(0), 
    parent_pattern_index(0) {}
   
  tppm_visitor(vertex_locator _vertex, 
    vertex_locator _target_vertex, 
    vertex_locator _parent, 
    size_t _itr_count, 
    size_t _max_itr_count, 
    size_t _source_index_pattern_indices, 
    size_t _parent_pattern_index, 
    bool _expect_target_vertex = true, 
    bool _do_pass_token = true, 
    bool _is_init_step = false) : 
    vertex(_vertex),
    target_vertex(_target_vertex),
    parent(_parent), 
    itr_count(_itr_count), 
    max_itr_count(_max_itr_count), 
    expect_target_vertex(_expect_target_vertex), 
    do_pass_token(_do_pass_token), 
    is_init_step(_is_init_step), 
    source_index_pattern_indices(_source_index_pattern_indices), 
    parent_pattern_index(_parent_pattern_index) {}  

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(!std::get<13>(alg_data)[vertex]) { // vertex_active
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);
    //auto& pattern_graph = std::get<4>(alg_data);
    auto g = std::get<11>(alg_data); // graph
    // std::get<12>(alg_data); // vertex_token_source_set   
    // std::get<13>(alg_data); // vertex_active
    // std::get<14>(alg_data)[vertex]; // token_source_edge_set
    // std::get<15>(alg_data); // pattern_seleted_edges

    int mpi_rank = havoqgt_env()->world_comm().rank();
    
    // verify if this vertex has already forwarded a token from the originating vertex
    if (!is_init_step && max_itr_count > itr_count) {
      auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g->locator_to_label(target_vertex));
      if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
        return false;
      }
    }
   
    if (!do_pass_token && is_init_step && itr_count == 0) { // probably it never gets here
      // TODO: this is probbaly wrong
//      if (vertex.is_delegate() && (g->master(vertex) != mpi_rank)) {
//        return false;
//      }

//      if (!(find_vertex->second.vertex_pattern_index == pattern_indices[0] && vertex_data == pattern[0])) {
      if(vertex_data == pattern[0]) { 
        return true;
      } else {
        return false; 
      } 
    } else if (!is_init_step) { // forward the token
       auto new_itr_count = itr_count + 1;
       auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index       
       auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      
       if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
         if (vertex_data == pattern[next_pattern_index]) {
           vertex_pattern_index = pattern_indices[next_pattern_index];
         } else {
           return false;
         } 
       } else {
         auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
         if (find_vertex == std::get<5>(alg_data).end()) {
           return false;
         }
         vertex_pattern_index = find_vertex->second.vertex_pattern_index;
       }

       if (vertex_data == pattern[next_pattern_index] &&
         vertex_pattern_index == pattern_indices[next_pattern_index]) {
         if (vertex_data == pattern[next_pattern_index] && 
           parent_pattern_index == pattern_indices[next_pattern_index - 1]) {

           if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
             return true;
           }  

           if (max_itr_count > itr_count) {
             // the token is good to be forwarded, now update vertex_token_source_set
             auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g->locator_to_label(target_vertex));
             if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
               auto insert_status = std::get<12>(alg_data)[vertex].insert(g->locator_to_label(target_vertex));
               if(!insert_status.second) {
                 std::cerr << "Error: failed to add an element to the set." << std::endl;
                 return false;
               }
               // std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(target_vertex)
	       //  << " to the vertex set" << std::endl; // Test     
	     } else {
               std::cerr << "Error: unexpected item in the set." << std::endl;
               return false;
             }
           } // if 

           //TODO: max_itr_count == itr_count // currently visit handle this 
             
           return true; // actual forwarding happen in the visitfunction  
         } else {
           return false;
         }    
       } else {
         return false;
       }  
    } else {
      return false;  
    }
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(!std::get<13>(alg_data)[vertex]) { // vertex_active
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();

    // TODO: do you need ths here? 
    // if vertex is a delegate (not a controller)   
    if (!is_init_step && vertex.is_delegate() && (g.master(vertex) != mpi_rank)) { 
      // verify if this vertex have already forwarded a token from the originating vertex
      if (!is_init_step && max_itr_count > itr_count) {
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
          return false; 
        } 		
      }
    }
/*
    auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
    if (find_vertex == std::get<5>(alg_data).end()) {
      return false;
    } //else { // Test
      //std::cout << find_vertex->first << " " << find_vertex->second.vertex_pattern_index 
      //<< " " << std::get<0>(alg_data)[vertex] << std::endl;
    //} // Test
*/    
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);   
    //auto& pattern_graph = std::get<4>(alg_data);

    auto pattern_cycle_length = std::get<7>(alg_data);   
    auto pattern_valid_cycle = std::get<8>(alg_data);
    //auto& pattern_found = std::get<9>(alg_data);
    //auto& edge_metadata = std::get<10>(alg_data); 
    // <11> // graph
    // std::get<12>(alg_data); // vertex_token_source_set
    // std::get<13>(alg_data); // vertex_active 
    // std::get<14>(alg_data); // token_source_edge_set
    // std::get<15>(alg_data); // pattern_seleted_edges

    if (!do_pass_token && is_init_step && itr_count == 0) {
      // only the source vertices create visitors
      
      //if (!(find_vertex->second.vertex_pattern_index == pattern_indices[0] && vertex_data == pattern[0])) {
      if (vertex_data != pattern[0]) { 
        return false;  
      }
   
      //std::cout << "Found source vertex " << g.locator_to_label(vertex) << " vertex_data " << vertex_data << std::endl; // Test
      
      // Important: only the token_source_map on the controller contains the source vertex (wrong I think)
       
      // if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
      // local, controller and delegates are added to the token_source_map    
      auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
      if (find_token_source == std::get<6>(alg_data).end()) {
        auto insert_status = std::get<6>(alg_data).insert({g.locator_to_label(vertex), false});
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the map." << std::endl;
          return false;   
        } 
        //std::cout << "Instrting " << vertex_data << " to token_source_map" << std::endl; // Test    
      }	
      // }

      // initiate token passing from the source vertex
      if (!std::get<15>(alg_data)) { // pattern_seleted_edges      
        
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        if (!std::get<10>(alg_data)[eitr]) {
          continue;
        }
 
        //std::get<10>(alg_data)[eitr] = 240; // Test
        //int neighbour_edge_data = std::get<10>(alg_data)[eitr]; // Test

        //std::cout << "Found source vertex " << g.locator_to_label(vertex) 
        //<< " vertex_data " << vertex_data 
        //<< " neighbour " << g.locator_to_label(neighbour) 
        //<< " edge data " << neighbour_edge_data << std::endl; // Test   
  
        tppm_visitor new_visitor(neighbour, vertex, vertex, 0, pattern_cycle_length, 0, pattern_indices[0], pattern_valid_cycle, true, false);
        vis_queue->queue_visitor(new_visitor);
      } 

      } else { // pattern_seleted_edges
        // send tokens only on selected edges
        for (auto& e : std::get<14>(alg_data)[vertex]) {  
          vertex_locator neighbour = g.label_to_locator(e);
          tppm_visitor new_visitor(neighbour, vertex, vertex, 0, pattern_cycle_length, 0, pattern_indices[0], pattern_valid_cycle, true, false);
          vis_queue->queue_visitor(new_visitor);
          std::cout << "MPI Rank " << mpi_rank << " " << g.locator_to_label(vertex) << " sending to " << e << std::endl; // Test 
        }

      } //pattern_seleted_edges

      return true;
 
    } else if (!is_init_step) { // else if      

      bool do_forward_token = false;
      auto new_itr_count = itr_count + 1;
      auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index  
      auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      //auto vertex_pattern_index = pattern_indices[source_index_pattern_indices + new_itr_count]; // TODO: read from the map

      if (vertex.is_delegate() && g.master(vertex) != mpi_rank) { // delegate but not the controller 
        if (vertex_data == pattern[next_pattern_index]) {
          vertex_pattern_index = pattern_indices[next_pattern_index];  
        } else {
          return false;  
        }  
      } else {
        auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
        if (find_vertex == std::get<5>(alg_data).end()) {
          return false;
        }
        vertex_pattern_index = find_vertex->second.vertex_pattern_index;
      }      
      
      // TODO: verify next_pattern_index < pattern_indices.size() before anythin else

      //int mpi_rank = havoqgt_env()->world_comm().rank(); // Test
      //if (mpi_rank == 51) { 
      //  std::cout << g.locator_to_label(vertex) << " vertex_pattern_index "
      //  << vertex_pattern_index << " received from parent_pattern_index " 
      //  << parent_pattern_index << " itr_count " << itr_count 
      //  << " max_itr_count " << max_itr_count << " next_pattern_index " << next_pattern_index
      //  << std::endl; // Test	 
      //} // Test  

      if (max_itr_count > itr_count) {

        // are vertex_data and vertex_pattern_index valid
        if (vertex_data == pattern[next_pattern_index] &&
          vertex_pattern_index == pattern_indices[next_pattern_index]) {
          // verify if received from a valid parent
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 
            do_forward_token = true; 
          }   
        } // if 
      } else if (max_itr_count == itr_count) {
        // are vertex_data and vertex_pattern_index valid
        bool match_found = false;

        if (vertex_data == pattern[next_pattern_index] &&        
          vertex_pattern_index == pattern_indices[next_pattern_index]) { 
          // verify if received from a valid parent
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) {
            match_found = true; 
          }      
        } // if

        // is this the target vertex
        if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
          && match_found && expect_target_vertex) {
          // found cycle 
          //std::cout << "found valid cycle - vertex " << g.locator_to_label(vertex) <<  " | parent_pattern_index " 
          //<<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex) 
          //<<  " vertex_pattern_index " << vertex_pattern_index << " itr " 
          //<< itr_count << std::endl; // Test

          //return false; // TODO: true ?	

          // TODO: this is probbaly wrong
          // Important: only the token_source_map on the controller contains the source vertex (wrong I think)
//          if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
//            std::get<9>(alg_data) = 1; // true; // pattern_found
//            return true;		        
//          } 
           
          // if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
          auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
          if (find_token_source == std::get<6>(alg_data).end()) {
            std::cerr << "Error: did not find the expected item in the map." << std::endl;
            //return true; // false ?           
            return false;
          }
          find_token_source->second = 1; // true;   
          // }
 
          std::get<9>(alg_data) = 1; // true; // pattern_found

	  // add edge to the token_source_edge_set
          auto find_token_source_edge = std::get<14>(alg_data)[vertex].find(g.locator_to_label(parent));
          if (find_token_source_edge == std::get<14>(alg_data)[vertex].end()) {
          auto insert_status = std::get<14>(alg_data)[vertex].insert(g.locator_to_label(parent));
          if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the set." << std::endl;
            return false;
          }
//          std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(parent)
//            << " to the vertex set" << std::endl; // Test
          } else {
//            std::cerr << "Error: unexpected item in the set." << g.locator_to_label(vertex) << " adding " << g.locator_to_label(parent)
//            << " to the vertex set" << std::endl; // Test
//            std::cerr << "Error: unexpected item in the set." << std::endl;
            return false;
          }

          return true; // Important : must return true to handle delegates

        } else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
          && match_found && !expect_target_vertex) {
          // found path
          // TODO: TBA 
          //std::cout << "found invalid cycle - vertex " <<  " | parent_pattern_index "
          //<<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex)
          //<<  " vertex_pattern_index " << vertex_pattern_index << " itr "
          //<< itr_count << std::endl; // Test  
          return true; 
        } else {
          // reached max iteration but did not find the target vertex or a loop
          //std::cout << "At " << g.locator_to_label(vertex) <<  ", did not find target " 
          //<< g.locator_to_label(target_vertex) <<  " after " << itr_count 
          //<< " iterations" <<std::endl; // Test
          return false;  
        }   
      } else {
        std::cerr << "Error: wrong code branch." << std::endl;    
        return false;
      }   

      if (!do_forward_token) {
        return false;
      } //else {
        //return false; // Test 
      //} 

      // all good, forward along the token
      
      // if vertex is a delegate (not handled in the previsit)
      if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
        // forwarded a token from a source, now update vertex_token_source_set 
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
          auto insert_status = std::get<12>(alg_data)[vertex].insert(g.locator_to_label(target_vertex));		
	  if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the set." << std::endl;
            return false;
          }
          //std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(target_vertex) 
          //  << " to the vertex set" << std::endl; // Test 
        } else {
          std::cerr << "Error: unexpected item in the set." << std::endl;
	  return false;
        }        
      } // if vertex is a delegate

      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        if (!std::get<10>(alg_data)[eitr]) {
          continue;
        }

        // do not forward the token to the parent the vertex received it from
        if (g.locator_to_label(neighbor) == g.locator_to_label(parent)) {
          continue;  
        }        
    
        tppm_visitor new_visitor(neighbor, target_vertex, vertex, new_itr_count, max_itr_count, 
          source_index_pattern_indices, vertex_pattern_index, expect_target_vertex); 
        // vertex_pattern_index = parent_pattern_index for the neighbours 
        vis_queue->queue_visitor(new_visitor);   
      }
	     
      return true;

      // else if
    } else {
      return false;
    }
    return false;		 
  }

  friend inline bool operator>(const tppm_visitor& v1, const tppm_visitor& v2) {
    return false;
    /*if (v1.itr_count > v2.itr_count) {
      return true;
    } else if (v1.itr_count < v2.itr_count) {
      return false; 
    }
    if (v1.vertex == v2.vertex) {
      return false;
    }
    return !(v1.vertex < v2.vertex);*/
    
    /*if (v1.itr_count <= v2.itr_count) {
      return true;
    } else {
      return false;
    }*/ 
  }

  //friend inline bool operator<(const tppm_visitor& v1, const tppm_visitor& v2) {
    //return false;
    //if (v1.itr_count < v2.itr_count) {
    //  return true;
    //} else if (v1.itr_count > v2.itr_count) {
    //  return false;
    //}
    //if (v1.vertex == v2.vertex) {
    //  return false;
    //}
    //return !(v1.vertex < v2.vertex);    
  //}

  vertex_locator vertex;
  vertex_locator target_vertex; // for path checking, this is the source vertex
  vertex_locator parent;
  size_t itr_count; // TODO: change type
  size_t max_itr_count; // equal to diameter - 1 of the pattern as itr_count is initialized to 0 // TODO: change type
  size_t source_index_pattern_indices; // index of the token source in the pattern_indices container // TODO: use bitmap to support repeating labels   
  size_t parent_pattern_index; // TODO: change to the same type as in the pattern_graph
  bool expect_target_vertex; // false for path checking, true for cycles
  bool do_pass_token;
  bool is_init_step; // TODO: could be removed, this one or do_pass_token
};

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename PatternGraph, 
  typename VertexStateMap, typename TokenSourceMap, typename EdgeMetaData, 
  typename VertexSetCollection, typename VertexActive>
void token_passing_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, 
  VertexRank& vertex_rank, PatternGraph& pattern_graph, VertexStateMap& vertex_state_map,
  TokenSourceMap& token_source_map, size_t pattern_cycle_length, bool pattern_valid_cycle, bool pattern_seleted_edges,
  std::vector<uint8_t>::reference pattern_found, EdgeMetaData& edge_metadata, 
  VertexSetCollection& vertex_token_source_set, VertexSetCollection& token_source_edge_set, VertexActive& vertex_active) { // TODO: bool& pattern_found does not work, why?
  //std::cout << "token_passing_pattern_matching_new.hpp" << std::endl;
 
  typedef tppm_visitor<TGraph> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank, 
    pattern_graph, vertex_state_map, token_source_map, pattern_cycle_length, pattern_valid_cycle, pattern_found, 
    edge_metadata, g, vertex_token_source_set, vertex_active, token_source_edge_set, pattern_seleted_edges);
  auto vq = create_visitor_queue<visitor_type, /*havoqgt::detail::visitor_priority_queue*/tppm_queue>(g, alg_data);
  vq.init_visitor_traversal_new();
  //vq.init_visitor_traversal_new_alt();
  MPI_Barrier(MPI_COMM_WORLD);
}

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_PATH_CHECKING_HPP_INCLUDED
