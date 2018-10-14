#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <havoqgt/approximate_pattern_matching/pattern_non_local_constraint.hpp>

template<typename Vertex, typename Edge, typename VertexData, 
  typename EdgeData, typename Boolean, typename BitSet, 
  typename Uint = uint64_t>
class pattern_graph_csr {
  public:
    /*pattern_graph_csr(std::string vertex_input_filename,  
      std::string edge_input_filename, 
      std::string vertex_data_input_filename, bool _directed) :
      //const bool _directed = true) :
      directed(_directed),
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0) {
      std::cout << "Building CSR graph ... " << std::endl;
      
      std::cout << "Reading vertex list ... " << std::endl;
      vertex_count = read_vertex_list(vertex_input_filename, vertices,
        vertex_degree);

      std::cout << "Reading vertex data list ... " << std::endl;
      read_vertex_data_list(vertex_data_input_filename, vertex_data);

      std::cout << "Reading edge list ... " << std::endl;
      edge_count = read_edge_list(edge_input_filename, edges);

      std::cout << "Completed building graph." << std::endl;
      output_graph_stat(); 
    }*/

    /*pattern_graph_csr(std::string edge_input_filename, 
      std::string vertex_input_filename,
      std::string vertex_data_input_filename, 
      std::string stat_input_filename,
      const bool _directed = true, const bool _mutable = false) :
      directed(_directed),
      vertex_count(0),
      edge_count(0),
      diameter(0),  
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0) {
//      std::cout << "Building CSR graph ... " << std::endl;
      
//      std::cout << "Reading edge list ... " << std::endl;
      edge_count = read_edge_list(edge_input_filename);

      //std::cout << "Generating vertex list ..." << std::endl;
      vertex_count = generate_vertex_list();       

//      std::cout << "Reading vertex list ... " << std::endl;
//      vertex_count = read_vertex_list(vertex_input_filename);

//      std::cout << "Reading vertex data list ... " << std::endl;
      read_vertex_data_list(vertex_data_input_filename);

      read_stat(stat_input_filename);

//      std::cout << "Completed building graph." << std::endl;
//      output_graph_stat(); 
    }*/

    /*pattern_graph_csr(std::string edge_input_filename, 
      std::string vertex_input_filename,
      std::string vertex_data_input_filename, 
      std::string edge_data_input_filename, 
      std::string stat_input_filename, //std::string tmp_input_filename,
      const bool _directed = false, const bool _mutable = false) :
      directed(_directed),
      vertex_count(0),
      edge_count(0),
      diameter(0),  
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0),
      edge_is_optional(0) { //, 
      //vertices_optional(0),
      //vertex_degree_optional(0),
      //vertex_data_optional(0),
      //edges_optional(0),
      //edge_list_optional(0) {
//      std::cout << "Building CSR graph ... " << std::endl;
      
//      std::cout << "Reading edge list ... " << std::endl;
      edge_count = read_edge_list(edge_input_filename);

      //std::cout << "Generating vertex list ..." << std::endl;
      vertex_count = generate_vertex_list();       
      //vertex_count = generate_vertex_list(edge_list, vertices, edges, vertex_degree);

//      std::cout << "Reading vertex list ... " << std::endl;
//      vertex_count = read_vertex_list(vertex_input_filename);

//      std::cout << "Reading vertex data list ... " << std::endl;
      read_vertex_data_list(vertex_data_input_filename);

//      std::cout << "Reading edge data list ... " << std::endl;
      //read_edge_data_list(edge_data_input_filename); 

      read_stat(stat_input_filename);

      generate_vertex_neighbor_data_count_map();

      generate_vertex_edges_bitset(); 

//      std::cout << "Comnleted building graph." << std::endl;
//      output_graph_stat(); 
    }*/

    pattern_graph_csr(std::string edge_input_filename, 
      std::string vertex_data_input_filename,        
      std::string stat_input_filename, 
      std::string edge_constraints_input_filename) :
      vertex_count(0),
      edge_count(0),
      diameter(0),  
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0),
      edge_is_optional(0), 
      vertex_edges_bitset(0), 
      vertex_edges_bitset_optional(0), 
      vertex_edges_bitset_all(0),
      vertex_min_optional_edge_count(0) {  

      edge_count = read_edge_list(edge_input_filename);

      vertex_count = generate_vertex_list();       

      read_vertex_data_list(vertex_data_input_filename);

      read_stat(stat_input_filename);

      generate_vertex_neighbor_data_count_map();

      generate_vertex_edges_bitset();

      read_vertex_local_constraints(edge_constraints_input_filename);
      //read_edge_constraints(edge_constraints_input_filename); 
    }

    /*pattern_graph_csr(size_t mpi_rank, std::string edge_input_filename, std::string vertex_input_filename,
      std::string vertex_data_input_filename,
      std::string edge_data_input_filename,
      std::string stat_input_filename,
      const bool _directed = false, const bool _mutable = false) :
      directed(_directed),
      vertex_count(0),
      edge_count(0),
      diameter(0),
      vertices(0),
      vertex_degree(0),
      vertex_data(0),
      edges(0),
      edge_list(0) {

      edge_count = read_edge_list(edge_input_filename);

      if (mpi_rank == 0) {
        output_graph_stat();
        for (auto e : edge_list) {
          std::cout << std::get<0>(e) << " " << 
            std::get<1>(e) << std::endl; 
        } 
        std::cout << "Max vertex: " << std::get<0>(edge_list[edge_list.size() - 1]) << std::endl;
      }

      //vertex_count = generate_vertex_list_2(mpi_rank);
      vertices.push_back(0);
      vertices.push_back(2);
      vertices.push_back(5);
      vertices.push_back(8);
      vertices.push_back(10);
      vertices.push_back(12);
      vertices.push_back(13);
      vertices.push_back(14);
      vertices.push_back(15);
      vertices.push_back(16);

      vertex_degree.push_back(2);
      vertex_degree.push_back(3); 
      vertex_degree.push_back(3); 
      vertex_degree.push_back(2);  
      vertex_degree.push_back(2);
      vertex_degree.push_back(1);
      vertex_degree.push_back(1);
      vertex_degree.push_back(1);
      vertex_degree.push_back(1);
      
      vertex_count = vertices.size() - 1;
 
      if (mpi_rank == 0) {
        output_graph_stat();
        std::cout << "Max vertex: " << std::get<0>(edge_list[edge_list.size() - 1]) << std::endl; 
      }

      read_vertex_data_list(vertex_data_input_filename);
      read_stat(stat_input_filename);
      
    }*/

    ~pattern_graph_csr() {
//      std::cout << "Disposing graph ... " << std::endl;
    } 

    //const bool directed; 
    bool directed;
    Vertex vertex_count;    
    Edge edge_count;
    Vertex max_vertex; 
    Edge diameter;  
    std::vector<Edge> vertices;
    std::vector<Edge> vertices_optional; 
    std::vector<Edge> vertex_degree;
    std::vector<Edge> vertex_degree_optional;
    std::vector<VertexData> vertex_data;
    std::vector<Vertex> edges;
    std::vector<Vertex> edges_optional;
    std::vector<Edge> edge_ID; 
    std::vector<EdgeData> edge_data;
    std::vector<bool> edge_is_optional; 
    std::vector<std::tuple<Vertex, Vertex, bool>> edge_list;
    std::vector<std::tuple<Vertex, Vertex>> edge_list_optional; 
    std::vector<std::unordered_map<VertexData, Uint>> vertex_neighbor_data_count_map; 
    std::vector<BitSet> vertex_edges_bitset;
    std::vector<BitSet> vertex_edges_bitset_optional;
    std::vector<BitSet> vertex_edges_bitset_all; 
    std::vector<int64_t> vertex_min_optional_edge_count;  

  private:
    Vertex read_vertex_list(std::string vertex_input_filename) { 
      std::ifstream vertex_input_file(vertex_input_filename, std::ifstream::in);
      std::string line;
      while (std::getline(vertex_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        Edge v_degree(0), v_offset(0);
        iss >> v_source >> v_degree >> v_offset;
        vertices.push_back(v_offset);
        vertex_degree.push_back(v_degree);
      }
      vertex_input_file.close();
      vertex_degree.erase(vertex_degree.end() - 1);
      return vertices.size() <= 0 ? 0 : vertices.size() - 1; // vertex count
    }

    void read_vertex_data_list(std::string vertex_data_input_filename) {
      std::ifstream vertex_data_input_file(vertex_data_input_filename, 
        std::ifstream::in);
      std::string line;
      while (std::getline(vertex_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        VertexData v_data(0);
        iss >> v_source >> v_data;
        vertex_data.push_back(v_data);
      }
      vertex_data_input_file.close();
    }

    void read_vertex_local_constraints(std::string edge_constraints_input_filename) {
      std::ifstream edge_constraints_input_file
        (edge_constraints_input_filename, std::ifstream::in);
      std::string line;
      while (std::getline(edge_constraints_input_file, line)) {
        std::istringstream iss(line);
        auto tokens = split(line, ':');
        assert(tokens.size() > 1);
        boost::trim(tokens[0]); // important  
        boost::trim(tokens[1]); // important         

        Vertex v_source(0);
        int64_t min_optional_edge_count(0);
        //iss >> v_source >> min_optional_edge_count; // -1 - vertex has no optional edge        
        
        try {
          v_source = boost::lexical_cast<Vertex>(tokens[0]); 
        } catch (boost::bad_lexical_cast& e) {
          std::cerr << "Exception caught : " << e.what() << std::endl;
          return; 
        }

        try {
          min_optional_edge_count = boost::lexical_cast<int64_t>(tokens[1]);
        } catch (boost::bad_lexical_cast& e) {
          std::cerr << "Exception caught : " << e.what() << std::endl;
          return;
        } 
        
        //std::cout << v_source << " " << min_optional_edge_count << std::endl; // Test        
        vertex_min_optional_edge_count.push_back(min_optional_edge_count);
      }
      edge_constraints_input_file.close();
    }  

    Edge read_edge_list(std::string edge_input_filename) {
      std::ifstream edge_input_file(edge_input_filename, std::ifstream::in);
      std::string line;
      while(std::getline(edge_input_file, line)) {
        std::istringstream iss(line);
        Vertex s(0), t(0); 
        bool is_optional(false);
        iss >> s >> t >> is_optional; // 0 - optional, 1 - mandatory        
        //if (!is_optional) {
          edge_list.push_back(std::forward_as_tuple(s, t, is_optional));
          edges.push_back(t);
          edge_is_optional.push_back(is_optional); 
        //} else {
        //  edges_optional.push_back(t);
        //  edge_list_optional.push_back(std::forward_as_tuple(s, t));
        //} 
      }
      edge_input_file.close();
      return edges.size(); // edge count
      //return edge_list.size();
    }

    void read_edge_constraints(std::string edge_constraints_input_filename) {
      std::ifstream edge_constraints_input_file
        (edge_constraints_input_filename, std::ifstream::in);
      std::string line;
      while (std::getline(edge_constraints_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        int64_t min_optional_edge_count(0);
        iss >> v_source >> min_optional_edge_count; // -1 - vertex has no optional edge
        //std::cout << v_source << " " << min_optional_edge_count << std::endl; // Test        
        vertex_min_optional_edge_count.push_back(min_optional_edge_count);
      }
      edge_constraints_input_file.close();
    }   

    void read_edge_data_list(std::string edge_data_input_filename) {
      std::ifstream edge_data_input_file(edge_data_input_filename, std::ifstream::in);
      std::string line;
      while(std::getline(edge_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex s(0), t(0);
        Edge e(0);
        EdgeData w(0); 
        iss >> s >> t >> e >> w;
        edge_ID.push_back(e); // TODO: edge IDs should be in a different file
        edge_data.push_back(w); 
      }
      edge_data_input_file.close();   
    }

    Vertex generate_vertex_list() { // assuming graph is undirected
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]); // TODO: update read_edge_list to get max_vertex
      Vertex l = 0; // edge list index
      Vertex degree = 0;
      Vertex current_vertex = vertex_count;
      Vertex source;
      //Vertex target;       
    
      do {
        auto edge = edge_list[l];
        source = std::get<0>(edge);
        //target = std::get<1>(edge);
        if (source == current_vertex) {
          degree++;
          l++;
        } else {
          //std::cout << current_vertex << std::endl;
          vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
          vertex_degree.push_back(degree);

          //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
          //vertex_data.push_back(v_data);

          // write vertex info to file
          //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";
   
          //vertex_data_file << current_vertex << " " << v_data << "\n";
      
          // update vertices array
          degree = 0;
          vertex_count++;
          current_vertex = vertex_count;
        }
      } while(current_vertex <= max_vertex); 

      // add the last dummy vertex to the vertex list
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();

      return vertex_count;  
    }

    // not used
    template<typename EdgeListVector, typename VertexVector, 
      typename EdgeVector, typename VertexDegreeVector>
    Vertex generate_vertex_list(EdgeListVector& edge_list, 
      VertexVector& vertices, EdgeVector& edges,
      VertexDegreeVector& vertex_degree, 
      const bool filter_optional_edges = false) { // assuming graph is undirected
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
      Vertex l = 0; // edge list index
      Vertex degree = 0;
      Vertex current_vertex = vertex_count;
      Vertex source;
      //Vertex target;       
    
      do {
        auto edge = edge_list[l];
        source = std::get<0>(edge);
        //target = std::get<1>(edge);
        auto is_optional = std::get<2>(edge);        
 
        if (source == current_vertex) {
          if (filter_optional_edges) {
            if (is_optional) { 
              degree++;
              l++;
              edges.push_back(edge);
            }
          } else {
            if (!is_optional) {
              degree++;
              l++;
              edges.push_back(edge); 
            }
          }   
        } else {
          //std::cout << current_vertex << std::endl;
          vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
          vertex_degree.push_back(degree);

          //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
          //vertex_data.push_back(v_data);

          // write vertex info to file
          //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";
   
          //vertex_data_file << current_vertex << " " << v_data << "\n";
      
          // update vertices array
          degree = 0;
          vertex_count++;
          current_vertex = vertex_count;
        }
      } while(current_vertex <= max_vertex); 

      // add the last dummy vertex to the vertex list
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();

      return vertex_count;  
    }

    Vertex generate_vertex_list_2(size_t mpi_rank) { // assuming graph is undirected
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
      Vertex l = 0; // edge list index
      Vertex degree = 0;
      Vertex current_vertex = vertex_count;
      Vertex source;
      //Vertex target;       
    
      /*do {
        auto edge = edge_list[l];
        source = std::get<0>(edge);
        //target = std::get<1>(edge);
        if (source == current_vertex) {
          degree++;
          l++;
        } else {
          //std::cout << current_vertex << std::endl;
//          vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
//                                                     vertices[vertices.size() - 1]));
//          vertex_degree.push_back(degree);

          //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
          //vertex_data.push_back(v_data);

          // write vertex info to file
          //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";
          if (mpi_rank == 0) {
            std::cout << current_vertex << " " << degree << " " << " " << "\n"; 
          }

          //vertex_data_file << current_vertex << " " << v_data << "\n";
      
          // update vertices array
          degree = 0;
          vertex_count++;
          current_vertex = vertex_count;
        }
      } while(current_vertex <= max_vertex); 

      // add the last dummy vertex to the vertex list
//      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
//                                                 vertices[vertices.size() - 1]));
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();*/
     
      vertices.resize(10);  
      vertices.push_back(0);
      vertices.push_back(2);
      vertices.push_back(8);
      vertices.push_back(10);
      vertices.push_back(12);
      vertices.push_back(13); 		 
      vertices.push_back(14);
      vertices.push_back(15);
      vertices.push_back(16);
 
      vertex_count = vertices.size() - 1; 

      return vertex_count;  
    }

    void read_stat(std::string stat_input_filename) {
      std::ifstream stat_input_file(stat_input_filename, 
        std::ifstream::in);
      std::string line;
      const char delim = ':';
      while(std::getline(stat_input_file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while(std::getline(iss, token, delim)) {
          tokens.push_back(token);
        }          
        //std::cout << tokens[0] << " " << tokens[1] << std::endl; 
        assert(tokens.size() > 1);
        boost::trim(tokens[0]);
        boost::trim(tokens[1]); 
        if (boost::iequals(tokens[0], "diameter")) {
          diameter = std::stoull(tokens[1]);         
        } 
      }
      stat_input_file.close();
    } 

    void generate_vertex_neighbor_data_count_map() {
      vertex_neighbor_data_count_map.resize(vertex_count); // vector of maps
      
      // make sure the maps are empty 
      for (auto v = 0; v < vertex_count; v++) {
        vertex_neighbor_data_count_map[v].empty();
      }

      for (auto v = 0; v < vertex_count; v++) {
        //std::cout << v << " vertex_data " << vertex_data[v]  
        //  << " vertex_degree " << vertex_degree[v] << std::endl;
        //  std::cout << " neighbours : ";
        for (auto e = vertices[v]; e < vertices[v + 1]; e++) {
          auto v_nbr = edges[e];    
          auto v_nbr_vertex_data = vertex_data[v_nbr];
          //std::cout << v_nbr << ", " << v_nbr_vertex_data << " | ";  

          auto find_nbr_vertex_data = vertex_neighbor_data_count_map[v].find(v_nbr_vertex_data);
          if (find_nbr_vertex_data == vertex_neighbor_data_count_map[v].end()) {
            auto insert_status = vertex_neighbor_data_count_map[v].insert({v_nbr_vertex_data, 1});    
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl; 
              return;
            }
          } else { 
            find_nbr_vertex_data->second++;
          }            
  
        }
        //std::cout << std::endl;
      }  
    }

    void generate_vertex_edges_bitset() {
      vertex_edges_bitset.resize(vertex_count); // vector of bitsets
      vertex_edges_bitset_optional.resize(vertex_count);
      vertex_edges_bitset_all.resize(vertex_count); 
      for (auto v = 0; v < vertex_count; v++) {
        for (auto e = vertices[v]; e < vertices[v + 1]; e++) {
          assert(edges[e] < vertex_edges_bitset[v].size());
          assert(edges[e] < vertex_edges_bitset_optional[v].size());
          auto v_nbr_is_mandatory = edge_is_optional[e]; // 0 - optional, 1 - mendatory
	  if (!v_nbr_is_mandatory) { // 0 - optional
            vertex_edges_bitset_optional[v].set(edges[e]);
          }
          else if (v_nbr_is_mandatory) { // 1 - mendatory
            vertex_edges_bitset[v].set(edges[e]);
          }
        } // for
        vertex_edges_bitset_all[v] = vertex_edges_bitset[v] | vertex_edges_bitset_optional[v];
        assert(!vertex_edges_bitset_all[v].none()); 
      } // for
    } 

    void output_graph_stat() {
      std::cout << "Number of vertices: " << vertex_count << std::endl;
      std::cout << "Number of edges: " << edge_count << std::endl;
      std::cout << "Number of vertex data: " << vertex_data.size() << std::endl;
    } 
};
