#include <fstream>
#include <iostream>
#include <sstream>
#include <string>  
#include <vector>
#include <tuple>

#include <boost/algorithm/string.hpp>

#include <havoqgt/util.hpp> 

template <typename PatternType, typename Vertex = uint64_t, typename Edge = uint64_t> 
class pattern_util {
  public:
    //ypedef graph<Vertex, Edge, PatternType> pattern_graph;

    pattern_util(std::string pattern_input_filename) : 
      input_pattern_path_length(0),
      input_pattern(0),
      all_patterns(0),
      one_path_patterns(0), 
      two_path_patterns(0) {
    
      std::cout << "Reading the input pattern from file ..." << std::endl; 
      read_pattern_input(pattern_input_filename);
      
      input_pattern_path_length = input_pattern.size() - 1;
    
      std::cout << "Generating pattrens to search ..." << std::endl;
      generate_pattern_combinations(); 
      output_pattern_combinations(); 
    }

    pattern_util(std::string pattern_input_filename, bool is_integral_type) :
      input_pattern_path_length(0),
      input_pattern(0),
      all_patterns(0),
      one_path_patterns(0),
      two_path_patterns(0) {

      //std::cout << "Reading the input pattern list ..." << std::endl;
      read_pattern_list(pattern_input_filename, is_integral_type);
   
      //std::cout << "Pattrens to search ..." << std::endl;
      //for (auto ip : input_patterns) {
        //std::cout << "Data: " << std::endl; 
        //output_pattern(std::get<0>(ip));
        //std::cout << "Indices: " << std::endl;
        //output_pattern(std::get<1>(ip));               
        //std::cout << std::endl; 
      //}
    }

    /*pattern_util(std::string pattern_edge_filename, std::string pattern_vertex_filename, std::string pattern_vertex_data_filename, bool is_integral_type)
    {
      pattern_graph g(pattern_edge_filename, pattern_vertex_filename, pattern_vertex_data_filename, true, true);
      for (Vertex v = 0; v < g.vertex_count; v++) {
        std::cout << g.vertices[v] << " " << g.vertex_data[v] << " " << g.vertex_degree[v] << std::endl;
        for (auto e = g.vertices[v]; e < g.vertices[v + 1]; e++) {
          auto v_nbr = g.edges[e];   
          std::cout << v_nbr << ", " << std::endl;
        }
      } 
    }*/ 

    pattern_util(std::string pattern_input_filename, 
      std::string join_vertex_input_filename, bool is_integral_type) :
      input_pattern_path_length(0),
      input_pattern(0),
      all_patterns(0),
      one_path_patterns(0),
      two_path_patterns(0), 
      join_vertices(0) {

      //std::cout << "Reading the input pattern list ..." << std::endl;
      read_pattern_list(pattern_input_filename, is_integral_type);
      read_join_vertex_list(join_vertex_input_filename, is_integral_type);  
   
      //std::cout << "Pattrens to search ..." << std::endl;
      //for (auto ip : input_patterns) {
        //std::cout << "Data: " << std::endl; 
        //output_pattern(std::get<0>(ip));
        //std::cout << "Indices: " << std::endl;
        //output_pattern(std::get<1>(ip));               
        //std::cout << std::endl; 
      //}
    }   

    pattern_util(std::string pattern_input_filename, 
      std::string pattern_enumeration_input_filename, bool is_integral_type, bool is_enumeration) :
      input_pattern_path_length(0),
      input_pattern(0),
      all_patterns(0),
      one_path_patterns(0),
      two_path_patterns(0),
      enumeration_patterns(0),
      aggregation_steps(0),
      join_vertices(0) {

      //std::cout << "Reading the input pattern list ..." << std::endl;
      read_pattern_list(pattern_input_filename, is_integral_type);
      if (is_enumeration) { 
        //read_pattern_enumeration_list(pattern_enumeration_input_filename, is_integral_type);  
        read_pattern_enumeration_list_2(pattern_enumeration_input_filename, is_integral_type); 
      }
   
      //std::cout << "Pattrens to search ..." << std::endl;
      //for (auto ip : input_patterns) {
        //std::cout << "Data: " << std::endl; 
        //output_pattern(std::get<0>(ip));
        //std::cout << "Indices: " << std::endl;
        //output_pattern(std::get<1>(ip));               
        //std::cout << std::endl; 
      //}
    } 

    pattern_util(std::string pattern_input_filename, 
      std::string pattern_enumeration_input_filename, 
      std::string pattern_aggregation_input_filename,
      bool is_integral_type, bool is_enumeration) :
      input_pattern_path_length(0),
      input_pattern(0),
      all_patterns(0),
      one_path_patterns(0),
      two_path_patterns(0),
      enumeration_patterns(0),
      aggregation_steps(0),
      join_vertices(0) {

      //std::cout << "Reading the input pattern list ..." << std::endl;
      read_pattern_list(pattern_input_filename, is_integral_type);
      if (is_enumeration) { 
        //read_pattern_enumeration_list(pattern_enumeration_input_filename, is_integral_type);  
        //read_pattern_enumeration_list_2(pattern_enumeration_input_filename, is_integral_type);
        read_pattern_enumeration_list_3(pattern_enumeration_input_filename, is_integral_type);
        read_pattern_aggregation_list_3(pattern_aggregation_input_filename, is_integral_type); 
      }
   
      //std::cout << "Pattrens to search ..." << std::endl;
      //for (auto ip : input_patterns) {
        //std::cout << "Data: " << std::endl; 
        //output_pattern(std::get<0>(ip));
        //std::cout << "Indices: " << std::endl;
        //output_pattern(std::get<1>(ip));               
        //std::cout << std::endl; 
      //}
    } 

    ~pattern_util() {
    }
  
    static void output_pattern(std::vector<PatternType> pattern) {
      for (auto p : pattern) {
        std::cout << std::to_string(p) << ", ";
      }
      std::cout << std::endl;
    }

    size_t input_pattern_path_length;
    std::vector<PatternType> input_pattern; 
    std::vector<std::vector<PatternType>> all_patterns;
    std::vector<std::vector<PatternType>> one_path_patterns;
    std::vector<std::vector<PatternType>> two_path_patterns; 

    std::vector< std::tuple< std::vector<PatternType>, std::vector<Vertex>, Edge, bool, bool, bool> > input_patterns;     std::vector<std::vector<Vertex>> enumeration_patterns;
    std::vector<std::vector<uint8_t>> aggregation_steps;
      
    std::vector<Vertex> join_vertices;
 
  private:

    void read_pattern_list(std::string pattern_input_filename, bool is_integral_type) {
      std::ifstream pattern_input_file(pattern_input_filename,
        std::ifstream::in);
      std::string line;
      while (std::getline(pattern_input_file, line)) {
        std::istringstream iss(line);
        //std::cout << line << std::endl;

        auto tokens = split(line, ':');
        assert(tokens.size() > 1);

        boost::trim(tokens[0]); // important  
        boost::trim(tokens[1]); // important      
        boost::trim(tokens[2]); // important
        boost::trim(tokens[3]); // important
        boost::trim(tokens[4]); // important
        boost::trim(tokens[5]); // important    
 
        if (is_integral_type) {
          input_patterns.push_back(
            std::forward_as_tuple(
              split<PatternType>(tokens[0], ' '), split<Vertex>(tokens[1], ' '), 
              std::stoull(tokens[2]), std::stoull(tokens[3]), 
              std::stoull(tokens[4]), std::stoull(tokens[5])
            )
          );
        } else { 
          input_patterns.push_back( 
            std::forward_as_tuple( 
              split_char<PatternType>(tokens[0], ' '), split<Vertex>(tokens[1], ' '), 
              std::stoull(tokens[2]), std::stoull(tokens[3]),
              std::stoull(tokens[4]), std::stoull(tokens[5]) 
            ) 
          );   
        }

      }
      pattern_input_file.close();  
    }   

    void read_join_vertex_list(std::string join_vertex_input_filename, 
      bool is_integral_type) {
      std::ifstream join_vertex_input_file(join_vertex_input_filename,
        std::ifstream::in);
      std::string line;
      while (std::getline(join_vertex_input_file, line)) {
        boost::trim(line); // important     

        //auto tokens = split<Vertex>(line, ' ');
        //assert(tokens.size() > 0);

        if (is_integral_type) {
          join_vertices.push_back(std::stoull(line));     
        }
      }  
      join_vertex_input_file.close();
    }

    void read_pattern_enumeration_list(std::string pattern_enumeration_input_filename,
      bool is_integral_type) {
      std::ifstream pattern_enumeration_input_file(pattern_enumeration_input_filename,
        std::ifstream::in);
      std::string line;  
      while (std::getline(pattern_enumeration_input_file, line)) {
        //std::cout << line << std::endl;   
        boost::trim(line); // important 

	std::istringstream iss(line);

        auto tokens = split(line, ':');
        assert(tokens.size() > 1); 	   

        //boost::trim(tokens[0]); // important  
        boost::trim(tokens[1]); // important 

        if (is_integral_type) {
          enumeration_patterns.push_back(split<Vertex>(tokens[1], ' '));  
        }    
      }
      pattern_enumeration_input_file.close();
    } 

    void read_pattern_enumeration_list_2(std::string pattern_enumeration_input_filename,
      bool is_integral_type) {
      std::ifstream pattern_enumeration_input_file(pattern_enumeration_input_filename,
        std::ifstream::in);
      std::string line;  
      while (std::getline(pattern_enumeration_input_file, line)) {
        //std::cout << line << std::endl;   
        boost::trim(line); // important 

	std::istringstream iss(line);

        auto tokens = split(line, ':');
        assert(tokens.size() > 2); 	   

        //boost::trim(tokens[0]); // important  
        boost::trim(tokens[1]); // important
        boost::trim(tokens[2]); // important 

        if (is_integral_type) {
          enumeration_patterns.push_back(split<Vertex>(tokens[1], ' ')); 
          aggregation_steps.push_back(split<uint8_t>(tokens[2], ' ')); 
        }    
      }
      pattern_enumeration_input_file.close();
    } 

    void read_pattern_enumeration_list_3(std::string pattern_enumeration_input_filename,
      bool is_integral_type) {
      std::ifstream pattern_enumeration_input_file(pattern_enumeration_input_filename,
        std::ifstream::in);
      std::string line;  
      while (std::getline(pattern_enumeration_input_file, line)) {
        //std::cout << line << std::endl;   
        boost::trim(line); // important 

	std::istringstream iss(line);

        auto tokens = split(line, ':');
        assert(tokens.size() > 2); 	   

        //boost::trim(tokens[0]); // important  
        boost::trim(tokens[1]); // important
        boost::trim(tokens[2]); // important 

        if (is_integral_type) {
          enumeration_patterns.push_back(split<Vertex>(tokens[1], ' ')); 
          //aggregation_steps.push_back(split<uint8_t>(tokens[2], ' ')); 
        }    
      }
      pattern_enumeration_input_file.close();
    } 

    void read_pattern_aggregation_list_3(std::string pattern_aggregation_input_filename, 
      bool is_integral_type) {
      std::ifstream pattern_aggregation_input_file(pattern_aggregation_input_filename,
        std::ifstream::in);
      std::string line;  
      while (std::getline(pattern_aggregation_input_file, line)) {
        //std::cout << line << std::endl;   
        boost::trim(line); // important 

	std::istringstream iss(line);        

        if (is_integral_type) {
          aggregation_steps.push_back(split<uint8_t>(line, ' ')); 
        }    
      }
      pattern_aggregation_input_file.close();
    } 

    void read_pattern_input(std::string pattern_input_filename) {
      std::ifstream pattern_input_file(pattern_input_filename, 
        std::ifstream::in);
      std::string line;
      while (std::getline(pattern_input_file, line)) {
        std::istringstream iss(line);
        std::cout << line << std::endl;
        PatternType p(0);    
        iss >> p;
        input_pattern.push_back(p); 
      }
      pattern_input_file.close(); 
    }     

    void generate_pattern_combinations() {
      size_t max_path_length = input_pattern.size() - 1;
      size_t path_length = 1;
      while (path_length <= max_path_length) {
        for (size_t i = 0, j = i + path_length; j < input_pattern.size(); i++, j++) {
          std::vector<PatternType> pattern(input_pattern.begin() + i, input_pattern.begin() + j + 1);
          all_patterns.push_back(pattern);
          if (path_length == 1) {
            one_path_patterns.push_back(pattern);
          } else if (path_length == 2) {
            two_path_patterns.push_back(pattern);
          } 
        }
        path_length++; 
      }  
    }
  
    void output_pattern_combinations() {
      for (auto pattern : all_patterns) {
        for (auto p : pattern) {
          std::cout << p << ", ";
        }
        std::cout << std::endl;
      } 
    }     
};
