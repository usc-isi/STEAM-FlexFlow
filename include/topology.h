/* Copyright 2021 Massachusetts Institue of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _FLEXFLOW_TOPOLOGY_H_
#define _FLEXFLOW_TOPOLOGY_H_

#include <vector>
#include <deque>
#include <map> 
#include <unordered_map>
#include <unordered_set>
#include <random>

#include "ffconst.h"
#include "config.h"

/**
 * class NetworkTopology
 * Defines the network topology by a connection matrix 
 */
class NetworkTopology {

public: 
  /**
   * construct a network topology without a provided starting one.
   * This will fill the matrix as a degree-constrained spanning tree,
   * and then allocate the rest of the interfaces randomly
   */
  NetworkTopology(int num_nodes, uint64_t bw, int if_cnt,  
                  std::unordered_map<size_t, Device*>& device_map);

  /**
   * construct a network topology as specified by topofile
   */
  NetworkTopology(int num_nodes, int if_cnt, uint64_t bw,
              const std::string& topofile, 
              std::unordered_map<size_t, Device*>& device_map);

  /* general getter functions */
  uint64_t get_bandwidth_bps(int i, int j);
  uint64_t get_bandwidth_Bps(int i, int j); 
  int get_if_in_use(int node);
  int get_port_in_use();

  /**
   * returns the shortest path route between node src and dst. Ties are broken
   * based on whatever C++ decides
   */
  std::deque<int>& get_shortest_route(int src, int dst);

  /* print the connection matrix in matrix form */
  void print_conn_matrix();

  /**
   * Return all the shortest path route in the system 
   */
  std::vector<std::pair<std::pair<unsigned, unsigned>, 
    std::deque<int> > > get_allroutes();

  /**
   * print all routes in to fname
   */
  void print_routes(const std::string & fname,
                  std::vector<std::pair<std::pair<unsigned, unsigned>, 
                  std::deque<int> > > & result);

  /* helpers to check if the current topology remains connected */
  bool check_connected();
  /* This one takes a stored mulatigraph topology */
  bool check_connected(const std::unordered_map<uint64_t, size_t> & topo);

  /* helpers to import and export topologies */
  std::unordered_map<uint64_t, size_t> get_topology();
  /* The first one sets a multi graph (full topology), the second one
   * sets a spanning tree based topology (t is the edge set)
   */
  void set_topology(std::unordered_map<uint64_t, size_t> & t);
  void set_topology(std::unordered_set<uint64_t> & t);

  /* helpers to initialize the topology as a ST or check if its a ST */
  void fill_spanning_tree(bool allocate_rest);
  bool check_spt(const std::unordered_set<uint64_t> & edgeset);

  // TODO: port these
  // bool can_add();
  // int add_link(int i, int j, bool check = true);
  // int remove_link(int i, int j, bool check = true);
  // void update(std::vector<std::pair<uint64_t, uint64_t> > & demand,
  //            std::vector<std::pair<uint64_t, uint64_t> > & cp_demand);
  // void dheuristic_plus(std::vector<std::pair<uint64_t, uint64_t> > & logical_demand, 
  //                             std::vector<std::pair<uint64_t, uint64_t> > & physical_demand,
  //                             bool & forced);
  // void connect_cc(std::unordered_map<uint64_t, uint64_t> & logical_id_to_demand);
  // void load_route(const TaskGraphProtoBuf::Topology& topopb);

  /* remove all cached routes */
  inline void clear_routes() {
    for (auto& item: routes) {
      if (item.second != nullptr) {
        delete item.second;
      }
    }
    routes.clear();
  }

  std::unordered_map<size_t, Device*>& ids_to_inter_node_comm_device;
  int num_nodes;      // number of nodes in the system
  int if_cnt;         // per device
  uint64_t bandwidth; // per link

  std::vector<std::vector<int> > conn_matrix;

#ifdef FF_USE_ECMP
  std::unordered_multimap<uint64_t, std::deque<int>*> routes;
#else 
  std::unordered_map<uint64_t, std::deque<int>*> routes;
#endif

  std::random_device rd;
  std::mt19937 gen;

};

#endif