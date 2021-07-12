#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <utility>
#include <unordered_set>

#include "simulator.h"
// #define EDGE(a, b, n) ((a) > (b) ? ((a) * (n) + (b)) : ((b) * (n) + (a)))
#define PRINT_EDGE(e, n) do {std::cout << "(" << e / n << ", " << e % n << ")";} while (0);

static std::random_device rd; 
static std::mt19937 gen = std::mt19937(rd()); 

ShortestPathNetworkRoutingStrategy::ShortestPathNetworkRoutingStrategy(
    const ConnectionMatrix & c, 
    const std::map<size_t, CommDevice*>& devmap,
    int total_devs) 
: conn(c), devmap(devmap), total_devs(total_devs)
{} 

EcmpRoutes ShortestPathNetworkRoutingStrategy::get_routes(int src_node, int dst_node) 
{
  int key = src_node * total_devs + dst_node;

  if (conn[key] > 0) {
    return std::make_pair(std::vector<float>({1}), std::vector<Route>({Route({devmap.at(key)})}));
  }

  // one-shortest path routing
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::priority_queue<std::pair<uint64_t, uint64_t>, 
                      std::vector<std::pair<uint64_t, uint64_t> >,
                      std::greater<std::pair<uint64_t, uint64_t> > > pq;
  pq.push(std::make_pair(dist[src_node], src_node));
  dist[src_node] = 0;

  /*
   * dijkstra implementation. Really BFS would work but this is easier for
   * future copy-pasting... 
   */
  while (!pq.empty()) {
    int min_node = pq.top().second;
    pq.pop();
    visited[min_node] = true;

    if (min_node == dst_node)
      break;

    for (int i = 0; i < total_devs; i++) {
      if (visited[i] || conn[min_node * total_devs + i] == 0) {
        continue;
      }
      double new_dist = dist[min_node] + 1; // numeric_limits<uint64_t>::max() / get_bandwidth_bps(min_node, i);
      if (new_dist < dist[i]) {
        dist[i] = new_dist;
        prev[i] = min_node;
        pq.push(std::make_pair(new_dist, i));
      }
    }
  }

  Route result = Route();
  int curr = dst_node;
  while (prev[curr] != -1) {
    result.insert(result.begin(), devmap.at(prev[curr] * total_devs + curr));
    curr = prev[curr];
  }

  return std::make_pair(std::vector<float>({1}), std::vector<Route>({result}));

}

FlatDegConstraintNetworkTopologyGenerator::FlatDegConstraintNetworkTopologyGenerator(int num_nodes, int degree) 
: num_nodes(num_nodes), degree(degree)
{}

ConnectionMatrix FlatDegConstraintNetworkTopologyGenerator::generate_topology() const
{
  ConnectionMatrix conn = std::vector<int>(num_nodes*num_nodes, 0);
  
  int allocated = 0;
  int curr_node = 0;
  std::unordered_set<int> visited_node;
  visited_node.insert(0);

  std::uniform_int_distribution<> distrib(0, num_nodes - 1);

  while ((long)visited_node.size() != num_nodes) {
    distrib(gen);
    int next_step = distrib(gen);
    if (next_step == curr_node) {
      continue;
    } 
    if (visited_node.find(next_step) == visited_node.end()) {
      if (conn[get_id(curr_node, next_step)] == degree) {
        continue;
      }
      conn[get_id(curr_node, next_step)]++;
      conn[get_id(next_step, curr_node)]++;
      visited_node.insert(next_step);
      curr_node = next_step;
      allocated += 2;
    }
  }

  assert(allocated == (num_nodes - 1) * 2);

  std::vector<std::pair<int, int> > node_with_avail_if;
  for (int i = 0; i < num_nodes; i++) {
    int if_inuse = get_if_in_use(i, conn);
    if (if_inuse < degree) {
      node_with_avail_if.emplace_back(i, degree - if_inuse);
    }
  }

  distrib = std::uniform_int_distribution<>(0, node_with_avail_if.size() - 1);
  int a = 0, b = 0;

  while (node_with_avail_if.size() > 1) {
    a = distrib(gen);
    while ((b = distrib(gen)) == a);

    assert(conn[get_id(node_with_avail_if[a].first, node_with_avail_if[b].first)] < degree);
    conn[get_id(node_with_avail_if[a].first, node_with_avail_if[b].first)]++;
    conn[get_id(node_with_avail_if[b].first, node_with_avail_if[a].first)]++;
    allocated += 2;

    bool changed = false;
    if (--node_with_avail_if[a].second == 0) {
      if (a < b) {
        b--;
      }
      node_with_avail_if.erase(node_with_avail_if.begin() + a);
      changed = true;
    }
    if (--node_with_avail_if[b].second == 0) {
      node_with_avail_if.erase(node_with_avail_if.begin() + b);
      changed = true;
    }
    if (changed) {
      distrib = std::uniform_int_distribution<>(0, node_with_avail_if.size() - 1);
    }
  }

  std::cout << "Topology generated: " << std::endl;
  NetworkTopologyGenerator::print_conn_matrix(conn, num_nodes, 0);
  return conn;
  
}

int FlatDegConstraintNetworkTopologyGenerator::get_id(int i, int j) const
{
  return i * num_nodes + j;
}

int FlatDegConstraintNetworkTopologyGenerator::get_if_in_use(int node, const ConnectionMatrix & conn) const
{
  int result = 0;
  for (int i = 0; i < num_nodes; i++) {
    result += conn[get_id(node, i)];
  }
  return result;
}

BigSwitchNetworkTopologyGenerator::BigSwitchNetworkTopologyGenerator(int num_nodes)
: num_nodes(num_nodes)
{}

ConnectionMatrix BigSwitchNetworkTopologyGenerator::generate_topology() const 
{
  ConnectionMatrix conn = std::vector<int>((num_nodes+1)*(num_nodes+1), 0);
  for (int i = 0; i < num_nodes; i++) {
    conn[i * (num_nodes+1) + num_nodes] = 1;
    conn[num_nodes * (num_nodes+1) + i] = 1;
  }
  return conn;
}


DemandHeuristicNetworkOptimizer::DemandHeuristicNetworkOptimizer(MachineModel* machine) 
: L1Optimizer(machine)
{
  alpha = 0.1;
}

void DemandHeuristicNetworkOptimizer::task_added(SimTask * task) 
{
  uint64_t key;
  CommDevice *commDev;
  NominalCommDevice *ncommDev;
  switch (task->type) {
  
  case SimTask::TASK_COMM:
    commDev = reinterpret_cast<CommDevice*>(task->device);
    if (physical_traffic_demand.find(commDev->device_id) != physical_traffic_demand.end()) 
      physical_traffic_demand[commDev->device_id] = task->xfer_size;
    else
      physical_traffic_demand[commDev->device_id] += task->xfer_size;
  case SimTask::TASK_BACKWARD: 
  case SimTask::TASK_FORWARD:
    key = reinterpret_cast<uint64_t>(task->device);
    if (dev_busy_time.find(key) == dev_busy_time.end())
      dev_busy_time[key] = task->run_time;
    else
      dev_busy_time[key] += task->run_time;
  break;

  case SimTask::TASK_NOMINAL_COMM:
    ncommDev = reinterpret_cast<NominalCommDevice*>(task->device);
    if (logical_traffic_demand.find(ncommDev->device_id) != logical_traffic_demand.end()) 
      logical_traffic_demand[ncommDev->device_id] = task->xfer_size;
    else
      logical_traffic_demand[ncommDev->device_id] += task->xfer_size;
  break;

  case SimTask::TASK_ALLREDUCE:

  break;

  case SimTask::TASK_BARRIER:

  break;

  case SimTask::TASK_UPDATE:

  break;
  }
}

size_t DemandHeuristicNetworkOptimizer::edge_id(int i, int j) 
{
  return i * machine->get_total_devs() + j;
}

size_t DemandHeuristicNetworkOptimizer::unordered_edge_id(int i, int j) 
{
  return i > j ? edge_id(i, j) : edge_id(j, i);
}

typedef std::pair<uint64_t, uint64_t> DemandToIdMap;

void DemandHeuristicNetworkOptimizer::optimize(int mcmc_iter, float sim_iter_time)
{
  if (sim_iter_time < best_sim_time) 
    best_sim_time = sim_iter_time;
  float diff = sim_iter_time - curr_sim_time;
  bool change = diff < 0 ? true :
    static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) < std::exp(-alpha * diff);
  if (change) {
    curr_sim_time = sim_iter_time;
  }
  else {
    num_iter_nochange++;
  }

  if (!change && num_iter_nochange < no_improvement_th)
    return;
  
  num_iter_nochange = 0;
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  
  // This only works for flat network at the moment. 
  // to extend this to a rack based design do this for the other part of
  // the connection matrix, but the demand need to be summed.
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
    
  std::unordered_map<size_t, uint64_t> max_of_bidir;
  for (int i = 0; i < nnode; i++) {
    for (int j = 0; j < nnode; j++) {
      size_t eid = edge_id(i, j);
      if (logical_traffic_demand.find(eid) == logical_traffic_demand.end()) {
        size_t ueid = unordered_edge_id(i, j);
        uint64_t traffic_amount = logical_traffic_demand[eid];
        if (max_of_bidir.find(ueid) == max_of_bidir.end() 
            || traffic_amount > max_of_bidir[ueid]) {
          max_of_bidir[ueid] = traffic_amount;
        }
      }
    }
  }
  std::set<DemandToIdMap, std::greater<DemandToIdMap>> pq;
  for (auto &item: max_of_bidir) {
    pq.insert(DemandToIdMap(item.second, item.first));
  }

  // TODO: copy machine-switch link?
  ConnectionMatrix conn = std::vector<int>(ndevs*ndevs, 0);

  std::unordered_map<size_t, size_t> node_if_allocated;

  while (pq.size() > 0) {

    DemandToIdMap target = *pq.begin();
    pq.erase(pq.begin());

    size_t node0 = target.second / ndevs;
    size_t node1 = target.second % ndevs;
    
    // conn[target.second]++;
    conn[edge_id(node0, node1)]++;
    conn[edge_id(node1, node0)]++;

    if (node_if_allocated.find(node0) == node_if_allocated.end()) {
      node_if_allocated[node0] = 1;
    }
    else {
      node_if_allocated[node0]++;
    }
    
    if (node_if_allocated.find(node1) == node_if_allocated.end()) {
      node_if_allocated[node1] = 1;
    }
    else {
      node_if_allocated[node1]++;
    }

    // std::cout << "first is " << target.first << std::endl;
    target.first /= 2; //*= (double)conn[target.second]/(conn[target.second] + 1);
    if (target.first > 0) {
      pq.insert(target);
    }

    if (node_if_allocated[node0] == if_cnt || node_if_allocated[node1] == if_cnt) {
      for (auto it = pq.begin(); it != pq.end(); ) {
        if (node_if_allocated[node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node0, ndevs)) {
          // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else if (node_if_allocated[node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node1, ndevs)) {
          // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else {
          ++it;
        }
      } 
    }
  }

  nm->set_topology(conn);
  // simulator->topo->set_topology(conn); 
  // simulator->print_conn_matrix();

  // set<size_t> used_nodes_set;
  std::set<size_t> linked_nodes;
  for (auto item: node_if_allocated) {
    linked_nodes.insert(item.first);
  }
  std::vector<size_t> unlinked_nodes;

  for (size_t i = 0; i < ndevs; i++) {
    if (linked_nodes.find(i) == linked_nodes.end())
      unlinked_nodes.push_back(i);
  }
  // add all un-used nodes to a CC
  // std::cout << "unused node: " << std::endl;
  // for (auto n: unlinked_nodes) {
    // std::cout << "\t" << n;
  // }
  // std::cout << std::endl;
  if (unlinked_nodes.size() > 1) {

    int allocated = 0;

    size_t num_nodes = unlinked_nodes.size();
    size_t curr_node = unlinked_nodes[0];
    std::unordered_set<size_t> visited_node;
    visited_node.insert(curr_node);

    std::uniform_int_distribution<> distrib(0, num_nodes - 1);

    // conn_matrix = vector<vector<int> >(num_nodes, vector<int>(num_nodes, 0));

    while ((long)visited_node.size() != num_nodes) {
      // distrib(gen);
      int next_step = unlinked_nodes[distrib(gen)];
      if (next_step == curr_node) {
        continue;
      } 
      if (visited_node.find(next_step) == visited_node.end()) {
        // uint64_t edge_id = unordered_edge_id(next_step, curr_node);
        conn[edge_id(next_step, curr_node)]++;
        conn[edge_id(curr_node, next_step)]++;

        if (node_if_allocated.find(next_step) == node_if_allocated.end()) {
          node_if_allocated[next_step] = 1;
        }
        else {
          node_if_allocated[next_step]++;
        }
        
        if (node_if_allocated.find(curr_node) == node_if_allocated.end()) {
          node_if_allocated[curr_node] = 1;
        }
        else {
          node_if_allocated[curr_node]++;
        }
        visited_node.insert(next_step);
        curr_node = next_step;
        allocated += 1;
      }
    }
    assert(allocated == (num_nodes - 1));

    std::vector<std::pair<size_t, size_t> > node_with_avail_if;
    for (size_t i = 0; i < ndevs; i++) {
      /*
      size_t if_inuse = node_if_allocated[unlinked_nodes[i]];
      if (if_inuse < if_cnt) {
        node_with_avail_if.emplace_back(unlinked_nodes[i], if_cnt - if_inuse);
      }
      */
      size_t if_inuse = node_if_allocated[i];
      if (if_inuse < if_cnt) {
        node_with_avail_if.emplace_back(i, if_cnt - if_inuse);
      }
    }

    distrib = std::uniform_int_distribution<>(0, node_with_avail_if.size() - 1);
    size_t a = 0, b = 0;

    std::unordered_set<size_t> unused_node_set = 
      std::unordered_set<size_t>(unlinked_nodes.begin(), unlinked_nodes.end());
    while (!maxed(node_if_allocated, if_cnt, ndevs)/*maxed(node_if_allocated, unused_node_set, if_cnt)*/) {
      a = distrib(gen);
      while ((b = distrib(gen)) == a);

      size_t node0 = node_with_avail_if[a].first;
      size_t node1 = node_with_avail_if[b].first;

        conn[edge_id(node0, node1)]++;
        conn[edge_id(node1, node0)]++;

      if (node_if_allocated.find(node0) == node_if_allocated.end()) {
        node_if_allocated[node0] = 1;
      }
      else {
        node_if_allocated[node0]++;
      }
      
      if (node_if_allocated.find(node1) == node_if_allocated.end()) {
        node_if_allocated[node1] = 1;
      }
      else {
        node_if_allocated[node1]++;
        }
      allocated += 1;

      bool changed = false;
      if (--node_with_avail_if[a].second == 0) {
        if (a < b) {
          b--;
        }
        node_with_avail_if.erase(node_with_avail_if.begin() + a);
        changed = true;
      }
      if (--node_with_avail_if[b].second == 0) {
        node_with_avail_if.erase(node_with_avail_if.begin() + b);
        changed = true;
      }
      if (changed) {
        distrib = std::uniform_int_distribution<>(0, node_with_avail_if.size() - 1);
      }
    }
    std::cerr << "finished allocating CC for unused nodes. Network:" << std::endl;
    //simulator->print_conn_matrix();
  }

  // Make all CC connected
  std::unordered_map<uint64_t, uint64_t> logical_id_to_demand;
  for (auto & item: max_of_bidir) {
    logical_id_to_demand[item.second] = item.first;
  }

  connect_cc(logical_id_to_demand, conn); 
  nm->set_topology(conn); 
  nm->update_route();
  // simulator->print_conn_matrix();

}

size_t DemandHeuristicNetworkOptimizer::get_if_in_use(size_t node, const ConnectionMatrix & conn) 
{
  size_t result = 0;
  for (int i = 0; i < machine->get_num_nodes(); i++) {
    result += conn[edge_id(node, i)];
  }
  return result;
}

bool DemandHeuristicNetworkOptimizer::add_link(size_t i, size_t j, ConnectionMatrix & conn) 
{
  assert(i != j);
  if (get_if_in_use(i, conn) >= if_cnt || get_if_in_use(j, conn) >= if_cnt) {
    return false;
  }
  conn[edge_id(i, j)]++;
  conn[edge_id(j, i)]++;
  return true;
}

void DemandHeuristicNetworkOptimizer::remove_link(size_t i, size_t j, ConnectionMatrix & conn) 
{
  assert(i != j);
  if (conn[edge_id(i, j)] > 0) {
    conn[edge_id(i, j)]--;
    conn[edge_id(j, i)]--;
  }
}

void DemandHeuristicNetworkOptimizer::connect_cc(
        std::unordered_map<uint64_t, uint64_t> & logical_id_to_demand, 
        ConnectionMatrix &conn) 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  size_t num_nodes = nm->num_nodes;
  size_t ndevs = nm->num_switches + num_nodes;

  // reconnect phase
  // find connected components 
  int n_cc = 0;
  std::vector<int> node_to_ccid = std::vector<int>(num_nodes, -1);
  std::vector<std::set<size_t> > ccs;
  // node_to_ccid[0] = 0;
  std::queue<size_t> search_q;

  for (size_t i = 0; i < num_nodes; i++) {
    if (node_to_ccid[i] == -1) {
      search_q.push(i);
      node_to_ccid[i] = n_cc++;
      ccs.emplace_back();
      ccs.back().insert(i);
      while (!search_q.empty()) {
        size_t curr = search_q.front();
        // node_to_ccid[curr] = n_cc;
        search_q.pop();
        for (size_t j = 0; j < num_nodes; j++) {
          if (curr != j && conn[edge_id(curr, j)] > 0 && node_to_ccid[j] == -1) {
            node_to_ccid[j] = node_to_ccid[curr];
            ccs.back().insert(j);
            search_q.push(j);
          }
        }
      }
      // n_cc++;
    }
    else {
      continue;
    }
  }

  // std::cout << "n_cc " << n_cc << std::endl;
  // std::cout << "node_to_ccid:" << std::endl;

  // for (size_t i = 0; i < node_to_ccid.size(); i++) {
  //   std::cout << "\t" << i << ", " << node_to_ccid[i] << std::endl;
  // }
  
  // for (size_t i = 0; i < ccs.size(); i++) {
  //   std::cout << "CC " << i << ": " << std::endl;
  //   for (size_t v: ccs[i]) {
  //     std::cout << "\t" << v;
  //   }
  //   std::cout << std::endl;
  // }
  
  assert(n_cc > 0);
  if (n_cc > 1) {

    // find the two lowest demanded line in the two CC and do a 2er
    // size_t cc0 = 0;
    // size_t cc1 = 1;

    //std::vector<pair<uint64_t, uint64_t> > cc0_d, cc1_d;
    int v00, v01, v10, v11;

    while (n_cc > 1) {

      if (ccs[0].size() == 1 && ccs[1].size() == 1) {
        bool success = add_link(*ccs[0].begin(), *ccs[1].begin(), conn);
        assert(success);
        success = add_link(*ccs[0].begin(), *ccs[1].begin(), conn);
        assert(success);
      }

      else if (ccs[0].size() == 1 || ccs[1].size() == 1) { // ccs[1].size > 1

        size_t singleton = ccs[0].size() == 1 ? 0 : 1;
        size_t group = singleton == 0 ? 1 : 0;

        uint64_t e_to_remove = 0;
        uint64_t min_demand = std::numeric_limits<uint64_t>::max();

        for (size_t i = 0; i < num_nodes; i++) {
          for (size_t j = i + 1; j < num_nodes; j++) {
            if (ccs[group].find(i) != ccs[group].end() && 
                ccs[group].find(j) != ccs[group].end() && 
                conn[edge_id(i, j)] > 0) {
              uint64_t ueid = unordered_edge_id(i, j);
              if (logical_id_to_demand.find(ueid) == logical_id_to_demand.end()) {
                e_to_remove = ueid;
                break;
              }
              else {
                if (logical_id_to_demand[ueid] < min_demand) {
                  min_demand = logical_id_to_demand[ueid];
                  e_to_remove = ueid;
                }
              }
            }
          }
        }
        assert(e_to_remove != 0);

        // std::cout << "1-n removing " << e_to_remove % ndevs << ", " <<  e_to_remove / ndevs << std::endl;
        remove_link(e_to_remove % ndevs, e_to_remove / ndevs, conn);
        bool success = add_link(*ccs[singleton].begin(), e_to_remove % ndevs, conn);
        assert(success);
        success = add_link(*ccs[singleton].begin(), e_to_remove / ndevs, conn);
        assert(success);

      }

      else {
        std::vector<uint64_t> new_links;
        for (size_t i = 0; i < num_nodes; i++) {
          for (size_t j = i + 1; j < num_nodes; j++) {
            if (conn[edge_id(i, j)] > 0) {
              new_links.emplace_back(unordered_edge_id(i, j));
            }
          }
        }

        sort(new_links.begin(), new_links.end(), [&] (uint64_t lhs, uint64_t rhs) {
          auto liter = logical_id_to_demand.find(lhs);
          uint64_t l = liter == logical_id_to_demand.end() ? 0 : liter->second;
          auto riter = logical_id_to_demand.find(rhs);
          uint64_t r = riter == logical_id_to_demand.end() ? 0 : riter->second;
          return l < r;
        });
        // cc0_d.clear();
        // cc1_d.clear();
        v00 = v01 = v10 = v11 = -1;

        for (auto & item: new_links) {
          size_t n0 = item % ndevs;
          size_t n1 = item / ndevs;
          if (v00 == -1 && 
              ccs[0].find(n0) != ccs[0].end() &&
              ccs[0].find(n1) != ccs[0].end()) {
            v00 = n0;
            v01 = n1;
          }
          else if (v10 == -1 && 
              ccs[1].find(n0) != ccs[1].end() &&
              ccs[1].find(n1) != ccs[1].end()) {
            v10 = n0;
            v11 = n1;
          }
          if (v00 != -1 && v10 != -1) {
            
            // std::cout << "swappig " << v00 << ", " << v01 << " and " << v10 << ", " << v11 << std::endl;
            remove_link(v00, v01, conn);
            remove_link(v10, v11, conn);
            bool success = add_link(v00, v11, conn);
            assert(success);
            success = add_link(v01, v10, conn);
            assert(success);

            break;
          }
        }
        assert(v00 != -1);
      }
      n_cc--;
      ccs[1].insert(ccs[0].begin(), ccs[0].end());
      ccs.erase(ccs.begin());
    }
  }
  // assert(check_connected());
}

void DemandHeuristicNetworkOptimizer::reset() 
{
  dev_busy_time.clear();
  physical_traffic_demand.clear();
  logical_traffic_demand.clear();

  best_sim_time = std::numeric_limits<float>::max();
  curr_sim_time = std::numeric_limits<float>::max();
  num_iter_nochange = 0;
  
}

// TODO
void* DemandHeuristicNetworkOptimizer::export_information()
{
  return nullptr;
}
