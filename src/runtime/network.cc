#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <utility>
#include <cmath>
#include <unordered_set>

#include "simulator.h"
// #define EDGE(a, b, n) ((a) > (b) ? ((a) * (n) + (b)) : ((b) * (n) + (a)))
#define PRINT_EDGE(e, n) do {std::cout << "(" << e / n << ", " << e % n << ")";} while (0);

#define INSERT_OR_ADD(_map, _key, _val) do {                                \
  if ((_map).find(_key) == (_map).end()) {                                  \
    (_map)[(_key)] = _val;                                                  \
  } else {                                                                  \
    (_map)[(_key)] += _val;                                                 \
  }                                                                         \
} while (0);                                                                \

static std::random_device rd; 
static std::mt19937 gen = std::mt19937(rd()); 
static std::uniform_real_distribution<double> unif(0, 1);

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
  assert(result.size() || src_node == dst_node);

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

#ifdef DEBUG_PRINT
  std::cout << "Topology generated: " << std::endl;
  NetworkTopologyGenerator::print_conn_matrix(conn, num_nodes, 0);
#endif
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
  no_improvement_th = 100;
  best_sim_time = std::numeric_limits<float>::max();
  curr_sim_time = std::numeric_limits<float>::max();
}

void DemandHeuristicNetworkOptimizer::task_added(SimTask * task) 
{
  uint64_t key;
  CommDevice *commDev;
  NominalCommDevice *ncommDev;
  switch (task->type) {
  
  case SimTask::TASK_COMM:
    commDev = reinterpret_cast<CommDevice*>(task->device);
    INSERT_OR_ADD(physical_traffic_demand, commDev->device_id, task->xfer_size);
  case SimTask::TASK_BACKWARD: 
  case SimTask::TASK_FORWARD:
    key = reinterpret_cast<uint64_t>(task->device);
    INSERT_OR_ADD(dev_busy_time, key, task->run_time);
  break;

  case SimTask::TASK_NOMINAL_COMM:
    ncommDev = reinterpret_cast<NominalCommDevice*>(task->device);
    INSERT_OR_ADD(logical_traffic_demand, ncommDev->device_id, task->xfer_size);
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
  if (sim_iter_time < best_sim_time) {
    best_sim_time = sim_iter_time;
  }
  float diff = sim_iter_time - curr_sim_time;
  std::cerr << "sim_iter_time: " << sim_iter_time << ", curr_sim_time: " << curr_sim_time 
            << ", best_iter_time: " << best_sim_time << std::endl;
  bool change = diff < 0 ? true :
    static_cast<float>(std::rand()) / static_cast<float>(static_cast<float>(RAND_MAX)) < std::exp(-alpha * diff);
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
  
  // TODO: copy machine-switch link?
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  ConnectionMatrix conn = std::vector<int>(ndevs*ndevs, 0);
  std::unordered_map<size_t, uint64_t> max_of_bidir;
  std::unordered_map<size_t, size_t> node_if_allocated;

  optimize_demand(conn, max_of_bidir, node_if_allocated);
#ifdef DEBUG_PRINT
  NetworkTopologyGenerator::print_conn_matrix(conn, nnode, 0);
#endif

  connect_unused_node(conn, node_if_allocated);

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

void DemandHeuristicNetworkOptimizer::optimize_demand(
  ConnectionMatrix &conn,
  std::unordered_map<size_t, uint64_t> &max_of_bidir,
  std::unordered_map<size_t, size_t> &node_if_allocated) 
{
  // This only works for flat network at the moment. 
  // to extend this to a rack based design do this for the other part of
  // the connection matrix, but the demand need to be summed.
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
    
  for (int i = 0; i < nnode; i++) {
    for (int j = 0; j < nnode; j++) {
      size_t eid = edge_id(i, j);
      if (logical_traffic_demand.find(eid) != logical_traffic_demand.end()) {
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
    // mod: pre-unscale the demand
    pq.insert(DemandToIdMap(item.second/(1 << conn[item.first]), item.first));
  }

  while (pq.size() > 0) {

    DemandToIdMap target = *pq.begin();
    pq.erase(pq.begin());

    size_t node0 = target.second / ndevs;
    size_t node1 = target.second % ndevs;
    
    // conn[target.second]++;
    conn[edge_id(node0, node1)]++;
    conn[edge_id(node1, node0)]++;

    INSERT_OR_ADD(node_if_allocated, node0, 1);
    INSERT_OR_ADD(node_if_allocated, node1, 1);

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
}

void DemandHeuristicNetworkOptimizer::connect_unused_node(
  ConnectionMatrix &conn,
  std::unordered_map<size_t, size_t> &node_if_allocated) 
{
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  // set<size_t> used_nodes_set;
  std::set<size_t> linked_nodes;
  for (auto item: node_if_allocated) {
    linked_nodes.insert(item.first);
  }
  std::vector<size_t> unlinked_nodes;

  for (size_t i = 0; i < nnode; i++) {
    if (linked_nodes.find(i) == linked_nodes.end())
      unlinked_nodes.push_back(i);
  }

  // add all un-used nodes to a CC
#ifdef DEBUG_PRINT
  std::cout << "unused node: " << std::endl;
  for (auto n: unlinked_nodes) {
    std::cout << "\t" << n;
  }
  std::cout << std::endl;
#endif

  if (unlinked_nodes.size() > 1) {

    int allocated = 0;

    size_t num_nodes = unlinked_nodes.size();
    size_t curr_node = unlinked_nodes[0];
    std::unordered_set<size_t> visited_node;
    visited_node.insert(curr_node);

    std::uniform_int_distribution<> distrib(0, num_nodes - 1);

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

        INSERT_OR_ADD(node_if_allocated, next_step, 1);
        INSERT_OR_ADD(node_if_allocated, curr_node, 1);
        
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

      INSERT_OR_ADD(node_if_allocated, node0, 1);
      INSERT_OR_ADD(node_if_allocated, node1, 1);
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
#ifdef DEBUG_PRINT
    std::cerr << "finished allocating CC for unused nodes. Network:" << std::endl;
    // simulator->print_conn_matrix();
    NetworkTopologyGenerator::print_conn_matrix(conn, num_nodes, 0);
#endif
  }
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

#ifdef DEBUG_PRINT
  std::cout << "n_cc " << n_cc << std::endl;
  std::cout << "node_to_ccid:" << std::endl;

  for (size_t i = 0; i < node_to_ccid.size(); i++) {
    std::cout << "\t" << i << ", " << node_to_ccid[i] << std::endl;
  }
  
  for (size_t i = 0; i < ccs.size(); i++) {
    std::cout << "CC " << i << ": " << std::endl;
    for (size_t v: ccs[i]) {
      std::cout << "\t" << v;
    }
    std::cout << std::endl;
  }
#endif
  
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

  // best_sim_time = std::numeric_limits<float>::max();
  // curr_sim_time = std::numeric_limits<float>::max();
  // num_iter_nochange = 0;
  
}

// TODO
void* DemandHeuristicNetworkOptimizer::export_information()
{
  return nullptr;
}

NSDI22Heuristic::NSDI22Heuristic(MachineModel* machine, size_t dp_deg, size_t mp_deg)
{
  this->net_machine = static_cast<NetworkedMachineModel*>(machine);
  this->dp_deg = dp_deg;
  this->mp_deg = mp_deg;
}


void NSDI22Heuristic::fw_task_added(SimTask* task) 
{
  return;
}

void NSDI22Heuristic::bw_task_added(SimTask* task) 
{
  return;
}

void NSDI22Heuristic::mp_lcomm_added(CommDevice* device, size_t xfer_size) 
{
  mp_tm_logical[device->device_id] += xfer_size;
}

void NSDI22Heuristic::dp_lcomm_added(SimTask* task) 
{
  assert(task->type == SimTask::TASK_NOMINAL_COMM);
  dp_tm_logical[task->device->device_id] += task->xfer_size;
}

void NSDI22Heuristic::mp_pcomm_added(CommDevice* device, size_t xfer_size, bool dir) 
{
  // assert(task->type == SimTask::TASK_COMM);
  if (dir)
    mp_tm_physical_dir[device->device_id] += xfer_size;
  else 
    mp_tm_physical_indir[device->device_id] += xfer_size;
}

void NSDI22Heuristic::dp_pcomm_added(SimTask* task, bool dir) 
{
  assert(task->type == SimTask::TASK_COMM);
  if (dir)
    dp_tm_physical_dir[task->device->device_id] += task->xfer_size;
  else 
    dp_tm_physical_indir[task->device->device_id] += task->xfer_size;
}

void NSDI22Heuristic::ar_task_added(SimTask* task) 
{
  assert(task->type == SimTask::TASK_ALLREDUCE);
  ar_tasks.push_back(task);
}

void NSDI22Heuristic::generate_dp_topology() 
{
  construct_dp_mat();
  DemandHeuristicNetworkOptimizer dheuristic{net_machine};
  dheuristic.if_cnt = dp_deg + mp_deg;
  dheuristic.logical_traffic_demand = dp_tm_logical;
  dheuristic.machine = net_machine;

  size_t nnode = net_machine->get_num_nodes();
  size_t ndevs = net_machine->get_total_devs();
  std::unordered_map<size_t, uint64_t> max_of_bidir;
  std::unordered_map<size_t, size_t> node_if_allocated;
  // ConnectionMatrix conn = std::vector<int>(ndevs*ndevs, 0);
  for (int i = 0; i < net_machine->num_nodes; i++) {
    for (int j = 0; j < i; j++) {
      if (net_machine->conn_matrix[edge_id(i, j)] > 0) {
        INSERT_OR_ADD(node_if_allocated, i, net_machine->conn_matrix[edge_id(i, j)]);
        INSERT_OR_ADD(node_if_allocated, j, net_machine->conn_matrix[edge_id(i, j)]);
      }
    }
  }

  dheuristic.optimize_demand(net_machine->conn_matrix, max_of_bidir, node_if_allocated);
  // for (int i = 0; i < net_machine->get_num_nodes(); i++) {
  //   for (int j = 0; net_machine->get_num_nodes(); j++) {
  //     net_machine->conn_matrix[edge_id(i, j)] += conn[edge_id(i, j)];
  //   }
  // }
  net_machine->set_topology(net_machine->conn_matrix);
  net_machine->update_route();
  // NetworkTopologyGenerator::print_conn_matrix(net_machine->conn_matrix, net_machine->num_nodes, 0);

}

void NSDI22Heuristic::simplify_mp_topology() 
{
  for (int i = 0; i < net_machine->get_num_nodes(); i++) {
    for (int j = 0; j < i; j++) {
      if (net_machine->conn_matrix[edge_id(i, j)] > 0 && 
          mp_tm_physical_dir[edge_id(i, j)] == 0 && 
          mp_tm_physical_indir[edge_id(i, j)] == 0 && 
          mp_tm_physical_dir[edge_id(j, i)] == 0 && 
          mp_tm_physical_indir[edge_id(j, i)] == 0)
      {
        net_machine->conn_matrix[edge_id(i, j)] = 0;
        net_machine->conn_matrix[edge_id(j, i)] = 0;
      }
    }
  }
  net_machine->set_topology(net_machine->conn_matrix);
}

void NSDI22Heuristic::optimize_indirection() 
{
  // std::map<uint64_t, uint64_t> edge_to_indirection;
  return;
}

void NSDI22Heuristic::construct_dp_mat() 
{
  for (SimTask * allreduce_task: ar_tasks) {
    int n_participants = allreduce_task->next_tasks.size();
  
    // recall that next_task stores node group in this case
    CompDevice * device = net_machine->get_gpu(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
    MemDevice * src_mem = net_machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
    MemDevice * dst_mem;
    float indiv_xfersize = (2.0 * (n_participants-1))/n_participants * allreduce_task->xfer_size;
    for (int i = 0; i < n_participants; i++) {
      dst_mem = net_machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[(i+1)%n_participants]));
      std::vector<CommDevice *> path = net_machine->get_comm_path(src_mem, dst_mem);
      if (path.size() == 0)
        continue;
      assert(path.size() == 1 && reinterpret_cast<CommDevice*>(path[0])->comm_type == CommDevice::NW_NOMINAL);
      dp_tm_logical[path[0]->device_id] += indiv_xfersize;
      src_mem = dst_mem;
      std::vector<CommDevice *> route = 
        static_cast<NominalCommDevice*>(path[0])->expand_to_physical(); 
      if (route.size() == 1)
        INSERT_OR_ADD(dp_tm_physical_dir, route[0]->device_id, indiv_xfersize);
      for (CommDevice * d: route) {
        INSERT_OR_ADD(dp_tm_physical_indir, d->device_id, indiv_xfersize);
      }
    }
  }
}


void NSDI22Heuristic::reset() 
{
  dp_tm_logical.clear();
  dp_tm_physical_dir.clear();
  dp_tm_physical_indir.clear();
  
  mp_tm_logical.clear();
  mp_tm_physical_dir.clear();
  mp_tm_physical_indir.clear();

  ar_tasks.clear();
}

size_t NSDI22Heuristic::edge_id(int i, int j) 
{
  return i * net_machine->get_total_devs() + j;
}

size_t NSDI22Heuristic::unordered_edge_id(int i, int j) 
{
  return i > j ? edge_id(i, j) : edge_id(j, i);
}

TwoDimTorusNetworkTopologyGenerator::TwoDimTorusNetworkTopologyGenerator(int num_nodes) 
{
  this->num_nodes = num_nodes;
  int edge = (int)std::sqrt(num_nodes) + 1;
  int extra = edge * edge - num_nodes;
  if (extra <= edge) {
    nrows = edge - 1;
    ncols = edge;
    nextras = edge - extra;
  }
  else {
    nrows = edge - 1;
    ncols = edge - 1;
    nextras = num_nodes - nrows * ncols;
  }
  assert(nrows * ncols + nextras == num_nodes);
}

ConnectionMatrix TwoDimTorusNetworkTopologyGenerator::generate_topology() const 
{
  ConnectionMatrix conn = ConnectionMatrix(num_nodes * num_nodes, 0);
  for (int row = 0; row < nrows + 1; row++) {
    for (int col = 0; col < ncols; col++) {
      if (row == nrows && col == nextras)
        break;
      int me = row * ncols + col;
      int my_right = row * ncols + ((col + 1) % (row == nrows ? nextras : ncols));
      int my_down = ((row + 1) % (col < nextras ? (nrows + 1) : nrows)) * ncols + col;
      // my right
      conn[get_id(me, my_right)] = 1;
      conn[get_id(my_right, me)] = 1;
      // my down
      conn[get_id(me, my_down)] = 1;
      conn[get_id(my_down, me)] = 1;
    }
  }
  // NetworkTopologyGenerator::print_conn_matrix(conn, num_nodes, 0);
  // assert(false);
  return conn;
}

void TwoDimTorusNetworkTopologyGenerator::test() 
{
  TwoDimTorusNetworkTopologyGenerator gen{144};
  assert(gen.nrows == 12);
  assert(gen.ncols == 12);
  assert(gen.nextras == 0);
  fprintf(stderr, "\n144:\n");
  ConnectionMatrix conn = gen.generate_topology();
  for (int i = 0; i < 144; i++) {
    for (int j = 0; j < 144; j++) {
      if (conn[gen.get_id(i, j)] > 0) {
        fprintf(stderr, "(%d, %d) - (%d, %d)\n", i%gen.ncols, i/gen.ncols, j%gen.ncols, j/gen.ncols);
      }
    }
  }
  gen = TwoDimTorusNetworkTopologyGenerator{138};
  assert(gen.nrows == 11);
  assert(gen.ncols == 12);
  assert(gen.nextras == 6);
  fprintf(stderr, "\n138:\n");
  conn = gen.generate_topology();
  for (int i = 0; i < 138; i++) {
    for (int j = 0; j < 138; j++) {
      if (conn[gen.get_id(i, j)] > 0) {
        fprintf(stderr, "(%d, %d) - (%d, %d)\n", i%gen.ncols, i/gen.ncols, j%gen.ncols, j/gen.ncols);
      }
    }
  }
  gen = TwoDimTorusNetworkTopologyGenerator{128};
  assert(gen.nrows == 11);
  assert(gen.ncols == 11);
  assert(gen.nextras == 7);
  fprintf(stderr, "\n128:\n");
  conn = gen.generate_topology();
  for (int i = 0; i < 128; i++) {
    for (int j = 0; j < 128; j++) {
      if (conn[gen.get_id(i, j)] > 0) {
        fprintf(stderr, "(%d, %d) - (%d, %d)\n", i%gen.ncols, i/gen.ncols, j%gen.ncols, j/gen.ncols);
      }
    }
  }
}

int TwoDimTorusNetworkTopologyGenerator::get_id(int i, int j) const 
{
  return i * num_nodes + j;
}

TwoDimTorusW2TURNRouting::TwoDimTorusW2TURNRouting(const ConnectionMatrix & c, 
                           const std::map<size_t, CommDevice*>& devmap, 
                           int nrows, int ncols, int nextras, int total_devs)
: conn(c), devmap(devmap), nrows(nrows), ncols(ncols), nextras(nextras), total_devs(total_devs)
{

}

void TwoDimTorusW2TURNRouting::clear()
{
  cached_routes.clear();
}

size_t TwoDimTorusW2TURNRouting::nodeid(int x, int y) const
{
  return x + y * ncols;
} 

int TwoDimTorusW2TURNRouting::getx(size_t nodeid) const
{
  return nodeid % ncols;
}

int TwoDimTorusW2TURNRouting::gety(size_t nodeid) const
{
  return nodeid / ncols;
}

bool TwoDimTorusW2TURNRouting::lastl(size_t nodeid) const
{
  return gety(nodeid) == nrows;
}

int TwoDimTorusW2TURNRouting::delta(int s, int d, int k) const
{
  return std::min(std::abs(s-d), k-std::abs(s-d));
}

uint64_t TwoDimTorusW2TURNRouting::edgeid(int src, int dst) const 
{
  return src * total_devs + dst;
}

#define MOD(a, b) ((a) % (b)) < 0 ? ((a) % (b)) + (b) : ((a) % (b))

int TwoDimTorusW2TURNRouting::mindir(int s, int d, int k) const
{
  int d1 = MOD((s - d), k);
  int d2 = MOD((d - s), k);
  if (d1 == d2) {
    return 0;
  }
  return d1 > d2 ? 1 : -1;
}

EcmpRoutes TwoDimTorusW2TURNRouting::get_routes(int src_node, int dst_node) 
{
  if (src_node == dst_node) {
    cached_routes[edgeid(src_node, dst_node)] = {std::vector<float>({1}), std::vector<Route>({})};
    return cached_routes[edgeid(src_node, dst_node)];
  }
  Route route;
  int xsrc = getx(src_node);
  int xdst = getx(dst_node);
  int ysrc = gety(src_node);
  int ydst = gety(dst_node);

  if (xsrc == xdst) {
    route_y_wrd(ysrc, ydst, xsrc, xsrc < nextras ? nrows+1 : nrows, route);
    cached_routes[edgeid(src_node, dst_node)] = {std::vector<float>({1}), std::vector<Route>({route})};
    return cached_routes[edgeid(src_node, dst_node)];
  }
  if (ysrc == ydst) {
    route_x_wrd(xsrc, xdst, ysrc, ysrc == nrows ? nextras : ncols, route);
    cached_routes[edgeid(src_node, dst_node)] = {std::vector<float>({1}), std::vector<Route>({route})};
    return cached_routes[edgeid(src_node, dst_node)];
  }
  int kxsrc = xsrc >= nextras ? ncols : nextras;
  int kysrc = xsrc >= nextras ? nrows : nrows + 1;
  int kxdst = xdst >= nextras ? ncols : nextras;
  int kydst = xdst >= nextras ? nrows : nrows + 1;
  int first_direction = (unif(gen) > 0.5) ? 0 : 1;
  if (first_direction) { // XYX
    std::uniform_int_distribution<> x_uni
        {0, std::min(kxsrc, kxdst) - 1};
    int xmid = x_uni(gen); 
    // first X
    if (kxsrc % 2 != 0) {
      if (route_in_mindir_end(xsrc, xmid, xdst, kxsrc) || 
          unif(gen) > (delta(xsrc, xdst, kxsrc) / (float)kxsrc)) {
        route_x(xsrc, xmid, ysrc, mindir(xsrc, xmid, kxsrc), kxsrc, route); 
      }
      else {
        route_x(xsrc, xmid, ysrc, -mindir(xsrc, xmid, kxsrc), kxsrc, route); 
      }
    } else {
      int mind = mindir(xsrc, xmid, kxsrc);
      if (mind == 0) {
        route_x(xsrc, xmid, ysrc, even_choose_dir(xsrc, xmid, xdst, kxsrc), kxsrc, route); 
      }
      else {
        route_x(xsrc, xmid, ysrc, mind, kxsrc, route); 
      }
    }
    // Y
    int kymid = xmid < nextras ? nrows + 1 : nrows;
    if (kymid % 2 != 0) {
      if (route_in_mindir_mid(xsrc, xdst, ysrc, ydst, xmid, kymid)) {
        route_y(ysrc, ydst, xmid, mindir(ysrc, ydst, kymid), kymid, route);
      }
      else {
        route_y_wrd(ysrc, ydst, xmid, kymid, route);
      }      
    } else {
      route_y_wrd(ysrc, ydst, xmid, kymid, route);
    }
    // last X
    if (kxdst % 2 != 0) {
      if (route_in_mindir_end(xmid, xdst, xsrc, kxdst) || 
          unif(gen) > (delta(xsrc, xdst, kxdst) / (float)kxdst)) {
        route_x(xmid, xdst, ydst, mindir(xmid, xdst, kxdst), kxdst, route); 
      }
      else {
        route_x(xmid, xdst, ydst, -mindir(xmid, xdst, kxdst), kxdst, route); 
      }
    } else {
      int mind = mindir(xmid, xdst, kxdst);
      if (mind == 0) {
        route_x(xmid, xdst, ydst, even_choose_dir(xmid, xdst, xsrc, kxdst), kxdst, route); 
      }
      else {
        route_x(xmid, xdst, ydst, mind, kxdst, route); 
      }
    }
  }
  else {
    std::uniform_int_distribution<> y_uni
      {0, std::min(kysrc, kydst) - 1};
    int ymid = y_uni(gen); 
    int dir;
    // first Y
    if (kysrc % 2 != 0) {
      if (route_in_mindir_end(ysrc, ymid, ydst, kysrc) || 
          unif(gen) > (delta(ysrc, ydst, kysrc) / (float)kysrc)) {
        route_y(ysrc, ymid, xsrc, mindir(ysrc, ymid, kysrc), kysrc, route); 
      }
      else {
        route_y(ysrc, ymid, xsrc, -mindir(ysrc, ymid, kysrc), kysrc, route); 
      }
    } else {
      int mind = mindir(ysrc, ymid, kysrc);
      if (mind == 0) {
        route_y(ysrc, ymid, xsrc, even_choose_dir(ysrc, ymid, ydst, kysrc), kysrc, route); 
      }
      else {
        route_y(ysrc, ymid, xsrc, mind, kysrc, route); 
      }
    }
    // X
    int kxmid = ymid == nrows ? nextras : ncols;
    if (kxmid % 2 != 0) {
      if (route_in_mindir_mid(ysrc, ydst, xsrc, xdst, ymid, kxmid)) {
        route_x(xsrc, xdst, ymid, mindir(xsrc, xdst, kxmid), kxmid, route);
      }
      else {
        route_x_wrd(xsrc, xdst, ymid, kxmid, route);
      }
    } else {
      route_x_wrd(xsrc, xdst, ymid, kxmid, route);
    }
    // LAST Y
    if (kydst %2 != 0) {
      if (route_in_mindir_end(ymid, ydst, ysrc, kydst) || 
          unif(gen) > (delta(ysrc, ydst, kydst) / (float)kydst)) {
        route_y(ymid, ydst, xdst, mindir(ymid, ydst, kydst), kydst, route); 
      }
      else {
        route_y(ymid, ydst, xdst, -mindir(ymid, ydst, kydst), kydst, route); 
      }
    } else {
      int mind = mindir(ymid, ydst, kydst);
      if (mind == 0) {
        route_y(ymid, ydst, xdst, even_choose_dir(ymid, ydst, ysrc, kydst), kydst, route); 
      }
      else {
        route_y(ymid, ydst, xdst, mind, kydst, route); 
      }
    }
  }
  cached_routes[edgeid(src_node, dst_node)] = {std::vector<float>{1}, std::vector<Route>{route}};
  return cached_routes[edgeid(src_node, dst_node)];
}

int TwoDimTorusW2TURNRouting::even_choose_dir(int s, int d, int m, int k) const
{
  int mind = mindir(s, d, k);
  if (mind != 0)
    return mind;
  if (m == d) {
    return unif(gen) > 0.5 ? 1 : -1;
  }
  if (MOD(d-s, k) > MOD(m-s, k)) {
    return -1;
  }
  return 1;
} 

bool TwoDimTorusW2TURNRouting::on_the_way(int s, int d, int m, int k) const
{
  if (s == d) 
    return false;
  int mind = mindir(s, d, k);
  assert(mind != 0);
  if (mind >= 0) {
    if (MOD(d-s, k) > MOD(m-s, k)) {
      return true;
    }
    return false;
  } 
  else {
    if (MOD(d-s, k) < MOD(m-s, k)) {
      return true;
    }
    return false;
  }
}

bool TwoDimTorusW2TURNRouting::route_in_mindir_end(int p1, int p2, int p3, int k) const
{
  if (delta(p1, p2, k) != k/2)
    return true;
  if (on_the_way(p1, p2, p3, k))
    return true;
  if (delta(p1, p3, k) == k/2)
    return true;
  return false;
}

bool TwoDimTorusW2TURNRouting::route_in_mindir_mid(int p1, int p2, int p3, int p4, int p5, int k) const
{
  if (p1 != p2 && delta(p3, p4, k) == k/2 && (p5 == p1 || p5 == p2))
    return true;
  return false;
}

void TwoDimTorusW2TURNRouting::route_x(int srcx, int dstx, int y, int dir, int k, Route& route)
{
  if (srcx == dstx)
    return;
  assert(dir != 0);
  int curr = srcx;
  do {
    int next = MOD(curr+dir,k);
    route.emplace_back(devmap.at(edgeid(nodeid(curr, y), nodeid(next, y))));
    curr = next;
  } while (curr != dstx);
}

void TwoDimTorusW2TURNRouting::route_y(int srcy, int dsty, int x, int dir, int k, Route& route)
{
  if (srcy == dsty)
    return;
  assert(dir != 0);
  int curr = srcy;
  do {
    int next = MOD(curr+dir, k);
    route.emplace_back(devmap.at(edgeid(nodeid(x, curr), nodeid(x, next))));
    curr = next;
  } while (curr != dsty);
}

void TwoDimTorusW2TURNRouting::route_x_wrd(int srcx, int dstx, int y, int k, Route& route)
{
  if (srcx == dstx)
    return;
  if (k % 2 != 0) {
    if (unif(gen) > (delta(srcx, dstx, k) / (float)k)) {
      route_x(srcx, dstx, y, mindir(srcx, dstx, k), k, route);
    } else {
      route_x(srcx, dstx, y, -mindir(srcx, dstx, k), k, route);
    }
  }
  else {
    if (k == 2) {
      route_x(srcx, dstx, y, 1, k, route);
    } else {
      int mind = mindir(srcx, dstx, k);
      if (mind == 0) mind = unif(gen) > 0.5 ? -1 : 1;
      if (unif(gen) > ((delta(srcx, dstx, k)-1) / (float)(k-2))) {
        route_x(srcx, dstx, y, mind, k, route);
      } else {
        route_x(srcx, dstx, y, -mind, k, route);
      }
    }
  }
}

void TwoDimTorusW2TURNRouting::route_y_wrd(int srcy, int dsty, int x, int k, Route& route)
{
  if (srcy == dsty)
    return;
  if (k % 2 != 0) {
    if (unif(gen) > (delta(srcy, dsty, k) / (float)k)) {
      route_y(srcy, dsty, x, mindir(srcy, dsty, k), k, route);
    } else {
      route_y(srcy, dsty, x, -mindir(srcy, dsty, k), k, route);
    }
  }
  else {
    if (k == 2) {
      route_y(srcy, dsty, x, 1, k, route);
    } else {
      int mind = mindir(srcy, dsty, k);
      if (mind == 0) mind = unif(gen) > 0.5 ? -1 : 1;
      if (unif(gen) > ((delta(srcy, dsty, k)-1) / (float)(k-2))) {
        route_y(srcy, dsty, x, mind, k, route);
      } else {
        route_y(srcy, dsty, x, -mind, k, route);
      }
    }
  }
}

void TwoDimTorusW2TURNRouting::test() 
{
  int total_devs = 432;
  TwoDimTorusNetworkTopologyGenerator gen{total_devs};
  ConnectionMatrix conn = gen.generate_topology();
  std::map<size_t, CommDevice*> devmap; 
  for (int i = 0; i < total_devs; i++) {
    for (int j = 0; j < total_devs; j++) {
      // if (conn_matrix[i * total_devs + j] > 0) {
        int device_id = i * total_devs + j;
        std::string link_name = "LINK " + std::to_string(i) + "-" + std::to_string(j);
        devmap[device_id] = new CommDevice(link_name, CommDevice::NW_COMM, 
          -1, -1, device_id, 0, conn[i * total_devs + j]);
      }
    // }
  }
  size_t total_hops = 0;
  int nlinks = 0;
  TwoDimTorusW2TURNRouting w2t{conn, devmap, gen.nrows, gen.ncols, gen.nextras, total_devs};
  std::map<size_t, EcmpRoutes> routes; 
  for (int i = 0; i < total_devs; i++) {
    for (int j = 0; j < total_devs; j++) {
      routes[i * total_devs + j] = w2t.get_routes(i, j);
      fprintf(stderr, "%d-%d::", i, j);
      if (routes[i * total_devs + j].second.size() > 0) {
        nlinks++;
        total_hops += routes[i * total_devs + j].second[0].size();
        for (CommDevice * d: routes[i * total_devs + j].second[0]) {
          fprintf(stderr, "%d-%d;", d->device_id/total_devs, d->device_id%total_devs);
        } 
      }
      fprintf(stderr, "\n");
    }
  }
  fprintf(stderr, "avg hop count: %f\n", (double)total_hops/nlinks);
}