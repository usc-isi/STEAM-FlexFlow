#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <utility>
#include <cmath>
#include <unordered_set>
#include <functional>
#include <numeric>

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

// all prime numbers below 2048. good enouggh.
const static uint16_t PRIMES[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039};

// for summing connections...
template <typename T>
static std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), 
                  std::back_inserter(result), std::plus<T>());
  return result;
}

WeightedShortestPathRoutingStrategy::WeightedShortestPathRoutingStrategy(
    const ConnectionMatrix & c, 
    const std::map<size_t, CommDevice*>& devmap,
    int total_devs) 
: conn(c), devmap(devmap), total_devs(total_devs)
{} 

EcmpRoutes WeightedShortestPathRoutingStrategy::get_routes(int src_node, int dst_node) 
{
  int key = src_node * total_devs + dst_node;

  if (conn[key] > 0) {
    return std::make_pair(std::vector<double>({1}), std::vector<Route>({Route({devmap.at(key)})}));
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
      if (new_dist < dist[i] || (new_dist == dist[i] && unif(gen) < 0.5)) {
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
  return std::make_pair(std::vector<double>{1}, std::vector<Route>{result});
}

void WeightedShortestPathRoutingStrategy::hop_count(int src_node, int dst_node, int & hop, int & narrowest)
{
  int key = src_node * total_devs + dst_node;

  if (conn[key] > 0) {
    hop = 0;
    narrowest = conn[key];
    return;
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
  hop = 0;
  narrowest = std::numeric_limits<int>::max();
  int curr = dst_node;
  while (prev[curr] != -1) {
    if (narrowest > conn[prev[curr] * total_devs + curr]) narrowest = conn[prev[curr] * total_devs + curr];
    hop++;
    curr = prev[curr];
  }
  assert(hop > 0 || src_node == dst_node);
}


std::vector<EcmpRoutes> WeightedShortestPathRoutingStrategy::get_routes_from_src(int src_node) 
{
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::priority_queue<std::pair<uint64_t, uint64_t>, 
                      std::vector<std::pair<uint64_t, uint64_t> >,
                      std::greater<std::pair<uint64_t, uint64_t> > > pq;
  pq.push(std::make_pair(dist[src_node], src_node));
  dist[src_node] = 0;
  while (!pq.empty()) {
    int min_node = pq.top().second;
    pq.pop();
    visited[min_node] = true;

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
  std::vector<EcmpRoutes> final_result;
  for (int i = 0; i < total_devs; i++) {
    if (i == src_node) {
      final_result.emplace_back(std::make_pair(std::vector<double>{}, std::vector<Route>{}));
      continue;
    }
    Route result = Route();
    int curr = i;
    while (prev[curr] != -1) {
      result.insert(result.begin(), devmap.at(prev[curr] * total_devs + curr));
      curr = prev[curr];
    }
    assert(result.size() > 0);
    final_result.emplace_back(std::make_pair(std::vector<double>{1}, std::vector<Route>{result}));
  }
  return final_result; 
}

std::vector<std::pair<int, int>> WeightedShortestPathRoutingStrategy::hop_count(int src_node) 
{
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::priority_queue<std::pair<uint64_t, uint64_t>, 
                      std::vector<std::pair<uint64_t, uint64_t> >,
                      std::greater<std::pair<uint64_t, uint64_t> > > pq;
  pq.push(std::make_pair(dist[src_node], src_node));
  dist[src_node] = 0;
  while (!pq.empty()) {
    int min_node = pq.top().second;
    pq.pop();
    visited[min_node] = true;

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

  std::vector<std::pair<int, int>> result;
  for (int i = 0; i < total_devs; i++) {
    if (i == src_node) {
      result.emplace_back(std::make_pair(-1, 0));
      continue;
    }
    int hop = -1;
    int narrowest = 0;
    int curr = i;
    while (prev[curr] != -1) {
      if (!narrowest || (narrowest > conn[prev[curr] * total_devs + curr])) 
        narrowest = conn[prev[curr] * total_devs + curr];
      hop++;
      curr = prev[curr];
    }
    result.emplace_back(std::make_pair(hop, narrowest));
  }
  return result; 
}

ShortestPathNetworkRoutingStrategy::ShortestPathNetworkRoutingStrategy(
    const ConnectionMatrix & c, 
    const std::map<size_t, CommDevice*>& devmap,
    int total_devs) 
: conn(c), devmap(devmap), total_devs(total_devs)
{} 

EcmpRoutes ShortestPathNetworkRoutingStrategy::get_routes(int src_node, int dst_node) 
{
  int key = src_node * total_devs + dst_node;
  // std::cerr << "routing " << src_node << ", " << dst_node << std::endl;

  if (conn[key] > 0) {
    return std::make_pair(std::vector<double>({1}), std::vector<Route>({Route({devmap.at(key)})}));
  }

  // one-shortest path routing
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::queue<uint64_t> q;
  q.push(src_node);
  dist[src_node] = 0;

  std::vector<int> rd_idx(total_devs);
  std::iota(std::begin(rd_idx), std::end(rd_idx), 0);

  // BFS
  while (!q.empty()) {
    int min_node = q.front();
    q.pop();
    visited[min_node] = true;

    if (min_node == dst_node)
      break;

    std::random_shuffle(rd_idx.begin(), rd_idx.end());
    for (int i: rd_idx) {
      if (visited[i] || conn[min_node * total_devs + i] == 0) {
        continue;
      }
      double new_dist = dist[min_node] + 1; 
      if (new_dist < dist[i] || (new_dist == dist[i] && unif(gen) < 0.5)) {
        dist[i] = new_dist;
        prev[i] = min_node;
        q.push(i);
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
  return std::make_pair(std::vector<double>{1}, std::vector<Route>{result});
}

std::vector<EcmpRoutes> ShortestPathNetworkRoutingStrategy::get_routes_from_src(int src_node) 
{
  
  std::vector<EcmpRoutes> final_result;
  for (int i = 0; i < total_devs; i++) {
    final_result.emplace_back(get_routes(src_node, i));
  }
  /*
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::queue<uint64_t> q;
  q.push(src_node);
  dist[src_node] = 0;

  std::vector<int> rd_idx(total_devs);
  std::iota(std::begin(rd_idx), std::end(rd_idx), 0);
  // BFS
  while (!q.empty()) {
    int min_node = q.front();
    q.pop();
    visited[min_node] = true;

    std::random_shuffle(rd_idx.begin(), rd_idx.end());
    for (int i: rd_idx) {
      if (visited[i] || conn[min_node * total_devs + i] == 0) {
        continue;
      }
      double new_dist = dist[min_node] + 1; 
      if (new_dist < dist[i] || (new_dist == dist[i] && unif(gen) < 0.5)) {
        dist[i] = new_dist;
        prev[i] = min_node;
        q.push(i);
      }
    }
  }

  std::vector<EcmpRoutes> final_result;
  for (int i = 0; i < total_devs; i++) {
    if (i == src_node) {
      final_result.emplace_back(std::make_pair(std::vector<double>{}, std::vector<Route>{}));
      continue;
    }
    Route result = Route();
    int curr = i;
    while (prev[curr] != -1) {
      result.insert(result.begin(), devmap.at(prev[curr] * total_devs + curr));
      curr = prev[curr];
    }
    // assert(result.size() > 0);
    final_result.emplace_back(std::make_pair(std::vector<double>{1}, std::vector<Route>{result}));
  }
  */
  return final_result; 
}

void ShortestPathNetworkRoutingStrategy::hop_count(int src_node, int dst_node, int & hop, int & narrowest)
{
  int key = src_node * total_devs + dst_node;

  if (conn[key] > 0) {
    hop = 0;
    narrowest = conn[key];
    return;
  }
  // one-shortest path routing
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::queue<uint64_t> q;
  q.push(src_node);
  dist[src_node] = 0;

  // BFS
  while (!q.empty()) {
    int min_node = q.front();
    q.pop();
    visited[min_node] = true;

    if (min_node == dst_node)
      break;

    for (int i = 0; i < total_devs; i++) {
      if (visited[i] || conn[min_node * total_devs + i] == 0) {
        continue;
      }
      double new_dist = dist[min_node] + 1; 
      if (new_dist < dist[i] || (new_dist == dist[i] && unif(gen) < 0.5)) {
        dist[i] = new_dist;
        prev[i] = min_node;
        q.push(i);
      }
    }
  }
  hop = 0;
  narrowest = std::numeric_limits<int>::max();
  int curr = dst_node;
  while (prev[curr] != -1) {
    if (narrowest > conn[prev[curr] * total_devs + curr]) narrowest = conn[prev[curr] * total_devs + curr];
    hop++;
    curr = prev[curr];
  }
  assert(hop > 0 || src_node == dst_node);
}

std::vector<std::pair<int, int>> ShortestPathNetworkRoutingStrategy::hop_count(int src_node) 
{
  std::vector<uint64_t> dist(total_devs, std::numeric_limits<uint64_t>::max());
  std::vector<int> prev(total_devs, -1);
  std::vector<bool> visited(total_devs, false);

  std::queue<uint64_t> q;
  q.push(src_node);
  dist[src_node] = 0;

  // BFS
  while (!q.empty()) {
    int min_node = q.front();
    q.pop();
    visited[min_node] = true;

    for (int i = 0; i < total_devs; i++) {
      if (visited[i] || conn[min_node * total_devs + i] == 0) {
        continue;
      }
      double new_dist = dist[min_node] + 1; 
      if (new_dist < dist[i] || (new_dist == dist[i] && unif(gen) < 0.5)) {
        dist[i] = new_dist;
        prev[i] = min_node;
        q.push(i);
      }
    }
  }

  std::vector<std::pair<int, int>> result;
  for (int i = 0; i < total_devs; i++) {
    if (i == src_node) {
      result.emplace_back(std::make_pair(-1, 0));
      continue;
    }
    int hop = -1;
    int narrowest = 0;
    int curr = i;
    while (prev[curr] != -1) {
      if (!narrowest || (narrowest > conn[prev[curr] * total_devs + curr])) 
        narrowest = conn[prev[curr] * total_devs + curr];
      hop++;
      curr = prev[curr];
    }
    result.emplace_back(std::make_pair(hop, narrowest));
  }
  return result; 
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
  alpha = 0.5;
  no_improvement_th = 50;
  best_sim_time = std::numeric_limits<double>::max();
  curr_sim_time = std::numeric_limits<double>::max();
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

size_t DemandHeuristicNetworkOptimizer::edge_id(int i, int j) const
{
  return i * machine->get_total_devs() + j;
}

size_t DemandHeuristicNetworkOptimizer::unordered_edge_id(int i, int j) const
{
  return i > j ? edge_id(i, j) : edge_id(j, i);
}

typedef std::pair<uint64_t, uint64_t> DemandToIdMap;

bool DemandHeuristicNetworkOptimizer::optimize(int mcmc_iter, double sim_iter_time, bool forced)
{
  double diff = sim_iter_time - curr_sim_time;
  std::cerr << "sim_iter_time: " << sim_iter_time << ", curr_sim_time: " << curr_sim_time 
            << ", best_iter_time: " << best_sim_time << std::endl;
  bool change = diff < 0 ? true : diff != 0 && unif(gen) < std::exp(-alpha * diff);
  if (sim_iter_time < best_sim_time) {
    best_sim_time = sim_iter_time;
    change = true;
  }
  if (change) {
    curr_sim_time = sim_iter_time;
  }
  else {
    num_iter_nochange++; 
  }

  if (!forced && !change && num_iter_nochange < no_improvement_th)
    return false;
  
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
  return true;
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

  // best_sim_time = std::numeric_limits<double>::max();
  // curr_sim_time = std::numeric_limits<double>::max();
  // num_iter_nochange = 0;
  
}

std::unique_ptr<L1OptimizerInformation> DemandHeuristicNetworkOptimizer::export_information()
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  return std::unique_ptr<L1OptimizerInformation>(new L1TopologyInformation{nm->conn_matrix});
}

void DemandHeuristicNetworkOptimizer::import_information(const std::unique_ptr<L1OptimizerInformation>& information) 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  nm->set_topology(static_cast<L1TopologyInformation*>(information.get())->conn);
  nm->update_route();
}

void DemandHeuristicNetworkOptimizer::delete_information(const std::unique_ptr<L1OptimizerInformation>& information) 
{
  // information.res
  // delete (ConnectionMatrix*)information;
}

void DemandHeuristicNetworkOptimizer::store_tm() const 
{

  std::ofstream output; 
  output.open("traffic_matrix.txt");
  output << "PHYSICAL_TM:" << std::endl;
  for (int i = 0; i < machine->get_total_devs(); i++) {
    for (int j = 0; j < machine->get_total_devs(); j++) {
      if (physical_traffic_demand.find(edge_id(i,j)) == physical_traffic_demand.end())
        output << "0" << "\t";
      else
        output << physical_traffic_demand.at(edge_id(i, j)) << "\t";
    }
    output << std::endl;
  }
  output << std::endl;
  output << "LOGICAL_TM:" << std::endl;
  for (int i = 0; i < machine->get_total_devs(); i++) {
    for (int j = 0; j < machine->get_total_devs(); j++) {
      if (logical_traffic_demand.find(edge_id(i,j)) == logical_traffic_demand.end())
        output << "0" << "\t";
      else
        output << logical_traffic_demand.at(edge_id(i, j)) << "\t";
    }
    output << std::endl;
  }
  output << std::endl;
  output << "INDIR_LIST:" << std::endl;
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  const ConnectionMatrix & conn = nm->conn_matrix;
  ShortestPathNetworkRoutingStrategy s{conn, nm->ids_to_nw_comm_device, nm->total_devs};
  std::unordered_map<size_t, std::pair<uint64_t, double>> result;
  for (auto& entry: logical_traffic_demand) {
    if (conn[entry.first] == 0 && entry.second > 0) {
      int hop_cnt, narrowest;
      s.hop_count(entry.first / nm->total_devs, entry.first % nm->total_devs, hop_cnt, narrowest);
      output << "(" << entry.first / nm->total_devs << ", " << entry.first % nm->total_devs 
             << "): " << entry.second << ", " << hop_cnt << ", " << narrowest << std::endl;
    }
  }
  output << std::endl;
}

DemandHeuristicNetworkOptimizerPlus::DemandHeuristicNetworkOptimizerPlus(MachineModel* machine)
: DemandHeuristicNetworkOptimizer(machine)
{}

void DemandHeuristicNetworkOptimizerPlus::connectivity_assign(
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
    
    // first stage: never assign parallel links across nodes
    if (conn[edge_id(node0, node1)] == 2) {
      for (auto it = pq.begin(); it != pq.end(); ) {
        if (DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node0, ndevs)) {
          // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else if (DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node1, ndevs)) {
          // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else {
          ++it;
        }
      }
      continue;
    }
    conn[edge_id(node0, node1)]++;
    conn[edge_id(node1, node0)]++;

    INSERT_OR_ADD(node_if_allocated, node0, 1);
    INSERT_OR_ADD(node_if_allocated, node1, 1);

    // std::cout << "first is " << target.first << std::endl;
    target.first /= (conn[target.second] + 1); //*= (double)conn[target.second]/(conn[target.second] + 1);
    if (target.first > 0) {
      pq.insert(target);
    }

    if (node_if_allocated[node0] == if_cnt/2 || node_if_allocated[node1] == if_cnt/2) {
      for (auto it = pq.begin(); it != pq.end(); ) {
        if (node_if_allocated[node0] == if_cnt/2 && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node0, ndevs)) {
          // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else if (node_if_allocated[node1] == if_cnt/2 && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node1, ndevs)) {
          // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = pq.erase(it);
        }
        else {
          ++it;
        }
      } 
    }
  }

  // for (int i = 0; i < nnode; i++) {
  //   conn[edge_id(i, (i + nnode / 5) % nnode)]++;
  //   conn[edge_id(i, 2 * (i + nnode / 5) % nnode)]++;
  //   conn[edge_id(i, 3 * (i + nnode / 5) % nnode)]++;
  //   conn[edge_id(i, 4 * (i + nnode / 5) % nnode)]++;
  // }
}

void DemandHeuristicNetworkOptimizerPlus::connect_topology(
        // const std::unordered_map<uint64_t, uint64_t> & logical_id_to_demand, 
        ConnectionMatrix &conn, 
        std::unordered_map<size_t, size_t> & node_if_allocated)
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
  while (n_cc > 1) {
    int64_t max_demand = std::numeric_limits<int64_t>::min();
    int node_from_cc0 = -1;
    int node_from_cc1 = -1;
    for (size_t n0: ccs[0]) {
      for (size_t n1: ccs[1]) {
        int64_t demand;
        if (logical_traffic_demand.find(edge_id(n0, n1)) == logical_traffic_demand.end())
          demand = 0;
        else 
          demand = logical_traffic_demand[edge_id(n0, n1)];
        if (demand > max_demand || (demand == max_demand && unif(gen) > 0.5)) {
          if (node_if_allocated[n0] == if_cnt || node_if_allocated[n1] == if_cnt) {
            continue;
          }
          node_from_cc0 = n0;
          node_from_cc1 = n1;
          max_demand = demand;
        }
      }
    }
    bool success = add_link(node_from_cc0, node_from_cc1, conn);
    node_if_allocated[node_from_cc0]++;
    node_if_allocated[node_from_cc1]++;
    assert(success);
    n_cc--;
    ccs[1].insert(ccs[0].begin(), ccs[0].end());
    ccs.erase(ccs.begin()); 
  }

}

void DemandHeuristicNetworkOptimizerPlus::utility_max_assign(
    ConnectionMatrix &conn,
    // const std::unordered_map<size_t, uint64_t> &max_of_bidir,
    std::unordered_map<size_t, size_t> &node_if_allocated)
{
  // link, <demand, discounted hopcount> 
  // std::unordered_map<size_t, std::pair<uint64_t, double>> indirect_traffic = construct_indir_traffic_list(conn);
  std::unordered_map<size_t, uint64_t> indirect_traffic = construct_bidir_negative_util(conn);
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  std::unordered_map<size_t, uint64_t> sum_of_bidir;
  for (int i = 0; i < nnode; i++) {
    for (int j = 0; j < nnode; j++) {
      size_t eid = edge_id(i, j);
      if (logical_traffic_demand.find(eid) != logical_traffic_demand.end()) {
        size_t ueid = unordered_edge_id(i, j);
        uint64_t traffic_amount = logical_traffic_demand[eid];
        INSERT_OR_ADD(sum_of_bidir, ueid, traffic_amount);
      }
    }
  }
  double utility = compute_utility(sum_of_bidir, indirect_traffic, conn);

  std::set<DemandToIdMap, std::greater<DemandToIdMap>> positive_pq;
  for (auto &item: sum_of_bidir) {
    // mod: pre-unscale the demand
    positive_pq.insert(DemandToIdMap(item.second/(1 << conn[item.first]), item.first));
  }

  std::set<DemandToIdMap, std::greater<DemandToIdMap>> negative_pq;
  for (auto &item: indirect_traffic) {
    negative_pq.insert(DemandToIdMap(item.second, item.first));
  }

  // clean up dead node
  for (auto & entry: node_if_allocated) {
    if (entry.second == if_cnt) {
      for (auto it = positive_pq.begin(); it != positive_pq.end(); ) {
        if (node_if_allocated[entry.first] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, entry.first, ndevs)) {
          // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = positive_pq.erase(it);
        }
        else {
          ++it;
        }
      }
      for (auto it = negative_pq.begin(); it != negative_pq.end(); ) {
        if (node_if_allocated[entry.first] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, entry.first, ndevs)) {
          // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
          it = negative_pq.erase(it);
        }
        else {
          ++it;
        }
      } 
    }
  }

  while (!positive_pq.empty() || !negative_pq.empty()) {
    if (positive_pq.empty()) {
      DemandToIdMap target = *negative_pq.begin();
      negative_pq.erase(negative_pq.begin());

      size_t node0 = target.second / ndevs;
      size_t node1 = target.second % ndevs;
      
      // conn[target.second]++;
      conn[edge_id(node0, node1)]++;
      conn[edge_id(node1, node0)]++;

      INSERT_OR_ADD(node_if_allocated, node0, 1);
      INSERT_OR_ADD(node_if_allocated, node1, 1);

      // std::cout << "first is " << target.first << std::endl;
      indirect_traffic = construct_bidir_negative_util(conn);
      negative_pq.clear();
      for (auto &item: indirect_traffic) {
        if (node_if_allocated[item.first/nnode] < if_cnt && node_if_allocated[item.first%nnode] < if_cnt)
          negative_pq.insert(DemandToIdMap(item.second, item.first));
      }
    }
    else if (negative_pq.empty()) {
      DemandToIdMap target = *positive_pq.begin();
      positive_pq.erase(positive_pq.begin());

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
        positive_pq.insert(target);
      }

      if (node_if_allocated[node0] == if_cnt || node_if_allocated[node1] == if_cnt) {
        for (auto it = positive_pq.begin(); it != positive_pq.end(); ) {
          if (node_if_allocated[node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node0, ndevs)) {
            // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
            it = positive_pq.erase(it);
          }
          else if (node_if_allocated[node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, node1, ndevs)) {
            // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
            it = positive_pq.erase(it);
          }
          else {
            ++it;
          }
        } 
      } 
    }
    else {
      DemandToIdMap n_target = *negative_pq.begin();
      ConnectionMatrix nconn = conn;
      size_t n_node0 = n_target.second / ndevs;
      size_t n_node1 = n_target.second % ndevs;
      nconn[edge_id(n_node0, n_node1)]++;
      nconn[edge_id(n_node1, n_node0)]++;
      // auto proposed_indir = construct_bidir_negative_util(nconn);
      double n_util = compute_utility(sum_of_bidir, indirect_traffic, nconn);
      std::unordered_map<size_t, uint64_t> proposed_indir;
      if (n_util < 0) {
        // proposed_indir = construct_bidir_negative_util(nconn);
        // n_util = compute_utility(sum_of_bidir, proposed_indir, nconn);
      }
      DemandToIdMap p_target = *positive_pq.begin();
      ConnectionMatrix pconn = conn;
      size_t p_node0 = p_target.second / ndevs;
      size_t p_node1 = p_target.second % ndevs;
      pconn[edge_id(p_node0, p_node1)]++;
      pconn[edge_id(p_node1, p_node0)]++;
      double p_util = compute_utility(sum_of_bidir, indirect_traffic, pconn);

      if (n_util > p_util || n_util < 0) {
        // use the negative link
        conn[edge_id(n_node0, n_node1)]++;
        conn[edge_id(n_node1, n_node0)]++;
        INSERT_OR_ADD(node_if_allocated, n_node0, 1);
        INSERT_OR_ADD(node_if_allocated, n_node1, 1);
        // std::cout << "first is " << target.first << std::endl;
        if (proposed_indir.empty()) {
          indirect_traffic.erase(negative_pq.begin()->second);
          negative_pq.erase(negative_pq.begin());
        }
        else {
          indirect_traffic = proposed_indir;
          negative_pq.clear();
          for (auto &item: indirect_traffic) {
            if (node_if_allocated[item.first/nnode] < if_cnt && node_if_allocated[item.first%nnode] < if_cnt)
              negative_pq.insert(DemandToIdMap(item.second, item.first));
          }
        }
        if (node_if_allocated[n_node0] == if_cnt || node_if_allocated[n_node1] == if_cnt) {
          for (auto it = positive_pq.begin(); it != positive_pq.end(); ) {
            if (node_if_allocated[n_node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, n_node0, ndevs)) {
              // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = positive_pq.erase(it);
            }
            else if (node_if_allocated[n_node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, n_node1, ndevs)) {
              // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = positive_pq.erase(it);
            }
            else {
              ++it;
            }
          }
          for (auto it = negative_pq.begin(); it != negative_pq.end(); ) {
            if (node_if_allocated[n_node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, n_node0, ndevs)) {
              // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = negative_pq.erase(it);
            }
            else if (node_if_allocated[n_node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, n_node1, ndevs)) {
              // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = negative_pq.erase(it);
            }
            else {
              ++it;
            }
          } 
        } 
      }
      else {
        positive_pq.erase(positive_pq.begin());
        conn[edge_id(p_node0, p_node1)]++;
        conn[edge_id(p_node1, p_node0)]++;

        INSERT_OR_ADD(node_if_allocated, p_node0, 1);
        INSERT_OR_ADD(node_if_allocated, p_node1, 1);

        // std::cout << "first is " << target.first << std::endl;
        p_target.first /= 2; //*= (double)conn[target.second]/(conn[target.second] + 1);
        if (p_target.first > 0) {
          positive_pq.insert(p_target);
        }
        if (node_if_allocated[p_node0] == if_cnt || node_if_allocated[p_node1] == if_cnt) {
          for (auto it = positive_pq.begin(); it != positive_pq.end(); ) {
            if (node_if_allocated[p_node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, p_node0, ndevs)) {
              // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = positive_pq.erase(it);
            }
            else if (node_if_allocated[p_node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, p_node1, ndevs)) {
              // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = positive_pq.erase(it);
            }
            else {
              ++it;
            }
          }
          for (auto it = negative_pq.begin(); it != negative_pq.end(); ) {
            if (node_if_allocated[p_node0] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, p_node0, ndevs)) {
              // std::cout << "node0 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = negative_pq.erase(it);
            }
            else if (node_if_allocated[p_node1] == if_cnt && DemandHeuristicNetworkOptimizer::has_endpoint(it->second, p_node1, ndevs)) {
              // std::cout << "node1 full, removing " << it->second /ndevs << ", " << it->second % ndevs << " with demand left " << it->first << std::endl; 
              it = negative_pq.erase(it);
            }
            else {
              ++it;
            }
          } 
        } 
      }
    }
  }
}

static double N_POWER2_MULFACTOR_LOOKUP[] = {0.000000000000000000000000000000,1.000000000000000000000000000000,1.500000000000000000000000000000,1.750000000000000000000000000000,1.875000000000000000000000000000,1.937500000000000000000000000000,1.968750000000000000000000000000,1.984375000000000000000000000000,1.992187500000000000000000000000,1.996093750000000000000000000000,1.998046875000000000000000000000,1.999023437500000000000000000000,1.999511718750000000000000000000,1.999755859375000000000000000000,1.999877929687500000000000000000,1.999938964843750000000000000000,1.999969482421875000000000000000,1.999984741210937500000000000000,1.999992370605468750000000000000,1.999996185302734375000000000000,1.999998092651367187500000000000,1.999999046325683593750000000000,1.999999523162841796875000000000,1.999999761581420898437500000000,1.999999880790710449218750000000,1.999999940395355224609375000000,1.999999970197677612304687500000,1.999999985098838806152343750000,1.999999992549419403076171875000,1.999999996274709701538085937500,1.999999998137354850769042968750,1.999999999068677425384521484375};

double DemandHeuristicNetworkOptimizerPlus::compute_utility(
  const std::unordered_map<size_t, std::pair<uint64_t, double>> &indirect_traffic,
  const ConnectionMatrix & conn)
{
  double result = 0;
  for (auto &entry: logical_traffic_demand) {
    if (entry.second > 0) {
      if (conn[entry.first] > 0) {
        result += entry.second * N_POWER2_MULFACTOR_LOOKUP[conn[entry.first]];
      }
      else {
        assert(indirect_traffic.find(entry.first) != indirect_traffic.end());
        result -= entry.second * indirect_traffic.at(entry.first).second;
      }
    }
  }
  return result;
}

double DemandHeuristicNetworkOptimizerPlus::compute_utility(
  const std::unordered_map<size_t,uint64_t> &sum_of_bidir,
  const std::unordered_map<size_t,uint64_t> &indirect_traffic,
  const ConnectionMatrix & conn)
{
  double result = 0;
  for (auto &entry: sum_of_bidir) {
    if (entry.second > 0) {
      if (conn[entry.first] > 0) {
        result += entry.second * N_POWER2_MULFACTOR_LOOKUP[conn[entry.first]];
      }
      else {
        assert(indirect_traffic.find(entry.first) != indirect_traffic.end());
        result -= indirect_traffic.at(entry.first);
      }
    }
  }
  std::cerr << "utility: " << result << std::endl;
  return result;
}

std::unordered_map<size_t, std::pair<uint64_t, double>>
DemandHeuristicNetworkOptimizerPlus::construct_indir_traffic_list(const ConnectionMatrix &conn) 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  ShortestPathNetworkRoutingStrategy s{conn, nm->ids_to_nw_comm_device, nm->total_devs};
  std::unordered_map<size_t, std::pair<uint64_t, double>> result;
  for (auto& entry: logical_traffic_demand) {
    if (conn[entry.first] == 0 && entry.second > 0) {
      int hop_cnt, narrowest;
      s.hop_count(entry.first / nm->total_devs, entry.first % nm->total_devs, hop_cnt, narrowest);
      double discounted_hop = (double) hop_cnt / narrowest;
      result[entry.first] = std::make_pair(entry.second, discounted_hop);
    }
  }
  return result;
}

std::unordered_map<size_t, size_t> 
DemandHeuristicNetworkOptimizerPlus::construct_bidir_negative_util(const ConnectionMatrix &conn)
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  ShortestPathNetworkRoutingStrategy s{conn, nm->ids_to_nw_comm_device, nm->total_devs};
  std::unordered_map<size_t, size_t> result;
  
  std::unordered_map<size_t, uint64_t> sum_of_bidir;
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  for (int i = 0; i < nnode; i++) {
    for (int j = 0; j < nnode; j++) {
      size_t eid = edge_id(i, j);
      if (logical_traffic_demand.find(eid) != logical_traffic_demand.end() && conn[eid] == 0) {
        size_t ueid = unordered_edge_id(i, j);
        uint64_t traffic_amount = logical_traffic_demand[eid];
        int hop_cnt, narrowest;
        s.hop_count(i, j, hop_cnt, narrowest);
        double discounted_hop = (double) hop_cnt / narrowest;
        // std::cerr << "from " << i << " to " << j << " discounted: " << discounted_hop 
        //   << "(hop cnt: " << hop_cnt << ", narrowest: " << narrowest << std::endl;
        INSERT_OR_ADD(result, ueid, traffic_amount * discounted_hop);
      }
    }
  }
  return result;
}

bool DemandHeuristicNetworkOptimizerPlus::optimize(int mcmc_iter, double sim_iter_time, bool forced)
{
  double diff = sim_iter_time - curr_sim_time;
  std::cerr << "sim_iter_time: " << sim_iter_time << ", curr_sim_time: " << curr_sim_time 
            << ", best_iter_time: " << best_sim_time << std::endl;
  bool change = diff < 0 ? true : diff != 0 && unif(gen) < std::exp(-alpha * diff);
  if (sim_iter_time < best_sim_time) {
    best_sim_time = sim_iter_time;
    change = true;
  }
  if (change) {
    curr_sim_time = sim_iter_time;
  }
  else {
    num_iter_nochange++; 
  }

  if (!forced && !change && num_iter_nochange < no_improvement_th)
    return false;
  
  num_iter_nochange = 0;
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  
  // TODO: copy machine-switch link?
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  ConnectionMatrix conn = std::vector<int>(ndevs*ndevs, 0);
  std::unordered_map<size_t, uint64_t> max_of_bidir;
  std::unordered_map<size_t, size_t> node_if_allocated;

  connectivity_assign(conn, max_of_bidir, node_if_allocated);
#ifdef DEBUG_PRINT
  std::cerr << "After connectivity_assign " << std::endl;
  NetworkTopologyGenerator::print_conn_matrix(conn, nnode, 0);
#endif
  // connect_topology(conn, node_if_allocated);
  // std::cerr << "After conn_topo " << std::endl;
  // NetworkTopologyGenerator::print_conn_matrix(conn, nnode, 0);
  // utility_max_assign(conn, node_if_allocated);
  // std::cerr << "After util_max " << std::endl;
  // NetworkTopologyGenerator::print_conn_matrix(conn, nnode, 0);

  // connect_unused_node(conn, node_if_allocated);

  // Make all CC connected
  // std::unordered_map<uint64_t, uint64_t> logical_id_to_demand;
  // for (auto & item: max_of_bidir) {
  //   logical_id_to_demand[item.second] = item.first;
  // }

  // connect_cc(logical_id_to_demand, conn); 
  nm->set_topology(conn); 
  nm->update_route();
  // simulator->print_conn_matrix();
  return true;

}

SpMulMat::SpMulMat(MachineModel * machine, int degree, bool bidir)
: DemandHeuristicNetworkOptimizer(machine), bidir(bidir)
{
  if_cnt = degree;
  constructed = false;
  construct_candidate_jumps();
}

void SpMulMat::task_added(SimTask * task) 
{
  NominalCommDevice *ncommDev;
  DPGroup dpg;
  switch (task->type) {
  
  case SimTask::TASK_COMM:
  case SimTask::TASK_BACKWARD: 
  case SimTask::TASK_FORWARD:
  break;

  // in this case this is guaranteed to be a MP traffic
  case SimTask::TASK_NOMINAL_COMM:
    ncommDev = reinterpret_cast<NominalCommDevice*>(task->device);
    INSERT_OR_ADD(mp_tm_logical, ncommDev->device_id, task->xfer_size);
  break;

  case SimTask::TASK_ALLREDUCE:
    dpg.group_size = task->next_tasks.size();
    dpg.starting_node = reinterpret_cast<uint64_t>(task->next_tasks[0]);
    dpg.xfer_size = task->xfer_size;
    dpgrps.emplace_back(dpg);
    INSERT_OR_ADD(dpgrpsz_xfersize, dpg.group_size, (double)task->xfer_size * (2.0 * (double)dpg.group_size - 1) / (double)dpg.group_size);
  break;

  case SimTask::TASK_BARRIER:

  break;

  case SimTask::TASK_UPDATE:

  break;
  }
}

// uint64_t SpMulMat::dpgrp_unique_key(const DPGroup & dpg) const
// {
//   return (((uint64_t)dpg.group_size) << 32) | (dpg.starting_node);
// }

// int SpMulMat::get_start_node(uint64_t id) const 
// {
//   return (int)(id & 0x00000000FFFFFFFFULL);
// }

// int SpMulMat::get_group_size(uint64_t id) const 
// {
//   return (int)((id & 0xFFFFFFFF00000000ULL) >> 32);
// }

std::vector<int> SpMulMat::negative(const std::vector<int>& v)
{
  std::vector<int> result{v};
  std::transform(result.cbegin(),result.cend(),result.begin(),
    [&](int x)->int{return machine->get_num_nodes() - x;});
  return result;
}

bool SpMulMat::segment_overlap(const std::vector<int>& a, const std::vector<int>& b) 
{
  assert(a.size() == b.size()); 
  for (int i = 0; i < a.size(); i++) {
    if (a.at(i) == b.at(i)) {
      return true;
    }
  } 
  return true;
}

std::vector<int> SpMulMat::choose_n(const std::vector<int>& cjs, int jmp, int n)
{
  if (n == 1) {
    return {jmp};
  }
  std::set<int> result{};
  for (int i = 0; i < n; i++) {
    if (jmp >= cjs.back()) {
      result.insert(cjs.front());
    }
    else {
      result.insert(*std::upper_bound(cjs.cbegin(), cjs.cend(), jmp));
    }
    // std::cerr << "Adding " << *std::upper_bound(cjs.cbegin(), cjs.cend(), jmp) 
    //           << ", jmp = " << jmp << std::endl;
    int next_jmp = jmp + machine->get_num_nodes()/n;
    jmp = MOD(next_jmp, machine->get_num_nodes());
  }
  return std::vector<int>{result.cbegin(), result.cend()};
}

std::vector<int> SpMulMat::choose_n_geo(const std::vector<int>& cjs, int n)
{
  if (n == 1) {
    return {cjs[0]};
  }
  double diff = (double)cjs[cjs.size() - 1] / cjs[0];
  double ratio = std::pow(diff, 1.0/(n-(bidir?0:1)));
  // std::cerr << "n: " << n << ",ratio: " << ratio << std::endl;

  std::set<int> result{};
  double curr = cjs[0];
  for (int i = 0; i < n; i++) {
    assert(curr <= cjs[cjs.size() - 1]);
    auto candidate = std::lower_bound(cjs.begin(), cjs.end(), curr);
    while (result.find(*candidate) != result.end()) {
      candidate++;
    }
    result.insert(*candidate);
    // std::cerr << "Adding " << *candidate << std::endl;
    curr *= ratio;
    // std::cerr << "curr " << curr << std::endl;
  }
  return std::vector<int>{result.cbegin(), result.cend()};
}

void SpMulMat::construct_candidate_jumps() {
  int total_nodes = machine->get_num_nodes();
  for (int i = 2; i <= total_nodes; i++) {
    if (total_nodes % i == 0) {
      int group_size = i;
      int base_jump = total_nodes / i;
      int nconfs = total_nodes % i;
      candidate_jumps.emplace(std::make_pair(group_size, std::vector<int>{}));
      for (int k = 1; ; k++) {
        if (k > group_size) {
          break;
        }
        if (Realm::gcd(group_size, k) == 1)
          candidate_jumps[group_size].emplace_back(base_jump * k);
      }
    }
  }
}

void SpMulMat::get_dp_mp_degree(int & dp_degree, int & mp_degree)
{

  // mp_degree = if_cnt - 2;
  // dp_degree = 2;
  // return; 

  double total_dp_traffic = 0;
  double total_mp_traffic = 0; 

  for (auto & entry: dpgrps) { 
    // ring allreduce
    total_dp_traffic += entry.xfer_size * (2 * entry.group_size - 1);
  }

  for (auto & entry: mp_tm_logical) {
    total_mp_traffic += entry.second;
  }

  dp_degree = lround(total_dp_traffic / (total_dp_traffic + total_mp_traffic) * if_cnt);
  mp_degree = lround(total_mp_traffic / (total_dp_traffic + total_mp_traffic) * if_cnt);

  if (bidir && dp_degree % 2 != 0) {
    mp_degree--;
    dp_degree++;
  }
  if (dp_degree == 0) {
    if (bidir) {
      dp_degree += 2;
      mp_degree -= 2;
    }
    else {
      dp_degree++;
      mp_degree--;
    }
    
  }
  assert(mp_degree + dp_degree == if_cnt);
#ifdef DEBUG_PRINT
  std::cerr << "dp deg: " << dp_degree << std::endl;
  std::cerr << "mp deg: " << mp_degree << std::endl;
#endif

}

std::vector<std::pair<uint64_t, int>> SpMulMat::generate_dp_topology(ConnectionMatrix & conn, int dp_degree) 
{
  size_t ndevs = machine->get_num_nodes();
  // reset conn
  conn = ConnectionMatrix(ndevs * ndevs, 0);

  // sort out the DP groups by weight
  std::set<std::pair<uint64_t, uint64_t>, std::greater<std::pair<uint64_t, uint64_t>>> sorted_dpgrp_id;
  uint64_t total_traffic = 0;

  for (auto & entry: dpgrpsz_xfersize) {
    if (entry.first != 1) {
      sorted_dpgrp_id.emplace(std::make_pair(entry.second, entry.first));
      total_traffic += entry.second;
    }
  }

  std::vector<std::pair<uint64_t, int>> n_parallel_rings;
  int n_rings_assigned = 0;
  bool all_ring = false;
  for (auto & entry: sorted_dpgrp_id) {
    int nlink = total_traffic * dp_degree / entry.first;
    if (nlink == 0) nlink = 1;
    if (bidir) {
      nlink *= 2;
    }
    if (nlink + n_rings_assigned > dp_degree) {
      nlink = dp_degree - n_rings_assigned;
      assert(!bidir || nlink % 2 == 0);
    }
    n_rings_assigned += nlink;
    n_parallel_rings.emplace_back(std::make_pair(entry.second, nlink));
    if (entry.second == machine->get_num_nodes()) all_ring = true;
    if (n_rings_assigned == dp_degree) 
      break;
  }
  // corner case: in DP only case or mp degree = 1, makes sure dp is connected
  // if (dp_degree >= if_cnt - 1 && !all_ring) {
  if (!all_ring) {
    n_parallel_rings.back().second -= bidir ? 2 : 1;
    if (n_parallel_rings.back().second == 0) {
      n_parallel_rings.pop_back();
    }
    n_parallel_rings.emplace_back(std::make_pair(machine->get_num_nodes(), bidir ? 2 : 1));
  }
    // }

  // rounding issue...
  while (n_rings_assigned < dp_degree) {
    for (auto & entry: n_parallel_rings) {
      entry.second += bidir ? 2 : 1;
      if (n_rings_assigned == dp_degree) 
        break;
    }
    assert(n_rings_assigned <= dp_degree);
  } 

  // this is for later solving the coin change problem. 
  // contains all the hops one can do on the ring
  std::set<int> coins;

 for (int i = 0; i < n_parallel_rings.size(); i++) {
    auto & curr_dpg = n_parallel_rings[i];
    double mp_satisfied = std::numeric_limits<double>::lowest();
    ConnectionMatrix best = conn;
    std::vector<int> best_jmps;
    selected_jumps.emplace(
      std::make_pair(curr_dpg.first, std::vector<std::vector<int>>{}));
    // for (int j = 0; j < curr_dpg.second; j += bidir ? 2 : 1) {
    // for (auto cj: candidate_jumps.at(curr_dpg.first)) {
      // std::vector<int> jmps = choose_n(candidate_jumps.at(curr_dpg.first), cj, curr_dpg.second / (bidir ? 2 : 1));
      std::vector<int> jmps = choose_n_geo(candidate_jumps.at(curr_dpg.first), curr_dpg.second / (bidir ? 2 : 1));
      ConnectionMatrix proposed = conn;
      for (int j = 0; j < jmps.size(); j++) {
        for (int k = 0; k < ndevs / curr_dpg.first; k++) {
          proposed = add_ring(proposed, k, jmps[j]);
        }
      }
      // ConnectionMatrix hopm = construct_hop_matrix(conn);
      double mp_this = compute_mp_satified(proposed);
#ifdef DEBUG_PRINT
      std::cerr << "MP satisfied: " << mp_this << std::endl;
#endif
      if (mp_this > mp_satisfied || (mp_this == mp_satisfied && unif(gen) > 0.5)) {
        mp_satisfied = mp_this;
        best = proposed;
        best_jmps = jmps;
      }
    // }
    // }
    assert(best_jmps.size() >= 1);
    conn = best;
    for (auto jmp: best_jmps) {
      selected_jumps[curr_dpg.first].emplace_back(std::vector<int>{jmp});
      coins.insert(jmp);
      if (bidir) {
        selected_jumps[curr_dpg.first].emplace_back(std::vector<int>{machine->get_num_nodes()-jmp});
        coins.insert(machine->get_num_nodes()-jmp);
      }
    }
  }

  // construct multihop rings
  // find where are the non-satisfied rings
  std::set<uint64_t> unsatisfied_rings;
  for (auto & entry: candidate_jumps) {
    unsatisfied_rings.insert(entry.first);
  }
  for (auto & entry: n_parallel_rings) {
    unsatisfied_rings.erase(unsatisfied_rings.find(entry.first));
  }

  std::function<bool(const std::vector<int>&, const std::vector<int>&)> vec_len_comp ( 
    [] (const std::vector<int> &lhs, const std::vector<int> &rhs) -> bool {
      return lhs.size() < rhs.size();
    });

  std::vector<int> all2all_coinchange = all_coin_change(coins);
  for (uint64_t entry: unsatisfied_rings) {
    std::set<std::vector<int>, std::function<bool(const std::vector<int>&, const std::vector<int>&)>>
      solutions(vec_len_comp);
    for (auto & cj: candidate_jumps.at(entry)) {
      auto sol = query_path(all2all_coinchange, cj);
      if (sol.size() > 0) 
        solutions.insert(sol);
    }
    std::vector<std::vector<int>> final_choice{};
    // TODO: what's the best solution here... 
    // backoff and add a +1 ring always?
    if (solutions.size() == 0) {
      selected_jumps.emplace(entry, final_choice);
      continue;
    }
    else {
      int shortest_path_len = solutions.begin()->size();
      for (auto &sol: solutions) {
        if (sol.size() == shortest_path_len) {
          for (auto & existing_sol: final_choice) {
            if (segment_overlap(existing_sol, sol)) {
              continue;
            }
          }
          final_choice.emplace_back(sol);
          if (bidir) {
            final_choice.emplace_back(negative(sol));
          }
        }
      }
    }
    selected_jumps.emplace(entry, final_choice);
  }

#ifdef DEBUG_PRINT
  print_all_rings();
#endif
  return n_parallel_rings;
}

void SpMulMat::generate_mp_matching(ConnectionMatrix & conn, int mp_degree) 
{
  // std::cerr << "MP degree: " << mp_degree << std::endl;
  assert(!bidir || mp_degree % 2 == 0);
  auto mp_tm = mp_tm_logical;
  for (int i = 0; i < mp_degree; i++) {
    generate_one_match(conn, mp_tm);
  }
}

void SpMulMat::generate_one_match(ConnectionMatrix & conn, std::unordered_map<uint64_t, uint64_t> & mp_tm)
{
  auto converted = convert_to_blsm_match_graph(mp_tm);
  blossom_match::Matching match(converted.first);
  auto solution = match.SolveMinimumCostPerfectMatching(converted.second);
  for(auto it = solution.first.begin(); it != solution.first.end(); it++){
		std::pair<int, int> e = converted.first.GetEdge( *it );
    conn[edge_id(e.first, e.second)] += 1;
    // std::cerr << "adding " << e.first << ", " << e.second << std::endl;
    if (mp_tm.find(edge_id(e.first, e.second)) != mp_tm.end()) {
      mp_tm[edge_id(e.first, e.second)] *= (double)conn[edge_id(e.first, e.second)]/(conn[edge_id(e.first, e.second)] + 1);
    }
    // if (bidir) {
      conn[edge_id(e.second, e.first)] += 1;
      if (mp_tm.find(edge_id(e.second, e.first)) != mp_tm.end()) {
        mp_tm[edge_id(e.second, e.first)] *= (double)conn[edge_id(e.second, e.first)]/(conn[edge_id(e.second, e.first)] + 1);
      }
    // }
	}
}

std::pair<blossom_match::Graph, std::vector<double>>
SpMulMat::convert_to_blsm_match_graph(std::unordered_map<uint64_t, uint64_t> & mp_tm)
{
  uint64_t max_entry = 0;
  for (auto & entry: mp_tm) {
    if (entry.second > max_entry) max_entry = entry.second;
  }
  size_t num_nodes = machine->get_num_nodes();
  std::vector<double> costs(num_nodes * (num_nodes - 1), 0);
  blossom_match::Graph graph(machine->get_num_nodes());
  for (int i = 0; i < machine->get_num_nodes(); i++) {
    for (int j = 0; j < machine->get_num_nodes(); j++) {
      if (i == j) continue;
      graph.AddEdge(i, j);
      costs[graph.GetEdgeIndex(i, j)] = 
        mp_tm.find(edge_id(i, j)) == mp_tm.end() ? max_entry : max_entry - mp_tm.at(edge_id(i, j));
    }
  }
  return std::make_pair(graph, costs);
}

ConnectionMatrix SpMulMat::add_ring(const ConnectionMatrix & conn, int start, int dist)
{
  ConnectionMatrix result{conn};
  int total_devs = machine->get_total_devs();
  int curr = start;
  int next = (start + dist) % total_devs;
  do {
    result[edge_id(curr, next)]++;
    if (bidir) result[edge_id(next, curr)]++;
    curr = next;
    next = (next + dist) % total_devs;
  } while (curr != start);
  return result;
}

double SpMulMat::compute_mp_satified(const ConnectionMatrix & conn) 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  ShortestPathNetworkRoutingStrategy s{conn, nm->ids_to_nw_comm_device, nm->total_devs};
  double result = 0;
  
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  for (int i = 0; i < nnode; i++) {
    // auto hopcnts = s.hop_count(i);
    // std::cerr << "from " << i << " to " << j << " discounted: " << discounted_hop 
    //   << "(hop cnt: " << hop_cnt << ", narrowest: " << narrowest << std::endl;
    for (int j = 0; j < nnode; j++) {
      if (j == i) { 
        continue; 
      }
      if (mp_tm_logical.find(edge_id(i, j)) != mp_tm_logical.end()) {
        if (conn[edge_id(i, j)] > 0) {
          result += mp_tm_logical.at(edge_id(i, j));
        } else {
          // if (hopcnts[j].first == -1) continue;
          // result += mp_tm_logical.at(edge_id(i, j)) / hopcnts[j].first;
        }
      }
    }
  }
  return result;
}

// std::vector<int> SpMulMat::coin_change(const std::set<int> & coins, int goal) 
// {
//   std::vector<int> coins_vec{coins.cbegin(), coins.cend()};
//   int min_coin = coins_vec[0];
//   if (min_coin < 0) {
//     std::transform(coins_vec.cbegin(), coins_vec.cend(), coins_vec.begin(),
//       [&] (int i) -> int {return i > 0 ? i : machine->get_num_nodes() + i;}
//     );
//     // goal = goal - min_coin + 1;
//     std::sort(coins_vec.begin(), coins_vec.end());
//     coins_vec.erase(std::unique(coins_vec.begin(), coins_vec.end()), coins_vec.end());
//   }

//   std::vector<int> dist(machine->get_num_nodes(), 0);
//   std::vector<int> back_idx(machine->get_num_nodes(), -1);
//   for (int c: coins_vec) {
//     dist[c - 1] = 1;
//     back_idx[c - 1] = c - 1;
//   }

//   for (int i = 0; i < goal; i++) {
//     if (dist[i] == 1) continue;
//     int index_candididate = -1;
//     int min_dist = 0;
//     for (int c: coins_vec) {
//       // if (i - (c - 1) >= 0) {
//         if (!min_dist || (dist[MOD(i-(c-1), machine->get_num_nodes())] + 1 < min_dist)) {
//           min_dist = dist[MOD(i-(c-1), machine->get_num_nodes())] + 1;
//           index_candididate = c-1;
//         }
//       // }
//     }
//   }

//   if (back_idx.back() == -1) {
//     return {};
//   }
//   std::vector<int> result;
//   int curr_coin = back_idx.back();
//   int curr_pos = goal;

//   while (curr_pos != curr_coin) {
//     result.push_back(coins.find(curr_coin) != coins.end() ? curr_coin : curr_coin - machine->get_num_nodes());
//     curr_pos = MOD(curr_pos - curr_coin, machine->get_num_nodes());
//     curr_coin = back_idx[curr_pos];
//   }
//   return result;
// }

std::vector<int> SpMulMat::all_coin_change(const std::set<int> & coins)
{
  // std::vector<int> coins_vec{coins.cbegin(), coins.cend()};
  // int min_coin = coins_vec[0];
  // if (min_coin < 0) {
  //   std::transform(coins_vec.cbegin(), coins_vec.cend(), coins_vec.begin(),
  //     [&] (int i) -> int {return i > 0 ? i : machine->get_num_nodes() + i;}
  //   );
  //   // goal = goal - min_coin + 1;
  //   std::sort(coins_vec.begin(), coins_vec.end());
  //   coins_vec.erase(std::unique(coins_vec.begin(), coins_vec.end()), coins_vec.end());
  // }
  std::vector<std::vector<int>> backtraces{};
  std::vector<std::vector<int>> dists{};
  std::vector<int> init_dists(machine->get_num_nodes(), -1);
  std::vector<int> init_backtrace(machine->get_num_nodes(), -1);
  for (int c: coins) {
    init_dists[c] = 0;
    init_backtrace[c] = c;
  }
  backtraces.emplace_back(init_backtrace);
  dists.emplace_back(init_dists);
  bool finished = std::find(dists.back().begin(), dists.back().end(), -1) == dists.back().end();
  while (!finished) {
  #ifdef DEBUG_PRINT
    std::cerr << "iter coin change: " << std::endl;
    for (int i : backtraces.back()) {
      std::cerr << i << ", ";
    }
    std::cerr << std::endl;
    for (int i : dists.back()) {
      std::cerr << i << ", ";
    }
    std::cerr << std::endl;
  #endif
    finished = true;
    std::vector<int> curr_dists(machine->get_num_nodes(), -1);
    std::vector<int> curr_backtrace(machine->get_num_nodes(), -1);
    for (int i = 0; i < machine->get_num_nodes(); i++) {
      if (dists.back()[i] != -1) {
        curr_dists[i] = dists.back()[i];
        curr_backtrace[i] = backtraces.back()[i];
      }
      else {
        for (int c: coins) {
          if (dists.back()[MOD(i-(c), machine->get_num_nodes())] != -1 && 
              (curr_dists[i] == -1 || 
              (dists.back()[MOD(i-(c), machine->get_num_nodes())] + 1 < curr_dists[i] || 
              (dists.back()[MOD(i-(c), machine->get_num_nodes())] + 1 == curr_dists[i] && unif(gen) > 0.5)))) {
            curr_dists[i] = dists.back()[MOD(i-(c), machine->get_num_nodes())] + 1;
            curr_backtrace[i] = c; 
          }
        }
        if (curr_dists[i] == -1) {
          finished = false;
        }
      }
    }
    backtraces.emplace_back(curr_backtrace);
    dists.emplace_back(curr_dists);
  }
#ifdef DEBUG_PRINT 
  std::cerr << "solution to coin change: " << std::endl;
  for (int i : backtraces.back()) {
    std::cerr << i << ", ";
  }
  std::cerr << std::endl;
#endif
  return backtraces.back();
}

std::vector<int> SpMulMat::query_path(const std::vector<int>& candidates, int jump)
{
  // jump;
  std::vector<int> result;
  int next_coin = candidates[jump];
  while (jump) {
    result.push_back(next_coin);
    jump = MOD(jump - (next_coin), machine->get_num_nodes());
    next_coin = candidates[jump];
  }
  // result.push_back(next_coin);
  return result;
}

bool SpMulMat::optimize(int mcmc_iter, double sim_iter_time, bool forced) 
{
  double diff = sim_iter_time - curr_sim_time;
  bool change = diff < 0 ? true : diff != 0 && unif(gen) < std::exp(-alpha * diff);
  if (sim_iter_time < best_sim_time) {
    best_sim_time = sim_iter_time;
    change = true;
  }
  std::cerr << "sim_iter_time: " << sim_iter_time << ", curr_sim_time: " << curr_sim_time 
          << ", best_iter_time: " << best_sim_time << ", change: " << change << std::endl;
  if (change) {
    curr_sim_time = sim_iter_time;
  }
  else {
    num_iter_nochange++; 
  }

  if (!forced && !change && num_iter_nochange < no_improvement_th)
    return false;

  std::cerr << "Changing... " << std::endl;
  selected_jumps.clear();
  for (auto & entry: dp_ncomms) {
    for (auto & c: entry.second) {
      if (c != nullptr) {
        delete c;
        c = nullptr;
      }
    }
  }
  dp_ncomms.clear();
  
  num_iter_nochange = 0;
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  
  // TODO: copy machine-switch link?
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  ConnectionMatrix dpconn = std::vector<int>(ndevs*ndevs, 0);
  ConnectionMatrix mpconn = std::vector<int>(ndevs*ndevs, 0);

  int dp_degree, mp_degree;
  get_dp_mp_degree(dp_degree, mp_degree);
  auto dp_rings = generate_dp_topology(dpconn, dp_degree);
  generate_mp_matching(mpconn, mp_degree);
  // NetworkTopologyGenerator::print_conn_matrix(mpconn, ndevs, 0);

  ConnectionMatrix final = mpconn + dpconn;
  // connect_topology(final, dp_rings, mp_degree);
  nm->set_topology(final);
  nm->update_route();
// #ifdef DEBUG_PRINT
  NetworkTopologyGenerator::print_conn_matrix(final, ndevs, 0);
// #endif
  return true;
}

void SpMulMat::construct_topology()
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  
  // TODO: copy machine-switch link?
  size_t nnode = machine->get_num_nodes();
  size_t ndevs = machine->get_total_devs();
  ConnectionMatrix dpconn = std::vector<int>(ndevs*ndevs, 0);
  ConnectionMatrix mpconn = std::vector<int>(ndevs*ndevs, 0);

  int dp_degree, mp_degree;
  get_dp_mp_degree(dp_degree, mp_degree);
  auto dp_rings = generate_dp_topology(dpconn, dp_degree);
  generate_mp_matching(mpconn, mp_degree);

  // NetworkTopologyGenerator::print_conn_matrix(mpconn, ndevs, 0);

  ConnectionMatrix final = mpconn + dpconn;
  // connect_topology(final, dp_rings, mp_degree);
  nm->set_topology(final);
  nm->update_route();
#ifdef DEBUG_PRINT
  NetworkTopologyGenerator::print_conn_matrix(final, ndevs, 0);
#endif
}

ConnectionMatrix SpMulMat::connect_topology(const ConnectionMatrix & conn, 
  ConnectionMatrix & mp_conn, ConnectionMatrix & dp_conn,
  const std::vector<std::pair<uint64_t, int>> & dp_rings, int mp_degree)
{
  // TODO: better way. This is to be implemented...
  return dp_conn+mp_conn;
#if 0
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  size_t num_nodes = nm->num_nodes;
  size_t ndevs = nm->num_switches + num_nodes;
  // find connected components 
  int n_cc = 0;
  std::vector<int> node_to_ccid = std::vector<int>(num_nodes, -1);
  std::vector<std::set<size_t> > ccs;
  std::queue<size_t> search_q;
  for (size_t i = 0; i < num_nodes; i++) {
    if (node_to_ccid[i] == -1) {
      search_q.push(i);
      node_to_ccid[i] = n_cc++;
      ccs.emplace_back();
      ccs.back().insert(i);
      while (!search_q.empty()) {
        size_t curr = search_q.front();
        search_q.pop();
        for (size_t j = 0; j < num_nodes; j++) {
          if (curr != j && conn[edge_id(curr, j)] > 0 && node_to_ccid[j] == -1) {
            node_to_ccid[j] = node_to_ccid[curr];
            ccs.back().insert(j);
            search_q.push(j);
          }
        }
      }
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
    if (mp_degree > 1 && bidir) {
      int v00, v01, v10, v11;
      while (n_cc > 1) {
        if (ccs[0].size() == 1 && ccs[1].size() == 1) {
          bool success = add_link(*ccs[0].begin(), *ccs[1].begin(), mp_conn);
          assert(success);
          success = add_link(*ccs[0].begin(), *ccs[1].begin(), mp_conn);
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
                  mp_conn[edge_id(i, j)] > 0) {
                uint64_t eid = edge_id(i, j);
                uint64_t ueid = unordered_edge_id(i, j);
                if (mp_tm_logical.find(eid) == mp_tm_logical.end()) {
                  e_to_remove = ueid;
                  break;
                }
                else {
                  if (mp_tm_logical[eid] < min_demand) {
                    min_demand = mp_tm_logical[eid];
                    e_to_remove = ueid;
                  }
                }
              }
            }
          }
          assert(e_to_remove != 0);

          // std::cout << "1-n removing " << e_to_remove % ndevs << ", " <<  e_to_remove / ndevs << std::endl;
          remove_link(e_to_remove % ndevs, e_to_remove / ndevs, mp_conn);
          bool success = add_link(*ccs[singleton].begin(), e_to_remove % ndevs, mp_conn);
          assert(success);
          success = add_link(*ccs[singleton].begin(), e_to_remove / ndevs, mp_conn);
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
            auto liter = mp_tm_logical.find(lhs);
            uint64_t l = liter == mp_tm_logical.end() ? 0 : liter->second;
            auto riter = mp_tm_logical.find(rhs);
            uint64_t r = riter == mp_tm_logical.end() ? 0 : riter->second;
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
              remove_link(v00, v01, mp_conn);
              remove_link(v10, v11, mp_conn);
              bool success = add_link(v00, v11, mp_conn);
              assert(success);
              success = add_link(v01, v10, mp_conn);
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
  }
#endif
}

// std::vector<std::vector<int>> SpMulMat::get_selected_jumps(int grp_sz) 
// {
//   return selected_jumps[grp_sz];
// }

void SpMulMat::reset() 
{
  mp_tm_logical.clear();
  dpgrpsz_xfersize.clear();
  dpgrps.clear();
  // selected_jumps.clear();
}

std::unique_ptr<L1OptimizerInformation> SpMulMat::export_information() 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  SpMulMatInformation * info = new SpMulMatInformation();
  info->conn = nm->conn_matrix;
  // info->dp_ncomms = dp_ncomms;
  info->selected_jumps = selected_jumps;
  return std::unique_ptr<L1OptimizerInformation>(info);
}


void SpMulMat::import_information(const std::unique_ptr<L1OptimizerInformation>& information) 
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  SpMulMatInformation * info = static_cast<SpMulMatInformation*>(information.get());
  // dp_ncomms = info->dp_ncomms;
  for (auto & entry: dp_ncomms) {
    for (auto & c: entry.second) {
      if (c != nullptr) {
        delete c;
        c = nullptr;
      }
    }
  }
  dp_ncomms.clear();
  selected_jumps = info->selected_jumps;
  nm->set_topology(info->conn);
  nm->update_route();
}

void SpMulMat::delete_information(const std::unique_ptr<L1OptimizerInformation>& information)
{
  // SpMulMatInformation * info = (SpMulMatInformation*)information;
  // // for (auto & entry: info->dp_ncomms) {
  // //   for (auto & c: entry.second) {
  // //     if (c != nullptr) {
  // //       delete c;
  // //       c = nullptr;
  // //     }
  // //   }
  // // }
  // delete info;
}

void SpMulMat::store_tm() const
{

}

const std::vector<NominalCommDevice*>& SpMulMat::get_dp_ncomms(int src, int grp_sz)
{
  NetworkedMachineModel * nm = static_cast<NetworkedMachineModel*>(this->machine);
  uint64_t key = ((uint64_t)src) << 32 | grp_sz;
  if (dp_ncomms.find(key) != dp_ncomms.end()) {
    return dp_ncomms[key];
  }
  int total_devs = machine->get_num_nodes();
  auto & paths = selected_jumps[grp_sz];
  dp_ncomms[key] = std::vector<NominalCommDevice*>{};
  for (int i = 0; i < paths.size(); i++) {
    int total_hops = 0;
    for (int h: paths[i]) total_hops += h;
    int dst = MOD(src + (total_hops), total_devs);
    std::string link_name = "NOMINAL_DP_" + std::to_string(src) + "-" + std::to_string(dst) + "-" + std::to_string(grp_sz) + "-" + std::to_string(i);
    int device_id = src * total_devs + dst;
    NominalCommDevice * ncomm = new NominalCommDevice(link_name, device_id, total_devs, nm->routing_strategy);
    Route r;
    int curr = src;
    int next;
    for (int j = 0; j < paths[i].size(); j++) {
      next = MOD(curr + paths[i][(j + 1) % paths[i].size()], total_devs);
      r.push_back(nm->ids_to_nw_comm_device.at(curr * total_devs + next));
      curr = next;
    }
    ncomm->set_physical_paths(std::make_pair(std::vector<double>{1}, std::vector<Route>{r}));
    dp_ncomms[key].push_back(ncomm);
  }
  return dp_ncomms[key];
}

void SpMulMat::print_all_rings() const 
{
  for (auto & grp_paths: selected_jumps) {
    std::cerr << "Group size: " << grp_paths.first << std::endl;
    for (auto & path: grp_paths.second) {
      std::cerr << "\t";
      for (int hop: path) {
        std::cerr << hop << ", ";
      }
      std::cerr << std::endl;
    }
  }
  std::cerr << std::endl;
}
