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

// all prime numbers below 10000. I know how to write a sieve but this is good enouggh.
const static uint16_t PRIMES[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, 8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831, 8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973};

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
  return std::make_pair(std::vector<float>{1}, std::vector<Route>{result});
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

size_t DemandHeuristicNetworkOptimizer::edge_id(int i, int j) const
{
  return i * machine->get_total_devs() + j;
}

size_t DemandHeuristicNetworkOptimizer::unordered_edge_id(int i, int j) const
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

  for (int i = 0; i < nnode; i++) {
    conn[edge_id(i, (i + nnode / 5) % nnode)]++;
    conn[edge_id(i, 2 * (i + nnode / 5) % nnode)]++;
    conn[edge_id(i, 3 * (i + nnode / 5) % nnode)]++;
    conn[edge_id(i, 4 * (i + nnode / 5) % nnode)]++;
  }
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

void DemandHeuristicNetworkOptimizerPlus::optimize(int mcmc_iter, float sim_iter_time)
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

  connectivity_assign(conn, max_of_bidir, node_if_allocated);
// #ifdef DEBUG_PRINT
  std::cerr << "After connectivity_assign " << std::endl;
  NetworkTopologyGenerator::print_conn_matrix(conn, nnode, 0);
// #endif
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

}

SpMulMat::SpMulMat(MachineModel * machine, int degree, bool bidir)
: L1Optimizer(machine), degree(degree), bidir(bidir)
{

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
    dpg.starting_node = reinterpret_cast<int>(task->next_tasks[0]);
    dpg.xfer_size = task->xfer_size;
    dpgrps.emplace_back(dpg);
    INSERT_OR_ADD(dpgrp_xfersize, dpgrp_unique_key(dpg), (double)task->xfer_size * (2.0 * (double)dpg.group_size - 1) / (double)dpg.group_size);
  break;

  case SimTask::TASK_BARRIER:

  break;

  case SimTask::TASK_UPDATE:

  break;
  }
}

size_t SpMulMat::edge_id(int i, int j) const
{
  return i * machine->get_total_devs() + j;
}

size_t SpMulMat::unordered_edge_id(int i, int j) const
{
  return i > j ? edge_id(i, j) : edge_id(j, i);
}

uint64_t SpMulMat::dpgrp_unique_key(const DPGroup & dpg) const
{
  return (((uint64_t)dpg.group_size) << 32) | (dpg.starting_node);
}

int SpMulMat::get_start_node(uint64_t id) const 
{
  return (int)(id & 0x00000000FFFFFFFFULL);
}

int SpMulMat::get_group_size(uint64_t id) const 
{
  return (int)((id & 0xFFFFFFFF00000000ULL) >> 32);

}

void SpMulMat::construct_candidate_jumps() {
  int total_nodes = machine->get_num_nodes();
  for (int i = 1; i <= total_nodes; i++) {
    if (total_nodes % i == 0) {
      int group_size = i;
      int base_jump = total_nodes / i;
      int nconfs = total_nodes % i;
      candidate_jumps.emplace(std::make_pair(group_size, std::vector<int>{}));
      for (int k = 0; ; k++) {
        if (base_jump * PRIMES[k] >= group_size) {
          break;
        }
        candidate_jumps[group_size].emplace_back(base_jump * PRIMES[k]);
      }
    }
  }
}

void SpMulMat::get_dp_mp_degree(int & dp_degree, int & mp_degree)
{
  double total_dp_traffic = 0;
  double total_mp_traffic = 0; 

  for (auto & entry: dpgrps) { 
    // ring allreduce
    total_dp_traffic += entry.xfer_size * (2 * entry.group_size - 1);
  }

  for (auto & entry: mp_tm_logical) {
    total_dp_traffic += entry.second;
  }

  dp_degree = lround(total_dp_traffic / (total_dp_traffic + total_mp_traffic) * degree);
  mp_degree = lround(total_mp_traffic / (total_dp_traffic + total_mp_traffic) * degree);

  if (bidir && dp_degree % 2 != 0) {
    mp_degree--;
    dp_degree++;
  }
  assert(mp_degree + dp_degree == degree);
}

void SpMulMat::generate_dp_topology(ConnectionMatrix & conn, int dp_degree) 
{
  size_t ndevs = machine->get_num_nodes();
  // reset conn
  conn = ConnectionMatrix(ndevs * ndevs, 0);

  // sort out the DP groups by weight
  std::set<std::pair<uint64_t, uint64_t>, std::greater<std::pair<uint64_t, uint64_t>>> sorted_dpgrp_id;
  uint64_t total_traffic = 0;

  for (auto & entry: dpgrp_xfersize) {
    sorted_dpgrp_id.emplace(std::make_pair(entry.second, entry.first));
    total_traffic += entry.second;
  }

  std::vector<std::pair<uint64_t, int>> n_parallel_rings;
  int n_rings_assigned = 0;
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
    if (n_rings_assigned == dp_degree) 
      break;
  }

  // rounding issue...
  while (n_rings_assigned < dp_degree) {
    for (auto & entry: n_parallel_rings) {
      entry.second += bidir ? 2 : 1;
      if (n_rings_assigned == dp_degree) 
        break;
    }
    assert(n_rings_assigned <= dp_degree);
  } 

  // O(d. N/ln N . N^2)
  for (int i = 0; i < n_parallel_rings.size(); i++) {
    auto & curr_dpg = n_parallel_rings[i];
    double mp_satisfied = std::numeric_limits<double>::lowest();
    ConnectionMatrix best = conn;
    int best_jump = -1;
    selected_jumps.emplace(
      std::make_pair(curr_dpg.first, std::vector<std::vector<int>>{}));
    for (int j = 0; j < curr_dpg.second; j += bidir ? 2 : 1) {
      for (int k = 0; k < candidate_jumps.at(curr_dpg.first).size(); k++) {
        ConnectionMatrix proposed = 
          add_ring(conn, get_start_node(curr_dpg.first), candidate_jumps.at(curr_dpg.first)[k]);
        ConnectionMatrix hopm = construct_hop_matrix(conn);
        double mp_this = compute_mp_satified(hopm);
        if (mp_this > mp_satisfied || (mp_this == mp_satisfied && unif(gen) > 0.5)) {
          mp_satisfied = mp_this;
          best = proposed;
          best_jump = candidate_jumps.at(curr_dpg.first)[k];
        }
      }
    }
    assert(best_jump != -1);
    conn = best;
    selected_jumps[curr_dpg.first].emplace_back(std::vector<int>{best_jump});
  }

}

void SpMulMat::generate_mp_matching(ConnectionMatrix & conn, int dp_degree) 
{

}

ConnectionMatrix SpMulMat::add_ring(const ConnectionMatrix & conn, int start, int dist)
{

}

ConnectionMatrix SpMulMat::construct_hop_matrix(const ConnectionMatrix & conn)
{

}

double SpMulMat::compute_mp_satified(const ConnectionMatrix & conn) 
{

}

