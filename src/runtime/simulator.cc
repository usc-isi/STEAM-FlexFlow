/* Copyright 2020 Stanford
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

#include <random>
#include <queue>
#include <stack>
#include "simulator.h"
#include "model.h"

// #include "flatbuffers/util.h"
#include "taskgraph_generated.h"

// #define DEBUG_PRINT

int ParallelConfig::num_parts() const
{
  int nparts = 1;
  for (int i = 0; i < nDims; i++)
    nparts *= dim[i];
  return nparts;
}

bool ParallelConfig::is_data_parallel() const
{
  int nparts = 1;
  for (int i = 0; i < nDims; i++) {
    nparts *= dim[i];
    if ((i < nDims-1) && (dim[i] > 1))
      return false;
  }
  for (int i = 0; i < nparts; i++)
    if (device_ids[i] != i)
      return false;
  return true;
}

Device::Device(std::string const &name, DeviceType type, int node_id, int socket_id, int device_id)
: name(name), type(type), node_id(node_id), socket_id(socket_id), device_id(device_id)
{}

CompDevice::CompDevice(std::string const &name, CompDevType comp_type, int node_id, int socket_id, int device_id)
: Device(name, Device::DEVICE_COMP, node_id, socket_id, device_id), comp_type(comp_type)
{}

MemDevice::MemDevice(std::string const &name, MemDevType mem_type, int node_id, int socket_id, int device_id, size_t capacity)
: Device(name, Device::DEVICE_MEM, node_id, socket_id, device_id), mem_type(mem_type), capacity(capacity)
{}

CommDevice::CommDevice(std::string const &name, CommDevType comm_type, int node_id, int socket_id, int device_id, float latency, float bandwidth)
: Device(name, Device::DEVICE_COMM, node_id, socket_id, device_id), comm_type(comm_type), latency(latency), bandwidth(bandwidth)
{}

/* I hate this but this makes sense here... */
static std::random_device rd; 
static std::mt19937 gen = std::mt19937(rd()); 
static std::uniform_real_distribution<> std_uniform = std::uniform_real_distribution<>(0.0, 1.0); 

NominalCommDevice::NominalCommDevice(std::string const &name, int device_id, const EcmpRoutes& routes) 
: CommDevice(name, CommDevice::NW_NOMINAL, -1, -1, device_id, 0, 0), routes(routes)
{}

NominalCommDevice::NominalCommDevice(std::string const &name, int device_id) 
: CommDevice(name, CommDevice::NW_NOMINAL, -1, -1, device_id, 0, 0)
{}
    
Route NominalCommDevice::expand_to_physical() const 
{
  int pick = 0;
  double choice = std_uniform(gen);
  for (int i = 0; i < routes.first.size(); i++) {
    if (choice > routes.first[i]) break;
    pick = i;
  }
  Route ret = Route(routes.second[pick].begin(), routes.second[pick].end());
  return ret;
}

void NominalCommDevice::set_physical_paths(const EcmpRoutes &rs) 
{
  routes = rs;
}

const EcmpRoutes & NominalCommDevice::get_all_routes() 
{
  return routes;
}

SimTask::SimTask()
{}

void SimTask::add_next_task(SimTask* task)
{
  next_tasks.push_back(task);
  task->counter ++;
}

std::string SimTask::get_type_str() const {
  switch (type) {
    case TASK_FORWARD:
      return "Forward";
    case TASK_BACKWARD:
      return "Backward";
    case TASK_COMM:
      return "Comm";
    case TASK_UPDATE:
      return "Update";
    case TASK_BARRIER:
      return "Barrier";
    default:
      assert(false && "Unknown task type");
  }
}

TaskManager::TaskManager(size_t _max_num_tasks)
: max_num_tasks(_max_num_tasks)
{
  tasks = (SimTask**) malloc(sizeof(SimTask*) * max_num_tasks);
  for (size_t i = 0; i < max_num_tasks; i++) {
    tasks[i] = new SimTask();
  }
}

void TaskManager::reset()
{
  global_task_id = 0;
  hash_to_forward_task.clear();
  hash_to_backward_task.clear();
}

SimTask* TaskManager::new_task()
{
  assert(global_task_id + 1 < max_num_tasks);
  SimTask* task = tasks[global_task_id++];
  task->ready_time = 0.0f;
  task->run_time = 0.0f;
  task->next_tasks.clear();
  task->counter = 0;
  task->device = NULL;
  task->mem = NULL;
  task->name.clear();

  // task->from_dev = -1;
  // task->to_dev = -1;
  task->xfer_size = 0;
  task->store = true;
  
  return task;
}

SimTask* TaskManager::new_update_task()
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_UPDATE;
  return task;
}

SimTask* TaskManager::new_barrier_task()
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_BARRIER;
  return task;
}

SimTask* TaskManager::new_comm_task()
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_COMM;
  return task;
}


SimTask* TaskManager::new_nominal_comm_task()
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_NOMINAL_COMM;
  return task;
}

SimTask* TaskManager::new_comm_task(std::string const &name, CommDevice *comm_device, size_t message_size)
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_COMM;
  task->name = name;
  task->device = comm_device;
  task->run_time = comm_device->latency + message_size / comm_device->bandwidth;
  return task;
}

SimTask* TaskManager::new_nominal_comm_task(std::string const &name, CommDevice *comm_device, size_t message_size)
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_NOMINAL_COMM;
  task->name = name;
  task->device = comm_device;
  task->run_time = comm_device->latency + message_size / comm_device->bandwidth;
  return task;
}

SimTask* TaskManager::new_forward_task(Op* op, int idx)
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_FORWARD;
  size_t hash = 17 * 31 + (size_t)(op);
  hash = hash * 31 + std::hash<int>()(idx);
  hash_to_forward_task[hash] = task;
  task->name = op->name;
  return task;
}

SimTask* TaskManager::new_backward_task(Op* op, int idx)
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_BACKWARD;
  size_t hash = 17 * 31 + (size_t)(op);
  hash = hash * 31 + std::hash<int>()(idx);
  hash_to_backward_task[hash] = task;
  task->name = op->name;
  return task;
}

SimTask* TaskManager::new_allreduce_task(Op *op, const std::vector<int> &node_ids, size_t message_size) 
{
  SimTask* task = new_task();
  task->type = SimTask::TASK_ALLREDUCE;
  // task->counter = node_ids[0];
  for (int i = 0; i < node_ids.size(); i++) {
    task->next_tasks.push_back(reinterpret_cast<SimTask*>(node_ids[i]));
  } 
  task->xfer_size = message_size;
  return task;
}

SimTask* TaskManager::get_forward_task(Op* op, int idx)
{
  size_t hash = 17 * 31 + (size_t)(op);
  hash = hash * 31 + std::hash<int>()(idx);
  assert(hash_to_forward_task.find(hash) != hash_to_forward_task.end());
  return hash_to_forward_task[hash];
}

SimTask* TaskManager::get_backward_task(Op* op, int idx)
{
  size_t hash = 17 * 31 + (size_t)(op);
  hash = hash * 31 + std::hash<int>()(idx);
  assert(hash_to_backward_task.find(hash) != hash_to_backward_task.end());
  return hash_to_backward_task[hash];
}

void Simulator::free_all()
{
  offset = 0;
}

size_t data_type_size(DataType type) {
  switch (type) {
    case DT_FLOAT:
      return sizeof(float);
    case DT_DOUBLE:
      return sizeof(double);
    case DT_INT32:
      return sizeof(int32_t);
    case DT_INT64:
      return sizeof(int64_t);
    case DT_BOOLEAN:
      return sizeof(bool);
    default:
      assert(false);
  }
}

void* Simulator::allocate(uint64_t num_elements, DataType type)
{
  uint64_t element_size = data_type_size(type);
  void* ret_ptr = base_ptr + offset;
  offset += element_size * num_elements;
  if (offset > capacity) {
    fprintf(stderr, "Simulator cannot measure some operators' performance."
        " Increate --simulator-workspace-size to at least %lu. element_size: %lu, num_elements: %lu, capacity: %lu\n", offset, element_size, num_elements, capacity);
    exit(0);
  }
  return ret_ptr;
}

void Simulator::add_task_dependencies_with_xfer(SimTask* src_task,
                                                SimTask* dst_task,
                                                size_t message_size)
{
  std::vector<CommDevice *> path = machine->get_comm_path(src_task->mem, dst_task->mem);
  // print the communication path
  // printf("Path from %s to %s is: ", src_task->mem->name.c_str(), dst_task->mem->name.c_str());
  // for (size_t i = 0; i < path.size(); i++) {
  //   printf("%s ", path[i]->name.c_str());
  // }
  // printf("\n");

  if (path.empty()) {
    src_task->add_next_task(dst_task);
    return;
  }
  assert(message_size > 0);
  std::vector<std::vector<SimTask *>> all_tasks;
  // Limit the max number of segments per message
  int seg_size = segment_size;
  int num_segment = message_size / seg_size;
  if (message_size % seg_size != 0) {
    num_segment += 1;
  }
  if (num_segment > max_num_segments) {
    num_segment = max_num_segments;
    seg_size = message_size / num_segment;
  }
  // Create all the comm tasks
  // Divide messages into segments
  for (size_t i = 0; i < path.size(); i++) {
    all_tasks.push_back({});
    for (int j = 0; j < num_segment; j++) {
      int cur_seg_size = seg_size;
      if (j == num_segment - 1) {
        cur_seg_size = message_size - (num_segment - 1) * seg_size;
      }
      std::string name = "seg " + std::to_string(j) + " from " + src_task->name + " to " + dst_task->name;
      SimTask *cur_task = task_manager->new_comm_task(name, path[i], cur_seg_size);
      all_tasks[i].push_back(cur_task);
    }
  }

  // Add dependencies among the comm tasks
  for (size_t i = 0; i < path.size(); i++) {
    for (int j = 0; j < num_segment; j++) {
      if (i == 0) {
        src_task->add_next_task(all_tasks[i][j]);
      }
      if (i == path.size() - 1) {
        all_tasks[i][j]->add_next_task(dst_task);
      }
      if (i > 0) {
        all_tasks[i-1][j]->add_next_task(all_tasks[i][j]);
      }
    }
  }

  // Add special dependencies for upi_ins, upi_outs, nic_ins, and nic_outs to prevent communication
  // overlap between upi_ins and upi_outs, and between nic_ins and nic_outs.
  if (num_segment > 1 and path.size() >= 2) {
    for (size_t i = 0; i < path.size(); i++) {
      for (int j = 1; j < num_segment; j++) {
        if (((CommDevice *)all_tasks[i][j]->device)->comm_type == CommDevice::NIC_OUT_COMM or
            ((CommDevice *)all_tasks[i][j]->device)->comm_type == CommDevice::UPI_OUT_COMM) {
          all_tasks[i+1][j-1]->add_next_task(all_tasks[i][j]);
        }
      }
    }
  }

  // call l1 optimizer's call back
  for (std::vector<SimTask *> & tv: all_tasks) {
    for (SimTask * t: tv) {
      if (l1optimizer) 
        l1optimizer->task_added(t);
    }
  }
}

void LogicalTaskgraphBasedSimulator::add_task_dependencies_with_xfer(
                                                SimTask* src_task,
                                                SimTask* dst_task,
                                                size_t message_size)
{
  std::vector<CommDevice *> path = machine->get_comm_path(src_task->mem, dst_task->mem);
#ifdef DEBUG_PRINT
  // print the communication path
  printf("Path from %s to %s is: ", src_task->mem->name.c_str(), dst_task->mem->name.c_str());
  for (size_t i = 0; i < path.size(); i++) {
    printf("%s ", path[i]->name.c_str());
  }
  printf("\n");
#endif

  if (path.empty()) {
    src_task->add_next_task(dst_task);
    return;
  }
  assert(message_size > 0);
  std::vector<SimTask*> final_tasks;
  for (CommDevice * d: path) {
    SimTask* task = task_manager->new_nominal_comm_task();
    task->device = d;
    task->run_time = 0;
    task->xfer_size = message_size;
    if (!final_tasks.empty()) {
      final_tasks.back()->add_next_task(task);
    }
    final_tasks.push_back(task);
    if (l1optimizer) 
      l1optimizer->task_added(task);
  }
  src_task->add_next_task(final_tasks[0]);
  final_tasks.back()->add_next_task(dst_task);
}

[[noreturn]] void handle_measure_operator_cost_unimplemented(Op const *op) {
    std::cerr << "measure_operator_cost not implemented for op "
              << op->name
              << " (type " << op->op_type << ")"
              << ". Please report this issue to the FlexFlow developers."
              << std::endl;
    std::abort();
}

CostMetrics Simulator::measure_operator_cost(Op* op, const ParallelConfig& config)
{
  size_t hash = 17 * 31 + (size_t)(op);
  hash = hash * 31 + std::hash<int>()(config.device_type);
  hash = hash * 31 + std::hash<int>()(config.nDims);
  for (int i = 0; i < config.nDims; i++)
    hash = hash * 31 + std::hash<int>()(config.dim[i]);
  std::map<size_t, CostMetrics>::const_iterator iter =
    hash_to_operator_cost.find(hash);
  if (iter == hash_to_operator_cost.end()) {
    CostMetrics cost_metrics;
    bool is_implemented = op->measure_operator_cost(this, config, cost_metrics);
    if (! is_implemented) {
      handle_measure_operator_cost_unimplemented(op);
    }
    hash_to_operator_cost[hash] = cost_metrics;
    return cost_metrics;
  } else {
    return iter->second;
  }
}

float Simulator::simulate_runtime(const FFModel* model,
                                  const std::map<Op*, ParallelConfig>& global,
                                  CompMode comp_mode)
{
  return this->simulate_runtime(model, global, comp_mode, "");
}

float Simulator::simulate_runtime(const FFModel* model,
                                  const std::map<Op*, ParallelConfig>& global,
                                  CompMode comp_mode,
                                  std::string const &export_file_name)
{
  // printf("%s\n", machine->to_string().c_str());
  task_manager->reset();
  if (l1optimizer)
    l1optimizer->reset();
  // Step 1: register forward and backward tasks
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    CostMetrics cost_metrics = measure_operator_cost(op, config);
    float forward_time = cost_metrics.forward_time;
    float backward_time = cost_metrics.backward_time;
    for (int j = 0; j < config.num_parts(); j++) {
      SimTask* task1 = task_manager->new_forward_task(op, j);
      task1->device = machine->get_gpu(config.device_ids[j]);
      task1->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
      task1->run_time = forward_time;
      if (l1optimizer) 
        l1optimizer->task_added(task1);
      if (comp_mode == COMP_MODE_TRAINING) {
        SimTask* task2 = task_manager->new_backward_task(op, j);
        task2->device = machine->get_gpu(config.device_ids[j]);
        task2->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
        task2->run_time = backward_time;
        task1->add_next_task(task2);
        if (l1optimizer) 
          l1optimizer->task_added(task2);
      }
    }
  }
  // Step 2: insert dependencies and comm. tasks before compute tasks
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    for (int j = 0; j < op->numInputs; j++) {
      Tensor t = op->inputs[j];
      Op* pre_op = t.owner_op;
      if (pre_op == NULL)
        continue;
      ParallelConfig pre_config = global.find(pre_op)->second;
      size_t element_size = data_type_size(t.data_type);
      for (int dstId = 0; dstId < config.num_parts(); dstId ++) {
        Domain dstR = op->get_input_tensor_shape(config, j, dstId);
        for (int srcId = 0; srcId < pre_config.num_parts(); srcId ++) {
          Domain srcR = pre_op->get_output_tensor_shape(pre_config, t.owner_idx, srcId);
          if (dstR.intersection(srcR).get_volume() > 0) {
            // Forward dependency
            {
              SimTask* dstT = task_manager->get_forward_task(op, dstId);
              SimTask* srcT = task_manager->get_forward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(srcT, dstT, dstR.intersection(srcR).get_volume() * element_size);
            }
            // Backward dependency
            if (comp_mode == COMP_MODE_TRAINING) {
              SimTask* dstT = task_manager->get_backward_task(op, dstId);
              SimTask* srcT = task_manager->get_backward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(dstT, srcT, dstR.intersection(srcR).get_volume() * element_size);
            }
          }
        }
      }
    }
  }
#ifdef FF_USE_NCCL
  // Do nothing since we will calculate NCCL cost at the end
  // fprintf(stderr, "USING NCCL\n");
#else
  // Step 2.5: add finals tasks for each compute device to capture the returning comm tasks
  // from parameter servers
  std::vector<SimTask*> finals;
  for (int d = 0; d < machine->get_num_gpus(); d++) {
    SimTask* t = task_manager->new_barrier_task();
    t->device = machine->get_gpu(d);
    t->mem = machine->get_gpu_fb_mem(d);
    t->run_time = 0;
    finals.push_back(t);
  }

  if (model->config.search_overlap_backward_update && comp_mode == COMP_MODE_TRAINING) {
    // Step 3a: consider backpropagation and weight update are overlapped
    for (int l = model->layers.size()-1; l >= 0; l--) {
      Op* op = model->layers[l];
      size_t element_size = data_type_size(DT_FLOAT); // assume all weights have float elements
      ParallelConfig pc = global.find(op)->second;
      for (int j = 0; j < op->numWeights; j++) {
        std::set<int> synched;
        for (int firstId = 0; firstId < pc.num_parts(); firstId++)
          if (synched.find(firstId) == synched.end()) {
            synched.insert(firstId);
            Domain firstR = op->get_weight_tensor_shape(pc, j, firstId);
            // Add a compute task for parameter update
            SimTask* updateT = task_manager->new_update_task();
            updateT->device = machine->get_gpu(pc.device_ids[firstId]);
            updateT->mem = machine->get_gpu_fb_mem(pc.device_ids[firstId]);
            // TODO add parameter synchronization time
            updateT->run_time = 0.0f; // Assume update task takes no time
            for (int nextId = firstId+1; nextId < pc.num_parts(); nextId++) {
              Domain nextR = op->get_weight_tensor_shape(pc, j, nextId);
              if (firstR.intersection(nextR).get_volume() > 0) {
                // Assert all or nothing:
                // The two weights must be fully overlapped or not at all
                assert(firstR == nextR);
                assert(synched.find(nextId) == synched.end());
                synched.insert(nextId);
                // Add comm. tasks from backT to updateT
                SimTask* backT = task_manager->get_backward_task(op, nextId);
                add_task_dependencies_with_xfer(backT, updateT, firstR.get_volume() * element_size);
                // Add comm. tasks from updateT to finalT
                SimTask* finalT = finals[backT->device->device_id];
                add_task_dependencies_with_xfer(updateT, finalT, firstR.get_volume() * element_size);
              }
            }
          }
      }
    }
  } else if (comp_mode == COMP_MODE_TRAINING) {
    // Step 3b: Bulk Synchronous Model
    // Add a per-device barrier before weight update
    std::vector<SimTask*> barriers;
    for (int d = 0; d < machine->get_num_gpus(); d++) {
      SimTask* t = task_manager->new_barrier_task();
      t->device = machine->get_gpu(d);
      t->mem = machine->get_gpu_fb_mem(d);
      t->run_time = 0;
      barriers.push_back(t);
    }
    for (size_t l = 0; l < model->layers.size(); l++) {
      Op* op = model->layers[l];
      ParallelConfig pc = global.find(op)->second;
      for (int j = 0; j < pc.num_parts(); j++) {
        SimTask* backT = task_manager->get_backward_task(op, j);
        backT->add_next_task(barriers[backT->device->device_id]);
      }
    }
    for (size_t l = 0; l < model->layers.size(); l++) {
      Op* op = model->layers[l];
      ParallelConfig pc = global.find(op)->second;
      size_t element_size = data_type_size(DT_FLOAT); // assume all weights have float elements
      for (int j = 0; j < op->numWeights; j++) {
        std::set<int> synched;
        for (int firstId = 0; firstId < pc.num_parts(); firstId++)
          if (synched.find(firstId) == synched.end()) {
            synched.insert(firstId);
            Domain firstR = op->get_weight_tensor_shape(pc, j, firstId);
            // Add a compute task for parameter update
            SimTask* updateT = task_manager->new_update_task();
            updateT->device = machine->get_gpu(pc.device_ids[firstId]);
            updateT->mem = machine->get_gpu_fb_mem(pc.device_ids[firstId]);
            updateT->run_time = 0.0f; // Assume update task takes no time
            barriers[updateT->device->device_id]->add_next_task(updateT);
            for (int nextId = firstId+1; nextId < pc.num_parts(); nextId++) {
              Domain nextR = op->get_weight_tensor_shape(pc, j, nextId);
              if (firstR.intersection(nextR).get_volume() > 0) {
                // Assert all or nothing:
                // The two weights must be fully overlapped or not at all
                assert(firstR == nextR);
                assert(synched.find(nextId) == synched.end());
                synched.insert(nextId);
                SimTask* backT = task_manager->get_backward_task(op, nextId);
                assert(backT->device->device_id == pc.device_ids[nextId]);
                SimTask* barrierT = barriers[backT->device->device_id];
                // Add comm. tasks from barrierT to updateT
                add_task_dependencies_with_xfer(barrierT, updateT, firstR.get_volume() * element_size);
                // Add comm. tasks from updateT to finalT
                SimTask* finalT = finals[backT->device->device_id];
                add_task_dependencies_with_xfer(updateT, finalT, firstR.get_volume() * element_size);
              }
            }
          }
      }
    }
  } else {
    assert(comp_mode == COMP_MODE_INFERENCE);
  }
#endif
  // Step 4: add ready tasks into ready_queue
  std::priority_queue<SimTask*, std::vector<SimTask*>, SimTaskCompare> ready_queue;
  for (size_t i = 0; i < task_manager->global_task_id; i++)
    if (task_manager->tasks[i]->counter == 0)
      ready_queue.push(task_manager->tasks[i]);
  // Step 5: perform simulation
  float sim_time = 0.0f;
  std::map<Device*, float> device_times;
  size_t idx = 0;
  DotFile<SimTask *> taskGraph;
  bool export_taskgraph = (export_file_name != "");
  if (export_taskgraph) {
    taskGraph.set_filename(export_file_name);
  }
  while (!ready_queue.empty()) {
    // Find the task with the earliest start time
    SimTask* cur_task = ready_queue.top();
    ready_queue.pop();
    float ready_time = 0;
    if (device_times.find(cur_task->device) != device_times.end()) {
      ready_time = device_times[cur_task->device];
    }
    float start_time = std::max(ready_time, cur_task->ready_time);
    float end_time = start_time + cur_task->run_time;
    device_times[cur_task->device] = end_time;
    if (export_taskgraph) {
      std::map<std::string, std::string> nodeAttrs;
      std::ostringstream label;
      label << "\"{ ";
      if (!(cur_task->name).empty()) {
        label << cur_task->name << " | ";
      }
      label << cur_task->get_type_str() << " | ";
      label << "{ " << start_time << " | " << end_time << " }";
      label << " }\"";
      nodeAttrs["label"] = label.str();
      nodeAttrs["shape"] = "record";
      taskGraph.add_node(cur_task, nodeAttrs);
    }
  #ifdef DEBUG_PRINT
    // printf("task[%lu] type(%d) run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
    //       idx, cur_task->type, cur_task->run_time, ready_time, start_time, (cur_task->device->name).c_str());
  #endif
    if (end_time > sim_time)
      sim_time = end_time;
    for (size_t i = 0; i < cur_task->next_tasks.size(); i++) {
      SimTask* next = cur_task->next_tasks[i];
      if (export_taskgraph) {
        taskGraph.add_edge(cur_task, next);
      }
      next->ready_time = std::max(next->ready_time, end_time);
      next->counter --;
      if (next->counter == 0) {
        ready_queue.push(next);
      }
    }
    idx++;
  }
  if (export_taskgraph) {
    taskGraph.close();
  }
  // Assert all tasks were processed
  assert(idx == task_manager->global_task_id);
#ifdef FF_USE_NCCL
  if (comp_mode == COMP_MODE_TRAINING) {
    for (size_t l = 0; l < model->layers.size(); l++) {
      Op* op = model->layers[l];
      size_t element_size = data_type_size(DT_FLOAT); // assume all weights have float elements
      ParallelConfig pc = global.find(op)->second;
      // Since all NCCL calls are blocking, we can add the NCCL cost
      // sequentially 
      for (int j = 0; j < op->numWeights; j++) {
        std::set<int> synched;
        for (int firstId = 0; firstId < pc.num_parts(); firstId++)
          if (synched.find(firstId) == synched.end()) {
            synched.insert(firstId);
            Domain firstR = op->get_weight_tensor_shape(pc, j, firstId);
            Device* firstDevice = machine->get_gpu(pc.device_ids[firstId]);
            float nccl_time = 0.0f;
            for (int nextId = firstId+1; nextId < pc.num_parts(); nextId++) {
              Domain nextR = op->get_weight_tensor_shape(pc, j, nextId);
              if (firstR.intersection(nextR).get_volume() > 0) {
                // Assert all or nothing:
                // The two weights must be fully overlapped or not at all
                assert(firstR == nextR);
                assert(synched.find(nextId) == synched.end());
                synched.insert(nextId);
                Device* nextDevice = machine->get_gpu(pc.device_ids[nextId]);
                // Compute the bandwidth between firstDevice/nextDevice
                float bandwidth = 0.0f;
                if (firstDevice->node_id == nextDevice->node_id) {
                  bandwidth = machine->get_intra_node_gpu_bandwidth();
                } else {
                  bandwidth = machine->get_inter_node_gpu_bandwidth();
                }
                nccl_time = std::max(nccl_time, (float)firstR.get_volume() * element_size / bandwidth);
              }
            }
            // Add ncclTime to sim_time given nccl calls are blocking
            sim_time += nccl_time;
          }
      }
    }
  } else {
    assert(comp_mode == COMP_MODE_INFERENCE);
  }
#endif
  // Step 6: add penalty to strategies that exceed the memory limits on devices
  std::vector<size_t> gpu_mem_usage(machine->get_num_gpus(), 0);
  float memory_penalty = 0.0f;
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    CostMetrics cost_metrics = measure_operator_cost(op, config);
    size_t memory_requirement = cost_metrics.memory_requirement;
    for (int j = 0; j < config.num_parts(); j++) {
      gpu_mem_usage[config.device_ids[j]] += memory_requirement;
    }
  }
  if (export_file_name != "") {  
    for (int i = 0; i < machine->get_num_gpus(); i++) {
        printf("Before penalty, dev id %d, usage %zu \n", i, gpu_mem_usage[i]); 
    }
  }
  // Penalize the total runtiem by 1ms if we exceed the memory budget by 1MB
  for (int i = 0; i < machine->get_num_gpus(); i++) {
    MemDevice* gpu_fb_mem = machine->get_gpu_fb_mem(i);
    if (gpu_mem_usage[i] > gpu_fb_mem->capacity and gpu_fb_mem->capacity >= 0)
      memory_penalty += (gpu_mem_usage[i] - gpu_fb_mem->capacity) * 1e-6;
  }
  //if (memory_penalty > 0.0f)
  //  printf("Memory penalty = %.4lf ms\n", memory_penalty);
  return sim_time + memory_penalty;
}


float LogicalTaskgraphBasedSimulator::simulate_runtime(
                                  const FFModel* model,
                                  const std::map<Op*, ParallelConfig>& global,
                                  CompMode comp_mode,
                                  std::string const &export_file_name) 
{
  // printf("%s\n", machine->to_string().c_str());
  task_manager->reset();
  if (l1optimizer)
    l1optimizer->reset();
  // Step 1: register forward and backward tasks
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    CostMetrics cost_metrics = measure_operator_cost(op, config);
    float forward_time = cost_metrics.forward_time;
    float backward_time = cost_metrics.backward_time;
    SimTask *ar_task = nullptr;
    for (int j = 0; j < config.num_parts(); j++) {
      SimTask* task1 = task_manager->new_forward_task(op, j);
      task1->device = machine->get_gpu(config.device_ids[j]);
      task1->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
      task1->run_time = forward_time;
      if (l1optimizer) 
        l1optimizer->task_added(task1);
      if (comp_mode == COMP_MODE_TRAINING) {
        SimTask* task2 = task_manager->new_backward_task(op, j);
        task2->device = machine->get_gpu(config.device_ids[j]);
        task2->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
        task2->run_time = backward_time;
        task1->add_next_task(task2);
        if (l1optimizer) 
          l1optimizer->task_added(task2);
        
      }
    }
  }

  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    size_t element_size = data_type_size(DT_FLOAT);
    // NER step: add allreduce task after backward propogation
    for (int j = 0; j < op->numWeights; j++) {
      std::set<int> synched;
      std::vector<int> node_ids;
      for (int firstId = 0; firstId < config.num_parts(); firstId++) {
        if (synched.find(firstId) == synched.end()) {
          synched.insert(firstId);
          Domain firstR = op->get_weight_tensor_shape(config, j, firstId);
          size_t xfer_size = firstR.get_volume() * element_size;
          node_ids.push_back(config.device_ids[firstId]);
          for (int nextId = firstId+1; nextId < config.num_parts(); nextId++) {
            Domain nextR = op->get_weight_tensor_shape(config, j, nextId);
            if (firstR.intersection(nextR).get_volume() > 0) {
              // Assert all or nothing:
              // The two weights must be fully overlapped or not at all
              assert(firstR == nextR);
              assert(synched.find(nextId) == synched.end());
              synched.insert(nextId);
              node_ids.push_back(config.device_ids[nextId]);
            }
          }
          
          SimTask* ar_task = task_manager->new_allreduce_task(op, node_ids, xfer_size);
          if (l1optimizer) 
            l1optimizer->task_added(ar_task);
          for (int dstId = 0; dstId < config.num_parts(); dstId ++) {
            task_manager->get_backward_task(op, dstId)->add_next_task(ar_task);
          }
        }
      
      }
    }
  }
        

  // Step 2: insert dependencies and comm. tasks before compute tasks
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    for (int j = 0; j < op->numInputs; j++) {
      Tensor t = op->inputs[j];
      Op* pre_op = t.owner_op;
      if (pre_op == NULL)
        continue;
      ParallelConfig pre_config = global.find(pre_op)->second;
      size_t element_size = data_type_size(t.data_type);
      for (int dstId = 0; dstId < config.num_parts(); dstId ++) {
        Domain dstR = op->get_input_tensor_shape(config, j, dstId);
        for (int srcId = 0; srcId < pre_config.num_parts(); srcId ++) {
          Domain srcR = pre_op->get_output_tensor_shape(pre_config, t.owner_idx, srcId);
          if (dstR.intersection(srcR).get_volume() > 0) {
            // Forward dependency
            {
              SimTask* dstT = task_manager->get_forward_task(op, dstId);
              SimTask* srcT = task_manager->get_forward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(srcT, dstT, dstR.intersection(srcR).get_volume() * element_size);
            }
            // Backward dependency
            if (comp_mode == COMP_MODE_TRAINING) {
              SimTask* dstT = task_manager->get_backward_task(op, dstId);
              SimTask* srcT = task_manager->get_backward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(dstT, srcT, dstR.intersection(srcR).get_volume() * element_size);
            }
          }
        }
      }
    }
  }
  
  // Step 4: add ready tasks into ready_queue
  std::priority_queue<SimTask*, std::vector<SimTask*>, SimTaskCompare> ready_queue;
  for (size_t i = 0; i < task_manager->global_task_id; i++)
    if (task_manager->tasks[i]->counter == 0)
      ready_queue.push(task_manager->tasks[i]);

  // Step 5: perform simulation

  float sim_time = 0.0f;
  std::map<Device*, float> device_times;
  // map<Device*, SimTask*> device_schedule;
  size_t idx = 0;
  while (!ready_queue.empty()) {
    // Find the task with the earliest start time
    SimTask* cur_task = ready_queue.top();
    ready_queue.pop();
    float ready_time = 0;
    float end_time;
    if (device_times.find(cur_task->device) != device_times.end()) {
      ready_time = device_times[cur_task->device];
    }
    float start_time = std::max(ready_time, cur_task->ready_time);
    if (cur_task->type == SimTask::TASK_NOMINAL_COMM) {
      end_time = route_transfer(cur_task, start_time, device_times);
    }
    else if (cur_task->type == SimTask::TASK_ALLREDUCE) {
      expand_allreduce(cur_task, start_time, ready_queue);
      idx++;
      continue;
    }
    else {
      end_time = start_time + cur_task->run_time;
      device_times[cur_task->device] = end_time;
    }

#ifdef DEBUG_PRINT
    printf("task[%lu/%lu] type(%d) run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          idx, task_manager->global_task_id, cur_task->type, cur_task->run_time, ready_time, start_time, (cur_task->device->name).c_str());
#endif

    if (end_time > sim_time) {
      sim_time = end_time;
    }

    for (size_t i = 0; i < cur_task->next_tasks.size(); i++) {
      SimTask* next = cur_task->next_tasks[i];
      // next->ready_time = max(next->ready_time, end_time);
      if (end_time > next->ready_time) {
        next->ready_time = end_time;
        // next->prev = t;
      }
      next->counter--;
      if (next->counter == 0) {
        ready_queue.push(next);
      }
    }
    idx++;
  }
  assert(idx == task_manager->global_task_id);
  
  // Step 6: add penalty to strategies that exceed the memory limits on devices
  // std::vector<size_t> gpu_mem_usage(machine->get_num_gpus(), 0);
  // float memory_penalty = 0.0f;
  // for (size_t l = 0; l < model->layers.size(); l++) {
  //   Op* op = model->layers[l];
  //   ParallelConfig config = global.find(op)->second;
  //   CostMetrics cost_metrics = measure_operator_cost(op, config);
  //   size_t memory_requirement = cost_metrics.memory_requirement;
  //   for (int j = 0; j < config.num_parts(); j++) {
  //     gpu_mem_usage[config.device_ids[j]] += memory_requirement;
  //   }
  // }
  // if (export_file_name != "") {  
  //   for (int i = 0; i < machine->get_num_gpus(); i++) {
  //       printf("Before penalty, dev id %d, usage %zu \n", i, gpu_mem_usage[i]); 
  //   }
  // }
  // // Penalize the total runtiem by 1ms if we exceed the memory budget by 1MB
  // for (int i = 0; i < machine->get_num_gpus(); i++) {
  //   MemDevice* gpu_fb_mem = machine->get_gpu_fb_mem(i);
  //   if (gpu_mem_usage[i] > gpu_fb_mem->capacity and gpu_fb_mem->capacity >= 0)
  //     memory_penalty += (gpu_mem_usage[i] - gpu_fb_mem->capacity) * 1e-6;
  // }
  //if (memory_penalty > 0.0f)
  //  printf("Memory penalty = %.4lf ms\n", memory_penalty);
  searlize_logical_taskgraph(model, "exp");
  return sim_time;//  + memory_penalty;
      
}

float LogicalTaskgraphBasedSimulator::simulate_runtime(const FFModel* model,
                                  const std::map<Op*, ParallelConfig>& global,
                                  CompMode comp_mode)
{
  return this->simulate_runtime(model, global, comp_mode, "");
}


float LogicalTaskgraphBasedSimulator::route_transfer(SimTask * transfer_task, 
                              float start_time,
                              std::map<Device*, float> &device_times) {
  std::vector<CommDevice *> route = 
    static_cast<NominalCommDevice*>(transfer_task->device)->expand_to_physical();

  float curr_task_start_time; 
  float curr_task_finish_time; 
  float curr_task_run_time = 0; 
  float curr_task_ready_time = transfer_task->ready_time; 
  float xfer_size = transfer_task->xfer_size;

  float final_start_time = 0;
  float final_finish_time = 0;

  SimTask * info_holder = new SimTask();
  info_holder->type = SimTask::TASK_COMM;

  for (unsigned int i = 0; i < route.size(); i++) {
    CommDevice * latency_task_device = route[i];
    float latency_task_run_time = machine->get_inter_node_gpu_latency();
    float latency_task_ready_time; 
    float latency_task_start_time; 
    if (i == 0) {
      latency_task_ready_time = curr_task_ready_time + curr_task_run_time;
      latency_task_start_time = std::max(device_times[latency_task_device], latency_task_ready_time);
      final_start_time = latency_task_start_time;
    }
    else {
      latency_task_ready_time = curr_task_finish_time;
      latency_task_start_time = std::max(device_times[latency_task_device], latency_task_ready_time);
    }
    float latency_task_finish_time = latency_task_start_time + latency_task_run_time;
    device_times[latency_task_device] = latency_task_finish_time;
    float dram_to_dram_run_time = xfer_size / latency_task_device->bandwidth;

    float dram_to_dram_start_time = latency_task_finish_time;
    float dram_to_dram_finish_time = dram_to_dram_start_time + dram_to_dram_run_time;
    device_times[latency_task_device] = dram_to_dram_finish_time;

    if (dram_to_dram_finish_time > final_finish_time) {
      final_finish_time = dram_to_dram_finish_time;
    }

    curr_task_ready_time = latency_task_ready_time;
    curr_task_start_time = latency_task_start_time;
    curr_task_finish_time = latency_task_finish_time;
    curr_task_run_time = latency_task_run_time;
    
#ifdef DEBUG_PRINT
    printf("\texpand: route[%u] run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          i, curr_task_run_time, curr_task_ready_time, curr_task_start_time, (latency_task_device->name).c_str());
    printf("\t\td2d: run_time(%.4lf) start_time(%.4lf) device(%s)\n",
          dram_to_dram_run_time, dram_to_dram_start_time, (latency_task_device->name).c_str());
#endif

    info_holder->device = latency_task_device;
    info_holder->run_time = dram_to_dram_run_time;
    info_holder->xfer_size = xfer_size;
    // info_holder->from_dev = CommDevice::get_from_dev(latency_task_device->device_id, mac);
    // info_holder->to_dev = 

    if (l1optimizer) 
      l1optimizer->task_added(info_holder);

  }

  delete info_holder;

  transfer_task->run_time = final_finish_time - final_start_time;
  return final_finish_time;
}

void LogicalTaskgraphBasedSimulator::expand_allreduce(SimTask * allreduce_task,
                                 float start_time,
                                 std::priority_queue<SimTask*, std::vector<SimTask*>, SimTaskCompare>& ready_queue) {

  int n_participants = allreduce_task->next_tasks.size();
  
  SimTask * final_task = new_update_task_unrecorded();

#ifdef FF_USE_NCCL
  // recall that next_task stores node group in this case
  final_task->device = machine->get_gpu(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
  MemDevice * src_mem = machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
  MemDevice * dst_mem;

  for (int i = 0; i < n_participants; i++) {
    dst_mem = machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[(i+1)%n_participants]));
    std::vector<CommDevice *> path = machine->get_comm_path(src_mem, dst_mem);
    for (CommDevice * d: path) {
      SimTask* task = new_comm_task_unrecorded();
      task->device = d;
      task->run_time = 0;
      task->ready_time = allreduce_task->ready_time;
      task->xfer_size = (2.0 * (n_participants-1))/n_participants * allreduce_task->xfer_size;
      task->add_next_task(final_task);
      ready_queue.push(task);
      if (l1optimizer)
        l1optimizer->task_added(task);
    }
    src_mem = dst_mem;
  }
  if (final_task->counter == 0) {
    final_task->ready_time = allreduce_task->ready_time;
    ready_queue.push(final_task);
  }
#else
  // assume parameter server in this case
  MemDevice * leader_mem = machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
  MemDevice * worker_mem;
  SimTask * ps_update_task = new_update_task_unrecorded();
  ps_update_task->device = machine->get_gpu(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
  final_task->device = machine->get_gpu(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[0]));
  ps_update_task->add_next_task(final_task);

  // ps gather
  for (int i = 0; i < n_participants; i++) {
    worker_mem = machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[i]));
    std::vector<CommDevice *> path = machine->get_comm_path(worker_mem, leader_mem);
    for (CommDevice * d: path) {
      SimTask* task = new_comm_task_unrecorded();
      task->device = d;
      task->run_time = 0;
      task->ready_time = allreduce_task->ready_time;
      task->xfer_size = allreduce_task->xfer_size;
      task->add_next_task(ps_update_task);
      ready_queue.push(task);
      if (l1optimizer)
        l1optimizer->task_added(task);
    }
  }

  // scatter
  for (int i = 0; i < n_participants; i++) {
    worker_mem = machine->get_gpu_fb_mem(reinterpret_cast<uint64_t>(allreduce_task->next_tasks[i]));
    std::vector<CommDevice *> path = machine->get_comm_path(leader_mem, worker_mem);
    for (CommDevice * d: path) {
      SimTask* task = new_comm_task_unrecorded();
      task->device = d;
      task->run_time = 0;
      task->ready_time = allreduce_task->ready_time;
      task->xfer_size = allreduce_task->xfer_size;
      ps_update_task->add_next_task(task);
      task->add_next_task(final_task);
      if (l1optimizer)
        l1optimizer->task_added(task);
    }
  }

  if (ps_update_task->counter == 0) {
    assert(final_task->counter == 1);
    ps_update_task->ready_time = allreduce_task->ready_time;
    ready_queue.push(ps_update_task);
  }

#endif

}

SimTask* LogicalTaskgraphBasedSimulator::new_comm_task_unrecorded() {
  SimTask* task = task_manager->new_task();
  task->type = SimTask::TASK_NOMINAL_COMM;
  task->store = false;
  return task;
}

SimTask* LogicalTaskgraphBasedSimulator::new_update_task_unrecorded() {
  SimTask* task = task_manager->new_task();
  task->type = SimTask::TASK_UPDATE;
  task->store = false;
  return task;
}

bool LogicalTaskgraphBasedSimulator::searlize_logical_taskgraph(const FFModel* model, std::string const &export_file_name) {
  flatbuffers::FlatBufferBuilder builder(262144);
  get_taskgraph_flatbuf(model, builder);
  std::ofstream ofs(export_file_name, std::ofstream::binary);
  if (!ofs.is_open()) return false;
  ofs.write((const char *) builder.GetBufferPointer(), (size_t)builder.GetSize());
  return !ofs.bad();
  // flatbuffers::SaveFile(export_file_name.c_str(),
  //                       (const char *) builder.GetBufferPointer(),
  //                       (size_t) builder.GetSize(), true);
  // return;
}

void LogicalTaskgraphBasedSimulator::get_taskgraph_flatbuf(const FFModel* model, flatbuffers::FlatBufferBuilder &builder) 
{
  builder.Clear();

  // Store topology
  // flatbuffers::FlatBufferBuilder builder = flatbuffers::FlatBufferBuilder();
  NetworkedMachineModel *nm = static_cast<NetworkedMachineModel*>(machine);
  size_t total_devs = nm->get_total_devs();
  std::vector<flatbuffers::Offset<FlatBufTaskGraph::Connection>> conns_v = 
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Connection>>();
  for (size_t i = 0; i < nm->get_total_devs(); i++) {
    for (size_t j = 0; j < i; j++) {
      size_t nlink;
      if ((nlink = nm->get_conn_matrix()[i * total_devs + j]) > 0) {
        conns_v.emplace_back(FlatBufTaskGraph::CreateConnection(builder, i, j, nlink));
      }
    }
  }
  auto conns = builder.CreateVector(conns_v);

  // store operators
  // builder.Clear();
  std::vector<flatbuffers::Offset<FlatBufTaskGraph::Operator>> op_v = 
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Operator>>();
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    auto opname = builder.CreateString(op->name);
    op_v.emplace_back(FlatBufTaskGraph::CreateOperator(builder, 
      reinterpret_cast<uint64_t>(op), (int)op->op_type, opname));
  }
  auto ops = builder.CreateVector(op_v);

  // store tasks
  // builder.Clear();
  std::vector<flatbuffers::Offset<FlatBufTaskGraph::Task>> task_v = 
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Task>>();
  // change: since there is no universal storage of device, creat a set of
  // all devices for the next entry
  std::unordered_set<Device *> devices;
  for (size_t i = 0; i < task_manager->global_task_id; i++) {
    SimTask * curr = task_manager->tasks[i];
    if (curr->store) {
      FlatBufTaskGraph::SimTaskType tasktype;
      uint64_t taskid = reinterpret_cast<uint64_t>(curr);
      std::vector<uint64_t> nexttasks = std::vector<uint64_t>();
      for (SimTask *t: curr->next_tasks) {
        nexttasks.push_back(reinterpret_cast<uint64_t>(t));
      }
      auto ntv = builder.CreateVector(nexttasks);
      switch (curr->type) {
      case SimTask::TASK_FORWARD:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_FORWARD;
      break;
      case SimTask::TASK_BACKWARD:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_BACKWARD;
      break;
      case SimTask::TASK_UPDATE:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_UPDATE;
      break;
      case SimTask::TASK_BARRIER:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_BARRIER;
      break;
      case SimTask::TASK_COMM:
        assert("Logical task graph shouldn't contain TASK_COMM!" && false);
      break;
      case SimTask::TASK_NOMINAL_COMM:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_NOMINAL_COMM;
      break;
      case SimTask::TASK_ALLREDUCE:
        tasktype = FlatBufTaskGraph::SimTaskType_TASK_ALLREDUCE;
      break;
      }
      task_v.emplace_back(FlatBufTaskGraph::CreateTask(
        builder,
        tasktype,
        taskid, 
        reinterpret_cast<uint64_t>(curr->device),
        curr->run_time,
        curr->xfer_size,
        ntv
      ));
    }
    if (curr->device)
      devices.insert(curr->device);
  }
  auto tasks = builder.CreateVector(task_v);

  // devices
  // builder.Clear();
  std::vector<flatbuffers::Offset<FlatBufTaskGraph::Device>> dev_v = 
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Device>>();
  for (Device *curr: devices) {
    FlatBufTaskGraph::DeviceType type;
    uint64_t deviceid = reinterpret_cast<uint64_t>(curr);
    CommDevice * comm_dev;
    switch (curr->type) {
    case Device::DEVICE_COMP: 
      dev_v.emplace_back(FlatBufTaskGraph::CreateDevice(
        builder, 
        reinterpret_cast<CompDevice*>(curr)->comp_type == CompDevice::LOC_PROC 
          ? FlatBufTaskGraph::DeviceType_DEVICE_COMP_CPU
          : FlatBufTaskGraph::DeviceType_DEVICE_COMP_GPU,
        deviceid, curr->node_id, curr->device_id, 0
      ));
    break;
    case Device::DEVICE_COMM: 
      comm_dev = reinterpret_cast<CommDevice*>(curr);
      switch (comm_dev->comm_type) {
      case CommDevice::MEMBUS_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_MEMBUS_COMM;
      break;
      case CommDevice::UPI_IN_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_UPI_IN_COMM;
      break;
      case CommDevice::UPI_OUT_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_UPI_OUT_COMM;
      break;
      case CommDevice::NIC_IN_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_NIC_IN_COMM;
      break;
      case CommDevice::NIC_OUT_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_NIC_OUT_COMM;
      break;
      case CommDevice::PCI_TO_HOST_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_PCI_TO_HOST_COMM;
      break;
      case CommDevice::PCI_TO_DEV_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_PCI_TO_DEV_COMM;
      break;
      case CommDevice::NVLINK_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_NVLINK_COMM;
      break;
      case CommDevice::NW_COMM:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_NW_COMM;
      break;
      case CommDevice::NW_NOMINAL:
        type = FlatBufTaskGraph::DeviceType_DEVICE_COMM_NW_NOMINAL;
      break;
      }
      dev_v.emplace_back(FlatBufTaskGraph::CreateDevice(
        builder, 
        type,
        deviceid, curr->node_id, curr->device_id, comm_dev->bandwidth
      ));
    break;
    case Device::DEVICE_MEM: 
      assert("Shouldn't store a memory device to taskgraph!" && false);
    }
  }
  auto devs = builder.CreateVector(dev_v);

  // routes
  // builder.Clear();
  std::vector<flatbuffers::Offset<FlatBufTaskGraph::Route>> route_v = 
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Route>>();
  for (auto ncd: nm->get_nomm_comm_devs()) {
    std::vector<flatbuffers::Offset<FlatBufTaskGraph::Path>> path_v = 
      std::vector<flatbuffers::Offset<FlatBufTaskGraph::Path>>();
    const EcmpRoutes& physical_routes = ncd.second->get_all_routes();
    for (size_t i = 0; i < physical_routes.first.size(); i++) {
      std::vector<uint32_t> hops_v = std::vector<uint32_t>();
      for (CommDevice * c: physical_routes.second[i]) {
        hops_v.push_back(c->node_id);
      }
      auto hops = builder.CreateVector(hops_v);
      auto path = FlatBufTaskGraph::CreatePath(builder, hops, physical_routes.first[i]);
      path_v.push_back(path);
    }
    auto paths = builder.CreateVector(path_v);
    route_v.push_back(FlatBufTaskGraph::CreateRoute(
      builder, 
      ncd.second->device_id / nm->get_total_devs(),
      ncd.second->device_id % nm->get_total_devs(),
      paths
    ));
  }

  FlatBufTaskGraph::TaskGraphBuilder tg_builder = FlatBufTaskGraph::TaskGraphBuilder(builder);

  tg_builder.add_ngpupernode(machine->get_num_gpus()/ machine->get_num_nodes());
  tg_builder.add_nnode(machine->get_num_nodes());
  tg_builder.add_nswitch(nm->get_num_switches());
  tg_builder.add_intergpubw(machine->get_intra_node_gpu_bandwidth());
  tg_builder.add_drambw(32 * 1024 * 1024.0f); // PCIE gen 4
  tg_builder.add_netbw(machine->get_inter_node_gpu_bandwidth());
  tg_builder.add_conn(conns);
  tg_builder.add_ops(ops);
  tg_builder.add_tasks(tasks);
  tg_builder.add_devices(devs);

  auto ftg = tg_builder.Finish();
  builder.Finish(ftg);
}

float DLSSchedulerBasedSimulator::simulate_runtime(const FFModel* model,
      const std::map<Op*, ParallelConfig>& global,
      CompMode comp_mode) 
{
  return this->simulate_runtime(model, global, comp_mode, "");
}

float DLSSchedulerBasedSimulator::simulate_runtime(const FFModel* model,
      const std::map<Op*, ParallelConfig>& global,
      CompMode comp_mode,
      std::string const &export_file_name)
{
  assert(comp_mode == CompMode::COMP_MODE_TRAINING);

  // printf("%s\n", machine->to_string().c_str());
  task_manager->reset();
  topofinder->reset();

  // Step 1: register forward and backward tasks. Add BP to a separate taskdag
  std::unordered_set<SimTask*> bp_tasks;
  DLSTaskDag bp_taskdag;
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    CostMetrics cost_metrics = measure_operator_cost(op, config);
    float forward_time = cost_metrics.forward_time;
    float backward_time = cost_metrics.backward_time;
    for (int j = 0; j < config.num_parts(); j++) {
      SimTask* task1 = task_manager->new_forward_task(op, j);
      // task1->device = machine->get_gpu(config.device_ids[j]);
      // task1->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
      task1->run_time = forward_time;
      // topofinder->fw_task_added(task1);
      if (comp_mode == COMP_MODE_TRAINING) {
        SimTask* task2 = task_manager->new_backward_task(op, j);
        // task2->device = machine->get_gpu(config.device_ids[j]);
        // task2->mem = machine->get_gpu_fb_mem(config.device_ids[j]);
        task2->run_time = backward_time;
        bp_tasks.insert(task2);
        task1->add_next_task(task2);
        // topofinder->bw_task_added(task2);
      }
    }
  }

  // acquire backward task dag
  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = global.find(op)->second;
    for (int j = 0; j < op->numInputs; j++) {
      Tensor t = op->inputs[j];
      Op* pre_op = t.owner_op;
      if (pre_op == NULL)
        continue;
      ParallelConfig pre_config = global.find(pre_op)->second;
      size_t element_size = data_type_size(t.data_type);
      for (int dstId = 0; dstId < config.num_parts(); dstId ++) {
        Domain dstR = op->get_input_tensor_shape(config, j, dstId);
        for (int srcId = 0; srcId < pre_config.num_parts(); srcId ++) {
          Domain srcR = pre_op->get_output_tensor_shape(pre_config, t.owner_idx, srcId);
          if (dstR.intersection(srcR).get_volume() > 0) {
            // Backward dependency
            SimTask* dstT = task_manager->get_backward_task(op, dstId);
            SimTask* srcT = task_manager->get_backward_task(pre_op, srcId);
            add_task_dependencies_with_xfer_sch(dstT, srcT, dstR.intersection(srcR).get_volume() * element_size, bp_taskdag);
          }
        }
      }
    }
  } 
  dls_schedule(bp_tasks, bp_taskdag);
  // now we should have device placement. get the new PC
  std::map<Op*, ParallelConfig> dev_placement = get_device_placements(bp_tasks, global);

  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = dev_placement.find(op)->second;
    // NER step: add allreduce task after backward propogation
    size_t element_size = data_type_size(DT_FLOAT); // assume all weights have float elements
    ParallelConfig pc = dev_placement.find(op)->second;
    for (int j = 0; j < op->numWeights; j++) {
      std::set<int> synched;
      std::vector<int> node_ids;
      for (int firstId = 0; firstId < pc.num_parts(); firstId++) {
        if (synched.find(firstId) == synched.end()) {
          synched.insert(firstId);
          Domain firstR = op->get_weight_tensor_shape(pc, j, firstId);
          size_t xfer_size = firstR.get_volume() * element_size;
          node_ids.push_back(pc.device_ids[firstId]);
          for (int nextId = firstId+1; nextId < pc.num_parts(); nextId++) {
            Domain nextR = op->get_weight_tensor_shape(pc, j, nextId);
            if (firstR.intersection(nextR).get_volume() > 0) {
              // Assert all or nothing:
              // The two weights must be fully overlapped or not at all
              assert(firstR == nextR);
              assert(synched.find(nextId) == synched.end());
              synched.insert(nextId);
              node_ids.push_back(pc.device_ids[nextId]);
            }
          }
          SimTask *ar_task = task_manager->new_allreduce_task(op, node_ids, xfer_size);
          for (int k = 0; k < config.num_parts(); k++) {
            task_manager->get_backward_task(op, k)->add_next_task(ar_task);
          }
          topofinder->ar_task_added(ar_task);
        }
      }
    }
  }

  for (size_t l = 0; l < model->layers.size(); l++) {
    Op* op = model->layers[l];
    ParallelConfig config = dev_placement.find(op)->second;
    for (int j = 0; j < op->numInputs; j++) {
      Tensor t = op->inputs[j];
      Op* pre_op = t.owner_op;
      if (pre_op == NULL)
        continue;
      ParallelConfig pre_config = dev_placement.find(pre_op)->second;
      size_t element_size = data_type_size(t.data_type);
      for (int dstId = 0; dstId < config.num_parts(); dstId ++) {
        Domain dstR = op->get_input_tensor_shape(config, j, dstId);
        for (int srcId = 0; srcId < pre_config.num_parts(); srcId ++) {
          Domain srcR = pre_op->get_output_tensor_shape(pre_config, t.owner_idx, srcId);
          if (dstR.intersection(srcR).get_volume() > 0) {
            // Forward dependency
            {
              SimTask* dstT = task_manager->get_forward_task(op, dstId);
              SimTask* srcT = task_manager->get_forward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(srcT, dstT, dstR.intersection(srcR).get_volume() * element_size);
            }
            // Backward dependency
            if (comp_mode == COMP_MODE_TRAINING) {
              SimTask* dstT = task_manager->get_backward_task(op, dstId);
              SimTask* srcT = task_manager->get_backward_task(pre_op, srcId);
              add_task_dependencies_with_xfer(dstT, srcT, dstR.intersection(srcR).get_volume() * element_size);
            }
          }
        }
      }
    }
  }

  topofinder->simplify_mp_topology();
  topofinder->generate_dp_topology();
  topofinder->optimize_indirection();

  // Step 4: add ready tasks into ready_queue
  std::priority_queue<SimTask*, std::vector<SimTask*>, SimTaskCompare> ready_queue;
  for (size_t i = 0; i < task_manager->global_task_id; i++)
    if (task_manager->tasks[i]->counter == 0)
      ready_queue.push(task_manager->tasks[i]);

  // Step 5: perform simulation

  float sim_time = 0.0f;
  std::map<Device*, float> device_times;
  // map<Device*, SimTask*> device_schedule;
  size_t idx = 0;
  while (!ready_queue.empty()) {
    // Find the task with the earliest start time
    SimTask* cur_task = ready_queue.top();
    ready_queue.pop();
    float ready_time = 0;
    float end_time;
    if (device_times.find(cur_task->device) != device_times.end()) {
      ready_time = device_times[cur_task->device];
    }
    float start_time = std::max(ready_time, cur_task->ready_time);
    if (cur_task->type == SimTask::TASK_NOMINAL_COMM) {
      end_time = route_transfer(cur_task, start_time, device_times);
    }
    else if (cur_task->type == SimTask::TASK_ALLREDUCE) {
      expand_allreduce(cur_task, start_time, ready_queue);
      idx++;
      continue;
    }
    else {
      end_time = start_time + cur_task->run_time;
      device_times[cur_task->device] = end_time;
    }

#ifdef DEBUG_PRINT
    printf("task[%lu/%lu] type(%d) run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          idx, task_manager->global_task_id, cur_task->type, cur_task->run_time, ready_time, start_time, (cur_task->device->name).c_str());
#endif

    if (end_time > sim_time) {
      sim_time = end_time;
    }

    for (size_t i = 0; i < cur_task->next_tasks.size(); i++) {
      SimTask* next = cur_task->next_tasks[i];
      // next->ready_time = max(next->ready_time, end_time);
      if (end_time > next->ready_time) {
        next->ready_time = end_time;
        // next->prev = t;
      }
      next->counter--;
      if (next->counter == 0) {
        ready_queue.push(next);
      }
    }
    idx++;
  }
  assert(idx == task_manager->global_task_id); 

  return sim_time;
}

// float DLSSchedulerBasedSimulator::compute_null_xfertime(SimTask* src, SimTask* dst, size_t message_size) 
// {
// }

void DLSSchedulerBasedSimulator::add_task_dependencies_with_xfer_sch(
    SimTask* src_task, SimTask* dst_task, size_t message_size,
    DLSTaskDag& bp_taskdag) 
{
  if (bp_taskdag.find(src_task) == bp_taskdag.end()) {
    bp_taskdag.emplace(std::make_pair(src_task, 
      std::vector<std::pair<SimTask*, size_t>>({std::make_pair(dst_task, message_size)})));
  }
  else {
    bp_taskdag[src_task].emplace_back(std::make_pair(dst_task, message_size));
  }
}

std::vector<SimTask*> DLSSchedulerBasedSimulator::rev_topological_sort(
    const std::unordered_set<SimTask*>& entry_nodes,
    const std::unordered_set<SimTask*>& bp_tasks,
    const DLSTaskDag& bp_taskdag) 
{
  std::vector<SimTask*> result;
  
  std::unordered_map<SimTask*, bool> visited;
  std::stack<std::pair<SimTask*, bool>> dfs_stack;

  for (SimTask* p: bp_tasks) {
    visited[p] = false;
  }
  for (SimTask *t: entry_nodes) {
    dfs_stack.push({t, false});
  }

  while (!dfs_stack.empty()) {
    auto node = dfs_stack.top();
    dfs_stack.pop();
    SimTask * task = node.first;

    if (node.second) {
      result.push_back(task);
      continue;
    }
    if(visited[task])
      continue;

    visited[task] = true;
    dfs_stack.push({task, true});
    
    if (bp_taskdag.find(task) == bp_taskdag.end())
      continue;

    for (auto & succ: bp_taskdag.at(task)) {
      if (!visited[succ.first]) {
        dfs_stack.push({succ.first, false});
      }
    }
  }
  assert(result.size() == bp_tasks.size());
  return result;
}

std::unordered_map<SimTask*, float> DLSSchedulerBasedSimulator::get_slevels(
    const std::vector<SimTask*>&& sorted_bp_tasks, 
    const DLSTaskDag& bp_taskdag) 
{
  std::unordered_map<SimTask*, float> slevels;
  for (SimTask * task: sorted_bp_tasks) {
    float maximum = 0;
    if (bp_taskdag.find(task) != bp_taskdag.end()) {
      for (auto & succ: bp_taskdag.at(task)) {
        if (slevels.find(succ.first) != slevels.end()) {
          if (slevels[succ.first] > maximum) {
            maximum = slevels[succ.first];
          }
        }
      }
    }
    float slevel = maximum + task->run_time;
    slevels[task] = slevel;
  }
  return slevels;
}

std::unordered_set<SimTask*> DLSSchedulerBasedSimulator::get_entry_nodes(
      const std::unordered_set<SimTask*>& bp_tasks,
      const DLSTaskDag& bp_taskdag)
{
  std::unordered_set<SimTask*> entryset;
  for (SimTask* p: bp_tasks) {
    entryset.insert(p);
  }
  for (auto & tp: bp_taskdag) {
    for (auto & p: tp.second) {
      if (entryset.find(p.first) != entryset.end()) {
        entryset.erase(entryset.find(p.first));
      }
    }
  }
  return entryset;
}

double DLSSchedulerBasedSimulator::dls_schedule(std::unordered_set<SimTask*>& bp_tasks,
    DLSTaskDag& bp_taskdag)  
{
  std::unordered_set<SimTask*> entry_nodes = get_entry_nodes(bp_tasks, bp_taskdag);
  DLSTaskDag predecessor_map = get_predecessors(bp_tasks, bp_taskdag);
  std::unordered_map<SimTask*, float> slevels = 
    get_slevels(rev_topological_sort(entry_nodes, bp_tasks, bp_taskdag), bp_taskdag);
  
  std::unordered_set<SimTask*> ready_pool;
  std::unordered_set<SimTask*> unsched_node;
  std::unordered_map<SimTask*, float> task_finish_time;
  std::unordered_map<SimTask*, size_t> num_predecessor_cleared;
  std::unordered_map<CompDevice*, float> processors;
  std::unordered_map<CommDevice*, float> links;
  
  // all nodes are not scheduled
  for (SimTask* tp: bp_tasks) {
    unsched_node.insert(tp);
  }
  // add entry nodes
  for (auto & task: entry_nodes) {
    ready_pool.insert(task);
  }
  // prepare processors
  for (auto & id_to_gpu: net_machine->id_to_gpu) {
    processors[id_to_gpu.second] = 0.0f;
  }
  // prepare links
  for (auto & l: net_machine->ids_to_nw_comm_device) {
    links[l.second] = 0.0f;
  }

  // TODOOO: initial assignemet
  // TODOOO: insersion based scheduling

  while (!unsched_node.empty()) {
    // std::set<std::pair<float, std::pair<CompDevice*, SimTask*>>> dlevels;
    float max_dl = std::numeric_limits<float>::lowest();
    CompDevice* final_gpu = nullptr;
    SimTask* final_task = nullptr;
    // for each node, for each processor, compute the DL
    for (auto & processor_time: processors) {
      for (SimTask* ready_task: ready_pool) {
        float da = compute_da(processor_time.first, ready_task, 
          predecessor_map, processors, links, task_finish_time);
        float dl = slevels[ready_task] - std::max(da, processor_time.second);
#ifdef TEST_DLSSCHEDULER
        fprintf(stderr, "Checking %p on GPU %d, da = %f, dl = %f (max: %f)\n", ready_task, processor_time.first->device_id, da, dl, max_dl);
#endif 
        if (dl > max_dl) {
          max_dl = dl;
          final_gpu = processor_time.first; 
          final_task = ready_task;
        }
        else if (dl == max_dl) {
          if ((float)std::rand() / RAND_MAX > 0.5) {
            max_dl = dl;
            final_gpu = processor_time.first; 
            final_task = ready_task;
          }
        }
      }
    }
    assert(final_gpu != nullptr && final_task != nullptr);
    schedule(final_gpu, final_task, predecessor_map, processors, links, task_finish_time);
    for (auto& s: bp_taskdag[final_task]) {
      if (num_predecessor_cleared.find(s.first) == num_predecessor_cleared.end()) {
        num_predecessor_cleared[s.first] = 1;
      } else {
        num_predecessor_cleared[s.first]++;
      }
      if (num_predecessor_cleared[s.first] == predecessor_map[s.first].size()) {
        // cleared, add this task to ready pool
        ready_pool.insert(s.first);
      }
    }
    ready_pool.erase(ready_pool.find(final_task));
    unsched_node.erase(unsched_node.find(final_task));
  }
  double final_makespan = std::numeric_limits<double>::min();
  for (auto &tf: task_finish_time) {
    if (tf.second > final_makespan)
      final_makespan = tf.second;
  }
  return final_makespan;
}

DLSTaskDag DLSSchedulerBasedSimulator::get_predecessors(
    const std::unordered_set<SimTask*>& bp_tasks, const DLSTaskDag& bp_taskdag) 
{
  DLSTaskDag rev_map;
  for (auto& s: bp_taskdag) {
    for (auto& succ: s.second)
    if (rev_map.find(succ.first) == rev_map.end()) {
      rev_map.emplace(std::make_pair(succ.first, 
        std::vector<std::pair<SimTask*, size_t>>({std::make_pair(s.first, succ.second)})));
    } else {
      rev_map[succ.first].emplace_back(std::pair<SimTask*, size_t>(std::make_pair(s.first, succ.second)));
    }
  }
  return rev_map;
}

inline float DLSSchedulerBasedSimulator::compute_da(
    const CompDevice* proc, SimTask* task, 
    const DLSTaskDag& pred, 
    const std::unordered_map<CompDevice*, float>& processors, 
    const std::unordered_map<CommDevice*, float>& links,
    const std::unordered_map<SimTask*, float>& task_finish_time)
{
  if (pred.find(task) == pred.end())
    return 0;
  // task->device = proc;
  task->mem = net_machine->get_gpu_fb_mem(proc->device_id);
  float final_da = std::numeric_limits<float>::lowest();
  for (auto & ps: pred.at(task)) {
    float da = sch_try_route_transfer(ps.first, task, ps.second, links, task_finish_time);
    if (da > final_da) {
      final_da = da;
    }
  }
  return final_da;
}

inline void DLSSchedulerBasedSimulator::schedule(CompDevice* proc, SimTask* task,
    const DLSTaskDag& pred,  
    std::unordered_map<CompDevice*, float>& processors, 
    std::unordered_map<CommDevice*, float>& links,
    std::unordered_map<SimTask*, float>& task_finish_time)
{
#ifdef TEST_DLSSCHEDULER
  fprintf(stderr, "Decided: %p on GPU %d\n", task, proc->device_id);
#endif
  task->device = proc;
  task->mem = net_machine->get_gpu_fb_mem(proc->device_id);
  float final_ready_time = 0;
  if (pred.find(task) != pred.end()) {
    for (auto & ps: pred.at(task)) {
      float da = sch_route_transfer(ps.first, task, ps.second, links, task_finish_time);
      if (da > final_ready_time) {
        final_ready_time = da;
      }
    }
  }
  processors[proc] = final_ready_time + task->run_time;
  task_finish_time[task] = final_ready_time + task->run_time;
}

    
std::map<Op*, ParallelConfig> DLSSchedulerBasedSimulator::get_device_placements(
    std::unordered_set<SimTask*>& bp_tasks,
    const std::map<Op*, ParallelConfig>& global)
{
  std::map<Op*, ParallelConfig> result = global;

  for (auto& op_pc: result) {
    Op* op = op_pc.first;
    ParallelConfig &config = result.find(op)->second;
    for (int j = 0; j < config.num_parts(); j++) {
      SimTask* bw = task_manager->get_backward_task(op, j);
      config.device_ids[j] = bw->device->device_id;
    }
  } 
  return result;
}

float DLSSchedulerBasedSimulator::sch_try_route_transfer(
    SimTask* src, SimTask* dst, size_t xfersize, 
    const std::unordered_map<CommDevice*, float>& links,
    const std::unordered_map<SimTask*, float>& task_finish_time) 
{
  std::vector<CommDevice *> path = net_machine->get_comm_path(src->mem, dst->mem);
  assert(path.size() <= 1);

  if (path.size() == 0)
    return task_finish_time.at(src);

  std::vector<CommDevice *> route = 
    static_cast<NominalCommDevice*>(path[0])->expand_to_physical();

  float curr_task_start_time; 
  float curr_task_finish_time; 
  float curr_task_run_time = 0; 
  float ready_time = task_finish_time.at(src);
  float curr_task_ready_time = ready_time; 
  float xfer_size = xfersize;
  
  std::unordered_map<CommDevice*, float> local_lat_device;

  float final_start_time = 0;
  float final_finish_time = 0;
#ifdef TEST_DLSSCHEDULER
  fprintf(stderr, "DLS_TRY_ROUTE: src: %p, dst: %p, xfersize: %lf\n", src, dst, xfer_size);
#endif
  for (unsigned int i = 0; i < route.size(); i++) {
    CommDevice * latency_task_device = route[i];
    float latency_task_run_time = net_machine->get_inter_node_gpu_latency();
    float latency_task_ready_time; 
    float latency_task_start_time; 
    if (i == 0) {
      latency_task_ready_time = curr_task_ready_time + curr_task_run_time;
      latency_task_start_time = std::max(links.at(latency_task_device), latency_task_ready_time);
      local_lat_device[latency_task_device] = latency_task_start_time;
      final_start_time = latency_task_start_time;
    }
    else {
      latency_task_ready_time = curr_task_finish_time;
      latency_task_start_time = std::max(links.at(latency_task_device), latency_task_ready_time);
      local_lat_device[latency_task_device] = latency_task_start_time;
    }
    float latency_task_finish_time = latency_task_start_time + latency_task_run_time;
    local_lat_device[latency_task_device] = latency_task_finish_time;
    float dram_to_dram_run_time = xfer_size / latency_task_device->bandwidth;

    float dram_to_dram_start_time = latency_task_finish_time;
    float dram_to_dram_finish_time = dram_to_dram_start_time + dram_to_dram_run_time;
    local_lat_device[latency_task_device] = dram_to_dram_finish_time;

    if (dram_to_dram_finish_time > final_finish_time) {
      final_finish_time = dram_to_dram_finish_time;
    }

    curr_task_ready_time = latency_task_ready_time;
    curr_task_start_time = latency_task_start_time;
    curr_task_finish_time = latency_task_finish_time;
    curr_task_run_time = latency_task_run_time;
    
#ifdef TEST_DLSSCHEDULER
    fprintf(stderr, "\tDLS_TRY_ROUTE: expand: route[%u] run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          i, curr_task_run_time, curr_task_ready_time, curr_task_start_time, (latency_task_device->name).c_str());
    fprintf(stderr, "\t\ttDLS_TRY_ROUTE: d2d: run_time(%.4lf) start_time(%.4lf) device(%s)\n",
          dram_to_dram_run_time, dram_to_dram_start_time, (latency_task_device->name).c_str());
#endif
  }
#ifdef TEST_DLSSCHEDULER
  fprintf(stderr, "DLS_TRY_ROUTE: fft: %f\n", final_finish_time);
#endif
  return final_finish_time;
}

float DLSSchedulerBasedSimulator::sch_route_transfer(
    SimTask* src, SimTask* dst, size_t xfersize, 
    std::unordered_map<CommDevice*, float>& links,
    std::unordered_map<SimTask*, float>& task_finish_time) 
{
  std::vector<CommDevice *> path = net_machine->get_comm_path(src->mem, dst->mem);
  assert(path.size() <= 1);
  if (path.size() == 0)
    return task_finish_time.at(src);

  std::vector<CommDevice *> route = 
    static_cast<NominalCommDevice*>(path[0])->expand_to_physical();

  float curr_task_start_time; 
  float curr_task_finish_time; 
  float curr_task_run_time = 0; 
  float ready_time = task_finish_time.at(src);
  float curr_task_ready_time = ready_time; 
  float xfer_size = xfersize;
  
  float final_start_time = 0;
  float final_finish_time = 0;
#ifdef TEST_DLSSCHEDULER
  fprintf(stderr, "DLS_TRY_ROUTE: src: %p, dst: %p, xfersize: %lf\n", src, dst, xfer_size);
#endif
  for (unsigned int i = 0; i < route.size(); i++) {
    CommDevice * latency_task_device = route[i];
    float latency_task_run_time = net_machine->get_inter_node_gpu_latency();
    float latency_task_ready_time; 
    float latency_task_start_time; 
    if (i == 0) {
      latency_task_ready_time = curr_task_ready_time + curr_task_run_time;
      latency_task_start_time = std::max(links.at(latency_task_device), latency_task_ready_time);
      links[latency_task_device] = latency_task_start_time;
      final_start_time = latency_task_start_time;
    }
    else {
      latency_task_ready_time = curr_task_finish_time;
      latency_task_start_time = std::max(links.at(latency_task_device), latency_task_ready_time);
      links[latency_task_device] = latency_task_start_time;
    }
    float latency_task_finish_time = latency_task_start_time + latency_task_run_time;
    links[latency_task_device] = latency_task_finish_time;
    float dram_to_dram_run_time = xfer_size / latency_task_device->bandwidth;

    float dram_to_dram_start_time = latency_task_finish_time;
    float dram_to_dram_finish_time = dram_to_dram_start_time + dram_to_dram_run_time;
    links[latency_task_device] = dram_to_dram_finish_time;

    if (dram_to_dram_finish_time > final_finish_time) {
      final_finish_time = dram_to_dram_finish_time;
    }

    curr_task_ready_time = latency_task_ready_time;
    curr_task_start_time = latency_task_start_time;
    curr_task_finish_time = latency_task_finish_time;
    curr_task_run_time = latency_task_run_time;
    
#ifdef TEST_DLSSCHEDULER
    fprintf(stderr, "\tDLS_TRY_ROUTE: expand: route[%u] run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          i, curr_task_run_time, curr_task_ready_time, curr_task_start_time, (latency_task_device->name).c_str());
    fprintf(stderr, "\t\ttDLS_TRY_ROUTE: d2d: run_time(%.4lf) start_time(%.4lf) device(%s)\n",
          dram_to_dram_run_time, dram_to_dram_start_time, (latency_task_device->name).c_str());
#endif
  }
#ifdef TEST_DLSSCHEDULER
  fprintf(stderr, "DLS_TRY_ROUTE: fft: %f\n", final_finish_time);
#endif
  return final_finish_time;
}

#ifdef TEST_DLSSCHEDULER
void DLSSchedulerBasedSimulator::test() 
{
  DLSTaskDag tg;
  std::unordered_set<SimTask*> tasks;
  std::unordered_map<SimTask*, std::string> printables;
  
  SimTask * t1 = new SimTask();
  t1->run_time = 2;
  tasks.insert(t1);
  printables[t1] = "t1";
  SimTask * t2 = new SimTask();
  t2->run_time = 3;
  tasks.insert(t2);
  printables[t2] = "t2";
  SimTask * t3 = new SimTask();
  t3->run_time = 3;
  tasks.insert(t3);
  printables[t3] = "t3";
  SimTask * t4 = new SimTask();
  t4->run_time = 4;
  tasks.insert(t4);
  printables[t4] = "t4";
  SimTask * t5 = new SimTask();
  t5->run_time = 5;
  tasks.insert(t5);
  printables[t5] = "t5";
  SimTask * t6 = new SimTask();
  t6->run_time = 4;
  tasks.insert(t6);
  printables[t6] = "t6";
  SimTask * t7 = new SimTask();
  t7->run_time = 4;
  tasks.insert(t7);
  printables[t7] = "t7";
  SimTask * t8 = new SimTask();
  t8->run_time = 4;
  tasks.insert(t8);
  printables[t8] = "t8";
  SimTask * t9 = new SimTask();
  t9->run_time = 1;
  tasks.insert(t9);
  printables[t9] = "t9";
  add_task_dependencies_with_xfer_sch(t1, t2, 4, tg);
  add_task_dependencies_with_xfer_sch(t1, t3, 1, tg);
  add_task_dependencies_with_xfer_sch(t1, t4, 1, tg);
  add_task_dependencies_with_xfer_sch(t1, t5, 1, tg);
  add_task_dependencies_with_xfer_sch(t1, t7, 10, tg);
  add_task_dependencies_with_xfer_sch(t2, t6, 1, tg);
  add_task_dependencies_with_xfer_sch(t2, t7, 1, tg);
  add_task_dependencies_with_xfer_sch(t3, t8, 1, tg);
  add_task_dependencies_with_xfer_sch(t4, t8, 1, tg);
  add_task_dependencies_with_xfer_sch(t6, t9, 5, tg);
  add_task_dependencies_with_xfer_sch(t7, t9, 6, tg);
  add_task_dependencies_with_xfer_sch(t8, t9, 5, tg);

  for (auto tn: printables) {
    fprintf(stderr, "Task %p is %s.\n", tn.first, tn.second.c_str());
  }

  auto entry_set = get_entry_nodes(tasks, tg);
  std::cerr << "Entry set: " << std::endl;
  for (auto t: entry_set) {
    std::cerr << printables[t] << std::endl;
  }
   
  auto toposorted = rev_topological_sort(entry_set, tasks, tg);
  std::cerr << "toposorted order: " << std::endl;
  for (auto t: toposorted) {
    std::cerr << printables[t] << " ";
  }
  std::cerr << std::endl;

  auto slevels = get_slevels(std::move(toposorted), tg);
  std::cerr << "slevels: " << std::endl;
  for (auto t: tasks) {
    std::cerr << printables[t] << ": " << slevels[t] << std::endl;
  }
  std::cerr << std::endl;

  auto preds = get_predecessors(tasks, tg);
  std::cerr << "preds: " << std::endl;
  for (auto t: tasks) {
    std::cerr << printables[t] << ": ";
    for (auto p: preds[t]) {
      std::cerr << printables[p.first] << "(" << p.second << ") ";
    }
  }
  std::cerr << std::endl;
  
  NetworkedMachineModel * nm = new NetworkedMachineModel();   
  nm->num_nodes = 4;
  nm->num_gpus_per_node = 1;
  nm->num_gpus = 4;
  nm->total_devs = 4;
  nm->num_switches = 0;
  nm->inter_gpu_bandwidth = 1;
  nm->gpu_dram_bandwidth = 1;
  nm->link_bandwidth = 1;
  nm->network_latency = 0;
  nm->set_pcie(false);
  nm->set_pipeline(true);
  for (int i = 0; i < 4; i++) {
    std::string gpu_name = "GPU " + std::to_string(i);
    nm->id_to_gpu[i] = new CompDevice(gpu_name, CompDevice::TOC_PROC, i, i, i);
    std::string gpu_mem_name = "GPU_FB_MEM " + std::to_string(i);
    nm->id_to_gpu_fb_mem[i] = new MemDevice(gpu_mem_name, MemDevice::GPU_FB_MEM, i, i, i, 100000000000ULL);
  }

  ConnectionMatrix conn = ConnectionMatrix(16, 0);
  auto addlink = [&](int i, int j) {
    conn[i * 4 + j] = 1;
    conn[j * 4 + i] = 1;
  };
  addlink(0, 1);
  addlink(1, 2);
  addlink(2, 3);
  addlink(3, 0);
  nm->conn_matrix = conn;
 
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      // if (conn_matrix[i * 4 + j] > 0) {
        int device_id = i * 4 + j;
        std::string link_name = "LINK " + std::to_string(i) + "-" + std::to_string(j);
        nm->ids_to_nw_comm_device[device_id] = new CommDevice(link_name, CommDevice::NW_COMM, 
          -1, -1, device_id, 0, conn[i * 4 + j]);
      }
    // }
  }

  auto get_routes = [&](int i, int j) -> EcmpRoutes {
    if ((i+1)%4 == j)
      return std::make_pair(std::vector<float>({1}), std::vector<Route>({Route({nm->ids_to_nw_comm_device[i*4+j]})}));
    if ((j+1)%4 == i)
      return std::make_pair(std::vector<float>({1}), std::vector<Route>({Route({nm->ids_to_nw_comm_device[j*4+i]})}));
    return std::make_pair(
      std::vector<float>({1}), 
      std::vector<Route>({Route{nm->ids_to_nw_comm_device[i*4+(i+1)%4], nm->ids_to_nw_comm_device[((i+1)%4)*4+j]}})
    );
  };

  for (int i = 0; i < nm->num_nodes; i++) {
    for (int j = 0; j < nm->num_nodes; j++) {
      int device_id = i * nm->total_devs + j;
      std::string link_name = "NOMINAL " + std::to_string(i) + "-" + std::to_string(j);
      nm->ids_to_nw_nominal_device[device_id] = new NominalCommDevice(link_name, device_id);
      nm->ids_to_nw_nominal_device[device_id]->set_physical_paths(get_routes(i, j));
    }
  }

  this->net_machine = nm;
  this->machine = nm;
  double sch_len = dls_schedule(tasks, tg);
  std::cerr << "final makespan: " << sch_len << std::endl;
  
}
#endif