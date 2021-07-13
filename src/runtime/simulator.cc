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
#include "simulator.h"
#include "model.h"
#include "queue"

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

void* Simulator::allocate(size_t num_elements, DataType type)
{
  off_t element_size = data_type_size(type);
  void* ret_ptr = base_ptr + offset;
  offset += element_size * num_elements;
  if (offset > capacity) {
    fprintf(stderr, "Simulator cannot measure some operators' performance."
        " Increate --simulator-workspace-size to at least %zd\n", offset);
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
    // printf("task[%lu] type(%d) run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
    //       idx, cur_task->type, cur_task->run_time, ready_time, start_time, (cur_task->device->name).c_str());
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
        // NER step: add allreduce task after backward propogation
        size_t element_size = data_type_size(DT_FLOAT); // assume all weights have float elements
        ParallelConfig pc = global.find(op)->second;
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
              task2->add_next_task(ar_task);
              if (l1optimizer) 
                l1optimizer->task_added(ar_task);
            }
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

    printf("task[%lu/%lu] type(%d) run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          idx, task_manager->global_task_id, cur_task->type, cur_task->run_time, ready_time, start_time, (cur_task->device->name).c_str());

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
    
    printf("\texpand: route[%u] run_time(%.4lf) ready_time(%.4lf) start_time(%.4lf) device(%s)\n",
          i, curr_task_run_time, curr_task_ready_time, curr_task_start_time, (latency_task_device->name).c_str());
    printf("\t\td2d: run_time(%.4lf) start_time(%.4lf) device(%s)\n",
          dram_to_dram_run_time, dram_to_dram_start_time, (latency_task_device->name).c_str());

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

// TODO
void LogicalTaskgraphBasedSimulator::searlize_logical_taskgraph(std::string const &export_file_name) {
  return;
}