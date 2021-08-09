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

#include "simulator.h"
#include "model.h"
//#include "realm/runtime_impl.h"
//#include "realm/cuda/cuda_module.h"
#include "cuda_helper.h"

typedef long long int coord_t;

typedef Realm::Point<1, coord_t> Point1;
typedef Realm::Rect<1, coord_t> Rect1;

Simulator::Simulator(const FFModel* model,
                     FFHandler _handler,
                     Memory _memory,
                     MachineModel *machine)
: memory(_memory), handler(_handler),
  offset(0), warmup_times(5), repeat_times(10),
  computationMode(model->config.computationMode)
{
  // Allocate simulator memory
  Rect1 bounds(Point1(0), Point1(0));
  std::vector<size_t> field_sizes;
  field_sizes.push_back(model->config.simulator_work_space_size);
  Realm::RegionInstance::create_instance(simulatorInst,
      memory, bounds, field_sizes, 0, Realm::ProfilingRequestSet()).wait();
  base_ptr = (char*)simulatorInst.pointer_untyped(0, sizeof(char));
  capacity = model->config.simulator_work_space_size;

  if (model->simonly) {
    // allocate memory for workspace
    Memory gpu_mem = _memory;
    Realm::Rect<1, coord_t> bounds(Realm::Point<1, coord_t>(0),
        Realm::Point<1, coord_t>(handler.workSpaceSize-1));
    std::vector<size_t> field_sizes;
    field_sizes.push_back(sizeof(char));
    Realm::RegionInstance workspaceInst;
    Realm::RegionInstance::create_instance(workspaceInst, gpu_mem, bounds,
        field_sizes, 0, Realm::ProfilingRequestSet()).wait();
    handler.workSpace = workspaceInst.pointer_untyped(0, sizeof(char));
  }

  size_t max_num_tasks = 256 * 1024 * 1024;

  cudaEventCreate(&start_event);
  cudaEventCreate(&end_event);
  conv2d_meta = new Conv2DMeta(handler);
  linear_meta = new LinearMeta(handler, 4096);
  pool2d_meta = new Pool2DMeta(handler);
  ele_unary_meta = new ElementUnaryMeta(handler);
  ele_binary_meta = new ElementBinaryMeta(handler);
  //softmax_meta = new SoftmaxMeta(handler);
  batch_matmul_meta = new BatchMatmulMeta(handler);
  concat_meta = new ConcatMeta(handler);
  //dropout_meta = new DropoutMeta(handler);
  transpose_meta = new TransposeMeta(handler);
  this->machine = machine;
  segment_size = model->config.simulator_segment_size;
  max_num_segments = model->config.simulator_max_num_segments;
  // Initialize task manager
  task_manager = new TaskManager(max_num_tasks);

  measurements = nullptr;
  l1optimizer = nullptr;
}

Simulator::~Simulator(void)
{
  simulatorInst.destroy();
}

__host__
void Simulator::strategy_search_task(const Task *task,
                                     const std::vector<PhysicalRegion> &regions,
                                     Context ctx, Runtime *runtime)
{
  const FFModel* model = *((FFModel**) task->args);
  Memory gpu_mem = Machine::MemoryQuery(Machine::get_machine())
         .only_kind(Memory::GPU_FB_MEM).best_affinity_to(task->target_proc).first();
  // Realm::MemoryImpl* memImpl =
  //     Realm::get_runtime()->get_memory_impl(gpu_mem);
  // Realm::Cuda::GPUFBMemory* memFBImpl = (Realm::Cuda::GPUFBMemory*) memImpl;
  // off_t offset = memFBImpl->alloc_bytes_local(model->config.simulator_work_space_size);
  // void* base_ptr = memFBImpl->get_direct_ptr(offset, 0);
  MachineModel *machine;
  if (model->config.machine_model_version == 0) {
    machine = (MachineModel *) new SimpleMachineModel(model->config.numNodes, model->config.workersPerNode, gpu_mem.capacity());
  }
  else if (model->config.machine_model_version == 1 and !model->config.machine_model_file.empty()) {
    machine = (MachineModel *) new EnhancedMachineModel(model->config.machine_model_file, gpu_mem.capacity());
  }
  else {
    assert(false && "machine model creation error: currently only support machine-model-version = 0 or 1. When machine-model-version = 1, machine-model-file should not be empty.");
  }
  // Assume this task is running on GPU0
  Simulator* simulator = new Simulator(model, model->handlers[0], gpu_mem, machine);
  // Set cublas/cudnn streams to allow Realm catch the events

  cudaStream_t stream;
  checkCUDA(get_legion_stream(&stream));
  checkCUDA(cublasSetStream(simulator->handler.blas, stream));
  checkCUDNN(cudnnSetStream(simulator->handler.dnn, stream));

  std::map<Op*, ParallelConfig> strategies;
  if (model->config.import_strategy_file.length() > 0) {
    // Load the strategy from config.strategies
    for (size_t l = 0; l < model->layers.size(); l++) {
      MappingTagID key = FFConfig::get_hash_id(std::string(model->layers[l]->name));
      std::map<MappingTagID, ParallelConfig>::const_iterator iter;
      iter = model->config.strategies.find(key);
      if (iter == model->config.strategies.end()) {
        fprintf(stderr, "ERROR: Cannot find strategy for operator %s in "
                "strategy file %s\n", model->layers[l]->name,
                model->config.import_strategy_file.c_str());
      }
      strategies[model->layers[l]] = iter->second;
    }
  } else {
    // Start from data parallel
    for (size_t l = 0; l < model->layers.size(); l++) {
      strategies[model->layers[l]] = model->layers[l]->get_data_parallel_config(*model);
    }
  }
  if (model->config.computationMode == COMP_MODE_TRAINING) {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(TRAINING)\n",
        model->config.search_budget, model->config.search_alpha);
  } else {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(INFERENCE)\n",
        model->config.search_budget, model->config.search_alpha);
  }
  model->optimize(simulator, strategies, model->config.search_budget,
      model->config.search_alpha, model->config.computationMode, model->config.enable_propagation);
  if (model->config.export_strategy_file.length() > 0) {
    fprintf(stderr, "Exporting the best discovered strategy to %s.\n",
        model->config.export_strategy_file.c_str());
    std::map<Op*, ParallelConfig>::const_iterator iter;
    std::map<std::string, ParallelConfig> strategy_output;
    for (iter = strategies.begin(); iter != strategies.end(); iter++) {
      strategy_output[iter->first->name] = iter->second;
    }
    save_strategies_to_file(model->config.export_strategy_file, strategy_output);
    fprintf(stderr, "To use the strategy for distributed training, restart"
        " FlexFlow and import the strategy (i.e., --import %s)\n",
        model->config.export_strategy_file.c_str());
    exit(0);
  }  else {
    fprintf(stderr, "The best discovered strategy is not exported.\n"
        "Please set a path to export the strategy using --export or --export-strategy.\n");
    exit(0);
  }
  // Start from data
  // memFBImpl->free_bytes_local(offset, model->config.simulator_work_space_size);
  delete(simulator);
  delete(machine);
}

__host__
void Simulator::simulation_task(const Task *task,
                                     const std::vector<PhysicalRegion> &regions,
                                     Context ctx, Runtime *runtime)
{
  FFModel* model = *((FFModel**) task->args);
  Memory gpu_mem = Machine::MemoryQuery(Machine::get_machine())
         .only_kind(Memory::GPU_FB_MEM).best_affinity_to(task->target_proc).first();
  // Realm::MemoryImpl* memImpl =
  //     Realm::get_runtime()->get_memory_impl(gpu_mem);
  // Realm::Cuda::GPUFBMemory* memFBImpl = (Realm::Cuda::GPUFBMemory*) memImpl;
  // off_t offset = memFBImpl->alloc_bytes_local(model->config.simulator_work_space_size);
  // void* base_ptr = memFBImpl->get_direct_ptr(offset, 0);
  BigSwitchNetworkTopologyGenerator topo_gen = BigSwitchNetworkTopologyGenerator(model->config.numNodes);
  
  // NetworkedMachineModel *nmachine = new NetworkedMachineModel(model->config.numNodes, 
  //   model->config.workersPerNode,  
  //   1,
  //   topo_gen.generate_topology(),
  //   gpu_mem.capacity(),
  //   20.0 * 1024 * 1024 / 8
  // );
  // nmachine->set_pcie(false);
  // nmachine->set_pipeline(false);

  SimpleMachineModel* nmachine = new SimpleMachineModel(model->config.numNodes, model->config.workersPerNode, gpu_mem.capacity());

  MachineModel *machine;
  machine = reinterpret_cast<MachineModel*>(nmachine);

  // Assume this task is running on GPU0
  Simulator* simulator = new Simulator(model, model->handlers[0], gpu_mem, machine);
  // Set cublas/cudnn streams to allow Realm catch the events

  if (model->config.mfile != "") {
    model->load_measurement(simulator, model->config.mfile);
  }

  cudaStream_t stream;
  checkCUDA(get_legion_stream(&stream));
  checkCUDA(cublasSetStream(simulator->handler.blas, stream));
  checkCUDNN(cudnnSetStream(simulator->handler.dnn, stream));

  std::map<Op*, ParallelConfig> strategies;
  if (model->config.import_strategy_file.length() > 0) {
    // Load the strategy from config.strategies
    for (size_t l = 0; l < model->layers.size(); l++) {
      MappingTagID key = FFConfig::get_hash_id(std::string(model->layers[l]->name));
      std::map<MappingTagID, ParallelConfig>::const_iterator iter;
      iter = model->config.strategies.find(key);
      if (iter == model->config.strategies.end()) {
        fprintf(stderr, "ERROR: Cannot find strategy for operator %s in "
                "strategy file %s\n", model->layers[l]->name,
                model->config.import_strategy_file.c_str());
      }
      strategies[model->layers[l]] = iter->second;
    }
  } else {
    // Start from data parallel
    for (size_t l = 0; l < model->layers.size(); l++) {
      strategies[model->layers[l]] = model->layers[l]->get_data_parallel_config(*model);
    }
  }
  if (model->config.computationMode == COMP_MODE_TRAINING) {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(TRAINING)\n",
        model->config.search_budget, model->config.search_alpha);
  } else {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(INFERENCE)\n",
        model->config.search_budget, model->config.search_alpha);
  }
  model->optimize(simulator, strategies, model->config.search_budget,
      model->config.search_alpha, model->config.computationMode, model->config.enable_propagation);
  if (model->config.export_strategy_file.length() > 0) {
    fprintf(stderr, "Exporting the best discovered strategy to %s.\n",
        model->config.export_strategy_file.c_str());
    std::map<Op*, ParallelConfig>::const_iterator iter;
    std::map<std::string, ParallelConfig> strategy_output;
    for (iter = strategies.begin(); iter != strategies.end(); iter++) {
      strategy_output[iter->first->name] = iter->second;
    }
    save_strategies_to_file(model->config.export_strategy_file, strategy_output);
    fprintf(stderr, "To use the strategy for distributed training, restart"
        " FlexFlow and import the strategy (i.e., --import %s)\n",
        model->config.export_strategy_file.c_str());
    exit(0);
  }  else {
    fprintf(stderr, "The best discovered strategy is not exported.\n"
        "Please set a path to export the strategy using --export or --export-strategy.\n");
    exit(0);
  }
  // Start from data
  // memFBImpl->free_bytes_local(offset, model->config.simulator_work_space_size);
  delete(simulator);
  delete(machine);
}


__host__
void Simulator::measurement_task(const Task *task,
                                     const std::vector<PhysicalRegion> &regions,
                                     Context ctx, Runtime *runtime)
{
  const FFModel* model = *((FFModel**) task->args);
  Memory gpu_mem = Machine::MemoryQuery(Machine::get_machine())
         .only_kind(Memory::GPU_FB_MEM).best_affinity_to(task->target_proc).first();

  SimpleMachineModel* nmachine = new SimpleMachineModel(model->config.numNodes, model->config.workersPerNode, gpu_mem.capacity());
  MachineModel *machine;
  machine = reinterpret_cast<MachineModel*>(nmachine);

  // Assume this task is running on GPU0
  Simulator* simulator = new Simulator(model, model->handlers[0], gpu_mem, machine);
  // Set cublas/cudnn streams to allow Realm catch the events

  cudaStream_t stream;
  checkCUDA(get_legion_stream(&stream));
  checkCUDA(cublasSetStream(simulator->handler.blas, stream));
  checkCUDNN(cudnnSetStream(simulator->handler.dnn, stream));

  const_cast<FFModel*>(model)->measure(simulator);
  
  delete(simulator);
  delete(machine);
}


LogicalTaskgraphBasedSimulator::LogicalTaskgraphBasedSimulator(const FFModel* model,
  FFHandler handler, Memory memory, MachineModel *machine)
: Simulator(model, handler, memory, machine)
{

}

__host__
void LogicalTaskgraphBasedSimulator::simulation_task(const Task *task,
                                     const std::vector<PhysicalRegion> &regions,
                                     Context ctx, Runtime *runtime)
{
  FFModel* model = *((FFModel**) task->args);
  Memory gpu_mem = Machine::MemoryQuery(Machine::get_machine())
         .only_kind(Memory::GPU_FB_MEM).best_affinity_to(task->target_proc).first();
  // Realm::MemoryImpl* memImpl =
  //     Realm::get_runtime()->get_memory_impl(gpu_mem);
  // Realm::Cuda::GPUFBMemory* memFBImpl = (Realm::Cuda::GPUFBMemory*) memImpl;
  // off_t offset = memFBImpl->alloc_bytes_local(model->config.simulator_work_space_size);
  // void* base_ptr = memFBImpl->get_direct_ptr(offset, 0);
  FlatDegConstraintNetworkTopologyGenerator topo_gen = FlatDegConstraintNetworkTopologyGenerator(model->config.numNodes, model->config.node_degree);
  // BigSwitchNetworkTopologyGenerator topo_gen = BigSwitchNetworkTopologyGenerator(model->config.numNodes);
  
  NetworkedMachineModel *nmachine = new NetworkedMachineModel(model->config.numNodes, 
    model->config.workersPerNode,  
    0,
    topo_gen.generate_topology(),
    gpu_mem.capacity(),
    model->config.iface_bandwidth
  );
  nmachine->set_pcie(false);
  nmachine->set_pipeline(true);

  // SimpleMachineModel* nmachine = new SimpleMachineModel(model->config.numNodes, model->config.workersPerNode, gpu_mem.capacity());

  MachineModel *machine;
  machine = reinterpret_cast<MachineModel*>(nmachine);

  // Assume this task is running on GPU0
  
  Simulator* simulator = new LogicalTaskgraphBasedSimulator(model, model->handlers[0], gpu_mem, machine);
  if (model->config.mfile != "") {
    model->load_measurement(simulator, model->config.mfile);
  }
  DemandHeuristicNetworkOptimizerPlus *dhopt = 
    new DemandHeuristicNetworkOptimizerPlus(machine);
  dhopt->if_cnt = model->config.node_degree;
  simulator->l1optimizer = dhopt;

  // Set cublas/cudnn streams to allow Realm catch the events

  cudaStream_t stream;
  checkCUDA(get_legion_stream(&stream));
  checkCUDA(cublasSetStream(simulator->handler.blas, stream));
  checkCUDNN(cudnnSetStream(simulator->handler.dnn, stream));

  std::map<Op*, ParallelConfig> strategies;
  if (model->config.import_strategy_file.length() > 0) {
    // Load the strategy from config.strategies
    for (size_t l = 0; l < model->layers.size(); l++) {
      MappingTagID key = FFConfig::get_hash_id(std::string(model->layers[l]->name));
      std::map<MappingTagID, ParallelConfig>::const_iterator iter;
      iter = model->config.strategies.find(key);
      if (iter == model->config.strategies.end()) {
        fprintf(stderr, "ERROR: Cannot find strategy for operator %s in "
                "strategy file %s\n", model->layers[l]->name,
                model->config.import_strategy_file.c_str());
      }
      strategies[model->layers[l]] = iter->second;
    }
  } else {
    // Start from data parallel
    for (size_t l = 0; l < model->layers.size(); l++) {
      // uint64_t opsz = std::numeric_limits<uint64_t>::max();
      // for (int i = 0; i < model->layers[l]->numWeights; i++) {
      //   if (model->layers[l]->weights[i].get_volume() < opsz)  {
      //     opsz = model->layers[l]->weights[i].get_volume() * sizeof(float);
      //   }
      // }
      // if (opsz * model->config.numNodes * model->config.workersPerNode < gpu_mem.capacity())
      if (model->layers[l]->op_type != OperatorType::OP_EMBEDDING)
        strategies[model->layers[l]] = model->layers[l]->get_data_parallel_config(*model);
      else
        strategies[model->layers[l]] = model->layers[l]->get_random_parallel_config(*model);
    }
  }
  if (model->config.computationMode == COMP_MODE_TRAINING) {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(TRAINING)\n",
        model->config.search_budget, model->config.search_alpha);
  } else {
    fprintf(stderr, "MCMC search configuration: budget(%zu) alpha(%.8lf) mode(INFERENCE)\n",
        model->config.search_budget, model->config.search_alpha);
  }
  model->optimize(simulator, strategies, model->config.search_budget,
      model->config.search_alpha, model->config.computationMode, model->config.enable_propagation);
  if (model->config.export_strategy_file.length() > 0) {
    fprintf(stderr, "Exporting the best discovered strategy to %s.\n",
        model->config.export_strategy_file.c_str());
    std::map<Op*, ParallelConfig>::const_iterator iter;
    std::map<std::string, ParallelConfig> strategy_output;
    for (iter = strategies.begin(); iter != strategies.end(); iter++) {
      strategy_output[iter->first->name] = iter->second;
    }
    save_strategies_to_file(model->config.export_strategy_file, strategy_output);
    fprintf(stderr, "To use the strategy for distributed training, restart"
        " FlexFlow and import the strategy (i.e., --import %s)\n",
        model->config.export_strategy_file.c_str());
    exit(0);
  }  else {
    fprintf(stderr, "The best discovered strategy is not exported.\n"
        "Please set a path to export the strategy using --export or --export-strategy.\n");
    exit(0);
  }
  // Start from data
  // memFBImpl->free_bytes_local(offset, model->config.simulator_work_space_size);
  delete(simulator);
  delete(machine);
}