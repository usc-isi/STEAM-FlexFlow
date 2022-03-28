/* Copyright 2020 Stanford, Facebook
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

#include "ncf.h"
#include "hdf5.h"
#include <sstream>

using namespace Legion;

LegionRuntime::Logger::Category log_app("NCF");

void parse_input_args(char **argv, int argc, NCFConfig& apConfig);

Tensor create_mlp(FFModel* model, const Tensor& input,
                  std::vector<int> ln, int sigmoid_layer)
{
  Tensor t = input;
  for (int i = 0; i < (int)(ln.size()-1); i++) {
    float std_dev = sqrt(2.0f / (ln[i+1] + ln[i]));
    Initializer* weight_init = new NormInitializer(std::rand(), 0, std_dev);
    std_dev = sqrt(2.0f / ln[i+1]);
    Initializer* bias_init = new NormInitializer(std::rand(), 0, std_dev);
    ActiMode activation = i == sigmoid_layer ? AC_MODE_SIGMOID : AC_MODE_RELU;
    t = model->dense(t, ln[i+1], activation, true/*bias*/, NULL/*weight_sharing*/, weight_init, bias_init);
  }
  return t;
}

Tensor create_emb(FFModel* model, const Tensor& input,
                  int input_dim, int output_dim, int idx)
{
  float range = sqrt(1.0f / input_dim);
  Initializer* embed_init = new UniformInitializer(std::rand(), -range, range);
  return model->embedding(input, input_dim, output_dim, AGGR_MODE_SUM, NULL/*weight_sharing*/, embed_init);
}

Tensor interact_features(FFModel* model,
                         const std::vector<Tensor>& ly,
                         std::string interaction)
{
  // Currently only support cat
  // TODO: implement dot attention
  if (interaction == "cat") {
    Tensor* inputs = (Tensor*) malloc(sizeof(Tensor) * (ly.size()));
    for (size_t i = 0; i < ly.size(); i++)
      inputs[i] = ly[i];
    return model->concat(ly.size(), inputs, 1/*axis*/);
    free(inputs);
  } else {
    assert(false);
  }
}

void print_vector(const std::string& name, const std::vector<int>& vector)
{
  std::ostringstream out;
  for (size_t i = 0; i < vector.size() - 1; i++)
    out << vector[i] << " ";
  if (vector.size() > 0)
    out << vector[vector.size() - 1];
  log_app.print("%s: %s", name.c_str(), out.str().c_str());
}

void top_level_task(const Task* task,
                    const std::vector<PhysicalRegion>& regions,
                    Context ctx, Runtime* runtime)
{
  FFConfig ffConfig;
  // Parse input arguments
  NCFConfig ncfConfig;
  {
    const InputArgs &command_args = HighLevelRuntime::get_input_args();
    char **argv = command_args.argv;
    int argc = command_args.argc;
    parse_input_args(argv, argc, ncfConfig);
    log_app.print("batchSize(%d) workersPerNodes(%d) numNodes(%d)",
        ffConfig.batchSize, ffConfig.workersPerNode, ffConfig.numNodes);
    log_app.print("EmbeddingBagSize(%d)", ncfConfig.embedding_bag_size);
    print_vector("MLP Top", ncfConfig.mlp_top);
  }

  ffConfig.numNodes = ncfConfig.nsimnode;
  FFModel ff(ffConfig, true);

  std::vector<Tensor> user_gmf_input;
  std::vector<Tensor> item_gmf_input;
  std::vector<Tensor> user_mlp_input;
  std::vector<Tensor> item_mlp_input;

  for (size_t i = 0; i < ncfConfig.ntables; i++) {
    const int dims[] = {ffConfig.batchSize, ncfConfig.embedding_bag_size};
    Tensor input = ff.create_tensor<2>(dims, DT_INT64);
    user_gmf_input.push_back(input);
    item_gmf_input.push_back(input);
    user_mlp_input.push_back(input);
    item_mlp_input.push_back(input);
  }

  std::vector<Tensor> ly_user_gmf;
  for (size_t i = 0; i < ncfConfig.ntables; i++) {
    int input_dim = ncfConfig.user_num_per_table;
    int output_dim = ncfConfig.factor_num;
    ly_user_gmf.push_back(create_emb(&ff, user_gmf_input[i], input_dim, output_dim, i));
  }
  std::vector<Tensor> ly_item_gmf;
  for (size_t i = 0; i < ncfConfig.ntables; i++) {
    int input_dim = ncfConfig.item_num_per_table;
    int output_dim = ncfConfig.factor_num;
    ly_item_gmf.push_back(create_emb(&ff, item_gmf_input[i], input_dim, output_dim, i));
  }
  std::vector<Tensor> ly_mlp;
  for (size_t i = 0; i < ncfConfig.ntables; i++) {
    int input_dim = ncfConfig.user_num_per_table;
    int output_dim = ncfConfig.embed_dim_per_table;
    ly_mlp.push_back(create_emb(&ff, user_mlp_input[i], input_dim, output_dim, i));
  }
  std::vector<Tensor> ly_item_mlp;
  for (size_t i = 0; i < ncfConfig.ntables; i++) {
    int input_dim = ncfConfig.item_num_per_table;
    int output_dim = ncfConfig.embed_dim_per_table;
    ly_mlp.push_back(create_emb(&ff, item_mlp_input[i], input_dim, output_dim, i));
  }
  
  Tensor z_mlp = interact_features(&ff, ly_mlp, ncfConfig.arch_interaction_op);
  Tensor z_user_gmf = interact_features(&ff, ly_user_gmf, ncfConfig.arch_interaction_op);
  Tensor z_item_gmf = interact_features(&ff, ly_item_gmf, ncfConfig.arch_interaction_op);

  Tensor elem_prod = ff.multiply(z_user_gmf, z_item_gmf, false, "gmf");
  Tensor activation_gmf;
  {
    float std_dev = sqrt(2.0f);
    Initializer* weight_init = new NormInitializer(std::rand(), 0, std_dev);
    activation_gmf = ff.dense(elem_prod, 1, 
      AC_MODE_SIGMOID, false/*bias*/, NULL/*weight_sharing*/, weight_init, NULL);
  }
  Tensor p = create_mlp(&ff, z_mlp, ncfConfig.mlp_top, ncfConfig.mlp_top.size() - 2);

  std::vector<Tensor> final_concat_input {activation_gmf, p};
  Tensor final_concat = interact_features(&ff, final_concat_input, 
    ncfConfig.arch_interaction_op);

  Tensor out;
  {
    float std_dev = sqrt(2.0f);
    Initializer* weight_init = new NormInitializer(std::rand(), 0, std_dev);
    out = ff.dense(final_concat, 1, 
      AC_MODE_SIGMOID, true/*bias*/, NULL/*weight_sharing*/, weight_init, weight_init);
  }
  
  if (ffConfig.measurement_only) {
    ff.run_measurement();
  }
  else {
    ff.simulate();
  }

}

void parse_input_args(char **argv, int argc, NCFConfig& config)
{
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--unum_per_table")) {
      config.user_num_per_table = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--inum_per_table")) {
      config.item_num_per_table = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--factor_num")) {
      config.factor_num = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--edim_per_table")) {
      config.factor_num = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--ntables")) {
      config.ntables = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--embedding-bag-size")) {
      config.embedding_bag_size = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--arch-mlp-top")) {
      std::stringstream ss{std::string(argv[++i])};
      std::string word;
      config.mlp_top.clear();
      while (std::getline(ss, word, '-')) {
        config.mlp_top.push_back(std::stoi(word));
      }
      continue;
    }
    if (!strcmp(argv[i], "--loss-threshold")) {
      config.loss_threshold = atof(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--sigmoid-top")) {
      config.sigmoid_top = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--sigmoid-bot")) {
      config.sigmoid_bot = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--arch-interaction-op")) {
      config.arch_interaction_op = std::string(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--dataset")) {
      config.dataset_path = std::string(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--data-size")) {
      config.data_size = atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--nsimnode")) {
      config.nsimnode = std::atoi(argv[++i]);
      continue;
    }
    if (!strcmp(argv[i], "--nsimgpu")) {
      config.nsimgpu = std::atoi(argv[++i]);
      continue;
    } 
  }
}

DataLoader::DataLoader(FFModel& ff,
                       const NCFConfig& ncf,
                       const std::vector<Tensor>& _sparse_inputs,
                       Tensor _dense_input, Tensor _label)
{
 
}

void DataLoader::load_entire_dataset(const Task *task,
                                     const std::vector<PhysicalRegion> &regions,
                                     Context ctx,
                                     Runtime* runtime)
{
 
}

void DataLoader::next_batch(FFModel& ff)
{
  
}

void DataLoader::shuffle()
{}

void DataLoader::reset()
{
  next_index = 0;
}

void DataLoader::load_sparse_input_cpu(const Task *task,
                                   const std::vector<PhysicalRegion> &regions,
                                   Context ctx,
                                   Runtime* runtime)
{
  std::cout << "load_sparse_input_cpu" << std::endl;
}

void register_custom_tasks()
{
  // Load entire dataset
  {
    TaskVariantRegistrar registrar(CUSTOM_CPU_TASK_ID_1, "Load Entire Dataset");
    registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<DataLoader::load_entire_dataset>(
        registrar, "Load Entire Dataset Task");
  }
  // Load Sparse Inputs
  {
    TaskVariantRegistrar registrar(CUSTOM_GPU_TASK_ID_1, "Load Sparse Inputs");
    registrar.add_constraint(ProcessorConstraint(Processor::TOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<DataLoader::load_sparse_input>(
        registrar, "Load Sparse Inputs Task");
  }
  // Load Dense Inputs
  {
    TaskVariantRegistrar registrar(CUSTOM_GPU_TASK_ID_2, "Load Dense Inputs");
    registrar.add_constraint(ProcessorConstraint(Processor::TOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<DataLoader::load_dense_input>(
        registrar, "Load Dense Inputs Task");
  }
  // Load Labels
  {
    TaskVariantRegistrar registrar(CUSTOM_GPU_TASK_ID_3, "Load Labels");
    registrar.add_constraint(ProcessorConstraint(Processor::TOC_PROC));
    registrar.set_leaf();
    Runtime::preregister_task_variant<DataLoader::load_label>(
        registrar, "Load Labels");
  }
}
