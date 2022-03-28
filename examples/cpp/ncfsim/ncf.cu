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
#include "cuda_helper.h"

void DataLoader::load_sparse_input(const Task *task,
                                   const std::vector<PhysicalRegion> &regions,
                                   Context ctx,
                                   Runtime* runtime)
{

}

void DataLoader::load_dense_input(const Task *task,
                                  const std::vector<PhysicalRegion> &regions,
                                  Context ctx,
                                  Runtime* runtime)
{

}

void DataLoader::load_label(const Task *task,
                            const std::vector<PhysicalRegion>& regions,
                            Context ctx,
                            Runtime* runtime)
{

}
