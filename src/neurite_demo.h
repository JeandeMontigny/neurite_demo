// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef NEURITE_DEMO_H_
#define NEURITE_DEMO_H_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

using namespace std;

struct NeuriteElongation : public BaseBiologyModule {
  NeuriteElongation() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    auto&& neurite = sim_object->template ReinterpretCast<
        experimental::neuroscience::NeuriteElement>();

    neurite->ElongateTerminalEnd(10, {0, 0, 1});
  }

  ClassDefNV(NeuriteElongation, 1);
};

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<NeuriteElongation>;
  using AtomicTypes =
      VariadicTypedef<experimental::neuroscience::NeuronSoma,
                      experimental::neuroscience::NeuriteElement>;
  using NeuronSoma = experimental::neuroscience::NeuronSoma;
  using NeuriteElement = experimental::neuroscience::NeuriteElement;
};

// define my cell creator
template <typename TSimulation = Simulation<>>
static void CellCreator(double min, double max, int num_cells) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto* random = sim->GetRandom();

  using Soma = experimental::neuroscience::NeuronSoma;
  auto* container = rm->template Get<Soma>();
  container->reserve(num_cells);

  for (int i = 0; i < num_cells; i++) {
    double x = random->Uniform(min, max);
    double y = random->Uniform(min, max);
    double z = random->Uniform(10);

    std::array<double, 3> position = {x, y, z};
    auto&& soma = rm->template New<Soma>(position);
    soma.SetDiameter(random->Uniform(7, 8));  // random diameter

    auto&& ne = soma.ExtendNewNeurite({0, 0, 1});
    ne->AddBiologyModule(NeuriteElongation());
  }
  container->Commit();

}  // end CellCreator

template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {
  Simulation<> simulation(argc, argv);
  auto* param = simulation.GetParam();
  auto* scheduler = simulation.GetScheduler();

  param->bound_space_ = true;
  param->min_bound_ = 0;
  param->max_bound_ = 50;

  CellCreator(param->min_bound_, param->max_bound_, 10);

  scheduler->Simulate(100);

  return 0;
}
}  // namespace bdm

#endif  // NEURITE_DEMO_H_
