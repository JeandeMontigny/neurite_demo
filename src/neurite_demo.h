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
#include "random.h"
#include "substance_initializers.h"
#include "model_initializer.h"

#include "neuroscience/compile_time_param.h"
#include "neuroscience/neuron_soma.h"
#include "neuroscience/neurite_element.h"

namespace bdm {

  using namespace std;

  template <typename TSimulation = Simulation<>>
  struct neurite_elongation : public BaseBiologyModule {
    neurite_elongation() : BaseBiologyModule(gAllBmEvents) {}

    template <typename T>
    void Run(T* sim_object) {
      auto&& neurite = sim_object->template ReinterpretCast<experimental::neuroscience::NeuriteElement>();

      neurite->ElongateTerminalEnd(10, {0, 0, 1});
    }

    ClassDefNV(neurite_elongation, 1);
  };


  // Define compile time parameter
  template <typename Backend>
  struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    using BiologyModules = Variant<neurite_elongation<>>;
    using AtomicTypes = VariadicTypedef<experimental::neuroscience::NeuronSoma, experimental::neuroscience::NeuriteElement>;
    using NeuronSoma = experimental::neuroscience::NeuronSoma;
    using NeuriteElement = experimental::neuroscience::NeuriteElement;
  };


  // define my cell creator
  template <typename Function, typename TSimulation = Simulation<>>
  static void CellCreator(double min, double max, int num_cells,
                          Function cell_builder) {
    auto* sim = TSimulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    // Determine simulation object type which is returned by the cell_builder
    using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

    auto container = rm->template Get<FunctionReturnType>();
    container->reserve(num_cells);

    for (int i = 0; i < num_cells; i++) {
      double x = random->Uniform(min, max);
      double y = random->Uniform(min, max);
      double z = random->Uniform(10);
      auto new_simulation_object = cell_builder({x, y, z});
      container->push_back(new_simulation_object);
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

    // Construct neurons
    auto construct_neurons = [](const array<double, 3>& position) {
      auto* simulation = TSimulation::GetActive();
      auto* random = simulation->GetRandom();
      experimental::neuroscience::NeuronSoma cell(position);
      cell.SetDiameter(random->Uniform(7, 8)); // random diameter

      auto ne = cell.ExtendNewNeurite({0, 0, 1});
      ne->GetSoPtr()->AddBiologyModule(neurite_elongation<>());

      return cell;
    };
    CellCreator(param->min_bound_, param->max_bound_, 10, construct_neurons);


    scheduler->Simulate(100);

    return 0;
    }
} // namespace bdm

#endif  // NEURITE_DEMO_H_
