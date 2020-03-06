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
#include "core/substance_initializers.h"
#include "biology-modules.h"
#include "extended-obj.h"
#include "util-methods.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {

  // number of cells in simulation
  int num_cells = 4;
  // size of simulation space
  int simulation_space = 500;
  auto set_param = [&](Param* param) {
    // Create an artificial bounds for the simulation space
    param->bound_space_ = true;
    param->min_bound_ = -simulation_space/2;
    param->max_bound_ = simulation_space/2;
  };

  // initialise neuroscience modlues
  experimental::neuroscience::InitModule();

  Simulation simulation(argc, argv, set_param);

  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();
  random->SetSeed(8794);

  // create chemical substances
  SubstanceCreator(param->max_bound_);

  // create cells
  CellCreator(param->min_bound_, param->max_bound_, num_cells);

  // simulate
  scheduler->Simulate(500);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;

} // end Simulate

}  // namespace bdm

#endif  // NEURITE_DEMO_H_
