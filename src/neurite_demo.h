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

struct NeuriteElongationBM : public BaseBiologyModule {
  NeuriteElongationBM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  NeuriteElongationBM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename TEvent, typename... TBms>
  void EventHandler(const TEvent&, TBms*...) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* ne) {
    // auto* sim = TSimulation::GetActive();
    ne->ElongateTerminalEnd(10, {0, 0, 1});
  }

  ClassDefNV(NeuriteElongationBM, 1);
};


// define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using SimObjectTypes = CTList<experimental::neuroscience::NeuronSoma,
    experimental::neuroscience::NeuriteElement>;

  BDM_CTPARAM_FOR(experimental::neuroscience, NeuriteElement) {
    using BiologyModules =
        CTList<NeuriteElongationBM>;
  };
};


inline int Simulate(int argc, const char** argv) {

  auto set_param = [&](auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 150;
  };

  Simulation<> simulation(argc, argv, set_param);

  experimental::neuroscience::NeuronSoma soma({75,75,10});
  auto&& ne = soma.ExtendNewNeurite({0, 0, 1});
  ne->AddBiologyModule(NeuriteElongationBM());

  simulation.GetScheduler()->Simulate(100);

  return 0;
}
}  // namespace bdm

#endif  // NEURITE_DEMO_H_
