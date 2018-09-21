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

//template <typename TSimulation = Simulation<>>
// struct NeuriteElongation : public BaseBiologyModule {
//   NeuriteElongation() : BaseBiologyModule(gAllEventIds) {}
//
//   /// Default event constructor
//   template <typename TEvent, typename TBm>
//     NeuriteElongation(const TEvent& event, TBm* other, uint64_t new_oid = 0) {
//   }
//
//   template <typename TEvent, typename... TBms>
//   void EventHandler(const TEvent&, TBms*...) {}
//
//   template <typename T, typename TSimulation = Simulation<>>
//   void Run(T* sim_object) {
//     auto* sim = TSimulation::GetActive();
//     auto&& ne = sim_object->template ReinterpretCast<experimental::neuroscience::NeuriteElement>();
//
//     int ownType = 0;
//     int otherType = 0;
//     auto countNeighbours = [&](auto&& neighbor, SoHandle neighbor_handle) {
//       // if neighbor is a NeuriteElement
//       if (neighbor->template IsSoType<experimental::neuroscience::NeuriteElement>()) {
//         auto&& neighbor_rc = neighbor->template
//           ReinterpretCast<experimental::neuroscience::NeuriteElement>();
//         auto n_soptr = neighbor_rc->GetSoPtr();
//         // if it is a direct relative
//         if (n_soptr->GetNeuronSomaOfNeurite() == ne->GetNeuronSomaOfNeurite()) {
//           ownType++;
//         }
//         else {
//           otherType++;
//         }
//       }
//     };
//
//     auto* grid = sim->GetGrid();
//     grid->ForEachNeighborWithinRadius(countNeighbours, ne, ne->GetSoHandle(), 4);
//     if (ownType >= otherType) {
//       ne->ElongateTerminalEnd(10, {0, 0, 1});
//     }
//   }
//
//   ClassDefNV(NeuriteElongation, 1);
// };

// define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  // using SimObjectTypes = CTList<NeuronSoma, NeuriteElement>;
  using NeuronSoma = experimental::neuroscience::NeuronSoma;
  // using NeuriteElement = experimental::neuroscience::NeuriteElement;

  BDM_CTPARAM_FOR(experimental::neuroscience, NeuriteElement) {
    // using BiologyModules =
    //     CTList<NeuriteElongation>;
  };
};

// define my cell creator
// template <typename TSimulation = Simulation<>>
// static void CellCreator(double min, double max, int num_cells) {
//   auto* sim = TSimulation::GetActive();
//   auto* rm = sim->GetResourceManager();
//   auto* random = sim->GetRandom();
//
//   using Soma = experimental::neuroscience::NeuronSoma;
//   auto* container = rm->template Get<Soma>();
//   container->reserve(num_cells);
//
//   for (int i = 0; i < num_cells; i++) {
//     double x = random->Uniform(min, max);
//     double y = random->Uniform(min, max);
//     double z = random->Uniform(10);
//
//     std::array<double, 3> position = {x, y, z};
//     auto&& soma = rm->template New<Soma>(position);
//     soma.SetDiameter(random->Uniform(7, 8));  // random diameter
//
//     auto&& ne = soma.ExtendNewNeurite({0, 0, 1});
//     ne->AddBiologyModule(NeuriteElongation());
//   }
//   container->Commit();
//
// }  // end CellCreator

template <typename TSimulation = Simulation<>>
inline int Simulate(int argc, const char** argv) {

  auto set_param = [&](auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 150;
  };

  Simulation<> simulation(argc, argv, set_param);
  auto* param = simulation.GetParam();
  auto* scheduler = simulation.GetScheduler();
  auto* rm = simulation.GetResourceManager();

//  CellCreator(param->min_bound_, param->max_bound_, 10);

  // auto* container = rm->template Get<experimental::neuroscience::NeuronSoma>();
  // container->reserve(1);
  // auto&& cell = rm->template
  //  New<experimental::neuroscience::NeuronSoma>({0,0,0});
  // container->Commit();

  experimental::neuroscience::NeuronSoma soma({0,0,0});
  // auto&& ne = soma.ExtendNewNeurite({0, 0, 1});
  // ne->AddBiologyModule(NeuriteElongation());

  simulation.GetScheduler()->Simulate(100);

  return 0;
}
}  // namespace bdm

#endif  // NEURITE_DEMO_H_
