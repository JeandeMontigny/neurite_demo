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
#include "substance_initializers.h"

using namespace std;

namespace bdm {

// enumerate substance in simulation
enum Substances { substance_apical, substance_basal };

// Define my custom neurite MyNeurite, which extends NeuriteElement
BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
  BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1, can_branch_);

 public:
  // creators
  MyNeuriteExt() {}
  MyNeuriteExt(const array<double, 3>& position) : Base(position) {}

  /// Default event constructor
  template <typename TEvent, typename TOther>
  MyNeuriteExt(const TEvent& event, TOther* other, uint64_t new_oid = 0)
      : Base(event, other, new_oid) {}

  /// Default event handler
  template <typename TEvent, typename... TOthers>
  void EventHandler(const TEvent& event, TOthers*... others) {
    Base::EventHandler(event, others...);
  }

  // define my custom functions for MyNeurite
  // set can_branch_ to specified bool
  void SetCanBranch(int b) { can_branch_[kIdx] = b; }
  // return can_branch_ value
  bool GetCanBranch() const { return can_branch_[kIdx]; }

 private:
  // new attribure for MyNeurite
  vec<bool> can_branch_;
};


// define BiologyModule. here, behaviour of apical dendrites
struct ApicalElongation_BM : public BaseBiologyModule {
  ApicalElongation_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  ApicalElongation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  // TODO: don't copy BM when split (elongate)

  // define Run() method, that will be called at each time step
  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* dendrite) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    // if diffusion grid hasn't been initialised
    if (!init_) {
      // initialised diffusion grid
      dg_guide_ = rm->GetDiffusionGrid("substance_apical");
      init_ = true;
    }

    // if the dendrite diameter is higher than 0.5
    if (dendrite->GetDiameter() > 0.5) {
      array<double, 3> gradient;
      // get gradient at dendrite position
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

      // define some parameters for next step direction
      double gradientWeight = 0.04;
      double randomnessWeight = 0.3;
      double oldDirectionWeight = 3;
      // define a random direction
      array<double, 3> random_axis = {random->Uniform(-1, 1),
                                      random->Uniform(-1, 1),
                                      random->Uniform(-1, 1)};
      auto oldDirection =
          Math::ScalarMult(oldDirectionWeight, dendrite->GetSpringAxis());
      auto gradDirection = Math::ScalarMult(gradientWeight, gradient);
      auto randomDirection = Math::ScalarMult(randomnessWeight, random_axis);

      // calculate step direction
      array<double, 3> newStepDirection =
          Math::Add(Math::Add(oldDirection, randomDirection), gradDirection);

      // actually elongate the dendrite in new step direction
      dendrite->ElongateTerminalEnd(25, newStepDirection);
      // reduce neurite diameter as it elongates
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.0002);

      // if the dendrite can branch and is a terminal segment
      if (dendrite->GetCanBranch() && dendrite->IsTerminal() &&
          random->Uniform() < 0.033) {
        // define a random direction
        auto rand_noise = random->template UniformArray<3>(-0.1, 0.1);
        // calculate new direction
        array<double, 3> branchDirection = Math::Add(
            Math::Perp3(Math::Add(dendrite->GetUnitaryAxisDirectionVector(),
                                  rand_noise),
                        random->Uniform(0, 1)),
            dendrite->GetSpringAxis());
        // branch neurite 1, and get its new daugher, dendrite_2
        auto dendrite_2 = dendrite->Branch(branchDirection);
        // set values for this new dendrite
        dendrite_2->SetCanBranch(false);
        // set diameter of daugher to diameter of it's mother
        dendrite_2->SetDiameter(dendrite->GetDiameter());
      }

    }  // end if diameter
  }    // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  ClassDefNV(ApicalElongation_BM, 1);
};  // end ApicalElongation_BM

// define BiologyModule. here, behaviour of basal dendrites
struct BasalElongation_BM : public BaseBiologyModule {
  BasalElongation_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  BasalElongation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  // TODO: don't copy BM when split (elongate)

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* dendrite) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid("substance_basal");
      init_ = true;
    }

    if (dendrite->IsTerminal() && dendrite->GetDiameter() > 0.75) {
      array<double, 3> gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

      double gradientWeight = 0.02;
      double randomnessWeight = 0.5;
      double oldDirectionWeight = 5;

      array<double, 3> random_axis = {random->Uniform(-1, 1),
                                      random->Uniform(-1, 1),
                                      random->Uniform(-1, 1)};
      auto oldDirection =
          Math::ScalarMult(oldDirectionWeight, dendrite->GetSpringAxis());
      auto gradDirection = Math::ScalarMult(gradientWeight, gradient);
      auto randomDirection = Math::ScalarMult(randomnessWeight, random_axis);

      array<double, 3> newStepDirection =
          Math::Add(Math::Add(oldDirection, randomDirection), gradDirection);

      dendrite->ElongateTerminalEnd(25, newStepDirection);
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.001);

      if (random->Uniform() < 0.008) {
        dendrite->SetDiameter(dendrite->GetDiameter() - 0.01);
        dendrite->Bifurcate();
      }

    }  // end if diameter

  }  // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  ClassDefNV(BasalElongation_BM, 1);
};  // end BasalElongation_BM


// Define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using NeuriteElement = MyNeurite;
  using SimObjectTypes = CTList<experimental::neuroscience::NeuronSoma, MyNeurite>;

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules = CTList<ApicalElongation_BM, BasalElongation_BM>;
  };

};  // end BDM_CTPARAM

// main function, called by neurite_demo.cc
inline int Simulate(int argc, const char** argv) {

  // define some parameters for simulation
  auto set_param = [&](auto* param) {
    // our simulation will have boundaries
    param->bound_space_ = true;
    // set boundaries values
    param->min_bound_ = -200;
    param->max_bound_ = 200;
    // max length for neurite segments
    param->neurite_max_length_ = 2.0;
  };

  // create and initialised the simulation
  Simulation<> simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  // set seed for random number generation
  random->SetSeed(random->Uniform(0, 1e8));

  // define a lambda function for neuron creation
  auto neuron_builder = [&rm](const std::array<double, 3>& position) {
    experimental::neuroscience::NeuronSoma soma(position);
    soma.SetDiameter(10);
    auto soma_soptr = soma.GetSoPtr();
    rm->push_back(soma);

    auto&& dendrite_apical = soma_soptr->ExtendNewNeurite({0, 0, 1});
    dendrite_apical->AddBiologyModule(ApicalElongation_BM());
    dendrite_apical->SetCanBranch(true);

    auto&& dendrite_basal1 = soma_soptr->ExtendNewNeurite({0, 0, -1});
    dendrite_basal1->AddBiologyModule(BasalElongation_BM());

    auto&& dendrite_basal2 = soma_soptr->ExtendNewNeurite({0, 0.6, -0.8});
    dendrite_basal2->AddBiologyModule(BasalElongation_BM());

    auto&& dendrite_basal3 = soma_soptr->ExtendNewNeurite({0.3, -0.6, -0.8});
    dendrite_basal3->AddBiologyModule(BasalElongation_BM());
  };

  // call function for neuron creation
  neuron_builder({0, 0, 0});

  // define substance for neurite attraction
  ModelInitializer::DefineSubstance(0, "substance_apical", 0, 0,
                                    param->max_bound_ / 2);
  ModelInitializer::DefineSubstance(1, "substance_basal", 0, 0,
                                    param->max_bound_ / 2);
  // create substance with gaussian distribution for neurite attraction
  ModelInitializer::InitializeSubstance(
      0, "substance_apical",
      GaussianBand(param->max_bound_, 200, Axis::kZAxis));
  ModelInitializer::InitializeSubstance(
      1, "substance_basal", GaussianBand(param->min_bound_, 200, Axis::kZAxis));

  auto* scheduler = simulation.GetScheduler();

  // simuate for 500 steps
  scheduler->Simulate(500);

  std::cout << "Simulation completed successfully!" << std::endl;

  return 0;
}  // end Simulate

} // namespace bdm

#endif  // NEURITE_DEMO_H_
