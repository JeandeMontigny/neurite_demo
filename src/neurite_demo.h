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

enum Substances { substance_apical, substance_basal };

// Define my custom cell MyCell extending NeuronSoma
BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
  BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, labelSWC_);

 public:
  MyCellExt() {}
  MyCellExt(const array<double, 3>& position) : Base(position) {}

  /// Default event constructor
  template <typename TEvent, typename TOther>
  MyCellExt(const TEvent& event, TOther* other, uint64_t new_oid = 0)
    : Base(event, other, new_oid) {}

  /// Default event handler
  template <typename TEvent, typename... TOthers>
  void EventHandler(const TEvent& event, TOthers*... others) {
    Base::EventHandler(event, others...);
  }

  inline void SetLabel(int label) { labelSWC_[kIdx] = label; }
  inline int GetLabel() const { return labelSWC_[kIdx]; }
  inline void IncreaseLabel() { labelSWC_[kIdx] = labelSWC_[kIdx] + 1; }

 private:
  vec<int> labelSWC_;
};

// Define my custom neurite MyNeurite, which extends NeuriteElement
BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
  BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1, can_branch_);

 public:
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

  void SetCanBranch(int b) { can_branch_[kIdx] = b; }
  bool GetCanBranch() const { return can_branch_[kIdx]; }

 private:
  vec<bool> can_branch_;
};


struct ApicalElongation_BM : public BaseBiologyModule {
  ApicalElongation_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  ApicalElongation_BM(const TEvent& event, TBm* other, uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  // TODO: don't copy BM when split (elongate)

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* dendrite) {
    auto* sim = TSimulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid("substance_apical");
      init_ = true;
    }

    if (dendrite->GetDiameter() > 0.5) {
      array<double, 3> gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);
      // double concentration = 0;

      double gradientWeight = 0.04;
      double randomnessWeight = 0.3;
      double oldDirectionWeight = 3;

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

      if (dendrite->GetCanBranch() && dendrite->IsTerminal() &&
          random->Uniform() < 0.033) {
        auto rand_noise = random->template UniformArray<3>(-0.1, 0.1);
        array<double, 3> branchDirection = Math::Add(
            Math::Perp3(Math::Add(dendrite->GetUnitaryAxisDirectionVector(),
                                  rand_noise),
                        random->Uniform(0, 1)),
            dendrite->GetSpringAxis());
        auto dendrite_2 = dendrite->Branch(branchDirection);
        dendrite_2->SetCanBranch(false);
        dendrite_2->SetDiameter(0.65);
      }

    }  // end if diameter
  }    // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  ClassDefNV(ApicalElongation_BM, 1);
};  // end ApicalElongation_BM

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
      // double concentration = 0;

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


template <typename TSimulation = Simulation<>>
inline void morpho_exporteur() {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto* param = sim->GetParam();
  int seed = sim->GetRandom()->GetSeed();

  rm->ApplyOnAllElements([&](auto&& so, SoHandle) {
    if (so->template IsSoType<MyCell>()) {
      auto&& cell = so.template ReinterpretCast<MyCell>();
      int thisCellType = cell.GetCellType();
      auto cellPosition = cell.GetPosition();
      ofstream swcFile;
      string swcFileName = Concat("./cell", cell.GetUid(),
                                  "_type", thisCellType, "_seed", seed, ".swc")
                               .c_str();
      swcFile.open(swcFileName);
      cell->SetLabel(1);
      // swcFile << labelSWC_ << " 1 " << cellPosition[0] << " "
      //         << cellPosition[1]  << " " << cellPosition[2] << " "
      //         << cell->GetDiameter()/2 << " -1";
      swcFile << cell->GetLabel() << " 1 0 0 0 " << cell->GetDiameter() / 2
              << " -1";

      for (auto& ne : cell->GetDaughters()) {
        swcFile << swc_neurites(ne, 1, cellPosition);
      }  // end for neurite in cell
      swcFile.close();
    }
  });  // end for cell in simulation
  std::cout << "swc export done" << std::endl;
}  // end morpho_exporteur


template <typename T>
inline string swc_neurites(T ne, int labelParent,
                           array<double, 3> somaPosition) {
  array<double, 3> nePosition = ne->GetPosition();
  nePosition[0] = nePosition[0] - somaPosition[0];
  nePosition[1] = nePosition[1] - somaPosition[1];
  nePosition[2] = nePosition[2] - somaPosition[2];
  string temps;

  ne->GetMySoma()->IncreaseLabel();
  // set explicitly the value of GetLabel() other wise it is not properly set
  int currentLabel = ne->GetMySoma()->GetLabel();

  // if branching point
  if (ne->GetDaughterRight() != nullptr) {
    // FIXME: segment indice should be 5, no 3. If set here,
    // it's not the actual branching point, but the following segment
    // need to run correction.py to correct file
    temps =
        Concat(temps, "\n", currentLabel, " 3 ", nePosition[0], " ",
               nePosition[1], " ", nePosition[2], " ", ne->GetDiameter() / 2,
               " ", labelParent,
               swc_neurites(ne->GetDaughterRight(), currentLabel, somaPosition))
            .c_str();
    ne->GetMySoma()->IncreaseLabel();
  }
  // if is straigh dendrite
  // need to update currentLabel
  currentLabel = ne->GetMySoma()->GetLabel();
  if (ne->GetDaughterLeft() != nullptr) {
    temps =
        Concat(temps, "\n", currentLabel, " 3 ", nePosition[0], " ",
               nePosition[1], " ", nePosition[2], " ", ne->GetDiameter() / 2,
               " ", labelParent,
               swc_neurites(ne->GetDaughterLeft(), currentLabel, somaPosition))
            .c_str();
  }
  // if ending point
  if (ne->GetDaughterLeft() == nullptr && ne->GetDaughterRight() == nullptr) {
    temps = Concat(temps, "\n", currentLabel, " 6 ", nePosition[0], " ",
                   nePosition[1], " ", nePosition[2], " ",
                   ne->GetDiameter() / 2, " ", labelParent)
                .c_str();
  }

  return temps;
} // end swc_neurites


// Define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;
  using SimObjectTypes = CTList<MyCell, MyNeurite>;

  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules = CTList<NullBiologyModule>;
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules = CTList<ApicalElongation_BM, BasalElongation_BM>;
  };

};  // end BDM_CTPARAM

inline int Simulate(int argc, const char** argv) {
  uint64_t num_cells = 1;
  auto set_param = [&](auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = -200;
    param->max_bound_ = 200;
    param->neurite_max_length_ = 2.0;
  };

  Simulation<> simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  random->SetSeed(8794);

  auto neuron_builder = [&rm](const std::array<double, 3>& position) {
    MyCell soma(position);
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

  // neuron_builder({0, 0, 0});
  // neuron_builder({20, 20, 0});
  uint64_t cells_per_dim = std::sqrt(num_cells);
  double space = 10;
  for (size_t x = 0; x < cells_per_dim; x++) {
      auto x_pos = x * space;
      for (size_t y = 0; y < cells_per_dim; y++) {
        auto y_pos = y * space;
        neuron_builder({x_pos, y_pos, 0});
      }
  }

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
  // substance initialization happens in first iteration - exclude it from
  // runtime measurments
  scheduler->Simulate(1);

  auto start = Timing::Timestamp();
  scheduler->Simulate(500);
  auto stop = Timing::Timestamp();
  std::cout << "RUNTIME " << (stop - start) << std::endl;

  std::cout << "Simulation completed successfully!" << std::endl;
  std::cout << "num sim objects: " << rm->GetNumSimObjects() << std::endl;
  return 0;
}  // end Simulate

} // namespace bdm

#endif  // NEURITE_DEMO_H_
