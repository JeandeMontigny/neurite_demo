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


  // Define my custom cell MyCell extending NeuronSoma
  BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
    BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, foo_);

  public:
    MyCellExt() {}

    MyCellExt(const array<double, 3>& position) : Base(position) {}

    /// Default event constructor
    template <typename TEvent, typename TOther>
    MyCellExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) {}

    /// Default event handler (exising biology module won't be modified on
    /// any event)
    template <typename TEvent, typename... TOthers>
    void EventHandler(const TEvent& event, TOthers*... others) {
      Base::EventHandler(event, others...);
    }

    void SetFoo(int f) { foo_[kIdx] = f; }
    int GetFoo() const { return foo_[kIdx]; }
    // This function is used by ParaView
    int* GetFooPtr() { return foo_.data(); }

  private:
    vec<int> foo_;
  };


  // Define my custom neurite MyNeurite, which extends NeuriteElement
  BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
    BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1, foo_);

  public:
    MyNeuriteExt() {}
    MyNeuriteExt(const array<double, 3>& position) : Base(position) {}

    using NeuronSoma = typename TCompileTimeParam::NeuronSoma;
    using NeuronSomaSoPtr = ToSoPtr<NeuronSoma>;

    /// Default event constructor
    template <typename TEvent, typename TOther>
    MyNeuriteExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {}

    /// Default event handler
    template <typename TEvent, typename... TOthers>
    void EventHandler(const TEvent& event, TOthers*... others) {
      Base::EventHandler(event, others...);
    }

    void SetFoo(int f) { foo_[kIdx] = f; }
    int GetFoo() const { return foo_[kIdx]; }
    // This function is used by ParaView
    int* GetFooPtr() { return foo_.data(); }

  private:
    vec<int> foo_;
  };


  // biology module for neurite growth
  struct NeuriteElongationBM : public BaseBiologyModule {
    NeuriteElongationBM() : BaseBiologyModule(
      // don't copy biology module for bifurcation event
      {experimental::neuroscience::NeuriteBifurcationEvent::kEventId}
    ) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    NeuriteElongationBM(const TEvent& event, TBm* other, uint64_t new_oid = 0) : BaseBiologyModule(event, other, new_oid) {}

    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent& event, TBms*... others) {
      BaseBiologyModule::EventHandler(event, others...);
    }

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* ne) {
      auto* sim = TSimulation::GetActive();
      auto* random = sim->GetRandom();
      ne->ElongateTerminalEnd(100, {random->Uniform(-1, 1), random->Uniform(-1, 1), 1});
      if (ne->IsTerminal() && random->Uniform() < 0.01) {
        ne->Bifurcate();
      }
    }

    ClassDefNV(NeuriteElongationBM, 1);
  };


  // biology module for neurite creation
  struct NeuriteCreationBM: public BaseBiologyModule {
    NeuriteCreationBM() : BaseBiologyModule(gNullEventId) {}

    /// Default event constructor
    template <typename TEvent, typename TBm>
    NeuriteCreationBM(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

    template <typename TEvent, typename... TBms>
    void EventHandler(const TEvent&, TBms*...) {}

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* soma) {

      if (!init_) {
        auto&& ne = soma->ExtendNewNeurite({0, 0, 1});
        ne->AddBiologyModule(NeuriteElongationBM());
        init_ = true;
      }
    }

  private:
    bool init_ = false;
    ClassDefNV(NeuriteCreationBM, 1);
  };


  // define compile time parameter
  BDM_CTPARAM(experimental::neuroscience) {
    BDM_CTPARAM_HEADER(experimental::neuroscience);

    using SimObjectTypes = CTList<MyCell, MyNeurite>;
    using NeuronSoma = MyCell;
    using NeuriteElement = MyNeurite;

    BDM_CTPARAM_FOR(bdm, MyCell) {
      using BiologyModules = CTList<NeuriteCreationBM>;
    };

    BDM_CTPARAM_FOR(bdm, MyNeurite) {
      using BiologyModules = CTList<NeuriteElongationBM>;
    };
  };


  inline int Simulate(int argc, const char** argv) {

    auto set_param = [&](auto* param) {
      param->neurite_max_length_ = 2.0;
    };

    Simulation<> simulation(argc, argv, set_param);

    auto* rm = simulation.GetResourceManager();
    std::array<double, 3> pos = {0,0,0};
    auto&& soma = rm->New<MyCell>(pos);
    soma.SetDiameter(10);
    soma.AddBiologyModule(NeuriteCreationBM());

    simulation.GetScheduler()->Simulate(200);

    return 0;
  }
}  // namespace bdm

#endif  // NEURITE_DEMO_H_
