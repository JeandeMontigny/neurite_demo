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
#ifndef EXTENDED_OBJ_H_
#define EXTENDED_OBJ_H_

#include "core/sim_object/sim_object.h"
#include "neuroscience/neurite_element.h"
#include "neuroscience/neuron_soma.h"

namespace bdm {

  // Define my custom cell object MyCell extending NeuronSoma
  BDM_SIM_OBJECT(MyCell, experimental::neuroscience::NeuronSoma) {
    BDM_SIM_OBJECT_HEADER(MyCell, experimental::neuroscience::NeuronSoma, 1, labelSWC_);

   public:
    // creators
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

    // define my custom functions for MyCell
    // set labelSWC_ to specified value
    inline void SetLabel(int label) { labelSWC_[kIdx] = label; }
    // return labelSWC_ value
    inline int GetLabel() const { return labelSWC_[kIdx]; }
    // increase by 1 the value of labelSWC_
    inline void IncreaseLabel() { labelSWC_[kIdx] = labelSWC_[kIdx] + 1; }

   private:
    // new attribure for MyCell
    vec<int> labelSWC_;
  };

  // Define my custom neurite MyNeurite, which extends NeuriteElement
  BDM_SIM_OBJECT(MyNeurite, experimental::neuroscience::NeuriteElement) {
    BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement, 1,
      can_branch_, its_soma_);

   public:
    // creators
    MyNeuriteExt() {}
    MyNeuriteExt(const array<double, 3>& position) : Base(position) {}

    using NeuronSoma = typename TCompileTimeParam::NeuronSoma;
    using NeuronSomaSoPtr = ToSoPtr<NeuronSoma>;

    /// Default event constructor
    template <typename TEvent, typename TOther>
    MyNeuriteExt(const TEvent& event, TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {
      its_soma_[kIdx] = other->its_soma_[other->kIdx];
    }

    template <typename TOther>
    MyNeuriteExt(const experimental::neuroscience::NewNeuriteExtensionEvent& event,
      TOther* other, uint64_t new_oid = 0) : Base(event, other, new_oid) {
      its_soma_[kIdx] = other->GetSoPtr();
  }

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


    void SetMySoma(NeuronSomaSoPtr soma) { its_soma_[kIdx] = soma; }
    NeuronSomaSoPtr GetMySoma() { return its_soma_[kIdx]; }

   private:
    // new attribure for MyNeurite
    vec<bool> can_branch_;
    vec<NeuronSomaSoPtr> its_soma_;
  };
}  // namespace bdm

#endif  // MY_NEURITE_H_
