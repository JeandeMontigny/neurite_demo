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

// Define my custom neurite MyNeurite, which extends NeuriteElement
class MyNeurite : public experimental::neuroscience::NeuriteElement {
  BDM_SIM_OBJECT_HEADER(MyNeurite, experimental::neuroscience::NeuriteElement,
                        1, can_branch_);

  public:
    MyNeurite() : Base() {}

    virtual ~MyNeurite() {}

    /// Default event constructor
    MyNeurite(const Event& event, SimObject* other, uint64_t new_oid = 0)
    : Base(event, other, new_oid) {}

    // Default event handler
    void EventHandler(const Event& event, SimObject* other1,
      SimObject* other2 = nullptr) {
        Base::EventHandler(event, other1, other2);
      }

    void SetCanBranch(int b) { can_branch_ = b; }
    bool GetCanBranch() const { return can_branch_; }

  private:
    bool can_branch_ = false;
};

}  // namespace bdm

#endif  // EXTENDED_OBJ_H_
