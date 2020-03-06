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
#ifndef BIOLOGY_MODULES_H_
#define BIOLOGY_MODULES_H_

#include "core/biology_module/biology_module.h"
#include "core/resource_manager.h"
#include "extended-obj.h"

namespace bdm {

// enumerate substances in simulation
enum Substances { kSubstanceApical, kSubstanceBasal };


// ApicalElongation_BM
struct ApicalElongation_BM : public BaseBiologyModule {
  ApicalElongation_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  ApicalElongation_BM(const Event& event, BaseBiologyModule* other,
                      uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new ApicalElongation_BM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new ApicalElongation_BM(*this);
  }

  void Run(SimObject* so) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid(kSubstanceApical);
      init_ = true;
    }

    auto* dendrite = bdm_static_cast<MyNeurite*>(so);
    if (dendrite->GetDiameter() > 0.5) {
      Double3 gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);
      // double concentration = 0;

      double gradient_weight = 0.06;
      double randomness_weight = 0.3;
      double old_direction_weight = 4;

      Double3 random_axis = {random->Uniform(-1, 1), random->Uniform(-1, 1),
                             random->Uniform(-1, 1)};
      auto old_direction = dendrite->GetSpringAxis() * old_direction_weight;
      auto grad_direction = gradient * gradient_weight;
      auto random_direction = random_axis * randomness_weight;

      Double3 new_step_direction =
          old_direction + random_direction + grad_direction;

      dendrite->ElongateTerminalEnd(100, new_step_direction);
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.00071);

      if (dendrite->GetCanBranch() && dendrite->IsTerminal() &&
          dendrite->GetDiameter() > 0.55 && random->Uniform() < 0.030) {
        auto rand_noise = random->template UniformArray<3>(-0.1, 0.1);
        Double3 branch_direction =
            Math::Perp3(dendrite->GetUnitaryAxisDirectionVector() + rand_noise,
                        random->Uniform(0, 1)) +
            dendrite->GetSpringAxis();
        auto* dendrite_2 =
            bdm_static_cast<MyNeurite*>(dendrite->Branch(branch_direction));
        dendrite_2->SetCanBranch(false);
        dendrite_2->SetDiameter(0.65);
      }

    }  // end if diameter
  }    // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  BDM_CLASS_DEF_OVERRIDE(ApicalElongation_BM, 1);
};  // end ApicalElongation_BM


// BasalElongation_BM
struct BasalElongation_BM : public BaseBiologyModule {
  BasalElongation_BM() : BaseBiologyModule(gAllEventIds) {}

  /// Default event constructor
  BasalElongation_BM(const Event& event, BaseBiologyModule* other,
                     uint64_t new_oid = 0)
      : BaseBiologyModule(event, other, new_oid) {}

  /// Create a new instance of this object using the default constructor.
  BaseBiologyModule* GetInstance(const Event& event, BaseBiologyModule* other,
                                 uint64_t new_oid = 0) const override {
    return new BasalElongation_BM(event, other, new_oid);
  }

  /// Create a copy of this biology module.
  BaseBiologyModule* GetCopy() const override {
    return new BasalElongation_BM(*this);
  }

  // TODO: don't copy BM when split (elongate)

  void Run(SimObject* so) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();

    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid(kSubstanceBasal);
      init_ = true;
    }

    auto* dendrite = bdm_static_cast<MyNeurite*>(so);
    if (dendrite->IsTerminal() && dendrite->GetDiameter() > 0.7) {
      Double3 gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);
      // double concentration = 0;

      double gradient_weight = 0.03;
      double randomness_weight = 0.4;
      double old_direction_weight = 6;

      Double3 random_axis = {random->Uniform(-1, 1), random->Uniform(-1, 1),
                             random->Uniform(-1, 1)};

      auto old_direction = dendrite->GetSpringAxis() * old_direction_weight;
      auto grad_direction = gradient * gradient_weight;
      auto random_direction = random_axis * randomness_weight;

      Double3 new_step_direction =
          old_direction + random_direction + grad_direction;

      dendrite->ElongateTerminalEnd(50, new_step_direction);
      dendrite->SetDiameter(dendrite->GetDiameter() - 0.00085);

      if (random->Uniform() < 0.006) {
        dendrite->Bifurcate();
      }
    }  // end if diameter
  }  // end run

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
  BDM_CLASS_DEF_OVERRIDE(BasalElongation_BM, 1);
}; // end BasalElongation_BM


}  // namespace bdm

#endif  // BIOLOGY_MODULES_H_
