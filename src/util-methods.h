#ifndef UTIL_NETHODS_
#define UTIL_NETHODS_

#include "biodynamo.h"
#include "biology-modules.h"
#include "extended-obj.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

  // define my substance creator
  inline void SubstanceCreator(double max_space) {
    ModelInitializer::DefineSubstance(kSubstanceApical, "substance_apical", 0, 0,
                                      max_space / 20);
    ModelInitializer::DefineSubstance(kSubstanceBasal, "substance_basal", 0, 0,
                                      max_space / 20);
    // initialise substance with gaussian distribution for neurite attraction
    ModelInitializer::InitializeSubstance(
        kSubstanceApical, "substance_apical",
        GaussianBand(max_space, 200, Axis::kZAxis));
    ModelInitializer::InitializeSubstance(
        kSubstanceBasal, "substance_basal",
        GaussianBand(max_space, 200, Axis::kZAxis));
  }

  // define my cell creator
  inline void CellCreator(double min_space, double max_space, int num_cells) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    for (int i = 0; i < num_cells; i++) {
      // +- 100 so cells are not too close from boundaries
      double x = random->Uniform(min_space + 100, max_space - 100);
      double y = random->Uniform(min_space + 100, max_space - 100);
      double z = 0;

      // create cell
      experimental::neuroscience::NeuronSoma* cell = new experimental::neuroscience::NeuronSoma({x, y, z});
      cell->SetDiameter(10);
      auto cell_soptr = cell->GetSoPtr<experimental::neuroscience::NeuronSoma>();
      rm->push_back(cell);

       // create dendrites
      MyNeurite my_neurite;
      auto* dendrite_apical = bdm_static_cast<MyNeurite*>(
          cell->ExtendNewNeurite({0, 0, 1}, &my_neurite));
      dendrite_apical->AddBiologyModule(new ApicalElongation_BM());
      dendrite_apical->SetCanBranch(true);

      auto* dendrite_basal1 = bdm_static_cast<MyNeurite*>(
          cell->ExtendNewNeurite({0, 0, -1}, &my_neurite));
      dendrite_basal1->AddBiologyModule(new BasalElongation_BM());
      dendrite_basal1->SetCanBranch(true);

      auto* dendrite_basal2 = bdm_static_cast<MyNeurite*>(
          cell->ExtendNewNeurite({0, 0.6, -0.8}, &my_neurite));
      dendrite_basal2->AddBiologyModule(new BasalElongation_BM());
      dendrite_basal2->SetCanBranch(true);

      auto* dendrite_basal3 = bdm_static_cast<MyNeurite*>(
          cell->ExtendNewNeurite({0.3, -0.6, -0.8}, &my_neurite));
      dendrite_basal3->AddBiologyModule(new BasalElongation_BM());
      dendrite_basal3->SetCanBranch(true);

    } // end for num_cells
  }  // end CellCreator

} // namespace bdm

#endif
