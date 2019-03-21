#ifndef UTIL_NETHODS_
#define UTIL_NETHODS

namespace bdm {

  // define morpho_exporteur for neuron morphology export in swc format
  template <typename TSimulation = Simulation<>>
  inline void morpho_exporteur() {
    auto* sim = TSimulation::GetActive();
    auto* rm = sim->GetResourceManager();
    int seed = sim->GetRandom()->GetSeed();

    int cell_nb = 0;

    rm->ApplyOnAllElements([&](auto&& so, SoHandle) {
      if (so->template IsSoType<MyCell>()) {
        auto&& cell = so.template ReinterpretCast<MyCell>();
        auto cellPosition = cell.GetPosition();
        ofstream swcFile;
        string swcFileName = Concat("./cell", cell_nb,"_seed", seed, ".swc").c_str();
        swcFile.open(swcFileName);
        cell->SetLabel(1);
        cell_nb++;
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

} // namespace bdm

#endif
