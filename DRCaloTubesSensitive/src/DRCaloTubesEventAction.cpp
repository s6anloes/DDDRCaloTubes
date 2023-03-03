#include "DRCaloTubesEventAction.h"

//#include "G4RunManager.hh"

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(): G4UserEventAction() {}

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            NofScinDetEA = 0;
            NofCherDetEA = 0;
        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDetEA = "<<NofScinDetEA<<std::endl;
            std::cout<<"NofCherDetEA = "<<NofCherDetEA<<std::endl;
        }
    } // namespace sim
} // namespace drc
