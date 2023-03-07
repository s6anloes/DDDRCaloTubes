#include "DRCaloTubesEventAction.h"

//#include "G4RunManager.hh"

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(DRCaloTubesRunAction* runAction): G4UserEventAction(),
            fRunAction(runAction) {}

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            NofScinDet = 0;
            NofCherDet = 0;
            fRunAction->Reset();
        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            fRunAction->Fill(NofCherDet, NofScinDet);
        }
    } // namespace sim
} // namespace drc
