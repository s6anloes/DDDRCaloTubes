#include "DRCaloTubesEventAction.h"

//#include "G4RunManager.hh"

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(DRCaloTubesRunAction* runAction): G4UserEventAction(),
            fRunAction(runAction), NofCherDet(0), NofScinDet(0){}

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            NofScinDet = 0;
            NofCherDet = 0;
            FibreSignalsCher.resize(2880);
            FibreSignalsScin.resize(2880);
            FibreSignalsCher.clear();
            FibreSignalsScin.clear();
            fRunAction->Reset();
        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            fRunAction->Fill(NofCherDet, NofScinDet, FibreSignalsCher, FibreSignalsScin);
        }
    } // namespace sim
} // namespace drc
