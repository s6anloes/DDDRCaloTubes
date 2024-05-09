#include "DRCaloTubesEventAction.h"
#include "G4SystemOfUnits.hh"

//#include "G4RunManager.hh"
#include "G4RootAnalysisManager.hh"

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(): G4UserEventAction(),
            NofCherDet(0), NofScinDet(0){
                //FibreSignalsCher.resize(nfibres);
                //FibreSignalsScin.resize(nfibres);
                //for (int i=0; i<nfibres; i++){
                //    FibreSignalsCher.push_back(0);
                //    FibreSignalsScin.push_back(0);
                //}
            }

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            NofScinDet = 0;
            NofCherDet = 0;
            // FibreSignalsCher.clear();
            // FibreSignalsScin.clear();

        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, NofCherDet);
            analysisManager->FillNtupleDColumn(1, NofScinDet);
            analysisManager->AddNtupleRow();
        }
    } // namespace sim
} // namespace drc
