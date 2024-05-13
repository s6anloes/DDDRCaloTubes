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
            fibre_signals_map.clear();
            tower_ids.clear();
            fibre_ids.clear();
            fibre_signals.clear();

        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            // Iterate through the fibre_signals_map 
            for (const auto& tower_entry : fibre_signals_map) {
                int tower_id = tower_entry.first;
                const auto& fibre_map = tower_entry.second;

                for (const auto& fibre_entry : fibre_map) {
                    unsigned int fibre_id = fibre_entry.first;
                    G4int signal_count = fibre_entry.second;

                    tower_ids.push_back(tower_id);
                    fibre_ids.push_back(fibre_id);
                    fibre_signals.push_back(signal_count);
                }
            }


            G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, NofCherDet);
            analysisManager->FillNtupleDColumn(1, NofScinDet);
            for (size_t i = 0; i < tower_ids.size(); ++i) {
                analysisManager->FillNtupleIColumn(2, tower_ids.at(i)); // Fill each tower ID separately
                analysisManager->FillNtupleIColumn(3, fibre_ids.at(i));
                analysisManager->FillNtupleIColumn(4, fibre_signals.at(i));
            }
            
            analysisManager->AddNtupleRow();
        }

        void DRCaloTubesEventAction::AddFibreSignal(int tower_id, unsigned int fibre_id, G4int n) {
            // Check if the tower_id exists in the map
            auto tower_it = fibre_signals_map.find(tower_id);
            
            // If tower_id doesn't exist, create a new map for it
            if (tower_it == fibre_signals_map.end()) {
                fibre_signals_map.insert({tower_id, std::unordered_map<unsigned int, G4int>()});
                tower_it = fibre_signals_map.find(tower_id);
            }
            
            auto& fibre_map = tower_it->second;
            auto fibre_it = fibre_map.find(fibre_id);
            
            // If fibre_id doesn't exist, create a new entry for it
            if (fibre_it == fibre_map.end()) {
                fibre_map.insert({fibre_id, 0});
                fibre_it = fibre_map.find(fibre_id);
            }
            
            fibre_it->second += n;
        }
    } // namespace sim
} // namespace drc
