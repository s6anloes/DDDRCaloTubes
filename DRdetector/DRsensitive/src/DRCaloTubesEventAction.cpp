#include "DRCaloTubesEventAction.h"
#include "G4SystemOfUnits.hh"

//#include "G4RunManager.hh"
#include "G4RootAnalysisManager.hh"
#include "G4ParticleDefinition.hh"

#include <limits>

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(): G4UserEventAction(),
            EnergyScin(0), EnergyCher(0), NofCherDet(0), NofScinDet(0), Leakage(0), NeutrinoLeakage(0){
                //FibreSignalsCher.resize(nfibres);
                //FibreSignalsScin.resize(nfibres);
                //for (int i=0; i<nfibres; i++){
                //    FibreSignalsCher.push_back(0);
                //    FibreSignalsScin.push_back(0);
                //}
            }

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            EnergyScin = 0;
            EnergyCher = 0;
            NofScinDet = 0;
            NofCherDet = 0;
            // FibreSignalsCher.clear();
            // FibreSignalsScin.clear();
            fibre_signals_map.clear();
            tower_ids.clear();
            fibre_ids.clear();
            fibre_signals.clear();

            PrimaryParticleEnergy = 0;
            PrimaryPDGID = 0;
            PrimaryX = std::numeric_limits<double>::quiet_NaN();
            PrimaryY = std::numeric_limits<double>::quiet_NaN();
            PrimaryZ = std::numeric_limits<double>::quiet_NaN();
            PrimaryPX = std::numeric_limits<double>::quiet_NaN();
            PrimaryPY = std::numeric_limits<double>::quiet_NaN();
            PrimaryPZ = std::numeric_limits<double>::quiet_NaN();

            Leakage = 0;
            NeutrinoLeakage = 0;


        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event* evt) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            // Iterate through the fibre_signals_map 
            for (const auto& tower_entry : fibre_signals_map) {
                int tower_id = tower_entry.first;
                const auto& fibre_map = tower_entry.second;

                for (const auto& fibre_entry : fibre_map) {
                    int fibre_id = fibre_entry.first;
                    G4int signal_count = fibre_entry.second;

                    tower_ids.push_back(tower_id);
                    fibre_ids.push_back(fibre_id);
                    fibre_signals.push_back(signal_count);
                }
            }

            G4PrimaryParticle* priPar = evt->GetPrimaryVertex()->GetPrimary();
            PrimaryPDGID = priPar->GetParticleDefinition()->GetPDGEncoding();
            PrimaryParticleEnergy = priPar->GetKineticEnergy();
            PrimaryX = evt->GetPrimaryVertex()->GetX0();
            PrimaryY = evt->GetPrimaryVertex()->GetY0();
            PrimaryZ = evt->GetPrimaryVertex()->GetZ0();
            PrimaryPX = priPar->GetMomentumDirection().x();
            PrimaryPY = priPar->GetMomentumDirection().y();
            PrimaryPZ = priPar->GetMomentumDirection().z();


            G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, EnergyScin);
            analysisManager->FillNtupleDColumn(1, EnergyCher);
            analysisManager->FillNtupleDColumn(2, NofCherDet);
            analysisManager->FillNtupleDColumn(3, NofScinDet);
            analysisManager->FillNtupleDColumn(7, PrimaryParticleEnergy);
            analysisManager->FillNtupleDColumn(8, PrimaryPDGID);
            analysisManager->FillNtupleDColumn(9, PrimaryX);
            analysisManager->FillNtupleDColumn(10, PrimaryY);
            analysisManager->FillNtupleDColumn(11, PrimaryZ);
            analysisManager->FillNtupleDColumn(12, PrimaryPX);
            analysisManager->FillNtupleDColumn(13, PrimaryPY);
            analysisManager->FillNtupleDColumn(14, PrimaryPZ);
            analysisManager->FillNtupleDColumn(15, Leakage);
            analysisManager->FillNtupleDColumn(16, NeutrinoLeakage);
            
            analysisManager->AddNtupleRow();
        }

        void DRCaloTubesEventAction::AddFibreSignal(int tower_id,  int fibre_id, G4int n) {
            // Check if the tower_id exists in the map
            auto tower_it = fibre_signals_map.find(tower_id);
            
            // If tower_id doesn't exist, create a new map for it
            if (tower_it == fibre_signals_map.end()) {
                fibre_signals_map.insert({tower_id, std::unordered_map<int, G4int>()});
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
