#ifndef DRCaloTubesEventAction_h
#define DRCaloTubesEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include <map>


namespace dd4hep {
    namespace sim{
        class DRCaloTubesEventAction : public G4UserEventAction {
            public:
                DRCaloTubesEventAction();
                virtual ~DRCaloTubesEventAction();

                virtual void BeginOfEventAction(const G4Event*) final;
                virtual void EndOfEventAction(const G4Event*) final;
                inline void AddCher(G4int n) { NofCherDet += n; }
                inline void AddScin(G4int n) { NofScinDet += n; }
                // void AddFibreCher(unsigned int fibre_id, G4int n) { FibreSignalsCher[fibre_id] += n; }
                // void AddFibreScin(unsigned int fibre_id, G4int n) { FibreSignalsScin[fibre_id] += n; }
                void AddFibreSignal(int tower_id, unsigned int fibre_id, G4int n);

                inline std::vector<double>& GetTowerIDs() { return tower_ids; }
                inline std::vector<double>& GetFibreIDs() { return fibre_ids; }
                inline std::vector<double>& GetFibreSignals() { return fibre_signals; }

                inline void SavePrimaryPDGID(G4int pdgid){ PrimaryPDGID = pdgid; }
                inline void SavePrimaryXYZ(G4double x, G4double y, G4double z){
                    PrimaryX = x;
                    PrimaryY = y;
                    PrimaryZ = z;
                }

                inline void SavePrimaryEnergy(G4double primaryparticleenergy){ PrimaryParticleEnergy = primaryparticleenergy; }
                void SavePrimaryMomentumDirection(G4double x, G4double y, G4double z){
                    PrimaryTheta = std::acos(z/std::sqrt(x*x + y*y + z*z));
                    PrimaryPhi = std::atan2(y, x);
                }

                std::vector<double> tower_ids;
                std::vector<double> fibre_ids;
                std::vector<double> fibre_signals;


            private:
                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected


                G4int    PrimaryPDGID; //PDGID of primary particle
                G4double PrimaryParticleEnergy; //Primary particle energy
                G4double PrimaryX; //Primary particle x position
                G4double PrimaryY; //Primary particle y position
                G4double PrimaryZ; //Primary particle z position
                G4double PrimaryTheta; //Primary particle theta
                G4double PrimaryPhi; //Primary particle phi
                
                std::unordered_map<int, std::unordered_map<unsigned int, G4int>> fibre_signals_map;
                

        };
    }
}

#endif