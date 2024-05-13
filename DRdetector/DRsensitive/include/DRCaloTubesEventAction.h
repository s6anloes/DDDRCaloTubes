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

                // inline std::map<unsigned int, G4int> GetFibreSignalsCher() { return FibreSignalsCher; }
                // inline std::map<unsigned int, G4int> GetFibreSignalsScin() { return FibreSignalsScin; }
                std::vector<double> tower_ids;
                std::vector<double> fibre_ids;
                std::vector<double> fibre_signals;


            private:
                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected
                // std::map<unsigned int, G4int> FibreSignalsCher;
                // std::map<unsigned int, G4int> FibreSignalsScin;
                std::unordered_map<int, std::unordered_map<unsigned int, G4int>> fibre_signals_map;
                //std::vector<G4int> FibreSignalsCher;
                //std::vector<G4int> FibreSignalsScin;

        };
    }
}

#endif