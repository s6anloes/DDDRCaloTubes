#ifndef DRCaloTubesEventAction_h
#define DRCaloTubesEventAction_h 1

#include "DRCaloTubesRunAction.h"

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include <map>


namespace dd4hep {
    namespace sim{
        class DRCaloTubesEventAction : public G4UserEventAction {
            public:
                DRCaloTubesEventAction(DRCaloTubesRunAction* runAction);
                virtual ~DRCaloTubesEventAction();

                virtual void BeginOfEventAction(const G4Event*) final;
                virtual void EndOfEventAction(const G4Event*) final;
                inline void AddCher(G4int n) { NofCherDet += n; }
                inline void AddScin(G4int n) { NofScinDet += n; }
                // void AddFibreCher(unsigned int fibre_id, G4int n) { FibreSignalsCher[fibre_id] += n; }
                // void AddFibreScin(unsigned int fibre_id, G4int n) { FibreSignalsScin[fibre_id] += n; }

                // inline std::map<unsigned int, G4int> GetFibreSignalsCher() { return FibreSignalsCher; }
                // inline std::map<unsigned int, G4int> GetFibreSignalsScin() { return FibreSignalsScin; }


            private:

                DRCaloTubesRunAction* fRunAction;

                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected
                // std::map<unsigned int, G4int> FibreSignalsCher;
                // std::map<unsigned int, G4int> FibreSignalsScin;
                
                //std::vector<G4int> FibreSignalsCher;
                //std::vector<G4int> FibreSignalsScin;

        };
    }
}

#endif