#ifndef DRCaloTubesEventAction_h
#define DRCaloTubesEventAction_h 1

#include "DRCaloTubesRunAction.h"

#include "G4UserEventAction.hh"
#include "G4Event.hh"


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
                inline void AddFibreCher(int fibre_id, G4int n) { FibreSignalsCher.at(fibre_id) += n; }
                inline void AddFibreScin(int fibre_id, G4int n) { FibreSignalsScin.at(fibre_id) += n; }

                inline std::vector<G4int> GetFibreSignalsCher() { return FibreSignalsCher; }
                inline std::vector<G4int> GetFibreSignalsScin() { return FibreSignalsScin; }



            private:

                DRCaloTubesRunAction* fRunAction;

                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected
                int nfibres = 48*60;
                std::vector<G4int> FibreSignalsCher;
                std::vector<G4int> FibreSignalsScin;

        };
    }
}

#endif