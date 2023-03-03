#ifndef DRCaloTubesEventAction_h
#define DRCaloTubesEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"


namespace dd4hep {
    namespace sim{
        class DRCaloTubesEventAction : public G4UserEventAction {
            public:
                DRCaloTubesEventAction();
                virtual ~DRCaloTubesEventAction();

                virtual void BeginOfEventAction(const G4Event*) final;
                virtual void EndOfEventAction(const G4Event*) final;
                inline void AddCher(G4int n) { NofCherDetEA += n; }
                inline void AddScin(G4int n) { NofScinDetEA += n; }


            private:


                G4int     NofCherDetEA; //Number of Cherenkov p.e. detected 
                G4int     NofScinDetEA; //Number of Scintillating p.e. detected

        };
    }
}

#endif