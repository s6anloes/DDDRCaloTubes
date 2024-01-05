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
                void AddFibreCher(unsigned int fibre_id, G4int n) { FibreSignalsCher[fibre_id] += n; }
                void AddFibreScin(unsigned int fibre_id, G4int n) { FibreSignalsScin[fibre_id] += n; }

                inline std::map<unsigned int, G4int> GetFibreSignalsCher() { return FibreSignalsCher; }
                inline std::map<unsigned int, G4int> GetFibreSignalsScin() { return FibreSignalsScin; }

                inline void SetX(G4double x) { x_truth = x; }
                inline void SetY(G4double y) { y_truth = y; }
                inline void SetZ(G4double z) { z_truth = z; }

                inline void SetFirstStepNumber(G4int n) { first_step_number = n; }
                inline unsigned int GetFirstStepNumber() { return first_step_number; }

                inline void AddLeakage(G4double l) { leakage += l; }



            private:

                DRCaloTubesRunAction* fRunAction;

                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected
                std::map<unsigned int, G4int> FibreSignalsCher;
                std::map<unsigned int, G4int> FibreSignalsScin;

                unsigned int first_step_number = -1;  // to make sure first position in truth volume is saved
                unsigned int last_step_number = 0;    // to make sure only last step in volume is considered
                G4double x_truth;
                G4double y_truth;
                G4double z_truth;
                G4double leakage;
                
                //std::vector<G4int> FibreSignalsCher;
                //std::vector<G4int> FibreSignalsScin;

        };
    }
}

#endif