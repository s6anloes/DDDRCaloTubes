#ifndef DRCaloTubesRunAction_h
#define DRCaloTubesRunAction_h 1

#include "DRCaloTubesEventAction.h"

//Includers from Geant4
//
#include "G4UserRunAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"

class G4Run;


namespace dd4hep {
    namespace sim {
        class DRCaloTubesRunAction : public G4UserRunAction {
            
            public:
                //Constructor
                //
                DRCaloTubesRunAction(DRCaloTubesEventAction* eventAction);
                //De-constructor
                //
                virtual ~DRCaloTubesRunAction();

                //Methods
                //
                virtual void BeginOfRunAction(const G4Run*);
                virtual void EndOfRunAction(const G4Run*);

            private:
                DRCaloTubesEventAction* fEventAction;

        };
    }
}

#endif

//**************************************************
