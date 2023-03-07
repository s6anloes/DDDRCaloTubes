#ifndef DRCaloTubesRunAction_h
#define DRCaloTubesRunAction_h 1

//#include "DRCaloTubesEventAction.h"

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
                DRCaloTubesRunAction();
                //De-constructor
                //
                virtual ~DRCaloTubesRunAction();

                //Methods
                //
                virtual void BeginOfRunAction(const G4Run*);
                virtual void EndOfRunAction(const G4Run*);
                inline void Reset() { NofCherDet=0; NofScinDet=0; }
                void Fill(const G4int cher, const G4int scin);

            private:
                //DRCaloTubesEventAction* fEventAction;
                TFile* fFile;
                TTree* fTree;



                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected




        };
    }
}

#endif

//**************************************************
