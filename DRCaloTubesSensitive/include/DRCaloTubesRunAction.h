#ifndef DRCaloTubesRunAction_h
#define DRCaloTubesRunAction_h 1

//#include "DRCaloTubesEventAction.h"

//Includers from Geant4
//
#include "G4UserRunAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"

#include <unordered_map>

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
                void Reset(); 
                void Fill(const G4int cher, const G4int scin, std::map<int, G4int> fibrecher, std::map<int, G4int> fibrescin);

            private:
                //DRCaloTubesEventAction* fEventAction;
                TFile* fFile;
                TTree* fTree;



                G4int     NofCherDet; //Number of Cherenkov p.e. detected 
                G4int     NofScinDet; //Number of Scintillating p.e. detected

                int nfibres = 48*60;
                std::vector<G4int> FibreSignalsCher;//{std::vector<G4int>(nfibres,0)};  // Cherenkov signal in each fibre
                std::vector<G4int> FibreSignalsScin;//{std::vector<G4int>(nfibres,0)};  // Scinitillation signal in each fibre



        };
    }
}

#endif

//**************************************************
