//Prevent including headers multiple times
//
#ifndef DRCaloTubesGeneratorAction_h
#define DRCaloTubesGeneratorAction_h 1

#include "DRCaloTubesEventAction.h"

//Includers from Geant4
//
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

//Forwards declarations from Geant4
//
class G4GeneralParticleSource;
//class G4ParticleGun; //in case user want to switch to G4ParticleGun
class G4Event;

//Forward declarations from project
//
// class DRCaloTubesEventAction;

namespace dd4hep {
    namespace sim {
        class DRCaloTubesGeneratorAction : public G4VUserPrimaryGeneratorAction {
        
            public:
                //Constructor()
                //
                DRCaloTubesGeneratorAction(DRCaloTubesEventAction* evtAction);    

                //De-constructor()
                //
                virtual ~DRCaloTubesGeneratorAction();

                virtual void GeneratePrimaries(G4Event* event);
        
            private:
                G4GeneralParticleSource* fGeneralParticleSource;
                //G4ParticleGun*  fParticleGun; // G4ParticleGun
                DRCaloTubesEventAction* fEventAction;

        };
    }
}

#endif

//**************************************************
