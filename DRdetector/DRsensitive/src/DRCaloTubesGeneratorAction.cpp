//Includers from project files
//
#include "DRCaloTubesGeneratorAction.h"
#include "DRCaloTubesEventAction.h"

//Includers from Geant4
//
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"

namespace dd4hep {
    namespace sim {

        //Constructor
        //
        DRCaloTubesGeneratorAction::DRCaloTubesGeneratorAction(DRCaloTubesEventAction* evtAction)
        : G4VUserPrimaryGeneratorAction(),
        fGeneralParticleSource( nullptr ),
        fEventAction(evtAction)
        /*fParticleGun( nullptr )*/ {

            //G4int nofParticles = 1;                           //for particle gun
            //fParticleGun = new G4ParticleGun(nofParticles);   //for particle gun
            fGeneralParticleSource = new G4GeneralParticleSource();

            //default G4GeneralParticleSource parameters (can be changed via UI)
            //
            G4ParticleDefinition* particleDefinition =
                G4ParticleTable::GetParticleTable()->FindParticle("e-");

            fGeneralParticleSource->SetParticleDefinition(particleDefinition);
            fGeneralParticleSource->SetParticlePosition( G4ThreeVector( 0.,0.,0. ) );

        }

        //De-constructor
        //
        DRCaloTubesGeneratorAction::~DRCaloTubesGeneratorAction() {
        
            delete fGeneralParticleSource;
            //delete fParticleGun;

        }

        //GeneratePrimaries() method
        //
        void DRCaloTubesGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
            
            fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
            //fParticleGun->GeneratePrimaryVertex(anEvent);

            //Save primary particle energy, PDGID and x-y position
            //
            // fEventAction->SavePrimaryEnergy(fGeneralParticleSource->GetParticleEnergy());
            // fEventAction->SavePrimaryPDGID(fGeneralParticleSource->GetParticleDefinition()->GetPDGEncoding());
            // fEventAction->SavePrimaryXYZ(fGeneralParticleSource->GetParticlePosition().x(),
            //                             fGeneralParticleSource->GetParticlePosition().y(),
            //                             fGeneralParticleSource->GetParticlePosition().z());
            // fEventAction->SavePrimaryMomentumDirection(fGeneralParticleSource->GetParticleMomentumDirection().x(),
            //                                         fGeneralParticleSource->GetParticleMomentumDirection().y(),
            //                                         fGeneralParticleSource->GetParticleMomentumDirection().z());

        }

    }
}
