#include "DRCaloTubesSteppingAction.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4Tubs.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"

#include "DD4hep/Detector.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DDG4/Geant4Mapping.h"

namespace dd4hep {

    namespace sim {

        DRCaloTubesSteppingAction::DRCaloTubesSteppingAction(DRCaloTubesEventAction* eventAction)
        : fEventAction(eventAction){
        }

        DRCaloTubesSteppingAction::~DRCaloTubesSteppingAction() {}

        /// Apply Birks Law
        G4double DRCaloTubesSteppingAction::ApplyBirks( const G4double& de, const G4double& steplength ) 
        {
            return (de/steplength) / ( 1+fk_B*(de/steplength) ) * steplength;
        }

        /// Smear Scintillation Signal
        G4int DRCaloTubesSteppingAction::SmearSSignal( const G4double& satde ) 
        {
            return G4Poisson(satde*9.5);
        }

        /// Smear Cherenkov Signal
        G4int DRCaloTubesSteppingAction::SmearCSignal( )
        {
            return G4Poisson(0.153);
        }

        //Define GetDistanceToSiPM() method
        //
        G4double DRCaloTubesSteppingAction::GetDistanceToSiPM(const G4Step* step) {

            // Get the pre-step point
            const G4StepPoint* preStepPoint = step->GetPreStepPoint();
            // Get the global position of the pre-step point
            G4ThreeVector globalPos = preStepPoint->GetPosition();
            // Get the local position of the pre-step point in the current volume's coordinate system
            G4ThreeVector localPos = preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);
            // G4cout << "Local Position (X,Y,Z): (" << localPos.x()/CLHEP::mm << ", " << localPos.y()/CLHEP::mm << ", " << localPos.z()/CLHEP::mm << ") mm" << G4endl;

            // Get the logical volume of the current step
            G4LogicalVolume* currentVolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
            // Get the solid associated with the logical volume
            G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
            // Get the dimensions of the solid (size of the volume)
            G4double size = solid->GetZHalfLength();

            G4double distance_to_sipm = size - localPos.z();
            return distance_to_sipm;

        }

        //Define AttenuateHelper() method
        G4int DRCaloTubesSteppingAction::AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length) {
            double probability_of_survival = exp(-distance/attenuation_length);

            G4int survived_photons = 0;
            for (int i=0; i<signal; i++)
            {
                // Simulate drawing between 0 and 1 with probability x of getting 1
                if (G4UniformRand() <= probability_of_survival) survived_photons++;
            }

            return survived_photons;

        }

        //Define AttenuateSSignal() method
        //
        G4int DRCaloTubesSteppingAction::AttenuateSSignal(const G4int& signal, const G4double& distance) {

            return AttenuateHelper(signal, distance, fSAttenuationLength);    

        }

        //Define AttenuateCSignal() method
        //
        G4int DRCaloTubesSteppingAction::AttenuateCSignal(const G4int& signal, const G4double& distance) {

            return AttenuateHelper(signal, distance, fCAttenuationLength);    
            
        }

        void DRCaloTubesSteppingAction::UserSteppingAction(const G4Step* step) {// Get step info
            //{{{
            //std::cout<<"SDAction::process() call"<<std::endl;
            //if (step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()){
            //    std::cout<<"TrackID = " << step->GetTrack()->GetTrackID()<<std::endl;
            //    if (step->IsFirstStepInVolume() && step->IsLastStepInVolume() && step->GetTrack()->GetCurrentStepNumber()){
            //        std::cout<<"Track only one STEP" << std::endl;
            //    }
            //}
            //std::cout<<"SteppingAction:: Starting step!!!!!!!!" <<std::endl;
            //std::cout<<"Length of vectors. Cher: "<< std::to_string(fEventAction->GetFibreSignalsCher().size()) <<std::endl;
            G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
            G4double edep = step->GetTotalEnergyDeposit();
            G4double steplength = step->GetStepLength();
            std::string volume_name = volume->GetName(); 

            if (volume_name.substr(0, 7) == "leakage")
            {
                if (step->IsLastStepInVolume())
                {
                    auto name = step->GetTrack()->GetDefinition()->GetParticleName();
                    if (name=="nu_mu" || name=="nu_e" || name=="anti_nu_e" || name=="anti_nu_mu"){
                        fEventAction->AddNeutrinoLeakage(step->GetTrack()->GetKineticEnergy());
                        step->GetTrack()->SetTrackStatus(fStopAndKill);
                    }
                    else{
                        fEventAction->AddLeakage(step->GetTrack()->GetKineticEnergy());
                        step->GetTrack()->SetTrackStatus(fStopAndKill);
                    }
                } else {
                    // return;
                }
            } 
            
            //--------------------------------------------------
            //Store information from Scintillation and Cherenkov
            //signals
            //--------------------------------------------------

            G4int signalhit = 0;
            
            if ( volume_name.substr(0, 10) == "scin_fibre" ) //scintillating fiber/tube
            { 
                fEventAction->AddEdepScin(edep);
                //G4VPhysicalVolume* step_vol  = step->GetTrack()->GetVolume();
                //std::cout<<"Step Volume in Geant4: " << step_vol->GetName() <<" : " << std::to_string(step_vol->GetCopyNo())<<std::endl; 
                if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) 
                {
                    step->GetTrack()->SetTrackStatus( fStopAndKill ); 
                }

                if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle

                G4double distance_to_sipm = DRCaloTubesSteppingAction::GetDistanceToSiPM(step);
                // std::cout<<"UserSteppingAction:: Distance to SiPM: " << distance_to_sipm/CLHEP::mm << " mm" <<std::endl;
                signalhit = DRCaloTubesSteppingAction::SmearSSignal( DRCaloTubesSteppingAction::ApplyBirks( edep, steplength ) );
                signalhit = DRCaloTubesSteppingAction::AttenuateSSignal(signalhit, distance_to_sipm);
                fEventAction->AddScin(signalhit);
                auto handle = step->GetPreStepPoint()->GetTouchableHandle();
                unsigned int fibre_id = handle->GetCopyNumber(2);
                short int layer_id = handle->GetCopyNumber(4);
                unsigned short int stave_id = handle->GetCopyNumber(5);
                int tower_id = (layer_id << 16) | stave_id;
                // G4cout << "Fibre ID: " << fibre_id << " Layer ID: " << layer_id << " Stave ID: " << stave_id << " Tower ID: " << tower_id << G4endl;
                fEventAction->AddFibreSignal(tower_id, fibre_id, signalhit);
                // fEventAction->AddFibreScin(fibre_id, signalhit); 
            } // end of scintillating fibre

            else if ( volume_name.substr(0, 10) == "cher_fibre" ) //Cherenkov fiber/tube
            {
                fEventAction->AddEdepCher(edep);
                if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
                {
                    G4OpBoundaryProcessStatus theStatus = Undefined;

                    G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

                    if (OpManager) 
                    {
                        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
                        G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

                        for ( G4int i=0; i<MAXofPostStepLoops; i++) 
                        {
                            G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
                            fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
                            if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
                        }
                    }

                    switch ( theStatus ){
                                            
                        case TotalInternalReflection: 
                        {
                            //std::cout<<"SteppingAction:: Total Internal Refelction"<<std::endl;
                            G4double distance_to_sipm = DRCaloTubesSteppingAction::GetDistanceToSiPM(step);
                            G4int c_signal = DRCaloTubesSteppingAction::SmearCSignal( );
                            signalhit = DRCaloTubesSteppingAction::AttenuateCSignal(c_signal, distance_to_sipm);
                            fEventAction->AddCher(signalhit);
                            auto handle = step->GetPreStepPoint()->GetTouchableHandle();
                            unsigned int fibre_id = handle->GetCopyNumber(2);
                            short int layer_id = handle->GetCopyNumber(4);
                            unsigned short int stave_id = handle->GetCopyNumber(5);
                            int tower_id = (layer_id << 16) | stave_id;
                            // G4cout << "Fibre ID: " << fibre_id << " Layer ID: " << layer_id << " Stave ID: " << stave_id << " Tower ID: " << tower_id << G4endl;
                            fEventAction->AddFibreSignal(tower_id, fibre_id, signalhit);
                            // unsigned int fibre_id = step->GetTrack()->GetVolume()->GetCopyNo();
                            // fEventAction->AddFibreCher(fibre_id, c_signal); 
                            step->GetTrack()->SetTrackStatus( fStopAndKill );
                            break;
                        }
                        default:
                            step->GetTrack()->SetTrackStatus( fStopAndKill );
                    } //end of swich cases

                } //end of optical photon

            } //end of Cherenkov fiber

            //std::cout<<"SteppingAction:: Finishing step!!!!!!!!" <<std::endl;

            return;
        }

    } // namespace sim

} // namespace dd4hep
