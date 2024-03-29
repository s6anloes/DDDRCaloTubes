#include "DRCaloTubesSteppingAction.h"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"

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
            const G4double k_B = 0.126; //Birks constant
            return (de/steplength) / ( 1+k_B*(de/steplength) ) * steplength;
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
            
            //--------------------------------------------------
            //Store information from Scintillation and Cherenkov
            //signals
            //--------------------------------------------------

            std::string volume_name;
            std::string scin_fibre_name = "scin_fibre_0";
            std::string scin_clad_name  = "scin_clad_0";
            std::string cher_fibre_name = "cher_fibre_0";
            std::string cher_clad_name  = "cher_clad_0";
            volume_name = volume->GetName(); 
            G4int signalhit = 0;
            
            // truth information from truth box
            if ( volume_name=="truth_volume_0" ) 
            {  
                unsigned int step_number = step->GetTrack()->GetCurrentStepNumber();
                if (step->IsFirstStepInVolume() && step->GetTrack()->GetParentID()==0 && step_number<fEventAction->GetFirstStepNumber())
                {
                    //std::cout<<"I WANT THE TRUTH"<<std::endl;
                    fEventAction->SetFirstStepNumber(step_number);
                    G4ThreeVector position  = step->GetPreStepPoint()->GetPosition();
                    G4ThreeVector direction = step->GetPreStepPoint()->GetMomentumDirection().unit();

                    G4double truth_margin = 100.0*CLHEP::um;
                    G4double z = direction.getZ(); 
                    G4double count = (z!=0) ? truth_margin/z : 0;
                    position += G4ThreeVector(count*direction.getX(), count*direction.getY(), truth_margin);

                    fEventAction->SetX(position.getX()/CLHEP::mm);
                    fEventAction->SetY(position.getY()/CLHEP::mm);
                    fEventAction->SetZ(position.getZ()/CLHEP::mm);

                } // end of first step in volume

                else if (step->IsLastStepInVolume() && step->GetTrack()->GetNextVolume()->GetName()=="world_volume_1")
                {
                    //std::cout<<"YOU CANT HANDLE THE TRUTH"<<std::endl;
                    G4double leakage = step->GetPreStepPoint()->GetTotalEnergy();
                    fEventAction->AddLeakage(leakage/CLHEP::MeV);

                } // end of last step in volume

            } // end of truth box            

            
            else if ( volume_name.substr(0, 4) == "scin" ) //scintillating fiber/tube
            { 
                //G4VPhysicalVolume* step_vol  = step->GetTrack()->GetVolume();
                //std::cout<<"Step Volume in Geant4: " << step_vol->GetName() <<" : " << std::to_string(step_vol->GetCopyNo())<<std::endl; 
                if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) 
                {
                    step->GetTrack()->SetTrackStatus( fStopAndKill ); 
                }

                if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle

                signalhit = DRCaloTubesSteppingAction::SmearSSignal( DRCaloTubesSteppingAction::ApplyBirks( edep, steplength ) );
                fEventAction->AddScin(signalhit);
                unsigned int fibre_id = step->GetTrack()->GetVolume()->GetCopyNo();
                fEventAction->AddFibreScin(fibre_id, signalhit); 
            } // end of scintillating fibre

            else if ( volume_name.substr(0, 4) == "cher" ) //Cherenkov fiber/tube
            {
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
                            G4int c_signal = DRCaloTubesSteppingAction::SmearCSignal( );
                            fEventAction->AddCher(c_signal);
                            unsigned int fibre_id = step->GetTrack()->GetVolume()->GetCopyNo();
                            fEventAction->AddFibreCher(fibre_id, c_signal); 
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
