#include "DD4hep/Version.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4Poisson.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"


#include "DRCaloTubesHit.h"


//Forward declarations from Geant4
//
class G4OpBoundaryProcess;
enum G4OpBoundaryProcessStatus;
class G4ProcessVector;
enum G4ProcessVectorTypeIndex;


#if DD4HEP_VERSION_GE(1, 21)
#define GEANT4_CONST_STEP const
#else
#define GEANT4_CONST_STEP
#endif



class DRCaloTubesSD {
    public:
        typedef DRCaloTubesHit Hit;
        int placeholer;
};


int NofCherDet, NofScinDet = 0;

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  
    /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
    namespace sim {


        /// Define collections created by this sensitivie action object
        template <> void Geant4SensitiveAction<DRCaloTubesSD>::defineCollections() {
            m_collectionID = declareReadoutFilteredCollection<DRCaloTubesSD::Hit>();
        }

        /// Apply Birks Law
        G4double ApplyBirks( const G4double& de, const G4double& steplength ) 
        {
            const G4double k_B = 0.126; //Birks constant
            return (de/steplength) / ( 1+k_B*(de/steplength) ) * steplength;
        }

        /// Smear Scintillation Signal
        G4int SmearSSignal( const G4double& satde ) 
        {
            return G4Poisson(satde*9.5);
        }

        /// Smear Cherenkov Signal
        G4int SmearCSignal( )
        {
            return G4Poisson(0.153);
        }

        /// Method for generating hit(s) using the information of G4Step object.
        template <> bool Geant4SensitiveAction<DRCaloTubesSD>::process(const G4Step* step, G4TouchableHistory* ) {

            // Get step info
            //{{{
            //std::cout<<"SDAction::process() call"<<std::endl;
            //if (step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()){
            //    std::cout<<"TrackID = " << step->GetTrack()->GetTrackID()<<std::endl;
            //    if (step->IsFirstStepInVolume() && step->IsLastStepInVolume() && step->GetTrack()->GetCurrentStepNumber()){
            //        std::cout<<"Track only one STEP" << std::endl;
            //    }
            //}
            G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
            G4double edep = step->GetTotalEnergyDeposit();
            G4double steplength = step->GetStepLength();
            
            //--------------------------------------------------
            //Store information from Scintillation and Cherenkov
            //signals
            //--------------------------------------------------
        
            std::string fibre_name;
            std::string scin_fibre_name = "scin_fibre_0";
            std::string scin_clad_name  = "scin_clad_0";
            std::string cher_fibre_name = "cher_fibre_0";
            std::string cher_clad_name  = "cher_clad_0";
            fibre_name = volume->GetName(); 
            //std::cout<<"SDAction::process: fibre_name = " <<fibre_name<<std::endl;
            G4int signalhit = 0;

            
            if ( fibre_name==scin_fibre_name || fibre_name==scin_clad_name ) //scintillating fiber/tube
            { 
                if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) 
                {
                    //std::cout<<"SCINTILLATION:: PostStep Position = "<<step->GetPostStepPoint()->GetPosition()<<std::endl;
                    step->GetTrack()->SetTrackStatus( fStopAndKill ); 
                }

                if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return true; } //not ionizing particle


                
                            
                signalhit = SmearSSignal( ApplyBirks( edep, steplength ) );
                NofScinDet += signalhit; 
                //std::cout<<"SteppingAction:: NofScinDet = "<<NofScinDet<<std::endl;
            }

            if ( fibre_name==cher_fibre_name || fibre_name==cher_clad_name ) //Cherenkov fiber/tube
            {
                //std::cout<<"SteppingAction:: Cherenkov Fibre"<<std::endl;
                //std::cout<<"SteppingAction:: Particle name = " <<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<std::endl;
                //std::cout<<"Corresponding TrackID = "<<step->GetTrack()->GetTrackID()<<std::endl;
                //std::cout<<"CHERENKOV:: PreStep Position = "<<step->GetPreStepPoint()->GetPosition()<<std::endl;    
                //std::cout<<"SteppingAction Optical Photon name = " <<G4OpticalPhoton::Definition()->GetParticleName()<<std::endl;
                if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
                {
                    //std::cout<<"SteppingAction:: Optical Photon"<<std::endl;
                    //std::cout<<"CHERENKOV:: PreStep Position = "<<step->GetPreStepPoint()->GetPosition()<<std::endl;
                    //std::cout<<"OpticalPhoton IsFirstStepInVolume = "<<step->IsFirstStepInVolume()<<std::endl;
                    //std::cout<<"OpticalPhoton IsLastStepInVolume  = "<<step->IsLastStepInVolume()<<std::endl;
                    //std::cout<<"Corresponding TrackID = "<<step->GetTrack()->GetTrackID()<<std::endl;
                    //std::cout<<"SDAction::process: fibre_name = " <<fibre_name<<std::endl;
                    //std::cout<<"Next Volume = " << step->GetTrack()->GetNextVolume()->GetName() <<std::endl;
                    //std::cout<<"Stepnumber = "<<step->GetTrack()->GetCurrentStepNumber()<<std::endl;   
                    G4OpBoundaryProcessStatus theStatus = Undefined;

                    G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

                    if (OpManager) 
                    {
                        //std::cout<<"SteppingAction:: OpManager"<<std::endl;
                        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
                        G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

                        for ( G4int i=0; i<MAXofPostStepLoops; i++) 
                        {
                            G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
                            G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
                            if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
                        }
                    }

                    //std::cout<<"SteppingAction:: theStatus = " <<theStatus<<std::endl;
                    

                    //G4int c_signal = SmearCSignal( );
                    //NofCherDet += c_signal;

                    switch ( theStatus ){
                                            
                        case TotalInternalReflection: 
                        {
                            //std::cout<<"SteppingAction:: Total Internal Refelction"<<std::endl;
                            G4int c_signal = SmearCSignal( );
                            NofCherDet += c_signal;
                            step->GetTrack()->SetTrackStatus( fStopAndKill );
                            //break;
                        }
                        default:
                            step->GetTrack()->SetTrackStatus( fStopAndKill );
                    } //end of swich cases

                } //end of optical photon

            } //end of Cherenkov fiber

            return true;

        }

    }
}


//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DRCaloTubesSD> DRCaloTubesSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DRCaloTubesSDAction)





/*
//Forward declarations from Geant4
//
class G4OpBoundaryProcess;
enum G4OpBoundaryProcessStatus;
class G4ProcessVector;
enum G4ProcessVectorTypeIndex; 

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  
    /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
    namespace sim   {

        struct DRCaloTubesSDData {

            Geant4Sensitive*  sensitive{};

            G4int     NofCherDet; //Number of Cherenkov p.e. detected 
            G4int     NofScinDet; //Number of Scintillating p.e. detected

            G4OpBoundaryProcess* fOpProcess;

            DRCaloTubesSDData() :
            NofCherDet(0),
            NofScinDet(0)
            {
            }


                /// Clear collected information and restart for new hit
            void clear()  {
                // nothing to clear
            }

            /// Apply Birks Law
            G4double ApplyBirks( const G4double& de, const G4double& steplength ) 
            {
                const G4double k_B = 0.126; //Birks constant
                return (de/steplength) / ( 1+k_B*(de/steplength) ) * steplength;
            }

            /// Smear Scintillation Signal
            G4int SmearSSignal( const G4double& satde ) 
            {
                return G4Poisson(satde*9.5);
            }

            /// Smear Cherenkov Signal
            G4int SmearCSignal( )
            {
                return G4Poisson(0.153);
            }


            /// Method for generating hit(s) using the information of G4Step object.
            G4bool process(G4Step GEANT4_CONST_STEP * step, G4TouchableHistory* ) 
            {
                // Get step info
                //
                //std::cout<<"SDAction::process() call"<<std::endl;
                G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
                G4double edep = step->GetTotalEnergyDeposit();
                G4double steplength = step->GetStepLength();
                
                //--------------------------------------------------
                //Store information from Scintillation and Cherenkov
                //signals
                //--------------------------------------------------
            
                std::string fibre_name;
                std::string scin_name = "scin_fibre_0";
                std::string cher_name = "cher_fibre_0";
                fibre_name = volume->GetName(); 
                //std::cout<<"SDAction::process: fibre_name = " <<fibre_name<<std::endl;
                G4int signalhit = 0;

                if ( fibre_name==scin_name ) //scintillating fiber/tube
                { 
                    if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) 
                    {
                        //std::cout<<"SCINTILLATION:: PostStep Position = "<<step->GetPostStepPoint()->GetPosition()<<std::endl;
                        step->GetTrack()->SetTrackStatus( fStopAndKill ); 
                    }

                    if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return true; } //not ionizing particle


                    
                                
                    signalhit = SmearSSignal( ApplyBirks( edep, steplength ) );
                    NofScinDet += signalhit; 
                }

                if ( fibre_name==cher_name ) //Cherenkov fiber/tube
                {
                    //std::cout<<"SteppingAction:: Cherenkov Fibre"<<std::endl;
                    //std::cout<<"SteppingAction:: Particle name = " <<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<std::endl;
                    //std::cout<<"SteppingAction Optical Photon name = " <<G4OpticalPhoton::Definition()->GetParticleName()<<std::endl;
                    if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
                    {
                        //std::cout<<"SteppingAction:: Optical Photon"<<std::endl;
                        std::cout<<"CHERENKOV:: PreStep Position = "<<step->GetPreStepPoint()->GetPosition()<<std::endl;
                        std::cout<<"OpticalPhoton IsLastStepInVolume = "<<step->IsLastStepInVolume()<<std::endl;
                        std::cout<<"Corresponding TrackID = "<<step->GetTrack()->GetTrackID()<<std::endl;    
                        G4OpBoundaryProcessStatus theStatus = Undefined;

                        G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

                        if (OpManager) 
                        {
                            //std::cout<<"SteppingAction:: OpManager"<<std::endl;
                            G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
                            G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

                            for ( G4int i=0; i<MAXofPostStepLoops; i++) 
                            {
                                G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
                                fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
                                if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
                            }
                        }

                        //std::cout<<"SteppingAction:: theStatus = " <<theStatus<<std::endl;
                        

                        //G4int c_signal = SmearCSignal( );
                        //NofCherDet += c_signal;

                        switch ( theStatus ){
                                                
                            case TotalInternalReflection: 
                            {
                                std::cout<<"SteppingAction:: Total Internal Refelction"<<std::endl;
                                G4int c_signal = SmearCSignal( );
                                NofCherDet += c_signal;
                                step->GetTrack()->SetTrackStatus( fStopAndKill );
                            }
                            default:
                                step->GetTrack()->SetTrackStatus( fStopAndKill );
                        } //end of swich cases

                    } //end of optical photon

                } //end of Cherenkov fiber


                return true;
            } // end of process()


            /// Pre-event action callback
            void beginEvent(const G4Event* event)
            {
                NofScinDet = 0;
                NofCherDet = 0;
            }

            /// Post-event action callback
            void endEvent(const G4Event* event)   {
                std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
                std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
                std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;
            }
    

        
        };


        /// Initialization overload for specialization
        template <> void Geant4SensitiveAction<DRCaloTubesSDData>::initialize() {
        eventAction().callAtBegin(&m_userData, &DRCaloTubesSDData::beginEvent);
        eventAction().callAtEnd(&m_userData,&DRCaloTubesSDData::endEvent);


        //runAction().callAtEnd();
        //steppingAction().

        m_userData.sensitive = this;

        }

        /// Define collections created by this sensitivie action object
        template <> void Geant4SensitiveAction<DRCaloTubesSDData>::defineCollections() {
        m_collectionID = defineCollection<Geant4Tracker::Hit>(m_sensitive.readout().name());
        }

        /// Method for clear
        template <> void Geant4SensitiveAction<DRCaloTubesSDData>::clear(G4HCofThisEvent*) {
        m_userData.clear();
        }

        /// Method for generating hit(s) using the information of G4Step object.
        template <> G4bool
        Geant4SensitiveAction<DRCaloTubesSDData>::process(G4Step GEANT4_CONST_STEP * step, G4TouchableHistory* history) {
        return m_userData.process(step, history);
        }

        typedef Geant4SensitiveAction<DRCaloTubesSDData>  DRCaloTubesSDAction;

    }
}


#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE( DRCaloTubesSDAction )
*/