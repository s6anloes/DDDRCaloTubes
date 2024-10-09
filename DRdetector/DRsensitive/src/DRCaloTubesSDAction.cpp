#include "DD4hep/Version.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4RunAction.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4Poisson.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

#include "DRCaloTubesRunAction.h"
#include "DRCaloTubesEventAction.h"
#include "DRCaloTubesSteppingAction.h"


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



//class dd4hep::sim::DRCaloTubesEventAction;

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  
    /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
    namespace sim   {

        

        struct DRCaloTubesSDData {

            Geant4Sensitive*  sensitive{};

            DRCaloTubesRunAction*       fRunAction;
            DRCaloTubesEventAction*     fEventAction;
            DRCaloTubesSteppingAction*  fSteppingAction;


            DRCaloTubesSDData()
            {
            }


                /// Clear collected information and restart for new hit
            void clear()  {
                // nothing to clear
            }


            /// Method for generating hit(s) using the information of G4Step object.
            G4bool process(G4Step GEANT4_CONST_STEP * step, G4TouchableHistory* ) 
            {
                fSteppingAction->UserSteppingAction(step);
                return true;
            } // end of process()


            void beginRun(const G4Run* run)
            {
                fEventAction = new DRCaloTubesEventAction();
                fRunAction = new DRCaloTubesRunAction(fEventAction);
                fSteppingAction = new DRCaloTubesSteppingAction(fEventAction);
                fRunAction->BeginOfRunAction(run);
            }

            void endRun(const G4Run* run)
            {
                fRunAction->EndOfRunAction(run);
            }

            /// Pre-event action callback
            void beginEvent(const G4Event* event)
            {
                
                fEventAction->BeginOfEventAction(event);
            }

            /// Post-event action callback
            void endEvent(const G4Event* event)   {

                fEventAction->EndOfEventAction(event);
            }
        
        };


        /// Initialization overload for specialization
        template <> void Geant4SensitiveAction<DRCaloTubesSDData>::initialize() 
        {
            
            eventAction().callAtBegin(&m_userData, &DRCaloTubesSDData::beginEvent);
            eventAction().callAtEnd(&m_userData,&DRCaloTubesSDData::endEvent);

            runAction().callAtBegin(&m_userData, &DRCaloTubesSDData::beginRun);
            runAction().callAtEnd(&m_userData, &DRCaloTubesSDData::endRun);

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
