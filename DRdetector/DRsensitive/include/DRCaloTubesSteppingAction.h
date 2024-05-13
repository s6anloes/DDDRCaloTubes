#ifndef DRCaloTubesSteppingAction_h
#define DRCaloTubesSteppingAction_h 1

#include "DRCaloTubesEventAction.h"

#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4Poisson.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

// Data model
//#include "edm4hep/MCParticleCollection.h"
//#include "edm4hep/SimCalorimeterHitCollection.h"

//Forward declarations from Geant4
//
class G4OpBoundaryProcess;
enum G4OpBoundaryProcessStatus;
class G4ProcessVector;
enum G4ProcessVectorTypeIndex;

namespace dd4hep {
    namespace sim {
        class DRCaloTubesSteppingAction : public G4UserSteppingAction {
            public:
            DRCaloTubesSteppingAction(DRCaloTubesEventAction* eventAction);
            virtual ~DRCaloTubesSteppingAction();

            virtual void UserSteppingAction(const G4Step*);

            G4double ApplyBirks(const G4double& de, const G4double& steplength);
            G4int SmearSSignal( const G4double& satde );
            G4int SmearCSignal( );

            G4double GetDistanceToSiPM(const G4Step* step);
            G4int AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length);
            G4int AttenuateSSignal(const G4int& signal, const G4double& distance);
            G4int AttenuateCSignal(const G4int& signal, const G4double& distance);

            DRCaloTubesEventAction*  fEventAction;

            //void setSegmentation(dd4hep::DDSegmentation::GridDRcalo* seg) { pSeg = seg; }
            //void setEdepsCollection(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps = data; }
            //void setEdeps3dCollection(edm4hep::SimCalorimeterHitCollection* data) { m_Edeps3d = data; }
            //void setLeakagesCollection(edm4hep::MCParticleCollection* data) { m_Leakages = data; }

            //void setThreshold(const double thres) { m_thres = thres; }

            private:
                //void accumulate(unsigned int &prev, dd4hep::DDSegmentation::CellID& id64, float edep);
                //bool checkId(edm4hep::SimCalorimeterHit edep, dd4hep::DDSegmentation::CellID& id64);

                //void saveLeakage(G4Track* track, G4StepPoint* pre);
            const G4double fk_B = 0.126; //Birks constant

            const G4double fSAttenuationLength = 191.6*CLHEP::cm; // from test beam data
            const G4double fCAttenuationLength = 388.9*CLHEP::cm; // from test beam data

                

            G4OpBoundaryProcess* fOpProcess;

            // collections owned by DRCaloTubesEventAction
            //edm4hep::SimCalorimeterHitCollection* m_Edeps;
            //edm4hep::SimCalorimeterHitCollection* m_Edeps3d;
            //edm4hep::MCParticleCollection* m_Leakages;
            //double m_thres;
        };
}
}

#endif