#include "DRCaloTubesEventAction.h"
#include "G4SystemOfUnits.hh"

//#include "G4RunManager.hh"

namespace dd4hep {
    namespace sim {
        DRCaloTubesEventAction::DRCaloTubesEventAction(DRCaloTubesRunAction* runAction): G4UserEventAction(),
            fRunAction(runAction), NofCherDet(0), NofScinDet(0){
                //FibreSignalsCher.resize(nfibres);
                //FibreSignalsScin.resize(nfibres);
                //for (int i=0; i<nfibres; i++){
                //    FibreSignalsCher.push_back(0);
                //    FibreSignalsScin.push_back(0);
                //}
            }

        DRCaloTubesEventAction::~DRCaloTubesEventAction() {}

        void DRCaloTubesEventAction::BeginOfEventAction(const G4Event*) {
            NofScinDet = 0;
            NofCherDet = 0;
            FibreSignalsCher.clear();
            FibreSignalsScin.clear();

            
            x_truth = -999999999*m;
            y_truth = -999999999*m;
            z_truth = -999999999*m;
            leakage = 0*eV; 
            first_step_number = -1;

            fRunAction->Reset();
        }

        void DRCaloTubesEventAction::EndOfEventAction(const G4Event*) {
            std::cout<<"**********************************END OF EVENT********************************************"<<std::endl;
            std::cout<<"NofScinDet = "<<NofScinDet<<std::endl;
            std::cout<<"NofCherDet = "<<NofCherDet<<std::endl;

            fRunAction->Fill(NofCherDet, NofScinDet, FibreSignalsCher, FibreSignalsScin, x_truth, y_truth, z_truth, leakage);
        }
    } // namespace sim
} // namespace drc
