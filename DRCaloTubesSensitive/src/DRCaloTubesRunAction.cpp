

//Includers from project files
//
#include "DRCaloTubesRunAction.h"

//Includers from Geant4
//
//#include "g4root.hh"
#include "G4RootAnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//Includers from C++
//
#include <string>

namespace dd4hep {
    namespace sim {

        //Define constructor
        //
        DRCaloTubesRunAction::DRCaloTubesRunAction()
            : G4UserRunAction(){ 
        
            //print event number per each event (default, can be overwritten with macro)
            //
            G4RunManager::GetRunManager()->SetPrintProgress(1);     

            //Instantiate analysis manager
            //
            auto analysisManager = G4RootAnalysisManager::Instance();
            analysisManager->SetVerboseLevel( 1 );
            analysisManager->SetNtupleMerging( 1 );

            //Using ROOT as analysisManager type, print it
            //
            G4cout << "DRCaloTubes-> Using " << analysisManager->GetType() << G4endl;

            //Define ntuple structure
            //
            analysisManager->CreateNtuple("DRCaloTubesout", "DRCaloTubesoutput");
            analysisManager->CreateNtupleDColumn("NofCherDet");                     //0
            analysisManager->CreateNtupleDColumn("NofScinDet");                     //1
            analysisManager->CreateNtupleDColumn("FibreIDsCher");
            analysisManager->CreateNtupleDColumn("FibreSignalsCher");
            analysisManager->CreateNtupleDColumn("FibreIDsScin");
            analysisManager->CreateNtupleDColumn("FibreSignalsScin");
            analysisManager->FinishNtuple();
            
        }

        //Define de-constructor
        //
        DRCaloTubesRunAction::~DRCaloTubesRunAction(){
        
            //Delete only instance of G4RootAnalysisManager
            //
            delete G4RootAnalysisManager::Instance();  

        }

        //Define BeginOfRunAction() and EndOfRunAction() methods
        //
        void DRCaloTubesRunAction::BeginOfRunAction( const G4Run* Run )  { 
            
            //Save random seeds (optional)
            //
            //G4RunManager::GetRunManager()->SetRandomNumberStore( true );
            
            //Open output file, one per Run
            //
            auto analysisManager = G4RootAnalysisManager::Instance();
            std::string runnumber = std::to_string( Run->GetRunID() );
            G4String outputfile = "DRCaloTubesout_Run"+runnumber;
            analysisManager->OpenFile( outputfile );

            std::string rootfile_name = "DRCaloTubesout_Run"+runnumber+".root";
            fFile = new TFile(rootfile_name.c_str(), "RECREATE");
            fTree = new TTree("DRCaloTubesData", "Tree with DRCaloTubes data");
            fTree->Branch("NofCherDet", &NofCherDet);
            fTree->Branch("NofScinDet", &NofScinDet);
            fTree->Branch("FibreIDsCher", &FibreIDsCher);
            fTree->Branch("FibreSignalsCher", &FibreSignalsCher);
            fTree->Branch("FibreIDsScin", &FibreIDsScin);
            fTree->Branch("FibreSignalsScin", &FibreSignalsScin);

        }

        void DRCaloTubesRunAction::EndOfRunAction( const G4Run* ) {
        
            auto analysisManager = G4RootAnalysisManager::Instance();

            analysisManager->Write();
            analysisManager->CloseFile();

            fTree->Write();
            fFile->Close();

        }

        void DRCaloTubesRunAction::Reset() {
            NofCherDet=0; 
            NofScinDet=0; 
            FibreSignalsCher.clear(); 
            FibreSignalsScin.clear();
        }

        void DRCaloTubesRunAction::Fill(const G4int cher, const G4int scin, 
                                        std::map<unsigned int, G4int> fibrecher, 
                                        std::map<unsigned int, G4int> fibrescin) {

            FibreSignalsCher.clear();
            FibreSignalsScin.clear();
            FibreIDsCher.clear();
            FibreIDsScin.clear();
            
            NofCherDet = cher;
            NofScinDet = scin;

            /* unsigned int cher_max = fibrecher.empty() ? 0 : fibrecher.rbegin()->first;
            unsigned int scin_max = fibrescin.empty() ? 0 : fibrescin.rbegin()->first;
            unsigned int nfibres = cher_max > scin_max ? cher_max : scin_max;
            for (unsigned int i=0; i<=nfibres; i++) 
            {
                if (fibrecher.count(i)) {
                    FibreSignalsCher.push_back(fibrecher.at(i));
                } else {
                    FibreSignalsCher.push_back(0);
                }

                if (fibrescin.count(i)) {
                    FibreSignalsScin.push_back(fibrescin.at(i));
                } else {
                    FibreSignalsScin.push_back(0);
                }
            } */
            for (auto const& [key, val]: fibrecher)
            {
                FibreIDsCher.push_back(key);
                FibreSignalsCher.push_back(val);
            }

            for (auto const& [key, val]: fibrescin)
            {
                FibreIDsScin.push_back(key);
                FibreSignalsScin.push_back(val);
            }

            fTree->Fill();

            FibreSignalsCher.clear();
            FibreSignalsScin.clear();
            FibreIDsCher.clear();
            FibreIDsScin.clear();
        }

    }
}
