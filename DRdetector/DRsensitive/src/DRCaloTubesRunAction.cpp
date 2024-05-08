

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

            // std::string rootfile_name = "DRCaloTubesout_Run"+runnumber+".root";
            // fFile = new TFile(rootfile_name.c_str(), "RECREATE");
            // fTree = new TTree("DRCaloTubesData", "Tree with DRCaloTubes data");
            // fTree->Branch("NofCherDet", &NofCherDet);
            // fTree->Branch("NofScinDet", &NofScinDet);
            // fTree->Branch("FibreIDsCher", &FibreIDsCher);
            // fTree->Branch("FibreSignalsCher", &FibreSignalsCher);
            // fTree->Branch("FibreIDsScin", &FibreIDsScin);
            // fTree->Branch("FibreSignalsScin", &FibreSignalsScin);
            // fTruth = new TTree("MCTruth", "Tree with MC truth information");
            // fTruth->Branch("x_truth", &x_truth);
            // fTruth->Branch("y_truth", &y_truth);
            // fTruth->Branch("z_truth", &z_truth);
            // fTruth->Branch("leakage", &leakage);

        }

        void DRCaloTubesRunAction::EndOfRunAction( const G4Run* ) {
        
            auto analysisManager = G4RootAnalysisManager::Instance();

            analysisManager->Write();
            analysisManager->CloseFile();

            // fTree->Write();
            // // fTruth->Write();
            // fFile->Close();

        }

        void DRCaloTubesRunAction::Reset() {
            NofCherDet=0; 
            NofScinDet=0; 
            // FibreSignalsCher.clear(); 
            // FibreSignalsScin.clear();
        }

        void DRCaloTubesRunAction::Fill(const G4int cher, const G4int scin) {

            // FibreSignalsCher.clear();
            // FibreSignalsScin.clear();
            // FibreIDsCher.clear();
            // FibreIDsScin.clear();
            
            NofCherDet = cher;
            NofScinDet = scin;

            // for (auto const& [key, val]: fibrecher)
            // {
            //     FibreIDsCher.push_back(key);
            //     FibreSignalsCher.push_back(val);
            // }

            // for (auto const& [key, val]: fibrescin)
            // {
            //     FibreIDsScin.push_back(key);
            //     FibreSignalsScin.push_back(val);
            // }

            G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, NofCherDet);
            analysisManager->FillNtupleDColumn(1, NofScinDet);
            analysisManager->AddNtupleRow();

            // fTree->Fill();

            // FibreSignalsCher.clear();
            // FibreSignalsScin.clear();
            // FibreIDsCher.clear();
            // FibreIDsScin.clear();
        }

    }
}
