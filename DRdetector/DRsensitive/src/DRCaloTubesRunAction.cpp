

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
        DRCaloTubesRunAction::DRCaloTubesRunAction(DRCaloTubesEventAction* eventAction)
            : G4UserRunAction(), fEventAction(eventAction){ 
        
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
            analysisManager->CreateNtupleDColumn("EnergyScin");                     //0
            analysisManager->CreateNtupleDColumn("EnergyCher");                     //1
            analysisManager->CreateNtupleDColumn("NofCherDet");                     //2
            analysisManager->CreateNtupleDColumn("NofScinDet");                     //3
            analysisManager->CreateNtupleDColumn("TowerID", fEventAction->GetTowerIDs());                        //4
            analysisManager->CreateNtupleDColumn("FibreID", fEventAction->GetFibreIDs());                        //5
            analysisManager->CreateNtupleDColumn("NofDet", fEventAction->GetFibreSignals());                         //6
            analysisManager->CreateNtupleDColumn("PrimaryEnergy");                  //7
            analysisManager->CreateNtupleDColumn("PrimaryPDGID");                   //8
            analysisManager->CreateNtupleDColumn("PrimaryX");                       //9
            analysisManager->CreateNtupleDColumn("PrimaryY");                       //10
            analysisManager->CreateNtupleDColumn("PrimaryZ");                       //11
            analysisManager->CreateNtupleDColumn("PrimaryPX");                      //12
            analysisManager->CreateNtupleDColumn("PrimaryPY");                      //13
            analysisManager->CreateNtupleDColumn("PrimaryPZ");                      //14
            analysisManager->CreateNtupleDColumn("Leakage");                        //15
            analysisManager->CreateNtupleDColumn("NeutrinoLeakage");                //16
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

        }

        void DRCaloTubesRunAction::EndOfRunAction( const G4Run* ) {
        
            auto analysisManager = G4RootAnalysisManager::Instance();

            analysisManager->Write();
            analysisManager->CloseFile();

        }

    }
}
