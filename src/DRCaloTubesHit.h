/// Framework include files
#include "DDG4/Geant4Data.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"


typedef ROOT::Math::XYZVector Position;
typedef ROOT::Math::XYZVector Direction;

class DRCaloTubesHit : public dd4hep::sim::Geant4Calorimeter::Hit {

    public:
        int NofCherDet, NofScinDet;

        /// Default constructor
        DRCaloTubesHit() = default;
        /// Initializing constructor
        DRCaloTubesHit(const Position& cell_pos) : dd4hep::sim::Geant4Calorimeter::Hit(cell_pos),NofCherDet(0),NofScinDet(0) {
        }

        /// Default destructor
        virtual ~DRCaloTubesHit() = default;

};


// CINT configuration
#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/// Define namespaces
#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::sim;
#pragma link C++ class     DRCaloTubesHit+;
#endif