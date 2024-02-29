#ifndef DRCaloTubesHit_h
#define DRCaloTubesHit_h 1

/// Framework include files
#include "DDG4/Geant4Data.h"
// #include "DD4hep/Objects.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"


// typedef ROOT::Math::XYZVector Position;
// typedef ROOT::Math::XYZVector Direction;

class DRCaloTubesHit : public dd4hep::sim::Geant4Calorimeter::Hit {

    public:
        int NofCherDet, NofScinDet;

        /// Default constructor
        DRCaloTubesHit() = default;
        /// Initializing constructor
        DRCaloTubesHit(const dd4hep::Position& cell_pos) : dd4hep::sim::Geant4Calorimeter::Hit(cell_pos),NofCherDet(0),NofScinDet(0) {
        }

        /// Default destructor
        virtual ~DRCaloTubesHit() = default;

};

#endif