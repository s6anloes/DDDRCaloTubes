#include "DRCaloTubesGrid.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
DRCaloTubesGrid::DRCaloTubesGrid(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "DRCaloTubesGrid";
  _description = "DRcalo segmentation based on the tower / (Cherenkov or Scintillation) fiber / SiPM hierarchy";

  // register all necessary parameters
  registerIdentifier("identifier_stave", "Cell ID identifier for stave", fStave, "stave");
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", fLayer, "layer");
  registerIdentifier("identifier_col", "Cell ID identifier for numCol", fCol, "col");
  registerIdentifier("identifier_row", "Cell ID identifier for numRow", fRow, "row");
  registerIdentifier("identifier_IsCherenkov", "Cell ID identifier for IsCherenkov", fIsCherenkov, "cherenkov");
  registerIdentifier("identifier_clad", "Cell ID identifier for clad", fClad, "clad");
  registerIdentifier("identifier_core", "Cell ID identifier for core", fCore, "core");

  registerParameter("grid_size", "Cell size", fGridSize, 2.0*mm, SegmentationParameter::LengthUnit);
  registerParameter("tower_theta", "Tower theta", fTheta, 1.0*deg, SegmentationParameter::AngleUnit);
  registerParameter("tower_phi", "Tower phi", fPhi, 1.0*deg, SegmentationParameter::AngleUnit);
}

DRCaloTubesGrid::DRCaloTubesGrid(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "DRCaloTubesGrid";
  _description = "DRcalo segmentation based on the tower / (Cherenkov or Scintillation) fiber / SiPM hierarchy";

  // register all necessary parameters
  registerIdentifier("identifier_stave", "Cell ID identifier for stave", fStave, "stave");
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", fLayer, "layer");
  registerIdentifier("identifier_col", "Cell ID identifier for numCol", fCol, "col");
  registerIdentifier("identifier_row", "Cell ID identifier for numRow", fRow, "row");
  registerIdentifier("identifier_IsCherenkov", "Cell ID identifier for IsCherenkov", fIsCherenkov, "cherenkov");
  registerIdentifier("identifier_clad", "Cell ID identifier for clad", fClad, "clad");
  registerIdentifier("identifier_core", "Cell ID identifier for core", fCore, "core");
  
  registerParameter("grid_size", "Cell size", fGridSize, 2.0*mm, SegmentationParameter::LengthUnit);
  registerParameter("tower_theta", "Tower theta", fTheta, 1.0*deg, SegmentationParameter::AngleUnit);
  registerParameter("tower_phi", "Tower phi", fPhi, 1.0*deg, SegmentationParameter::AngleUnit);
}

DRCaloTubesGrid::~DRCaloTubesGrid() {
}

Vector3D DRCaloTubesGrid::position(const CellID& cID) const {
  
  int col = Col(cID);
  int row = Row(cID);

  double x = fGridSize*col + (row&1)*fGridSize/2 + fGridSize/2;
  double y = fGridSize*sqrt(3)/2*row + fGridSize/2;
  
  return Vector3D(x, y, 0.0*mm);

}

Vector3D DRCaloTubesGrid::localPosition(const CellID& cID) const {
  int col = Col(cID);
  int row = Row(cID);

  double x = fGridSize*col + (row&1)*fGridSize/2 + fGridSize/2;
  double y = fGridSize*sqrt(3)/2*row + fGridSize/2;
  
  return Vector3D(x, y, 0.0*mm);
}


/// determine the cell ID based on the position
CellID DRCaloTubesGrid::cellID(const Vector3D& localPosition, const Vector3D& /* globalPosition */, const VolumeID& vID) const {
  // Following https://www.redblobgames.com/grids/hexagons/more-pixel-to-hex.html#charles-chambers
  double x = localPosition.X;
  double y = localPosition.Y;
  //std::cout<<"DRCaloTubeGrd::cellID: (x, y) = ("<<x/mm<<", "<<y/mm<<") mm"<<std::endl; 

  double q_float = 1/fGridSize * (x + y/sqrt(3.0));
  double r_float = 1/fGridSize * 2*y/sqrt(3.0);
  //std::cout<<"DRCaloTubeGrd::cellID: fGridSize = "<<fGridSize/mm<<" mm"<<std::endl; 
  //std::cout<<"DRCaloTubeGrd::cellID: (r, q)_float = ("<<r_float<<", "<<q_float<<")"<<std::endl; 

  std::pair<int, int> qr_reco = axial_round(q_float, r_float);
  int q_reco = qr_reco.first;
  int r_reco = qr_reco.second;
  //std::cout<<"DRCaloTubeGrd::cellID: (r, q)_reco = ("<<r_reco<<", "<<q_reco<<")"<<std::endl; 

  std::pair<int, int> col_row = axial2offset(q_reco, r_reco);
  int col = col_row.first;
  int row = col_row.second;

  int stave = Stave(vID);
  int layer = Layer(vID);

  return cellID(stave, layer, col, row);
}

CellID DRCaloTubesGrid::cellID(int stave, int layer, int col, int row) const {
  CellID stave_cID = static_cast<CellID>(stave);
  CellID layer_cID = static_cast<CellID>(layer);
  CellID col_cId = static_cast<CellID>(col);
  CellID row_cId = static_cast<CellID>(row);
  CellID cID = 0;

  _decoder->set(cID, fStave, stave_cID);
  _decoder->set(cID, fLayer, layer_cID);
  _decoder->set(cID, fCol, col_cId);
  _decoder->set(cID, fRow, row_cId);

  // CellID module = 1; // Fiber, SiPM, etc.
  // _decoder->set(cID, fModule, module);

  CellID isCeren = IsCherenkov(col, row) ? 1 : 0;
  _decoder->set(cID, fIsCherenkov, isCeren);

  return cID;
}


bool DRCaloTubesGrid::IsCherenkov(const CellID& aCellID) const {
  VolumeID isCeren = static_cast<VolumeID>(_decoder->get(aCellID, fIsCherenkov));
  return static_cast<bool>(isCeren);
}


bool DRCaloTubesGrid::IsCherenkov(int col, int row) const {
  return (row&1) ? true : false;
}

// int DRCaloTubesGrid::getLast32bits(const CellID& aCellID) const {
//   CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
//   int aId32 = (int)aId64;

//   return aId32;
// }

// CellID DRCaloTubesGrid::convertLast32to64(const int aId32) const {
//   CellID aId64 = (CellID)aId32;
//   aId64 <<= sizeof(int)*CHAR_BIT;

//   return aId64;
// }

int DRCaloTubesGrid::Row(const CellID& aCellID) const { 
  VolumeID row = static_cast<VolumeID>(_decoder->get(aCellID, fRow));
  return static_cast<int>(row);
}

int DRCaloTubesGrid::Col(const CellID& aCellID) const { 
  VolumeID col = static_cast<VolumeID>(_decoder->get(aCellID, fCol));
  return static_cast<int>(col);
}

int DRCaloTubesGrid::Stave(const VolumeID& volID) const { 
  VolumeID stave = static_cast<VolumeID>(_decoder->get(volID, fStave));
  return static_cast<int>(stave);
}

int DRCaloTubesGrid::Layer(const VolumeID& volID) const { 
  VolumeID layer = static_cast<VolumeID>(_decoder->get(volID, fLayer));
  return static_cast<int>(layer);
}

std::pair<int, int> DRCaloTubesGrid::axial_round(double q_float, double r_float) const {
  // Rounding to nearest hex according to https://www.redblobgames.com/grids/hexagons/#rounding

  double s_float = -q_float - r_float;

  int q_reco = round(q_float);
  int r_reco = round(r_float);
  int s_reco = round(s_float);

  double q_diff = abs(q_reco - q_float);
  double r_diff = abs(r_reco - r_float);
  double s_diff = abs(s_reco - s_float);

  if      (q_diff>r_diff && q_diff>s_diff)    q_reco = -r_reco-s_reco;
  else if (r_diff > s_diff)                   r_reco = -q_reco-s_reco;
/* else                                        s_reco = -q_reco-r_reco;    // Not needed because only q and r are returned */

  return std::make_pair(q_reco, r_reco);
}

std::pair<int, int> DRCaloTubesGrid::axial2offset(int q, int r) const {
  int row = -r;
  int col = q - (row + (row&1))/2;
  return std::make_pair(col, row);
}

} // namespace DDSegmentation
} // namespace dd4hep
