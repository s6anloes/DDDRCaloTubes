#include "DDDRCaloTubes/DRCaloTubeGrid.h"

#include <climits>
#include <cmath>
#include <stdexcept>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
DRCaloTubeGrid::DRCaloTubeGrid(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "DRCaloTubeGrid";
  _description = "DRcalo segmentation based on the tower / (Cherenkov or Scintillation) fiber / SiPM hierarchy";

  // register all necessary parameters
  registerIdentifier("identifier_r", "Cell ID identifier for numRow", fR, "r");
  registerIdentifier("identifier_q", "Cell ID identifier for numCol", fQ, "q");
  registerIdentifier("identifier_IsCherenkov", "Cell ID identifier for IsCherenkov", fIsCherenkov, "c");
  registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");
  registerParameter("grid_size", "Cell size", fGridSize, 2.0*mm, SegmentationParameter::LengthUnit);
}

DRCaloTubeGrid::DRCaloTubeGrid(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "DRCaloTubeGrid";
  _description = "DRcalo segmentation based on the tower / (Cherenkov or Scintillation) fiber / SiPM hierarchy";

  // register all necessary parameters
  registerIdentifier("identifier_r", "Cell ID identifier for numRow", fR, "r");
  registerIdentifier("identifier_q", "Cell ID identifier for numCol", fQ, "q");
  registerIdentifier("identifier_IsCherenkov", "Cell ID identifier for IsCherenkov", fIsCherenkov, "cherenkov");
  registerIdentifier("identifier_module", "Cell ID identifier for module", fModule, "module");
  registerParameter("grid_size", "Cell size", fGridSize, 2.0*mm, SegmentationParameter::LengthUnit);
}

DRCaloTubeGrid::~DRCaloTubeGrid() {
}

Vector3D DRCaloTubeGrid::position(const CellID& cID) const {
  // Following https://www.redblobgames.com/grids/hexagons/#hex-to-pixel-axial with different definitions for q and
  int q = Q(cID);
  int r = R(cID);

  // q-basevector: fGridSize*(   1,         0) (changes q without changing r)
  // r-basevector: fGridSize*(-1/2, sqrt(3)/2) (changes r without changing q)
  double x = fGridSize * (1*q - 1/2*r);
  double y = fGridSize * sqrt(3)/2*r;

  return Position(x, y, 0.0*mm);

}


/// determine the cell ID based on the position
CellID DRCaloTubeGrid::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition, const VolumeID& vID) const {
  // Following https://www.redblobgames.com/grids/hexagons/more-pixel-to-hex.html#charles-chambers
  double x = globalPosition.X;
  double y = globalPosition.Y;
  //std::cout<<"DRCaloTubeGrd::cellID: (x, y) = ("<<x/mm<<", "<<y/mm<<") mm"<<std::endl; 

  double q_float = 1/fGridSize * (x + y/sqrt(3.0));
  double r_float = 1/fGridSize * 2*y/sqrt(3.0);
  //std::cout<<"DRCaloTubeGrd::cellID: fGridSize = "<<fGridSize/mm<<" mm"<<std::endl; 
  //std::cout<<"DRCaloTubeGrd::cellID: (r, q)_float = ("<<r_float<<", "<<q_float<<")"<<std::endl; 

  std::pair<int, int> rq_reco = axial_round(r_float, q_float);
  int r_reco = rq_reco.first;
  int q_reco = rq_reco.second;
  //std::cout<<"DRCaloTubeGrd::cellID: (r, q)_reco = ("<<r_reco<<", "<<q_reco<<")"<<std::endl; 

  return cellID(r_reco, q_reco);
}

CellID DRCaloTubeGrid::cellID(int r, int q) const {
  CellID r_cId = static_cast<CellID>(r);
  CellID q_cId = static_cast<CellID>(q);
  CellID cID = 0;
  _decoder->set(cID, fR, r_cId);
  _decoder->set(cID, fQ, q_cId);

  CellID module = 1; // Fiber, SiPM, etc.
  _decoder->set(cID, fModule, module);

  CellID isCeren = IsCherenkov(r, q) ? 1 : 0;
  _decoder->set(cID, fIsCherenkov, isCeren);

  return cID;
}


bool DRCaloTubeGrid::IsCherenkov(const CellID& aCellID) const {
  VolumeID isCeren = static_cast<VolumeID>(_decoder->get(aCellID, fIsCherenkov));
  return static_cast<bool>(isCeren);
}


bool DRCaloTubeGrid::IsCherenkov(int r, int q) const {
  return (r&1) ? true : false;
}

int DRCaloTubeGrid::getLast32bits(const CellID& aCellID) const {
  CellID aId64 = aCellID >> sizeof(int)*CHAR_BIT;
  int aId32 = (int)aId64;

  return aId32;
}

CellID DRCaloTubeGrid::convertLast32to64(const int aId32) const {
  CellID aId64 = (CellID)aId32;
  aId64 <<= sizeof(int)*CHAR_BIT;

  return aId64;
}

int DRCaloTubeGrid::R(const CellID& aCellID) const { 
  VolumeID R = static_cast<VolumeID>(_decoder->get(aCellID, fR));
  return static_cast<int>(R);
}

int DRCaloTubeGrid::Q(const CellID& aCellID) const { 
  VolumeID Q = static_cast<VolumeID>(_decoder->get(aCellID, fQ));
  return static_cast<int>(Q);
}

std::pair<int, int> DRCaloTubeGrid::axial_round(double r_float, double q_float) const {
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

  return std::make_pair(r_reco, q_reco);
}

} // namespace DDSegmentation
} // namespace dd4hep
