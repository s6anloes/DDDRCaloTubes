#ifndef GridDRcalo_h
#define GridDRcalo_h 1

#include "DDSegmentation/Segmentation.h"
#include <DD4hep/DD4hepUnits.h>

namespace dd4hep {
namespace DDSegmentation {
class DRCaloTubeGrid : public Segmentation {
public:
  /// default constructor using an arbitrary type
  DRCaloTubeGrid(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  DRCaloTubeGrid(const BitFieldCoder* decoder);
  /// destructor
  virtual ~DRCaloTubeGrid() override;

  //  Determine the global(local) position based on the cell ID.
  virtual Vector3D position(const CellID& aCellID) const;

  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;

  CellID cellID(int r, int q) const;

  void setGridSize(double grid) { fGridSize = grid; }

  bool IsCherenkov(const CellID& aCellID) const;
  bool IsCherenkov(int r, int q) const;

  int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
  int getLast32bits(const CellID& aCellID) const;
  CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
  CellID convertLast32to64(const int aId32) const;

  bool IsCherenkov(const int& aId32) const { return IsCherenkov( convertLast32to64(aId32) ); }

  int R(const CellID& aCellID) const;
  int Q(const CellID& aCellID) const;

  std::pair<int, int> axial_round(double r_float, double q_float) const;

  inline const std::string& fieldNameR() const { return fR; }
  inline const std::string& fieldNameQ() const { return fQ; }
  inline const std::string& fieldNameIsCherenkov() const { return fIsCherenkov; }
  inline const std::string& fieldNameModule() const { return fModule; }

  inline void setFieldNameR(const std::string& fieldName) { fR = fieldName; }
  inline void setFieldNameQ(const std::string& fieldName) { fQ = fieldName; }
  inline void setFieldNameIsCherenkov(const std::string& fieldName) { fIsCherenkov = fieldName; }
  inline void setFieldNameModule(const std::string& fieldName) { fModule = fieldName; }


protected:
  std::string fR;
  std::string fQ;
  std::string fIsCherenkov;
  std::string fModule;

  double fGridSize;

};
}
}

#endif
