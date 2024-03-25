#ifndef GridDRcalo_h
#define GridDRcalo_h 1

#include "DDSegmentation/Segmentation.h"
#include <DD4hep/DD4hepUnits.h>

namespace dd4hep {
namespace DDSegmentation {
class DRCaloTubesGrid : public Segmentation {
public:
  /// default constructor using an arbitrary type
  DRCaloTubesGrid(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  DRCaloTubesGrid(const BitFieldCoder* decoder);
  /// destructor
  virtual ~DRCaloTubesGrid() override;

  //  Determine the global(local) position based on the cell ID.
  virtual Vector3D position(const CellID& aCellID) const;
  Vector3D localPosition(const CellID& aCellID) const;

  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;
  CellID cellID(int stave, int layer, int col, int row) const;

  void setGridSize(double grid) { fGridSize = grid; }

  bool IsCherenkov(const CellID& aCellID) const;
  bool IsCherenkov(int col, int row) const;

  // int getFirst32bits(const CellID& aCellID) const { return (int)aCellID; }
  // int getLast32bits(const CellID& aCellID) const;
  // CellID convertFirst32to64(const int aId32) const { return (CellID)aId32; }
  // CellID convertLast32to64(const int aId32) const;

  // bool IsCherenkov(const int& aId32) const { return IsCherenkov( convertLast32to64(aId32) ); }

  int Row(const CellID& aCellID) const;
  int Col(const CellID& aCellID) const;
  int Stave(const VolumeID& aCellID) const;
  int Layer(const VolumeID& aCellID) const;

  std::pair<int, int> axial_round(double q_float, double r_float) const;
  std::pair<int, int> axial2offset(int r, int q) const;

  inline const std::string& fieldNameStave() const { return fStave; }
  inline const std::string& fieldNameLayer() const { return fLayer; }
  inline const std::string& fieldNameAir() const { return fAir; }
  inline const std::string& fieldNameRow() const { return fRow; }
  inline const std::string& fieldNameCol() const { return fCol; }
  inline const std::string& fieldNameIsCherenkov() const { return fIsCherenkov; }
  inline const std::string& fieldNameClad() const { return fClad; }
  inline const std::string& fieldNameCore() const { return fCore; }


  inline void setFieldNameStave(const std::string& fieldName) { fStave = fieldName; }
  inline void setFieldNameLayer(const std::string& fieldName) { fLayer = fieldName; }
  inline void setFieldNameAir(const std::string& fieldName) { fAir = fieldName; }
  inline void setFieldNameRow(const std::string& fieldName) { fRow = fieldName; }
  inline void setFieldNameCol(const std::string& fieldName) { fCol = fieldName; }
  inline void setFieldNameIsCherenkov(const std::string& fieldName) { fIsCherenkov = fieldName; }
  inline void setFieldNameClad(const std::string& fieldName) { fClad = fieldName; }
  inline void setFieldNameCore(const std::string& fieldName) { fCore = fieldName; }



protected:
  std::string fStave;
  std::string fLayer;
  std::string fAir;

  std::string fRow;
  std::string fCol;
  std::string fIsCherenkov;
  std::string fClad;
  std::string fCore;

  double fGridSize;
  double fTheta;
  double fPhi;

};
}
}

#endif
