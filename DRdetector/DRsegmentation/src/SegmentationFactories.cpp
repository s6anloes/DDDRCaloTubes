#include "DD4hep/Factories.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace {
template <typename T>
dd4hep::SegmentationObject* create_segmentation(const dd4hep::BitFieldCoder* decoder) {
  return new dd4hep::SegmentationWrapper<T>(decoder);
}
}

#include "DRCaloTubesGrid.h"
DECLARE_SEGMENTATION(DRCaloTubesGrid, create_segmentation<dd4hep::DDSegmentation::DRCaloTubesGrid>)
