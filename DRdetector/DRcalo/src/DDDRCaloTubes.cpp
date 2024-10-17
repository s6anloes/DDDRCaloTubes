#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include "DRutils.h"
#include "DRconstructor.h"

using namespace dd4hep;


static Ref_t create_detector(Detector& description,
                             xml_h entities,
                             SensitiveDetector sens) 
{
    xml_det_t   x_det       = entities;    
    int         det_id      = x_det.id();
    std::string det_name    = x_det.nameStr();

    Material    air         = description.air();

    sens.setType("calorimeter");

    xml_dim_t   x_dim                  = x_det.dimensions();
    double      calo_inner_r           = x_dim.inner_radius();
    double      calo_outer_r           = x_dim.outer_radius();
    double      calo_inner_half_length = x_dim.z_length();

    double tower_length = calo_outer_r - calo_inner_r;
    if (tower_length<=0*mm) throw std::runtime_error("Outer calorimeter radius needs to be larger than inner radius");


    // Solid of barrel volume is created later (in the DRconstructor::construct_calorimeter) as union of the stave volumes
    // This prevents air gaps between the straight edges of the staves and the envelope (when using a tube as envelope for example)
    Volume barrel_volume("calorimeter_barrel");
    barrel_volume.setMaterial(air);
    barrel_volume.setVisAttributes(description, "scin_core_vis");


    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);


    DDDRCaloTubes::DRconstructor constructor(&description, entities, &sens);
    constructor.construct_calorimeter(barrel_volume);


    PlacedVolume barrel_placed = mother_volume.placeVolume(barrel_volume);
    barrel_placed.addPhysVolID("system", det_id);
    s_detElement.setPlacement(barrel_placed);



    return s_detElement;
}


DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)