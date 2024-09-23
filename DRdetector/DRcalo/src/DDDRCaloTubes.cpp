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

    // Cylinder encompassing entire calorimeter
    xml_dim_t   x_dim                  = x_det.dimensions();
    double      calo_inner_r           = x_dim.inner_radius();
    double      calo_outer_r           = x_dim.outer_radius();
    double      calo_inner_half_length = x_dim.z_length();

    double tower_length = calo_outer_r - calo_inner_r;
    if (tower_length<=0*mm) throw std::runtime_error("Outer calorimeter radius needs to be larger than inner radius");

    Tube        barrel_tube_solid(calo_inner_r, calo_outer_r, calo_inner_half_length+tower_length);
    Cone        place_for_endcaps((calo_inner_half_length+tower_length)/2, 0, 0, 0, calo_outer_r);
    SubtractionSolid front_endcap_subtraced(barrel_tube_solid, place_for_endcaps, Position(0, 0, (calo_inner_half_length+tower_length)/2));
    Transform3D back_endcap_transform(RotationZYX(0, 0, 180*deg), Position(0, 0, -(calo_inner_half_length+tower_length)/2));
    SubtractionSolid both_endcaps_subtraced(front_endcap_subtraced, place_for_endcaps, back_endcap_transform);

    Volume      barrel_volume("calorimeter_barrel", both_endcaps_subtraced, air);
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