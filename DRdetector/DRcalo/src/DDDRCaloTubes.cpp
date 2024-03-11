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

    Tube        calorimeter_solid(calo_inner_r, calo_inner_r+2*tower_length, calo_inner_half_length+2*tower_length); // Per design a "square" cylinder
    // Tube        calorimeter_solid(0, calo_inner_r+2*tower_length, 0); // Per design a "square" cylinder
    Volume      calorimeter_volume("calorimeter", calorimeter_solid, air);
    calorimeter_volume.setSolid(calorimeter_solid);
    calorimeter_volume.setMaterial(air);
    calorimeter_volume.setVisAttributes(description, "MyVis");

    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);


    DDDRCaloTubes::DRconstructor constructor(&description, entities, &sens);
    constructor.construct_calorimeter(calorimeter_volume);


    Transform3D calorimeter_tr(RotationZYX(0, 0, 0), Position(0, 0, 0));
    PlacedVolume calorimeter_placed = mother_volume.placeVolume(calorimeter_volume, calorimeter_tr);
    calorimeter_placed.addPhysVolID("system", det_id);
    s_detElement.setPlacement(calorimeter_placed);


    return s_detElement;
}


DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)