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
    Volume calorimeter_volume("calorimeter");
    DDDRCaloTubes::DRconstructor constructor(&description, entities, &sens, calorimeter_volume);
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
    double      tower_phi              = x_dim.deltaphi();

    double tower_length = calo_outer_r - calo_inner_r;
    if (tower_length<=0*mm) throw std::runtime_error("Outer calorimeter radius needs to be larger than inner radius");

    xml_comp_t  x_tube      = x_det.child(_Unicode(tube));
    xml_comp_t  x_capillary          = x_tube.child(_Unicode(capillary));
    Material    capillary_material   = description.material(x_capillary.materialStr()); 
    double      capillary_outer_r    = x_capillary.outer_r();

    double barrel_endcap_angle = std::atan2(calo_inner_half_length, calo_inner_r);
    Tube        calorimeter_solid(calo_inner_r, calo_inner_r+2*tower_length, calo_inner_half_length+2*tower_length); // Per design a "square" cylinder
    // Tube        calorimeter_solid(0, calo_inner_r+2*tower_length, 0); // Per design a "square" cylinder
    // std::string calorimeter_name = "calorimeter";
    // Volume      calorimeter_volume(calorimeter_name, calorimeter_solid, air);
    calorimeter_volume.setSolid(calorimeter_solid);
    calorimeter_volume.setMaterial(air);
    calorimeter_volume.setVisAttributes(description, "MyVis");

    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);

    /*****************************************
     *Calculation of required phi parameters *
     *****************************************/
    // ---------------------------------------------------------------------------------------------------
    int num_cols;
    double tan_phi = std::tan(tower_phi);
    unsigned int num_phi_towers;
    double num_phi_towers_d = 360.0*deg/tower_phi;
    // Check if num_phi_towers is a whole number
    if (DDDRCaloTubes::check_for_integer(num_phi_towers_d)) num_phi_towers = static_cast<unsigned int>(num_phi_towers_d);
    else throw std::runtime_error("Not an integer number of towers in phi direction");


    // ---------------------------------------------------------------------------------------------------
    unsigned int layer = 0;
    for (double covered_theta=0*deg; covered_theta<barrel_endcap_angle; layer++) 
    {
        // constructor.calculate_theta_parameters();
        // double theta = 90*deg - covered_theta;
        double delta_theta;
        Assembly tower_volume("tower");
        constructor.construct_tower(tower_volume, delta_theta);
        double phi = 0*deg;
        for (unsigned int stave=1; stave<=1; stave++, phi+=tower_phi)
        {
            unsigned int tower_id = stave + layer*num_phi_towers;
            constructor.place_tower(calorimeter_volume, tower_volume, stave, layer, tower_id, phi);
        }

        covered_theta += delta_theta;
        constructor.increase_covered_theta(delta_theta);
        
        // if (layer >= 1) break;
    }


    Transform3D calorimeter_tr(RotationZYX(0, 0, 0), Position(0, 0, 0));
    PlacedVolume calorimeter_placed = mother_volume.placeVolume(calorimeter_volume, calorimeter_tr);
    calorimeter_placed.addPhysVolID("system", det_id);
    s_detElement.setPlacement(calorimeter_placed);


    return s_detElement;
}


DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)