#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include "DRutils.h"
#include "DRconstructor.h"

using namespace dd4hep;

Assembly construct_barrel_tower(Detector& description,
                         xml_h& entities,
                         SensitiveDetector& sens,
                         Volume& calorimeter_volume,
                         double covered_theta,
                         double& delta_theta,
                         int num_cols,
                         double phi_back_shift,
                         Position& tower_position,
                         DDDRCaloTubes::DRconstructor& constructor);

void place_barrel_tower(Volume& calorimeter_volume,
                 Assembly& tower_volume,
                 unsigned int stave, 
                 unsigned int layer,
                 unsigned int tower_id,
                 Position tower_position,
                 double covered_theta,
                 double phi);

Assembly construct_endcap_tower(Detector& description,
                         xml_h& entities,
                         SensitiveDetector& sens,
                         Volume& calorimeter_volume,
                         double covered_theta,
                         double& delta_theta,
                         int num_cols,
                         double phi_back_shift,
                         Position& tower_position);

// int fast_floor(double x);
// int fast_ceil(double x);
// bool check_for_integer(double x);

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

    /* (Minimal) tower width with given inner calorimeter radius and tower_phi
       Might need to increase width to ensure integer number of tubes in phi direction
       at the cost of having to place the tower further back (calo_inner_r + epsilon) */
    double tower_min_frontface_width = calo_inner_r*tan_phi;
    int num_front_cols = DDDRCaloTubes::fast_ceil(tower_min_frontface_width/(2*capillary_outer_r));
    double tower_frontface_width = num_front_cols*2*capillary_outer_r;
    
    double phi_back_shift = tower_frontface_width/tan_phi - calo_inner_r;
    calo_inner_r = tower_frontface_width/tan_phi;

    // Calculate how many tubes there are in the back face
    double tower_outer_r_phi = calo_inner_r + tower_length; 
    double tower_max_backface_width = tower_outer_r_phi * tan_phi;
    int num_back_cols = DDDRCaloTubes::fast_floor(tower_max_backface_width/(2*capillary_outer_r));

    num_cols = num_back_cols;

    // ---------------------------------------------------------------------------------------------------
    unsigned int layer = 0;
    for (double covered_theta=0*deg; covered_theta<barrel_endcap_angle; layer++) 
    {
        constructor.calculate_theta_parameters();
        // double theta = 90*deg - covered_theta;
        double delta_theta;
        Position tower_position;
        // Assembly tower_volume = construct_barrel_tower(description, entities, sens, calorimeter_volume, covered_theta, delta_theta, num_cols, phi_back_shift, tower_position, constructor);
        Assembly tower_volume("tower");
        constructor.construct_tower(tower_volume, delta_theta, tower_position);
        double phi = 0*deg;
        for (unsigned int stave=1; stave<=1; stave++, phi+=tower_phi)
        {
            unsigned int tower_id = stave + layer*num_phi_towers;
            // place_barrel_tower(calorimeter_volume, tower_volume, stave, layer, tower_id, tower_position, covered_theta, phi);
            constructor.place_tower(calorimeter_volume, tower_volume, stave, layer, tower_id, tower_position, covered_theta, phi);
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

Assembly construct_barrel_tower(Detector& description,
                     xml_h& entities,
                     SensitiveDetector& sens,
                     Volume& calorimeter_volume,
                     double covered_theta,
                     double& delta_theta,
                     int num_cols,
                     double phi_back_shift,
                     Position& tower_position,
                     DDDRCaloTubes::DRconstructor& constructor) 
{
    xml_det_t   x_det       = entities;    
    std::string det_name    = x_det.nameStr();

    xml_dim_t   x_dim                  = x_det.dimensions();
    double      calo_inner_r           = x_dim.inner_radius();
    double      calo_outer_r           = x_dim.outer_radius();
    double      calo_inner_half_length = x_dim.z_length();
    double      tower_theta            = x_dim.deltatheta();
    double      tower_phi              = x_dim.deltaphi();

    double z_half = (calo_outer_r - calo_inner_r)/2;

    double tan_phi = std::tan(tower_phi);

    double barrel_endcap_angle = std::atan2(calo_inner_half_length, calo_inner_r);

    
    Assembly      tower_volume("tower");

    // Get parameters for tube construction
    xml_comp_t  x_tube      = x_det.child(_Unicode(tube));

    xml_comp_t  x_capillary          = x_tube.child(_Unicode(capillary));
    Material    capillary_material   = description.material(x_capillary.materialStr()); 
    double      capillary_outer_r    = x_capillary.outer_r();

    xml_comp_t  x_scin_clad          = x_tube.child(_Unicode(scin_clad));
    Material    scin_clad_material   = description.material(x_scin_clad.materialStr()); 
    double      scin_clad_outer_r    = x_scin_clad.outer_r();

    xml_comp_t  x_scin_core         = x_tube.child(_Unicode(scin_core));
    Material    scin_core_material  = description.material(x_scin_core.materialStr()); 
    double      scin_core_outer_r   = x_scin_core.outer_r();

    xml_comp_t  x_cher_clad          = x_tube.child(_Unicode(cher_clad));
    Material    cher_clad_material   = description.material(x_cher_clad.materialStr()); 
    double      cher_clad_outer_r    = x_cher_clad.outer_r();

    xml_comp_t  x_cher_core         = x_tube.child(_Unicode(cher_core));
    Material    cher_core_material  = description.material(x_cher_core.materialStr()); 
    double      cher_core_outer_r   = x_cher_core.outer_r();


    // Constants used through the function
    const double D = 4.0*capillary_outer_r/sqrt(3.0); // Long diagonal of hexagaon with capillary_outer_r as inradius
    const double V = 3.0*D/4.0;                       // Vertical spacing for pointy top oriented tubes

    /* (Minimal) tower width with given inner calorimeter radius and tower_phi
       Might need to increase width to ensure integer number of tubes in phi direction
       at the cost of having to place the tower further back (calo_inner_r + epsilon) */
    double tower_min_frontface_width = calo_inner_r*tan_phi;
    int num_front_cols = DDDRCaloTubes::fast_ceil(tower_min_frontface_width/(2*capillary_outer_r));
    double tower_frontface_width = num_front_cols*2*capillary_outer_r;
    
    calo_inner_r += phi_back_shift;

    // Calculate how many tubes there are in the back face
    double tower_outer_r_phi = calo_inner_r + 2*z_half; 
    double tower_max_backface_width = tower_outer_r_phi * tan_phi;
    int num_back_cols = DDDRCaloTubes::fast_floor(tower_max_backface_width/(2*capillary_outer_r));

    num_cols = num_back_cols;


    // Overlap between tubes cause by placing the tube in the gap of the previous row
    double overlap = 2*capillary_outer_r - V;


    // Calculate tower dimensions
    double covered_z = std::tan(covered_theta)*calo_inner_r;

    bool last_tower = (covered_theta+tower_theta>barrel_endcap_angle) ? true : false ;
    double tower_max_theta = covered_theta+tower_theta;
    double tower_max_z = std::tan(tower_max_theta)*calo_inner_r - covered_z; // Max distance the front face of this tower covers in z (not regarding how many cores actually fit)
    double tower_max_frontface_height = std::cos(covered_theta)*tower_max_z; // Tower height (in theta direction) without regarding how many tubes actually fit

    if (tower_max_frontface_height < 2*capillary_outer_r)
    {
        throw std::runtime_error("Can't construct tower with given tower_theta and calo_inner_radius");
    }
    
    // Calculate how many tubes fit at the front face for the given tower theta coverage.
    // This number will serve as the new covered theta since it is important to not have any gaps in the front face
    int num_front_rows = 1 + DDDRCaloTubes::fast_floor((tower_max_frontface_height-2*capillary_outer_r) / V);
    if (num_front_rows&1 ) num_front_rows++; // Make sure that front face ends on row with offset (i.e. even number of rows)

    double tower_frontface_height;
    double back_shift;
    double rad_distance;
    double this_tower_theta;

    tower_frontface_height = 2*capillary_outer_r + (num_front_rows-1)*V;
    
    // Distance by which straight edge of this tower is shifted backwards to ensure inner radius of calorimeter
    back_shift = std::tan(covered_theta)*tower_frontface_height;                      

    // Radial distance to exceed 2.5m inner radius of calorimeter for this tower
    rad_distance = calo_inner_r/std::cos(covered_theta);

    this_tower_theta = std::atan2(tower_frontface_height-overlap, rad_distance+back_shift);
    double tan_theta = std::tan(this_tower_theta);
    // double missing_theta = tower_theta - this_tower_theta;

    // Distance the front face of this tower covers in z
    double this_tower_z = std::tan(covered_theta+this_tower_theta)*calo_inner_r - covered_z; 

    // Calculate how many tubes there are in the back face
    double tower_outer_r = rad_distance + back_shift + 2*z_half; 
    double tower_max_backface_height = tower_outer_r * tan_theta;
    int num_back_rows = num_front_rows + DDDRCaloTubes::fast_floor((tower_max_backface_height-(tower_frontface_height-overlap))/V);

    int num_rows = num_back_rows;


    constructor.assemble_tower(tower_volume);
    
    
    double y_shift = std::cos(covered_theta)*(z_half+back_shift); // Where the tower reference point y coordinate is for this tower (not regarding inner calo radius)
    double z_shift = std::sin(covered_theta)*(z_half+back_shift); // How much the tower volume reference points moves in z wrt to previous tower

    double tower_x = 0*cm;
    double tower_y = y_shift + calo_inner_r;
    double tower_z = -(z_shift + covered_z);


    tower_position = Position(tower_x, tower_y, tower_z);
    delta_theta = (last_tower) ? barrel_endcap_angle : this_tower_theta ; // For last tower just return barrel_endcap_angle to make sure we cross boundary

    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "this_tower_theta = " << this_tower_theta/deg << std::endl;
    // std::cout << "this_tower_z     = " << this_tower_z     << std::endl;
    // std::cout << "tower_max_z      = " << tower_max_z      << std::endl;
    std::cout << "y_shift          = " << y_shift/mm          << std::endl;
    std::cout << "z_shift          = " << z_shift/mm          << std::endl;
    std::cout << "covered_theta    = " << covered_theta/deg    << std::endl;
    std::cout << "covered_z        = " << covered_z /mm       << std::endl;
    std::cout << "tower_x          = " << tower_x/mm         << std::endl;
    std::cout << "tower_y          = " << tower_y/mm         << std::endl;
    std::cout << "tower_z          = " << tower_z/mm         << std::endl;
    std::cout << "num_front_rows   = " << num_front_rows   << std::endl;
    std::cout << "num_front_cols   = " << num_front_cols   << std::endl;
    // std::cout << "phi_back_shift   = " << phi_back_shift/mm<< std::endl;
    std::cout << "num_back_cols    = " << num_back_cols    << std::endl;
    std::cout << "num_back_rows    = " << num_back_rows    << std::endl;
    // std::cout << "num_bad_rows     = " << num_bad_rows     << std::endl; 
    // std::cout << "weightA          = " << tower_volume->WeightA() << std::endl;
    // std::cout << "weight           = " << tower_volume->Weight() << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;



    return tower_volume;
}

void place_barrel_tower(Volume& calorimeter_volume,
                 Assembly& tower_volume,
                 unsigned int stave, 
                 unsigned int layer,
                 unsigned int tower_id,
                 Position tower_position,
                 double covered_theta,
                 double phi)
{

    double tower_x = std::sin(phi)*tower_position.Y();
    double tower_y = std::cos(phi)*tower_position.Y();
    double tower_z = tower_position.Z();

    /* std::cout<<"tower_id = "<<tower_id<<std::endl;
    std::cout<<"stave    = "<<stave<<std::endl;
    std::cout<<"layer    = "<<layer<<std::endl;
    std::cout<<"phi      = "<<phi/deg<<std::endl;
    std::cout<<"cov_theta= "<<covered_theta/deg<<std::endl;
    std::cout<<"deg      = "<<thetaDegrees<<std::endl;
    std::cout<<"dec      = "<<thetaDecimal<<std::endl;
    std::cout<<"----------------------------------------" << std::endl; */

    // Backward barrel region
    Transform3D tower_bwd_tr(RotationZYX(0, phi, -90*deg-covered_theta), Position(tower_x, tower_y, tower_z));
    PlacedVolume tower_bwd_placed = calorimeter_volume.placeVolume(tower_volume, -tower_id, tower_bwd_tr);
    tower_bwd_placed.addPhysVolID("stave", -stave).addPhysVolID("layer", -layer);
    
    // Forward barrel region
    Transform3D tower_fwd_tr(RotationZYX(180*deg, phi, -90*deg+covered_theta), Position(tower_x, tower_y, -tower_z));
    PlacedVolume tower_fwd_placed = calorimeter_volume.placeVolume(tower_volume, tower_id, tower_fwd_tr);
    tower_fwd_placed.addPhysVolID("stave", stave).addPhysVolID("layer", layer);

}

// int fast_floor(double x)
// {
//     return (int) x - (x < (int) x);
// }

// int fast_ceil(double x)
// {
//     return (int) x + (x > (int) x);
// }

// bool check_for_integer(double x)
// {
//     return (std::abs(x - std::round(x)) < 1e-10);
// }


DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)