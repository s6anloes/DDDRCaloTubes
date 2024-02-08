#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

double construct_tower(Detector& description,
                     xml_h& entities,
                     SensitiveDetector& sens,
                     Volume& calorimeter_volume,
                     const unsigned int& module_id,
                     double& covered_theta,
                     double& covered_phi);

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
    double      calo_inner_half_length = x_dim.z_length();
    double      tower_length           = 2*x_dim.zhalf();

    double barrel_endcap_angle = std::atan2(calo_inner_r, calo_inner_half_length);
    // Tube        calorimeter_solid(calo_inner_r, calo_inner_r+tower_length, calo_inner_half_length+tower_length); // Per design a "square" cylinder
    Tube        calorimeter_solid(0, calo_inner_r+tower_length, 0); // Per design a "square" cylinder
    std::string calorimeter_name = "calorimeter";
    Volume      calorimeter_volume(calorimeter_name, calorimeter_solid, air);
    calorimeter_volume.setVisAttributes(description, "MyVis");

    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);

    double covered_phi   = 0*deg;
    unsigned int tower_id = 0;

    for (double covered_theta=0*deg; covered_theta<barrel_endcap_angle; ) {
        double theta = 90*deg - covered_theta;
        double delta_theta = construct_tower(description, entities, sens, calorimeter_volume, tower_id, covered_theta, covered_phi);
        covered_theta += delta_theta;
        tower_id++;
    }

    Transform3D calorimeter_tr(RotationZYX(0, 0, 0), Position(0, 0, 0));
    PlacedVolume calorimeter_placed = mother_volume.placeVolume(calorimeter_volume, calorimeter_tr);
    calorimeter_placed.addPhysVolID("system", det_id);

    s_detElement.setPlacement(calorimeter_placed);

    return s_detElement;
}

double construct_tower(Detector& description,
                     xml_h& entities,
                     SensitiveDetector& sens,
                     Volume& calorimeter_volume,
                     const unsigned int& module_id,
                     double& covered_theta,
                     double& covered_phi) 
{
    xml_det_t   x_det       = entities;    
    std::string det_name    = x_det.nameStr();

    xml_dim_t   x_dim                  = x_det.dimensions();
    double      z_half                 = x_dim.zhalf();
    int         num_rows               = x_dim.number();
    int         num_cols               = x_dim.count();
    double      calo_inner_r           = x_dim.inner_radius();
    double      calo_inner_half_length = x_dim.z_length();
    double      tower_theta            = x_dim.deltatheta();
    double      tower_phi              = x_dim.deltaphi();


    double barrel_endcap_angle = std::atan2(calo_inner_r, calo_inner_half_length);

    
    Assembly      module_volume(det_name+"_module");

    // Get parameters for tube construction
    xml_comp_t  x_tube      = x_det.child(_Unicode(tube));

    xml_comp_t  x_capillary          = x_tube.child(_Unicode(capillary));
    Material    capillary_material   = description.material(x_capillary.materialStr()); 
    double      capillary_outer_r    = x_capillary.outer_r();

    xml_comp_t  x_scin_clad          = x_tube.child(_Unicode(scin_clad));
    Material    scin_clad_material   = description.material(x_scin_clad.materialStr()); 
    double      scin_clad_outer_r    = x_scin_clad.outer_r();

    xml_comp_t  x_scin_fibre         = x_tube.child(_Unicode(scin_fibre));
    Material    scin_fibre_material  = description.material(x_scin_fibre.materialStr()); 
    double      scin_fibre_outer_r   = x_scin_fibre.outer_r();

    xml_comp_t  x_cher_clad          = x_tube.child(_Unicode(cher_clad));
    Material    cher_clad_material   = description.material(x_cher_clad.materialStr()); 
    double      cher_clad_outer_r    = x_cher_clad.outer_r();

    xml_comp_t  x_cher_fibre         = x_tube.child(_Unicode(cher_fibre));
    Material    cher_fibre_material  = description.material(x_cher_fibre.materialStr()); 
    double      cher_fibre_outer_r   = x_cher_fibre.outer_r();


    // Constants used through the function
    const double D = 4.0*capillary_outer_r/sqrt(3.0); // Long diagonal of hexagaon with capillary_outer_r as inradius
    const double V = 3.0*D/4.0;                       // Vertical spacing for pointy top oriented tubes

    // Overlap between tubes cause by placing the tube in the gap of the previous row
    double overlap = 2*capillary_outer_r - V;


    // Calculate tower dimensions
    double covered_z = std::tan(covered_theta)*calo_inner_r;



    bool last_tower = (covered_theta+tower_theta>barrel_endcap_angle) ? true : false ;
    double tower_max_theta = (last_tower) ? barrel_endcap_angle : covered_theta+tower_theta;
    double tower_max_z = std::tan(tower_max_theta)*calo_inner_r - covered_z; // Max distance the front face of this tower covers in z (not regarding how many fibres actually fit)
    if (last_tower) tower_max_z += overlap/2; // to account for fine shift of the first tower
    double tower_max_frontface_height = std::cos(covered_theta)*tower_max_z; // Tower height (in theta direction) without regarding how many tubes actually fit

    if (tower_max_frontface_height < 2*capillary_outer_r)
    {
        throw std::runtime_error("Can't construct tower with given tower_theta and calo_inner_radius");
    }
    
    // Calculate how many tubes fit at the front face for the given tower theta coverage.
    // This number will serve as the new covered theta since it is important to not have any gaps in the front face
    int num_front_rows = 1 + floor((tower_max_frontface_height-2*capillary_outer_r) / V);
    if (num_front_rows&1 && !last_tower) num_front_rows++; // Make sure that front face ends on row with offset (i.e. even number of rows)

    double tower_frontface_height;
    double back_shift;
    double rad_distance;
    double this_tower_theta;
    // GOTO label
REDUCTION_REQUIRED:

    tower_frontface_height = 2*capillary_outer_r + (num_front_rows-1)*V;
    
    // Distance by which straight edge of this tower is shifted backwards to ensure inner radius of calorimeter
    back_shift = std::tan(covered_theta)*tower_frontface_height;                      

    // Radial distance to exceed 2.5m inner radius of calorimeter for this tower
    rad_distance = calo_inner_r/std::cos(covered_theta);

    this_tower_theta = std::atan2(tower_frontface_height-overlap, rad_distance+back_shift);
    // This catches cases where covered_theta+tower_theta does not exceed barrel_endcap_angle, but covered_theta+this_tower_theta does
    // since it can happen that this_tower_theta is larger than tower_theta
    if (covered_theta+this_tower_theta > barrel_endcap_angle) {
        num_front_rows--;
        goto REDUCTION_REQUIRED;
    }

    double tan_theta = std::tan(this_tower_theta);
    // double missing_theta = tower_theta - this_tower_theta;

    // Distance the front face of this tower covers in z
    double this_tower_z = std::tan(covered_theta+this_tower_theta)*calo_inner_r - covered_z; 

    // Calculate how many tubes there are in the back face
    double tower_outer_r = rad_distance + back_shift + 2*z_half; 
    double tower_max_backface_height = tower_outer_r * tan_theta;
    int num_back_rows = num_front_rows + floor((tower_max_backface_height-(tower_frontface_height-overlap))/V);

    num_rows = num_back_rows;

    // Placement and shortening of tube depends on whether the rows have an offset or not
    // The below method of calculating assumes that the final row at the front face before the staggering begins has an offset 
    // (even number of rows in front face) such that the next tower can be placed smoothly into the last row
    double tube_shortening_odd_stagger  = 2*capillary_outer_r/tan_theta;
    double tube_shortening_even_stagger = 2*V/tan_theta - tube_shortening_odd_stagger;


    for (int row=0; row<num_rows; row++)
    {
        for (int col=0; col<num_cols; col++)
        {
            // Definition of Air Volume which is created newly in each loop such that daughters can be 
            // repeatedly placed with identical copy numbers
            // TODO: Check what objects can be moved outside of loop (string _name, Tube _solid, etc.)

            // Configuration for placing the tube
            double offset = (row & 1) ? capillary_outer_r : 0.0*mm;
            double x = col*2*capillary_outer_r + offset;
            double y = row*V;                       // Vertical spacing for hexagonal grid (pointy-top)

            int staggered_row = (row >= num_front_rows) ? row-num_front_rows+1 : 0 ;
            int even_staggers = (staggered_row & 1) ? (staggered_row-1)/2 : staggered_row/2;
            int odd_staggers  = (staggered_row & 1) ? (staggered_row+1)/2 : staggered_row/2;

            double z = (even_staggers*tube_shortening_even_stagger + odd_staggers*tube_shortening_odd_stagger)/2;

            // Adding tube radius to x and y such that volume reference point is at the lower "corner" of the tower instead of in the middle of a fibre
            auto position = Position(x+capillary_outer_r, y+capillary_outer_r, z);

            // Axial coordinate conversion following https://www.redblobgames.com/grids/hexagons/#conversions-offset
            // Slighty changed to fit my q and r directions
            unsigned short int q = col + (row - (row&1))/2;
            unsigned short int r = row;
            //std::cout<<"(row, col) -> (r, q) : (" <<row<<", "<<col<<") -> (" << r<<", " <<q<<")" <<std::endl;

            // TubeID composed of q in first 16 bits, r in last 16 bits
            unsigned int tube_id = (q << 16) | r;
            
            //std::cout<<"(row, col) -> (r, q) -> (tubeID) : (" <<row<<", "<<col<<") -> (" <<r<<", " <<q<<") -> (" << tube_id << ")" <<std::endl; 

            double tube_half_length = z_half - z;
            // Capillary tube
            Tube        capillary_solid(0.0*mm, capillary_outer_r, tube_half_length);
            std::string capillary_name = "capillary";
            Volume      capillary_volume(capillary_name, capillary_solid, capillary_material);
            if (x_capillary.isSensitive()) capillary_volume.setSensitiveDetector(sens);
            capillary_volume.setVisAttributes(description, x_capillary.visStr()); 

            if (row & 1) // Cherenkov row
            {
                // Cherenkov cladding
                Tube        cher_clad_solid(0.0*mm, cher_clad_outer_r, tube_half_length);
                std::string cher_clad_name = "cher_clad";
                Volume      cher_clad_volume(cher_clad_name, cher_clad_solid, cher_clad_material);
                if (x_cher_clad.isSensitive()) cher_clad_volume.setSensitiveDetector(sens);
                PlacedVolume cher_clad_placed = capillary_volume.placeVolume(cher_clad_volume, tube_id);
                cher_clad_volume.setVisAttributes(description, x_cher_clad.visStr());
                cher_clad_placed.addPhysVolID("clad", 1);

                // Chrerenkov fibre
                Tube        cher_fibre_solid(0.0*mm, cher_fibre_outer_r, tube_half_length);
                std::string cher_fibre_name = "cher_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
                Volume      cher_fibre_volume(cher_fibre_name, cher_fibre_solid, cher_fibre_material);
                if (x_cher_fibre.isSensitive()) cher_fibre_volume.setSensitiveDetector(sens);
                PlacedVolume    cher_fibre_placed = cher_clad_volume.placeVolume(cher_fibre_volume, tube_id);
                cher_fibre_volume.setVisAttributes(description, x_cher_fibre.visStr());
                cher_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 1);
            }
            else // Scintillation row
            {
                // Scintillation cladding
                Tube        scin_clad_solid(0.0*mm, scin_clad_outer_r, tube_half_length);
                std::string scin_clad_name = "scin_clad";
                Volume      scin_clad_volume(scin_clad_name, scin_clad_solid, scin_clad_material);
                if (x_scin_clad.isSensitive()) scin_clad_volume.setSensitiveDetector(sens);
                PlacedVolume scin_clad_placed = capillary_volume.placeVolume(scin_clad_volume, tube_id);
                scin_clad_volume.setVisAttributes(description, x_scin_clad.visStr());
                scin_clad_placed.addPhysVolID("clad", 1);

                // Scintillation fibre
                Tube        scin_fibre_solid(0.0*mm, scin_fibre_outer_r, tube_half_length);
                std::string scin_fibre_name = "scin_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
                Volume      scin_fibre_volume(scin_fibre_name, scin_fibre_solid, scin_fibre_material);
                if (x_scin_fibre.isSensitive()) scin_fibre_volume.setSensitiveDetector(sens);
                PlacedVolume    scin_fibre_placed = scin_clad_volume.placeVolume(scin_fibre_volume, tube_id);
                scin_fibre_volume.setVisAttributes(description, x_scin_fibre.visStr());
                scin_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 0);
            }
            //auto tube_to_be_placed = (row & 1) ? &cher_air_volume : &scin_air_volume;

            // auto sensitive_to_be_placed = (row & 1) ? Volume(cher_tube_volume) : Volume(scin_tube_volume);

            // sensitive_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 1);


            PlacedVolume    tube_placed = module_volume.placeVolume(capillary_volume, tube_id, position);
            
            tube_placed.addPhysVolID("fibre", 0).addPhysVolID("q", q).addPhysVolID("r", r);

            
        }
    }
    
    double y_shift = std::cos(covered_theta)*(z_half+back_shift); // Where the tower reference point y coordinate is for this tower (not regarding inner calo radius)
    double z_shift = std::sin(covered_theta)*(z_half+back_shift); // How much the tower volume reference points moves in z wrt to previous tower

    double module_x = 0*cm;
    double module_y = y_shift + calo_inner_r;
    double module_z = -(z_shift + covered_z - overlap/2);
    // module_z = 0*cm;

    Transform3D module_tr(RotationZYX(0, 0, -90*deg-covered_theta), Position(module_x, module_y, module_z));
    // // Transform3D module_tr(RotationZYX(0, 0, 0), Position(0, 0, 0));
    PlacedVolume module_placed = calorimeter_volume.placeVolume(module_volume, module_tr);
    module_placed.addPhysVolID("module", module_id);

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "this_tower_theta = " << this_tower_theta/deg << std::endl;
    std::cout << "this_tower_z     = " << this_tower_z     << std::endl;
    std::cout << std::endl;
    std::cout << "tower_max_z      = " << tower_max_z      << std::endl;
    std::cout << "y_shift          = " << y_shift          << std::endl;
    std::cout << "z_shift          = " << z_shift          << std::endl;
    std::cout << "covered_theta    = " << covered_theta/deg    << std::endl;
    std::cout << "module_x         = " << module_x         << std::endl;
    std::cout << "module_y         = " << module_y         << std::endl;
    std::cout << "module_z         = " << module_z         << std::endl;
    std::cout << "num_front_rows   = " << num_front_rows   << std::endl;
    std::cout << "tower_id         = " << module_id        << std::endl;

    double return_value = (last_tower) ? barrel_endcap_angle : this_tower_theta ; // For last tower just return barrel_endcap_angle to make sure we cross boundary

    return return_value;
}

DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)