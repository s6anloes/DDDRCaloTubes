#include "DRconstructor.h"

#include <TMatrixD.h>

using namespace dd4hep;
using namespace DDDRCaloTubes;

DDDRCaloTubes::DRconstructor::DRconstructor(Detector* description,
                                            xml_h& entities,
                                            SensitiveDetector* sens):
                                            m_entities(entities)
{
    m_description = description;
    m_sens = sens;

    xml_dim_t x_dim = ((xml_det_t) entities).dimensions();

    // Calorimeter parameters
    m_calo_inner_r      = x_dim.inner_radius();
    m_calo_outer_r      = x_dim.outer_radius();
    m_calo_inner_half_z = x_dim.z_length();

    
    // Trap parameters
    xml_comp_t  x_trap = entities.child(_Unicode(trap));
    xml_comp_t  x_trap_support = x_trap.child(_Unicode(support));
    m_trap_wall_thickness_front = x_trap_support.depth();
    m_trap_wall_thickness_sides = x_trap_support.width();
    m_trap_wall_thickness_back  = x_trap_support.z2();


    // Tube parameters
    xml_comp_t  x_tube = entities.child(_Unicode(tube));

    xml_comp_t  x_capillary = x_tube.child(_Unicode(capillary));
    m_capillary_material    = m_description->material(x_capillary.materialStr()); 
    m_capillary_outer_r     = x_capillary.outer_r();
    m_capillary_visString   = x_capillary.visStr();
    m_capillary_isSensitive = x_capillary.isSensitive();

    xml_comp_t  x_scin_clad = x_tube.child(_Unicode(scin_clad));
    m_scin_clad_material    = m_description->material(x_scin_clad.materialStr()); 
    m_scin_clad_outer_r     = x_scin_clad.outer_r();
    m_scin_clad_visString   = x_scin_clad.visStr();
    m_scin_clad_isSensitive = x_scin_clad.isSensitive();

    xml_comp_t  x_scin_core = x_tube.child(_Unicode(scin_core));
    m_scin_core_material    = m_description->material(x_scin_core.materialStr()); 
    m_scin_core_outer_r     = x_scin_core.outer_r();
    m_scin_core_visString   = x_scin_core.visStr();
    m_scin_core_isSensitive = x_scin_core.isSensitive();

    xml_comp_t  x_cher_clad = x_tube.child(_Unicode(cher_clad));
    m_cher_clad_material    = m_description->material(x_cher_clad.materialStr()); 
    m_cher_clad_outer_r     = x_cher_clad.outer_r();
    m_cher_clad_visString   = x_cher_clad.visStr();
    m_cher_clad_isSensitive = x_cher_clad.isSensitive();

    xml_comp_t  x_cher_core = x_tube.child(_Unicode(cher_core));
    m_cher_core_material    = m_description->material(x_cher_core.materialStr()); 
    m_cher_core_outer_r     = x_cher_core.outer_r();
    m_cher_core_visString   = x_cher_core.visStr();
    m_cher_core_isSensitive = x_cher_core.isSensitive();


    // Tower parameters
    m_tower_theta = x_dim.deltatheta();
    m_tower_phi   = x_dim.deltaphi();


    // Construction parameters
    m_covered_theta = 0.0*deg;
    m_back_shift = 0.0*mm;
    // m_tower_volume = nullptr;

    this->calculate_tower_parameters();
    this->calculate_phi_parameters();

    m_num_rows = 0;
    m_num_front_rows = 0;
    m_this_tower_theta = 0.0*deg;
    m_tower_tan_theta = 0.0;

    m_air = m_description->material("Air");

    xml_comp_t  x_air = x_trap.child(_Unicode(air));
    m_air    = m_description->material(x_air.materialStr()); 
    m_air_visString     = x_air.visStr();

}

// Function to calculate all tower parameters which are derived from user given values
void DDDRCaloTubes::DRconstructor::calculate_tower_parameters()
{
    // Angle where endcap and barrel meet
    m_barrel_endcap_angle = std::atan2(m_calo_inner_half_z, m_calo_inner_r);

    m_capillary_diameter = 2*m_capillary_outer_r;

    // Constants used through the function
    double D = 4.0*m_capillary_outer_r/sqrt(3.0); // Long diagonal of hexagaon with capillary_outer_r as inradius
    m_V = 3.0*D/4.0;                        // Vertical spacing for pointy top oriented tubes
    m_overlap = m_capillary_diameter - m_V; // Overlap between tubes


    m_tower_tan_phi = std::tan(m_tower_phi); // Needed several times, calculate once here
    m_trap_half_length  = (m_calo_outer_r - m_calo_inner_r)/2; // Trapezoid half length
    
    m_tower_half_length = m_trap_half_length - m_trap_wall_thickness_front/2.0 - m_trap_wall_thickness_back/2.0; // Tower half length

}

// Function to calculate tower parameters specifically for phi direction
void DDDRCaloTubes::DRconstructor::calculate_phi_parameters()
{
    double num_phi_towers_d = 360.0*deg/m_tower_phi;
    // Check if num_phi_towers is a whole number
    if (check_for_integer(num_phi_towers_d)) m_num_phi_towers = static_cast<unsigned int>(std::round(num_phi_towers_d));
    else throw std::runtime_error("Not an integer number of towers in phi direction");

    
}

void DDDRCaloTubes::DRconstructor::calculate_theta_parameters()
{

    

    // Approximate tower position in z direction
    // In case we're already at endcap region, refers to approximate height at which tower is placed
    double covered_z;
    if (m_covered_theta < m_barrel_endcap_angle) {
        covered_z = std::tan(m_covered_theta)*m_calo_inner_r;
    } else {
        covered_z = std::tan(90*deg-m_covered_theta)*m_calo_inner_half_z;
    }




    double tower_max_theta = m_covered_theta+m_tower_theta;
    double tower_max_z = std::tan(tower_max_theta)*m_calo_inner_r - covered_z; // Max distance the front face of this tower covers in z (not regarding how many fibres actually fit)
    double tower_max_frontface_y = std::cos(m_covered_theta)*tower_max_z; // Tower height (in theta direction) without regarding how many tubes actually fit

    if (tower_max_frontface_y < m_capillary_diameter)
    {
        throw std::runtime_error("Can't construct tower with given tower_theta and calo_inner_radius");
    }

    m_tower_tan_theta = std::tan(m_tower_theta);
    m_trap_frontface_y = tower_max_frontface_y;
    m_trap_backface_y  = tower_max_frontface_y + 2*m_trap_half_length*m_tower_tan_theta;

    m_tower_frontface_y = m_trap_frontface_y - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));
    m_tower_backface_y  = m_trap_backface_y  - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));

    // Calculate how many tubes fit at the front face for the given tower theta coverage.
    // This number will serve as the new covered theta since it is important to not have any gaps in the front face
    m_num_front_rows = 1 + fast_floor((tower_max_frontface_y-m_capillary_diameter) / m_V);
    if (m_num_front_rows&1 ) m_num_front_rows++; // Make sure that front face ends on row with offset (i.e. even number of rows)

    double tower_frontface_y = m_capillary_diameter + (m_num_front_rows-1)*m_V;
    
    // Distance by which straight edge of this tower is shifted backwards to ensure inner radius of calorimeter
    m_back_shift = std::tan(m_covered_theta)*m_trap_frontface_y;                      

    // Radial distance to exceed 2.5m inner radius of calorimeter for this tower
    double rad_distance = m_calo_inner_r/std::cos(m_covered_theta);

    m_this_tower_theta = std::atan2(tower_frontface_y-m_overlap, rad_distance+m_back_shift);
    // m_tower_tan_theta = std::tan(m_tower_theta);
    // double missing_theta = m_tower_theta - this_tower_theta;

    // Calculate how many tubes there are in the back face
    double tower_outer_r = rad_distance + m_back_shift + 2*m_tower_half_length; 
    double tower_max_backface_y = tower_outer_r * m_tower_tan_theta;
    int num_back_rows = m_num_front_rows + fast_floor((tower_max_backface_y-(tower_frontface_y-m_overlap))/m_V);

    m_num_rows = num_back_rows;


    // Frontface width of trapezoid support
    // m_trap_frontface_x = m_calo_inner_r*m_tower_tan_phi;
    m_trap_frontface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*m_back_shift)*2*std::tan(m_tower_phi/2.0);
    m_trap_frontface_thetaangleedge_x = m_calo_inner_r*2*std::tan(m_tower_phi/2.0);

    // Right Trapezoids 
    // m_trap_frontface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*m_back_shift)*std::tan(m_tower_phi);
    // m_trap_frontface_thetaangleedge_x = m_calo_inner_r*std::tan(m_tower_phi);

    // m_tower_frontface_x = m_trap_frontface_x + m_tower_tan_phi*m_trap_wall_thickness_front - (1+1/std::cos(m_tower_phi))*m_trap_wall_thickness_sides; 

    // double num_front_cols_d = m_tower_frontface_x/m_capillary_diameter;
    // m_num_front_cols = fast_floor(num_front_cols_d);
    m_num_cols = 1;
    // double remainder = num_front_cols_d - m_num_front_cols;
    // if (remainder < 0.5)
    // tower_frontface_x = m_num_front_cols*m_capillary_diameter;


    // m_effective_inner_r = m_calo_inner_r; // Shifting of tower, change in calo radius
    double tower_backface_phi_increase_rightangleedge = 2*m_trap_half_length*std::cos(m_covered_theta) * 2*std::tan(m_tower_phi/2.0); // by how much the backface wides for both tower and trap
    double tower_backface_phi_increase_thetaangleedge = 2*m_trap_half_length*std::cos(m_covered_theta+m_tower_theta) * 2*std::tan(m_tower_phi/2.0); // by how much the backface wides for both tower and trap
    // m_trap_backface_x  = m_trap_frontface_x + tower_backface_phi_increase;
    m_trap_backface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*(m_back_shift+2*m_trap_half_length))*2*std::tan(m_tower_phi/2.0);
    m_trap_backface_thetaangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta+m_tower_theta)*2*m_trap_half_length/std::cos(m_tower_theta))*2*std::tan(m_tower_phi/2.0);


    // Right Trapezoids
    // m_trap_backface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*(m_back_shift+2*m_trap_half_length))*std::tan(m_tower_phi);
    // m_trap_backface_thetaangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta+m_tower_theta)*2*m_trap_half_length/std::cos(m_tower_theta))*std::tan(m_tower_phi);


    m_tower_frontface_x = m_trap_frontface_rightangleedge_x - (1+1/std::cos(m_tower_phi))*m_trap_wall_thickness_sides;
    m_tower_backface_x  = m_trap_backface_rightangleedge_x  - (1+1/std::cos(m_tower_phi))*m_trap_wall_thickness_sides;

    // Calculate how many tubes there are in the back face
    int num_back_cols = fast_floor(m_tower_backface_x/m_capillary_diameter);

    m_num_cols = num_back_cols;

}

void DDDRCaloTubes::DRconstructor::assemble_tower(Volume& tower_air_volume)
{

    // Placement and shortening of tube depends on whether the rows have an offset or not
    // The below method of calculating assumes that the final row at the front face before the staggering begins has an offset 
    // (even number of rows in front face) such that the next tower can be placed smoothly into the last row
    // double tube_shortening_odd_stagger  = m_capillary_diameter/m_tower_tan_theta;
    // double tube_shortening_even_stagger = 2*m_V/m_tower_tan_theta - tube_shortening_odd_stagger;
    double tower_centre_r = m_tower_half_length/std::cos(m_tower_polar_angle);
    double tower_centre_half_x = tower_centre_r*std::sin(m_tower_polar_angle)*std::cos(m_tower_azimuthal_angle) + m_tower_frontface_x/2.0;
    double tower_centre_half_y = tower_centre_r*std::sin(m_tower_polar_angle)*std::sin(m_tower_azimuthal_angle) + m_tower_frontface_y/2.0;
    double tower_centre_half_z = m_tower_half_length;


    int num_bad_rows = 0;
    int num_bad_cols = 0;

    double covered_tower_y = m_capillary_diameter;

    int row = 0;
    // for (unsigned int row=0; row<m_num_rows; row++)
    while(covered_tower_y < m_tower_backface_y)
    {
        // std::cout << "covered_tower_y = " << covered_tower_y/mm << std::endl;
        // std::cout << "m_tower_backface_y = " << m_tower_backface_y/mm << std::endl;
        // int staggered_row = (row >= m_num_front_rows) ? row-m_num_front_rows+1 : 0 ;
        // int even_staggers = (staggered_row & 1) ? (staggered_row-1)/2 : staggered_row/2;
        // int odd_staggers  = (staggered_row & 1) ? (staggered_row+1)/2 : staggered_row/2;

        // double row_staggered_z = (even_staggers*tube_shortening_even_stagger + odd_staggers*tube_shortening_odd_stagger)/2;
        double row_staggered_z = 0.0*mm;
        if (covered_tower_y > m_tower_frontface_y) row_staggered_z = (covered_tower_y-m_tower_frontface_y)/m_tower_tan_theta/2.0;

        // In row staggering it can happen in rare cases that the staggering becomes too large in final row
        // Has to do with fact that caclulation of num_rows does not take different staggering into account
        if (row_staggered_z > m_tower_half_length) {
            num_bad_rows++;
            std::cout << "Bad row by " << (row_staggered_z-m_tower_half_length)/um << " um" << std::endl;
            covered_tower_y += m_V;
            continue;
        }


        double offset = (row & 1) ? m_capillary_outer_r : 0.0*mm;

        double covered_tower_x = m_capillary_diameter+offset;
        int col = 0;
        // for (unsigned int col=0; col<m_num_cols; col++)
        while(covered_tower_x < m_tower_backface_x)
        {
            // TODO: Check what objects can be moved outside of loop (string _name, Tube _solid, etc.)

            // int staggered_col = (col >= m_num_front_cols) ? col-m_num_front_cols+1 : 0 ;
            double col_staggered_z = 0.0*mm;
            if (covered_tower_x > m_tower_frontface_x) col_staggered_z = (covered_tower_x-m_tower_frontface_x)/m_tower_tan_phi/2.0;
            // col_staggered_z = staggered_col*m_capillary_diameter/m_tower_tan_phi/2;
            if (col_staggered_z > m_tower_half_length) {
                num_bad_cols++;
                covered_tower_x += m_capillary_diameter;
                continue;
            }


            // Configuration for placing the tube
            double x = col*m_capillary_diameter + offset + m_capillary_outer_r;
            double y = row*m_V + m_capillary_outer_r;                       // Vertical spacing for hexagonal grid (pointy-top)

            double z = (row_staggered_z > col_staggered_z) ? row_staggered_z : col_staggered_z ;
            // std::cout << "row_staggered_z = " << row_staggered_z/mm << std::endl;
            // std::cout << "col_staggered_z = " << col_staggered_z/mm << std::endl;
            // std::cout << "m_tower_half_length = " << m_tower_half_length/mm << std::endl;

            // Adding tube radius to x and y such that volume reference point is at the lower "corner" of the tower instead of in the middle of a fibre
            // auto position = Position(x+m_capillary_outer_r, y+m_capillary_outer_r, z);
            auto position = Position(x-tower_centre_half_x, y-tower_centre_half_y, z);

            // Offset coordinates following https://www.redblobgames.com/grids/hexagons/#coordinates-offset
            unsigned short int q = col;
            unsigned short int r = row;

            // TubeID composed of q in first 16 bits, r in last 16 bits
            unsigned int tube_id = (q << 16) | r;
            // unsigned int tube_id = q + r*m_num_cols;
            
            // std::cout<<"(row, col) -> (r, q) -> (tubeID) : (" <<row<<", "<<col<<") -> (" <<r<<", " <<q<<") -> (" << tube_id << ")" <<std::endl; 

            double tube_half_length = m_tower_half_length - z;
            // std::cout << "tube_half_length = " << tube_half_length/mm << std::endl;

            // Capillary tube
            Tube        capillary_solid(0.0*mm, m_capillary_outer_r, tube_half_length);
            Volume      capillary_volume("capillary", capillary_solid, m_capillary_material);
            if (m_capillary_isSensitive) capillary_volume.setSensitiveDetector(*m_sens);
            capillary_volume.setVisAttributes(*m_description, m_capillary_visString); 

            if (row & 1) // Cherenkov row
            {
                // Cherenkov cladding
                Tube        cher_clad_solid(0.0*mm, m_cher_clad_outer_r, tube_half_length);
                Volume      cher_clad_volume("cher_clad", cher_clad_solid, m_cher_clad_material);
                if (m_cher_clad_isSensitive) cher_clad_volume.setSensitiveDetector(*m_sens);
                PlacedVolume cher_clad_placed = capillary_volume.placeVolume(cher_clad_volume, tube_id);
                cher_clad_volume.setVisAttributes(*m_description, m_cher_clad_visString);
                cher_clad_placed.addPhysVolID("clad", 1).addPhysVolID("cherenkov", 1);

                // Chrerenkov fibre
                Tube        cher_fibre_solid(0.0*mm, m_cher_core_outer_r, tube_half_length);
                Volume      cher_fibre_volume("cher_fibre", cher_fibre_solid, m_cher_core_material);
                if (m_cher_core_isSensitive) cher_fibre_volume.setSensitiveDetector(*m_sens);
                PlacedVolume    cher_fibre_placed = cher_clad_volume.placeVolume(cher_fibre_volume, tube_id);
                cher_fibre_volume.setVisAttributes(*m_description, m_cher_core_visString);
                cher_fibre_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);
            }
            else // Scintillation row
            {
                // Scintillation cladding
                Tube        scin_clad_solid(0.0*mm, m_scin_clad_outer_r, tube_half_length);
                Volume      scin_clad_volume("scin_clad", scin_clad_solid, m_scin_clad_material);
                if (m_scin_clad_isSensitive) scin_clad_volume.setSensitiveDetector(*m_sens);
                PlacedVolume scin_clad_placed = capillary_volume.placeVolume(scin_clad_volume, tube_id);
                scin_clad_volume.setVisAttributes(*m_description, m_scin_clad_visString);
                scin_clad_placed.addPhysVolID("clad", 1).addPhysVolID("cherenkov", 0);

                // Scintillation fibre
                Tube        scin_fibre_solid(0.0*mm, m_scin_core_outer_r, tube_half_length);
                Volume      scin_fibre_volume("scin_fibre", scin_fibre_solid, m_scin_core_material);
                if (m_scin_core_isSensitive) scin_fibre_volume.setSensitiveDetector(*m_sens);
                PlacedVolume    scin_fibre_placed = scin_clad_volume.placeVolume(scin_fibre_volume, tube_id);
                scin_fibre_volume.setVisAttributes(*m_description, m_scin_core_visString);
                scin_fibre_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);
            }

            PlacedVolume    tube_placed = tower_air_volume.placeVolume(capillary_volume, tube_id, position);
            tube_placed.addPhysVolID("air",0).addPhysVolID("col", col).addPhysVolID("row", row);
            // tube_placed.addPhysVolID("clad", 0).addPhysVolID("core", 0).addPhysVolID("q", q).addPhysVolID("r", r);

            covered_tower_x += m_capillary_diameter;
            col++;
            // std::cout << "col = " << col << std::endl;
            
        }

        covered_tower_y += m_V;
        row++;
        // std::cout << "row = " << row << std::endl;
    }

    // std::cout << "finishito" << std::endl;
    int check_num_back_rows = 1 + fast_floor((m_tower_backface_y-m_capillary_diameter) / m_V);
    int check_num_back_cols = fast_floor(m_tower_backface_x/m_capillary_diameter);

    // std::cout << "check_num_back_rows = " << check_num_back_rows << std::endl;
    // std::cout << "check_num_back_cols = " << check_num_back_cols << std::endl;
    // std::cout << "num_bad_rows = " << num_bad_rows << std::endl;
    // std::cout << "num_bad_cols = " << num_bad_cols << std::endl;

}


// Function to calculate the position of the tower in stave (for fixed phi = 0*deg)
void DDDRCaloTubes::DRconstructor::calculate_tower_position(double phi)
{
    double covered_z = std::tan(m_covered_theta)*m_calo_inner_r;
    // double x_shift = std::cos(m_covered_theta)*(m_tower_half_length+m_back_shift); // Where the tower reference point y coordinate is for this tower (not regarding inner calo radius)
    // double z_shift = std::sin(m_covered_theta)*(m_tower_half_length+m_back_shift); // How much the tower volume reference points moves in z wrt to previous tower

    // double tower_x = 0*cm;
    // double tower_y = y_shift + m_effective_inner_r;
    // double tower_z = -(z_shift + covered_z);
    double trap_centre_r = m_trap_half_length/std::cos(m_trap_polar_angle);

    // double trap_centre_half_x = m_trap_frontface_x/2.0 + std::tan(m_trap_azimuthal_angle)*m_tower_half_length;
    // // double trap_centre_half_y = m_trap_frontface_y/2.0 + std::tan(m_trap_polar_angle)*m_tower_half_length;
    // double trap_centre_half_y = m_trap_backface_y/2.0 - m_trap_frontface_y/2.0;

    double trap_centre_half_x = trap_centre_r*std::sin(m_trap_polar_angle)*std::cos(m_trap_azimuthal_angle) + m_trap_frontface_rightangleedge_x/2.0;
    double trap_centre_half_y = trap_centre_r*std::sin(m_trap_polar_angle)*std::sin(m_trap_azimuthal_angle) + m_trap_frontface_y/2.0;
    double trap_centre_half_z = m_trap_half_length;
    double trap_rad_centre = m_calo_inner_r/std::cos(m_covered_theta) + m_back_shift + m_trap_half_length;

    double stave_x = std::cos(m_covered_theta)*trap_rad_centre - std::sin(m_covered_theta)*trap_centre_half_y;
    double stave_y = trap_centre_half_x*0;
    double stave_z = std::sin(m_covered_theta)*trap_rad_centre + std::cos(m_covered_theta)*trap_centre_half_y;

    
    double tower_x = std::cos(phi)*stave_x - std::sin(phi)*stave_y;
    double tower_y = std::sin(phi)*stave_x + std::cos(phi)*stave_y;
    double tower_z = stave_z;

    // std::cout<<"***********************************************************"<<std::endl;
    // std::cout<< "trap_cntr_half_y = " << trap_centre_half_y/mm << std::endl;
    // std::cout<< "trap_rad_centre = " << trap_rad_centre/mm << std::endl;
    // std::cout<< "stave_x = " << stave_x/mm << std::endl;
    // std::cout<< "stave_y = " << stave_y/mm << std::endl;
    // std::cout<< "stave_z = " << stave_z/mm << std::endl;
    // std::cout<< "tower_x = " << tower_x/mm << std::endl;
    // std::cout<< "tower_y = " << tower_y/mm << std::endl;
    // std::cout<< "tower_z = " << tower_z/mm << std::endl;
    // std::cout<<"***********************************************************"<<std::endl;

    m_tower_position = dd4hep::Position(tower_x, tower_y, tower_z);
    // m_tower_position = Position(trap_centre_half_x, trap_centre_half_y, trap_centre_half_z);
}


// Function to construct the trapezoidal supoprt structure for the tower in which fibres are placed
void DDDRCaloTubes::DRconstructor::construct_tower_trapezoid(Volume& trap_volume)
{

        
        // polar coordinate conversion
        double delta_x = 0*(m_trap_backface_rightangleedge_x - m_trap_frontface_rightangleedge_x)/2.0;
        double delta_y = (m_trap_backface_y - m_trap_frontface_y)/2.0;
        double delta_z = 2.0*m_trap_half_length;
        m_trap_polar_angle = std::acos(delta_z/std::sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z));
        // m_trap_azimuthal_angle = std::atan2(delta_y, delta_x);   
        m_trap_azimuthal_angle = 90.0*deg;

        Trap trap_solid("trap_solid", m_trap_half_length, m_trap_polar_angle, m_trap_azimuthal_angle, 
                                      m_trap_frontface_y/2.0, m_trap_frontface_rightangleedge_x/2.0, m_trap_frontface_thetaangleedge_x/2.0, 0.,
                                      m_trap_backface_y/2.0,  m_trap_backface_rightangleedge_x/2.0,  m_trap_backface_thetaangleedge_x/2.0,  0.);

        // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
        // std::cout<< "m_trap_half_length = " << m_trap_half_length/mm << std::endl;
        // std::cout<< "m_trap_frontface_rightangleedge_x = " << m_trap_frontface_rightangleedge_x/mm << std::endl;
        // std::cout<< "m_trap_frontface_y = " << m_trap_frontface_y/mm << std::endl;
        // std::cout<< "m_trap_backface_rightangleedge_x  = " << m_trap_backface_rightangleedge_x/mm << std::endl;
        // std::cout<< "m_trap_backface_y  = " << m_trap_backface_y/mm << std::endl;
        // std::cout<< "m_trap_azimuthal_angle = " << m_trap_azimuthal_angle/deg << std::endl;
        // std::cout<< "m_trap_polar_angle = " << m_trap_polar_angle/deg << std::endl;
        // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;

        // Air volume in which fibres are placed
        double delta_x_air = 0*(m_tower_backface_x - m_tower_frontface_x)/2.0;
        double delta_y_air = (m_tower_backface_y - m_tower_frontface_y)/2.0;
        double delta_z_air = 2.0*m_tower_half_length;
        m_tower_polar_angle = std::acos(delta_z_air/std::sqrt(delta_x_air*delta_x_air + delta_y_air*delta_y_air + delta_z_air*delta_z_air));
        // m_tower_azimuthal_angle = std::atan2(delta_y_air, delta_x_air);
        m_tower_azimuthal_angle = 90.0*deg;

        Trap tower_air_solid("tower_solid", m_tower_half_length, m_tower_polar_angle, m_tower_azimuthal_angle, 
                                      m_tower_frontface_y/2.0, m_tower_frontface_x/2.0, m_tower_frontface_x/2.0, 0.,
                                      m_tower_backface_y/2.0,  m_tower_backface_x/2.0,  m_tower_backface_x/2.0,  0.);
                                
        Position tower_air_pos = Position((1.0-1.0/std::cos(m_tower_phi))*m_trap_wall_thickness_sides,
                                          (1.0-1.0/std::cos(m_tower_theta))*m_trap_wall_thickness_sides,
                                          (m_trap_wall_thickness_front-m_trap_wall_thickness_back)/2.0);

        // SubtractionSolid solid = SubtractionSolid("trap_final", trap_solid, tower_air_solid, tower_air_pos);
        Volume tower_air_volume("tower_air_volume", tower_air_solid, m_air);
        tower_air_volume.setVisAttributes(*m_description, m_air_visString);

        trap_volume.setSolid(trap_solid);
        trap_volume.setVisAttributes(*m_description, m_capillary_visString);

        // PlacedVolume tower_air_placed = trap_volume.placeVolume(tower_air_volume, tower_air_pos);
        // tower_air_placed.addPhysVolID("air", 1);

        // this->assemble_tower(tower_air_volume);
        // std::cout<< "ASDEMBLED TOWER" << std::endl;
    
}


void DDDRCaloTubes::DRconstructor::construct_tower(Volume& trap_volume,
                                                   double& delta_theta)
{
    
    this->calculate_theta_parameters();

    this->construct_tower_trapezoid(trap_volume);


    delta_theta = m_this_tower_theta;

    /* std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "this_tower_theta = " << m_this_tower_theta/deg << std::endl;
    std::cout << "y_shift          = " << y_shift/mm          << std::endl;
    std::cout << "z_shift          = " << z_shift/mm          << std::endl;
    std::cout << "covered_theta    = " << m_covered_theta/deg    << std::endl;
    std::cout << "covered_z        = " << covered_z/mm        << std::endl;
    std::cout << "tower_x          = " << tower_x/mm         << std::endl;
    std::cout << "tower_y          = " << tower_y/mm         << std::endl;
    std::cout << "tower_z          = " << tower_z/mm         << std::endl;
    std::cout << "num_front_rows   = " << m_num_front_rows   << std::endl;
    std::cout << "num_front_cols   = " << m_num_front_cols   << std::endl;
    std::cout << "num_cols         = " << m_num_cols    << std::endl;
    std::cout << "num_rows         = " << m_num_rows    << std::endl;
    // std::cout << "weightA          = " << tower_volume->WeightA() << std::endl;
    // std::cout << "weight           = " << tower_volume->Weight() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl; */

}


void DDDRCaloTubes::DRconstructor::reset_tower_parameters()
{
    // m_tower_volume = nullptr;
    m_num_rows = 0;
    m_num_front_rows = 0;
    m_this_tower_theta = 0.0*deg;
    m_tower_tan_theta = 0.0;
    // m_back_shift = 0.0*mm;
}

void DDDRCaloTubes::DRconstructor::place_tower(Volume& calorimeter_volume,
                 Volume& tower_volume,
                 unsigned int stave, 
                 unsigned int layer,
                 unsigned int tower_id,
                 double phi)
{
    double goal_x = std::sin(90*deg-m_covered_theta)*std::cos(phi);
    double goal_y = std::sin(90*deg-m_covered_theta)*std::sin(phi);
    double goal_z = std::cos(90*deg-m_covered_theta);
    Direction goal_direction = Direction(goal_x, goal_y, goal_z);
    Direction start_direction = Direction(0, 0, 1);

    Direction v = start_direction.Cross(goal_direction);
    double s = std::sqrt(v.Mag2());
    double c = start_direction.Dot(goal_direction);

    double unit_matrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    TMatrixD unit = TMatrixD(3, 3, unit_matrix);
    double v_x_matrix[9] = {0, -v.Z(), v.Y(), v.Z(), 0, -v.X(), -v.Y(), v.X(), 0}; 
    TMatrixD v_x = TMatrixD(3, 3, v_x_matrix);

    double factor = (1-c)/(s*s);
    TMatrixD rot_matrix = unit + v_x + v_x*v_x*factor;
    Rotation3D rot = Rotation3D(rot_matrix(0,0), rot_matrix(0,1), rot_matrix(0,2),
                                 rot_matrix(1,0), rot_matrix(1,1), rot_matrix(1,2),
                                 rot_matrix(2,0), rot_matrix(2,1), rot_matrix(2,2));

    /* std::cout<<"tower_id = "<<tower_id<<std::endl;
    std::cout<<"stave    = "<<stave<<std::endl;
    std::cout<<"layer    = "<<layer<<std::endl;
    std::cout<<"phi      = "<<phi/deg<<std::endl;
    std::cout<<"cov_theta= "<<covered_theta/deg<<std::endl;
    std::cout<<"deg      = "<<thetaDegrees<<std::endl;
    std::cout<<"dec      = "<<thetaDecimal<<std::endl;
    std::cout<<"----------------------------------------" << std::endl; */

    // // Backward barrel region
    // Transform3D tower_bwd_tr(RotationZYX(0, phi, -90*deg-m_covered_theta), Position(tower_x, tower_y, tower_z));
    // PlacedVolume tower_bwd_placed = calorimeter_volume.placeVolume(tower_volume, -tower_id, tower_bwd_tr);
    // tower_bwd_placed.addPhysVolID("stave", -stave).addPhysVolID("layer", -layer);

    // PlacedVolume trap_towerVol_placed = calorimeter_volume.placeVolume(trap_towerVol, -tower_id, tower_bwd_tr);
    // trap_towerVol_placed.addPhysVolID("stave", -stave).addPhysVolID("layer", -layer);
    
    // Forward barrel region
    // Transform3D tower_fwd_tr(RotationZYX(180*deg, phi, -90*deg+m_covered_theta), Position(tower_x, tower_y, -tower_z));
    // Transform3D tower_fwd_tr(RotationZYX(90*deg, 90*deg, 0*deg), Position(std::cos(phi)*2*m, std::sin(phi)*2*m, 0));
    // Transform3D tower_fwd_tr(RotationZYX(0*deg, 0*deg, 0), Position(0, 0, 0));
    RotationZ rot_first = RotationZ(90*deg);
    RotationY rot_second = RotationY(90*deg-m_covered_theta);
    RotationY rot_third = RotationY(-m_covered_theta);
    RotationZ rot_fourth = RotationZ(phi);
    // RotationZYX rot_a = RotationZYX(90*deg, 90*deg-m_covered_theta, 0*deg);
    EulerAngles layer_rot = EulerAngles(0*deg, -m_covered_theta, 0*deg);

    // std::cout << "ROtation phi = " << (phi*1.0)/deg << std::endl;
    // std::cout << "phi = " << phi/deg << std::endl;
    // Rotation3D rot = rot_fourth*rot_second*rot_first; 
    // rot_a = rot_fourth*rot_a;
    Transform3D tower_fwd_tr(rot*rot_fourth*rot_first, m_tower_position);
    // Transform3D tower_fwd_tr(RotationZYX(0,0,0), m_tower_position);
    PlacedVolume tower_fwd_placed = calorimeter_volume.placeVolume(tower_volume, tower_id, tower_fwd_tr);
    tower_fwd_placed.addPhysVolID("stave", stave).addPhysVolID("layer", layer);

}


void DDDRCaloTubes::DRconstructor::construct_calorimeter(Volume& calorimeter_volume)
{
    unsigned int layer = 0;
    while (m_covered_theta<m_barrel_endcap_angle) 
    {
        // constructor.calculate_theta_parameters();
        // double theta = 90*deg - covered_theta;
        double delta_theta;
        // Assembly tower_volume("tower");
        Volume trap_volume("tower");
        trap_volume.setMaterial(m_capillary_material);
        this->construct_tower(trap_volume, delta_theta);


        double phi = 0*deg;
        // if (m_covered_theta<m_barrel_endcap_angle/2.0)
        // if (layer==0)
        // {
        std::cout << "layer = " << layer << std::endl;
        for (unsigned int stave=1; stave<=m_num_phi_towers; stave++, phi+=m_tower_phi)
        {
            // if (layer != 1) continue;
            this->calculate_tower_position(phi);
            unsigned int tower_id = stave + layer*m_num_phi_towers;
            this->place_tower(calorimeter_volume, trap_volume, stave, layer, tower_id, phi);
        }
        // }

        // covered_theta += delta_theta;
        this->increase_covered_theta(m_tower_theta);
        
        // if (layer >= 1) break;
        layer++;
    }
}