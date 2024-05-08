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
    m_trap_material             = m_description->material(x_trap_support.materialStr());
    m_trap_visString            = x_trap_support.visStr();


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


    m_tower_half_phi = m_tower_phi/2.0; // Half of the tower phi angle
    m_tower_tan_half_phi = std::tan(m_tower_half_phi); // Needed several times, calculate once here
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

    m_tower_frontface_y = m_trap_frontface_y + m_trap_wall_thickness_front*m_tower_tan_theta - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));
    m_tower_backface_y  = m_trap_backface_y  - m_trap_wall_thickness_back*m_tower_tan_theta  - m_trap_wall_thickness_sides*(1+1/std::cos(m_tower_theta));
    
    // Distance by which straight edge of this tower is shifted backwards to ensure inner radius of calorimeter
    m_back_shift = std::tan(m_covered_theta)*m_trap_frontface_y;                      

    // Frontface width of trapezoid support
    m_trap_frontface_rightangleedge_x = (m_calo_inner_r+std::cos(m_covered_theta)*m_back_shift)*2*m_tower_tan_half_phi;
    m_trap_frontface_thetaangleedge_x = m_calo_inner_r*2*m_tower_tan_half_phi;


    double tower_backface_phi_increase_rightangleedge = 2*m_trap_half_length*std::cos(m_covered_theta) * 2*m_tower_tan_half_phi; // by how much the backface wides for both tower and trap
    double tower_backface_phi_increase_thetaangleedge = 2*m_trap_half_length/std::cos(m_tower_theta)*std::cos(m_covered_theta+m_tower_theta) * 2*m_tower_tan_half_phi; // by how much the backface wides for both tower and trap

    m_trap_backface_rightangleedge_x = m_trap_frontface_rightangleedge_x + tower_backface_phi_increase_rightangleedge;
    m_trap_backface_thetaangleedge_x = m_trap_frontface_thetaangleedge_x + tower_backface_phi_increase_thetaangleedge;

    double effective_side_wall_thickness_x = m_trap_wall_thickness_sides/std::cos(m_tower_phi);

    m_angle_edges_x = std::atan2((m_tower_backface_rightangleedge_x-m_tower_backface_thetaangleedge_x)/2.0, m_tower_backface_y);

    m_tower_frontface_rightangleedge_x = calculate_trap_width(m_trap_wall_thickness_sides, m_trap_wall_thickness_front, false) - 2.0*effective_side_wall_thickness_x;

    m_tower_backface_rightangleedge_x = calculate_trap_width(m_trap_wall_thickness_sides, m_trap_wall_thickness_back, true) - 2.0*effective_side_wall_thickness_x;

    m_tower_frontface_thetaangleedge_x = calculate_trap_width(m_tower_frontface_y+m_trap_wall_thickness_sides, m_trap_wall_thickness_front, false) - 2*effective_side_wall_thickness_x;

    m_tower_backface_thetaangleedge_x = calculate_trap_width(m_tower_backface_y+m_trap_wall_thickness_sides, m_trap_wall_thickness_back, true) - 2*effective_side_wall_thickness_x;

}


double DDDRCaloTubes::DRconstructor::calculate_trap_width(double given_y, double given_z, bool backface)
{
    // Calculate width (x_direction) of trapezoid at given y
    // Assuming y=0 corresponds to m_trap_frontface_rightangleedge_x
    if (given_z>2*m_trap_half_length) throw std::runtime_error("calculate_trap_width: Given z is larger than length of trapezoid");
    double max_y_for_given_z;
    if (backface) max_y_for_given_z = m_trap_backface_y - given_z*m_tower_tan_theta;
    else          max_y_for_given_z = m_trap_frontface_y + given_z*m_tower_tan_theta;
    if (given_y > max_y_for_given_z) throw std::runtime_error("calculate_trap_width: Given y is larger than maximum y for given z");

    double delta_y = -2.0*given_y*std::tan(m_angle_edges_x);
    double delta_z;
    if (backface) delta_z = -2.0*given_z*std::cos(m_covered_theta)*m_tower_tan_half_phi;
    else          delta_z =  2.0*given_z*std::cos(m_covered_theta)*m_tower_tan_half_phi;

    double trap_x;
    if (backface) trap_x = m_trap_backface_rightangleedge_x  + delta_y + delta_z;
    else          trap_x = m_trap_frontface_rightangleedge_x + delta_y + delta_z;

    return trap_x;

}

double DDDRCaloTubes::DRconstructor::calculate_tower_width(int given_row, bool backface)
{
    // Calculate width (x_direction) of tower at given row
    // Assuming row 0 is at the right angle edge

    // y distance where tube hits side of the wall with given angle between backfaces (y=0 corresponds to m_tower_backface_rightangleedge_x)
    double y = m_capillary_outer_r + given_row*m_V + std::cos(90*deg-m_angle_edges_x)*m_capillary_outer_r;


    double tower_x;
    if (backface) tower_x = m_tower_backface_rightangleedge_x  - 2.0*y*std::tan(m_angle_edges_x);
    else          tower_x = m_tower_frontface_rightangleedge_x - 2.0*y*std::tan(m_angle_edges_x);

    return tower_x;
}

void DDDRCaloTubes::DRconstructor::assemble_tower(Volume& tower_air_volume)
{

    // Placement and shortening of tube depends on whether the rows have an offset or not
    
    // Some coordinates used throughout the functions
    double tower_centre_r = m_tower_half_length/std::cos(m_tower_polar_angle);
    double tower_centre_half_y = tower_centre_r*std::sin(m_tower_polar_angle)*std::sin(m_tower_azimuthal_angle) + m_tower_frontface_y/2.0;


    Position back_upper_left_corner  = Position(-m_tower_backface_rightangleedge_x/2.0,  -tower_centre_half_y, m_tower_half_length);
    Position back_lower_left_corner  = Position(-m_tower_backface_thetaangleedge_x/2.0,  tower_centre_half_y,  m_tower_half_length);
    Position front_upper_left_corner = Position(-m_tower_frontface_rightangleedge_x/2.0, -tower_centre_half_y, -m_tower_half_length);

    Position back_upper_right_corner  = Position(m_tower_backface_rightangleedge_x/2.0,  -tower_centre_half_y, m_tower_half_length);
    Position back_lower_right_corner  = Position(m_tower_backface_thetaangleedge_x/2.0,  tower_centre_half_y,  m_tower_half_length);
    Position front_upper_right_corner = Position(m_tower_frontface_rightangleedge_x/2.0, -tower_centre_half_y, -m_tower_half_length);


    std::vector<double> plane_left_coefficients = get_plane_equation(back_upper_left_corner, back_lower_left_corner, front_upper_left_corner);
    std::vector<double> plane_right_coefficients = get_plane_equation(back_upper_right_corner, back_lower_right_corner, front_upper_right_corner);
    Direction line_direction = Direction(0, 0, -1);


    /* Position front_lower_left_corner  = Position(-m_tower_frontface_thetaangleedge_x/2.0, tower_centre_half_y, -m_tower_half_length);
    double distance_left = distance_from_plane(plane_left_coefficients, front_lower_left_corner);
    if (distance_left > 1e-6*mm) std::cout << "############################# FRONT LOWER LEFT CORNER NOT ON PLANE BY << "<< distance_left/mm << " mm #############################" << std::endl; */

    double tower_x = calculate_tower_width(0);
    unsigned int num_back_cols_rightangleedge = 1 + fast_floor((tower_x-2*std::sin(90*deg-m_angle_edges_x)*m_capillary_outer_r)/m_capillary_diameter); 
    unsigned int half_num_back_cols = num_back_cols_rightangleedge/2;
    double offset_col_x = (num_back_cols_rightangleedge&1) ? m_capillary_outer_r : 0.0*mm;

    double first_tube_x = half_num_back_cols*m_capillary_diameter + offset_col_x;
    double gap_to_wall = m_tower_backface_rightangleedge_x/2.0 - first_tube_x;


    int num_bad_rows = 0;
    // int num_bad_cols = 0;

    double covered_tower_y = m_capillary_diameter;

    // int row = 0;
    // while(covered_tower_y < m_tower_backface_y)
    unsigned int num_rows = fast_floor((m_tower_backface_y-m_capillary_diameter)/m_V) + 1; 
    std::cout << "TOTAL ROWS = " << num_rows << " COLS = " << num_back_cols_rightangleedge << std::endl;
    for (unsigned int row = 0; row < num_rows; row++, covered_tower_y+=m_V)
    {

        // Staggering of tubes at the lower edge (theta edge/face)
        double row_staggered_z = 0.0*mm;
        if (covered_tower_y > m_tower_frontface_y) row_staggered_z = (covered_tower_y-m_tower_frontface_y)/m_tower_tan_theta/2.0;
        // else std::cout << "Row: " << row << std::endl;

        // Failsafe for tubes which would have 'negative length'
        // Should not happen, but if it does, following rows will also be too short, so can skip the rest
        if (row_staggered_z > m_tower_half_length) {
            num_bad_rows = num_rows-row;
            std::cout << "Encountered bad row at row " << row << std::endl;
            std::cout << "Number of leftover bad rows: " << num_bad_rows << std::endl;
            break;
        }


        double offset = (row & 1) ? m_capillary_outer_r : 0.0*mm;

        double covered_tower_x = gap_to_wall+m_capillary_diameter+offset;

        tower_x = calculate_tower_width(row);

        double tower_front_x = calculate_tower_width(row, false);

        for (unsigned int col = 0; col < num_back_cols_rightangleedge; col++, covered_tower_x+=m_capillary_diameter)
        {
            // TODO: Check what objects can be moved outside of loop (string _name, Tube _solid, etc.)

    

            // Configuration for placing the tube
            double x = col*m_capillary_diameter + offset + m_capillary_outer_r;
            double y = row*m_V + m_capillary_outer_r;                       // Vertical spacing for hexagonal grid (pointy-top)

            double x_pos = x - first_tube_x;
            int sign = (x_pos > 0) ? 1 : ((x_pos < 0) ? -1 : 0);
            if (std::abs(x_pos+sign*m_capillary_outer_r*std::cos(m_angle_edges_x)) > tower_x/2.0)
            {
                // std::cout << "Skipping tube because of non-rectangular tower" << std::endl;
                continue;
            }

            double col_staggered_z = 0.0*mm;
            if (tower_x/2.0-(covered_tower_x-m_capillary_diameter) + m_capillary_outer_r*(1.0-std::cos(m_angle_edges_x)) > tower_front_x/2.0)
            {
                Position line_point = Position(x-first_tube_x-m_capillary_outer_r*std::cos(m_angle_edges_x), y-tower_centre_half_y-m_capillary_outer_r*std::sin(m_angle_edges_x), m_tower_half_length);
                Position intersection = get_intersection(plane_left_coefficients, line_point, line_direction);

                // Position intersection = get_intersection(Direction(plane_left_coefficients[0], plane_left_coefficients[1], plane_left_coefficients[2]), back_upper_left_corner, line_point, line_direction);

                
                col_staggered_z = (m_tower_half_length + intersection.z())/2.0;

            } else if (covered_tower_x-tower_x/2.0 > tower_front_x/2.0)
            {
                Position line_point = Position(x-first_tube_x+m_capillary_outer_r*std::cos(m_angle_edges_x), y-tower_centre_half_y-m_capillary_outer_r*std::sin(m_angle_edges_x), m_tower_half_length);
                Position intersection = get_intersection(plane_right_coefficients, line_point, line_direction);

                col_staggered_z = (m_tower_half_length + intersection.z())/2.0;
            }

            // if (col_staggered_z > m_tower_half_length) {
            //     num_bad_cols++;
            //     std::cout << "m_tower_half_length = " << m_tower_half_length/mm << " mm" << std::endl;
            //     continue;
            // }


            double z = (row_staggered_z > col_staggered_z) ? row_staggered_z : col_staggered_z ;

            // Reference point for tube placement is trapezoid centre
            auto position = Position(x_pos, y-tower_centre_half_y, z);

            // Offset coordinates following https://www.redblobgames.com/grids/hexagons/#coordinates-offset
            unsigned short int q = col;
            unsigned short int r = row;

            // TubeID composed of q in first 16 bits, r in last 16 bits
            unsigned int tube_id = (q << 16) | r;
            
            // std::cout<<"(row, col) -> (r, q) -> (tubeID) : (" <<row<<", "<<col<<") -> (" <<r<<", " <<q<<") -> (" << tube_id << ")" <<std::endl; 

            double tube_half_length = m_tower_half_length - z;


            if (tube_half_length < 0*mm) 
            {
                // covered_tower_x += m_capillary_diameter;
                // col++;
                std::cout << "m_tower_half_length = " << m_tower_half_length/um << " um" << std::endl;
                continue;
            }
            // std::cout << "Tube length row_"<<row<<"_col_"<<col<< " " << 2*tube_half_length/mm << " mm" << std::endl;

            

            // Capillary tube
            // std::string capillary_name = "capillary_" + std::to_string(row) + "_" + std::to_string(col);
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
                Tube        cher_core_solid(0.0*mm, m_cher_core_outer_r, tube_half_length);
                Volume      cher_core_volume("cher_fibre", cher_core_solid, m_cher_core_material);
                if (m_cher_core_isSensitive) cher_core_volume.setSensitiveDetector(*m_sens);
                PlacedVolume    cher_core_placed = cher_clad_volume.placeVolume(cher_core_volume, tube_id);
                cher_core_volume.setVisAttributes(*m_description, m_cher_core_visString);
                cher_core_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);
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
                Tube        scin_core_solid(0.0*mm, m_scin_core_outer_r, tube_half_length);
                Volume      scin_core_volume("scin_fibre", scin_core_solid, m_scin_core_material);
                if (m_scin_core_isSensitive) scin_core_volume.setSensitiveDetector(*m_sens);
                PlacedVolume    scin_core_placed = scin_clad_volume.placeVolume(scin_core_volume, tube_id);
                scin_core_volume.setVisAttributes(*m_description, m_scin_core_visString);
                scin_core_placed.addPhysVolID("core", 1).addPhysVolID("clad", 0);
            }

            

            PlacedVolume    tube_placed = tower_air_volume.placeVolume(capillary_volume, tube_id, position);
            tube_placed.addPhysVolID("air",0).addPhysVolID("col", col).addPhysVolID("row", row);
            // std::cout << "col = " << col << std::endl;
            
        }

    }

}


// Function to calculate the position of the tower in stave (for fixed phi = 0*deg)
void DDDRCaloTubes::DRconstructor::calculate_tower_position(double phi)
{
    
    double trap_centre_r = m_trap_half_length/std::cos(m_trap_polar_angle);

    // double trap_centre_half_x = trap_centre_r*std::sin(m_trap_polar_angle)*std::cos(m_trap_azimuthal_angle) + m_trap_frontface_rightangleedge_x/2.0;
    double trap_centre_half_y = trap_centre_r*std::sin(m_trap_polar_angle)*std::sin(m_trap_azimuthal_angle) + m_trap_frontface_y/2.0;
    double trap_rad_centre = m_calo_inner_r/std::cos(m_covered_theta) + m_back_shift + m_trap_half_length;

    double stave_x = std::cos(m_covered_theta)*trap_rad_centre - std::sin(m_covered_theta)*trap_centre_half_y;
    // double stave_y = 0;
    double stave_z = std::sin(m_covered_theta)*trap_rad_centre + std::cos(m_covered_theta)*trap_centre_half_y;

    
    double tower_x = std::cos(phi)*stave_x;
    double tower_y = std::sin(phi)*stave_x;
    double tower_z = stave_z;

    m_tower_position = dd4hep::Position(tower_x, tower_y, tower_z);
}


// Function to construct the trapezoidal supoprt structure for the tower in which fibres are placed
void DDDRCaloTubes::DRconstructor::construct_tower_trapezoid(Volume& trap_volume)
{

        
        // polar coordinate conversion
        // double delta_x = 0;
        double delta_y = (m_trap_backface_y - m_trap_frontface_y)/2.0;
        double delta_z = 2.0*m_trap_half_length;
        m_trap_polar_angle = std::acos(delta_z/std::sqrt(delta_y*delta_y + delta_z*delta_z));
        m_trap_azimuthal_angle = 90.0*deg;

        Trap trap_solid("trap_solid", m_trap_half_length, m_trap_polar_angle, m_trap_azimuthal_angle, 
                                      m_trap_frontface_y/2.0, m_trap_frontface_rightangleedge_x/2.0, m_trap_frontface_thetaangleedge_x/2.0, 0.,
                                      m_trap_backface_y/2.0,  m_trap_backface_rightangleedge_x/2.0,  m_trap_backface_thetaangleedge_x/2.0,  0.);


        // Air volume in which fibres are placed
        // double delta_x_air = 0;
        double delta_y_air = (m_tower_backface_y - m_tower_frontface_y)/2.0;
        double delta_z_air = 2.0*m_tower_half_length;
        m_tower_polar_angle = std::acos(delta_z_air/std::sqrt(delta_y_air*delta_y_air + delta_z_air*delta_z_air));
        m_tower_azimuthal_angle = 90.0*deg;

        Trap tower_air_solid("tower_solid", m_tower_half_length, m_tower_polar_angle, m_tower_azimuthal_angle, 
                                      m_tower_frontface_y/2.0, m_tower_frontface_rightangleedge_x/2.0, m_tower_frontface_thetaangleedge_x/2.0, 0.,
                                      m_tower_backface_y/2.0,  m_tower_backface_rightangleedge_x/2.0,  m_tower_backface_thetaangleedge_x/2.0,  0.);
                                
        Position tower_air_pos = Position(0,
                                         (1.0-1.0/std::cos(m_tower_theta))*m_trap_wall_thickness_sides,
                                         (m_trap_wall_thickness_front-m_trap_wall_thickness_back)/2.0/* +10*nm */);


        // Subtraction solid used sometimes for easier visualisation. NOT TO BE USED IN FINAL GEOMETRY
        // SubtractionSolid solid = SubtractionSolid("trap_final", trap_solid, tower_air_solid, tower_air_pos);
        Volume tower_air_volume("tower_air_volume", tower_air_solid, m_air);
        tower_air_volume.setVisAttributes(*m_description, m_air_visString);

        trap_volume.setSolid(trap_solid);
        trap_volume.setVisAttributes(*m_description, m_trap_visString);

        PlacedVolume tower_air_placed = trap_volume.placeVolume(tower_air_volume, tower_air_pos);
        tower_air_placed.addPhysVolID("air", 1);

        this->assemble_tower(tower_air_volume);
    
}


void DDDRCaloTubes::DRconstructor::construct_tower(Volume& trap_volume)
{
    
    this->calculate_theta_parameters();

    this->construct_tower_trapezoid(trap_volume);

}


void DDDRCaloTubes::DRconstructor::reset_tower_parameters()
{
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
        std::cout << "layer = " << layer << std::endl;
        Volume trap_volume("tower");
        trap_volume.setMaterial(m_trap_material);
        this->construct_tower(trap_volume);


        double phi = 0*deg;
        for (unsigned int stave=1; stave<=m_num_phi_towers; stave++, phi+=m_tower_phi)
        {
            // if (layer != 43) continue;
            this->calculate_tower_position(phi);
            unsigned int tower_id = stave + layer*m_num_phi_towers;
            this->place_tower(calorimeter_volume, trap_volume, stave, layer, tower_id, phi);
        }

        this->increase_covered_theta(m_tower_theta);
        
        // if (layer >= 1) break;
        layer++;
    }
}