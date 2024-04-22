#ifndef DRconstructor_H
#define DRconstructor_H 1

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include "DRutils.h"

using namespace dd4hep;

namespace DDDRCaloTubes {

class DRconstructor {
public:
    // Constructor
    DRconstructor(Detector* description,
                  xml_h& entities,
                  SensitiveDetector* sens);

    // Destructor
    ~DRconstructor() {}

    void calculate_tower_parameters();
    void calculate_phi_parameters();
    void calculate_theta_parameters();
    double calculate_trap_width(double given_y, double given_z, bool backface = false);
    double calculate_tower_width(int given_row, int& n_tubes, bool backface = true);
    void assemble_tower(Volume& tower_air_volume);
    void construct_tower_trapezoid(Volume& trap_volume);
    void calculate_tower_position(double phi);
    void construct_tower(Volume& trap_volume);
    void increase_covered_theta(const double& delta_theta) {m_covered_theta += delta_theta;}
    void reset_tower_parameters();
    void place_tower(Volume& calorimeter_volume,
                     Volume& tower_volume,
                     unsigned int stave, 
                     unsigned int layer,
                     unsigned int tower_id,
                     double phi);

    void construct_calorimeter(Volume& calorimeter_volume);

private:
    Detector* m_description;
    xml_h m_entities;
    SensitiveDetector* m_sens;

    // Calorimeter parameters
    double m_calo_inner_r;
    double m_calo_outer_r;
    double m_calo_inner_half_z;

    double m_barrel_endcap_angle; // calculated from m_calo_inner_half_z and m_calo_inner_r
    double m_effective_inner_r;  // inner radius of the calorimeter after the front face is shifted

    // Tube parameters
    double    m_capillary_outer_r;
    double    m_scin_clad_outer_r;
    double    m_scin_core_outer_r;
    double    m_cher_clad_outer_r;
    double    m_cher_core_outer_r;
    Material m_capillary_material;
    Material m_scin_clad_material;
    Material m_scin_core_material;
    Material m_cher_clad_material;
    Material m_cher_core_material;
    std::string m_capillary_visString;
    std::string m_scin_clad_visString;
    std::string m_scin_core_visString;
    std::string m_cher_clad_visString;
    std::string m_cher_core_visString;
    bool m_capillary_isSensitive;
    bool m_scin_clad_isSensitive;
    bool m_scin_core_isSensitive;
    bool m_cher_clad_isSensitive;
    bool m_cher_core_isSensitive;


    double m_capillary_diameter; // calculated from m_capillary_outer_r

    // Constants used through the function (calculated from other parameters)
    // double m_D; // Long diagonal of hexagaon with capillary_outer_r as inradius
    double m_V; // Vertical spacing for pointy top oriented tubes
    double m_overlap; // Overlap between tubes

    // Tower parameters
    double m_tower_theta;
    double m_tower_phi;
    
    double m_tower_half_phi;
    double m_tower_tan_half_phi; // calculated from m_tower_phi
    double m_tower_half_length; // calculated from m_calo_inner_r and m_calo_outer_r and m_trap_half_length

    // Tower Phi parameters
    unsigned int m_num_cols;             // number of fibres in phi direction (back face) 
    unsigned int m_num_front_cols;       // number of fibres in front face in phi direction
    unsigned int m_num_phi_towers;       // number of towers in phi direction
    double m_tower_frontface_rightangleedge_x;
    double m_tower_frontface_thetaangleedge_x;
    double m_tower_backface_rightangleedge_x;
    double m_tower_backface_thetaangleedge_x;

    // Tower Theta parameters
    unsigned int m_num_rows;             // number of fibres in theta direction (back face)
    unsigned int m_num_front_rows;       // number of fibres in front face in theta direction
    double m_this_tower_theta;
    double m_tower_tan_theta;
    double m_tower_frontface_y;
    double m_tower_backface_y;
    double m_tower_polar_angle;
    double m_tower_azimuthal_angle;

    // Trapezoid support parameters
    double m_trap_wall_thickness_sides;
    double m_trap_wall_thickness_front;
    double m_trap_wall_thickness_back;
    double m_trap_frontface_rightangleedge_x;       // width for frontface
    double m_trap_frontface_thetaangleedge_x;
    double m_trap_backface_rightangleedge_x;        // width for backface
    double m_trap_backface_thetaangleedge_x;
    double m_trap_frontface_y;      // height for frontface
    double m_trap_backface_y;       // height for backface
    double m_trap_azimuthal_angle;  // azimuthal angle for the trapezoid
    double m_trap_polar_angle;      // polar angle for the trapezoid
    double m_trap_half_length;      // half length for the trapezoid
    Material m_trap_material;
    std::string m_trap_visString;


    // Construction parameters
    double m_covered_theta;
    double m_back_shift;
    Position m_tower_position;
    // Assembly* m_tower_volume;

    Material m_air;
    std::string m_air_visString;

};

} // namespace DDDRCaloTubes

#endif // DRCONSTRUCTOR_H
