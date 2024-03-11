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
    void assemble_tower(Assembly& tower_volume);
    void calculate_tower_position();
    void construct_tower(Assembly& tower_volume, double& delta_theta);
    void increase_covered_theta(const float& delta_theta) {m_covered_theta += delta_theta;}
    void reset_tower_parameters();
    void place_tower(Volume& calorimeter_volume,
                     Assembly& tower_volume,
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
    float m_calo_inner_r;
    float m_calo_outer_r;
    float m_calo_inner_half_z;

    float m_barrel_endcap_angle; // calculated from m_calo_inner_half_z and m_calo_inner_r
    float m_effective_inner_r;  // inner radius of the calorimeter after the front face is shifted

    // Tube parameters
    float    m_capillary_outer_r;
    float    m_scin_clad_outer_r;
    float    m_scin_core_outer_r;
    float    m_cher_clad_outer_r;
    float    m_cher_core_outer_r;
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


    float m_capillary_diameter; // calculated from m_capillary_outer_r

    // Constants used through the function (calculated from other parameters)
    // float m_D; // Long diagonal of hexagaon with capillary_outer_r as inradius
    float m_V; // Vertical spacing for pointy top oriented tubes
    float m_overlap; // Overlap between tubes

    // Tower parameters
    float m_tower_theta;
    float m_tower_phi;

    float m_tower_tan_phi; // calculated from m_tower_phi
    float m_tower_half_length; // calculated from m_calo_inner_r and m_calo_outer_r

    // Tower Phi parameters
    unsigned int m_num_cols;             // number of fibres in phi direction (back face) 
    unsigned int m_num_front_cols;       // number of fibres in front face in phi direction
    unsigned int m_num_phi_towers;       // number of towers in phi direction

    // Tower Theta parameters
    unsigned int m_num_rows;             // number of fibres in theta direction (back face)
    unsigned int m_num_front_rows;       // number of fibres in front face in theta direction
    float m_this_tower_theta;
    float m_tower_tan_theta;

    // Construction parameters
    float m_covered_theta;
    float m_back_shift;
    Position m_tower_position;
    // Assembly* m_tower_volume;



};

} // namespace DDDRCaloTubes

#endif // DRCONSTRUCTOR_H
