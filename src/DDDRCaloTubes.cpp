#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

static Ref_t create_detector(Detector& description,
                             xml_h entities,
                             SensitiveDetector sens) 
{

    xml_det_t   x_det       = entities;    
    int         det_id      = x_det.id();
    std::string det_name    = x_det.nameStr();

    Material    air         = description.air();
    Material    brass       = description.material("Brass");

    xml_dim_t   x_dim       = x_det.dimensions();
    double      rmin        = x_dim.rmin();
    double      rmax        = x_dim.rmax();
    double      z_half      = x_dim.zmax();

    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);
    Assembly      module_volume(det_name+"_module");

    sens.setType("calorimeter");


    xml_comp_t  x_tube      = x_det.child(_Unicode(tube));
    int         num_rows    = x_tube.number();
    int         num_cols    = x_tube.count();

    xml_comp_t  x_capillary         = x_tube.child(_Unicode(capillary));
    Material    capillary_material  = description.material(x_capillary.materialStr()); 
    double      capillary_inner_r   = x_capillary.inner_r();
    double      capillary_outer_r   = x_capillary.outer_r();

    xml_comp_t  x_scin_fibre         = x_tube.child(_Unicode(scin_fibre));
    Material    scin_fibre_material  = description.material(x_scin_fibre.materialStr()); 
    double      scin_fibre_inner_r   = x_scin_fibre.inner_r();
    double      scin_fibre_outer_r   = x_scin_fibre.outer_r();

    xml_comp_t  x_cher_fibre         = x_tube.child(_Unicode(cher_fibre));
    Material    cher_fibre_material  = description.material(x_cher_fibre.materialStr()); 
    double      cher_fibre_inner_r   = x_cher_fibre.inner_r();
    double      cher_fibre_outer_r   = x_cher_fibre.outer_r();

    Tube            capillary_solid(0.0*mm, capillary_outer_r, z_half);
    std::string     capillary_name = "capillary";//_"+std::to_string(row)+"_"+std::to_string(col);
    // Need two volumes for scintillation and Cherenkov channels
    Volume          scin_tube_volume(capillary_name, capillary_solid, capillary_material);
    Volume          cher_tube_volume(capillary_name, capillary_solid, capillary_material);
    if (x_capillary.isSensitive())
    {
        scin_tube_volume.setSensitiveDetector(sens);
        cher_tube_volume.setSensitiveDetector(sens);
    }
    scin_tube_volume.setAttributes(description, x_capillary.regionStr(), x_capillary.limitsStr(), x_capillary.visStr());
    cher_tube_volume.setAttributes(description, x_capillary.regionStr(), x_capillary.limitsStr(), x_capillary.visStr());    

    Tube            scin_fibre_solid(0.0*mm, scin_fibre_outer_r, z_half);
    std::string     scin_fibre_name = "scin_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
    Volume          scin_fibre_volume(scin_fibre_name, scin_fibre_solid, scin_fibre_material);
    if (x_scin_fibre.isSensitive())
    {
        scin_fibre_volume.setSensitiveDetector(sens);
    }
    PlacedVolume    scin_fibre_placed = scin_tube_volume.placeVolume(scin_fibre_volume);
    scin_fibre_volume.setAttributes(description, x_scin_fibre.regionStr(), x_scin_fibre.limitsStr(), x_scin_fibre.visStr());

    Tube            cher_fibre_solid(0.0*mm, cher_fibre_outer_r, z_half);
    std::string     cher_fibre_name = "cher_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
    Volume          cher_fibre_volume(cher_fibre_name, cher_fibre_solid, cher_fibre_material);
    if (x_cher_fibre.isSensitive())
    {
        cher_fibre_volume.setSensitiveDetector(sens);
    }
    PlacedVolume    cher_fibre_placed = cher_tube_volume.placeVolume(cher_fibre_volume);
    cher_fibre_volume.setAttributes(description, x_cher_fibre.regionStr(), x_cher_fibre.limitsStr(), x_cher_fibre.visStr());


    int tube_id = 0;

    for (int row=0; row<num_rows; row++)
    {
        for (int col=0; col<num_cols; col++)
        {
            double offset = (row & 1) ? capillary_outer_r : 0.0*mm;
            double x = col*2*capillary_outer_r + offset;
            double D = 4.0*capillary_outer_r/sqrt(3.0); // Long diagonal of hexagaon with capillary_outer_r as inradius
            double y = row*D*3.0/4.0;                     // Vertical spacing for hexagonal grid (pointy-top)
            auto position = Position(x, y, 0.0*mm);

            auto tube_to_be_placed = (row & 1) ? &cher_tube_volume : &scin_tube_volume;

            PlacedVolume    capillary_placed = module_volume.placeVolume(*tube_to_be_placed, tube_id, position);
            capillary_placed.addPhysVolID("module", tube_id);

            
            scin_fibre_placed.addPhysVolID("layer", tube_id);

            tube_id++;
        }
    }

    PlacedVolume module_placed = mother_volume.placeVolume(module_volume);
    module_placed.addPhysVolID("system", det_id);
    s_detElement.setPlacement(module_placed);

/*    
    //Make a Cylinder
    Tube envelope(rmin, rmax, z_half);
    Volume envelopeVol(det_name+"_envelope", envelope, air);
    PlacedVolume physvol = description.pickMotherVolume(s_detElement).placeVolume(envelopeVol);

    // add system ID and identify as barrel (as opposed to endcap +/-1)
    physvol.addPhysVolID("system", s_detElement.id()).addPhysVolID(_U(side),0);
    s_detElement.setPlacement(physvol);

    double currentInnerRadius = 0.0*mm; // running inner radius
    Layering layering(x_det); // convenience class
    int layerNum = 0;

    for(xml_coll_t c(x_det,_U(layer)); c; ++c, ++layerNum) {
        xml_comp_t x_layer = c;
        const Layer* lay = layering.layer(layerNum); // Get the layer from the layering engine.
        const double layerThickness = lay->thickness();

        //loop over the number of repetitions
        for(int i=0, repeat=x_layer.repeat(); i<repeat; ++i, ++layerNum) {
            std::string layerName = det_name + _toString(layerNum,"_layer%d");
            
            //make a volume for the layer
            Tube layerTube(currentInnerRadius, currentInnerRadius + layerThickness, z_half);
            Volume layerVol(layerName, layerTube, air);
            DetElement layerElement(s_detElement, layerName, layerNum);
            PlacedVolume layerVolPlaced = envelopeVol.placeVolume(layerVol);
            layerVolPlaced.addPhysVolID("layer",layerNum);

            //loop over slices
            int sliceNum = 0;
            for(xml_coll_t slice(x_layer,_U(slice)); slice; ++slice, ++sliceNum) {
                xml_comp_t x_slice = slice;
                double sliceThickness = x_slice.thickness();
                Material sliceMat = description.material(x_slice.materialStr());
                std::string sliceName = layerName + _toString(sliceNum,"slice%d");

                Tube sliceTube(currentInnerRadius,
                currentInnerRadius + sliceThickness, z_half);
                Volume sliceVol(sliceName, sliceTube, sliceMat);
                
                if ( x_slice.isSensitive() ) {
                    sliceVol.setSensitiveDetector(sens);
                }

                //place the slice in the layer
                layerVol.placeVolume(sliceVol);
                currentInnerRadius += sliceThickness;

            } // slices

        } //repetitions

    }//layers
 */    
   


    return s_detElement;
}

DECLARE_DETELEMENT(DDDRCaloTubes,create_detector)