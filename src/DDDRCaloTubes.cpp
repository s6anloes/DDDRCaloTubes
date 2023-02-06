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

    xml_comp_t  x_scinfibre         = x_tube.child(_Unicode(scinfibre));
    Material    scinfibre_material  = description.material(x_scinfibre.materialStr()); 
    double      scinfibre_inner_r   = x_scinfibre.inner_r();
    double      scinfibre_outer_r   = x_scinfibre.outer_r();

    
    for (int row=0; row<num_rows; row++)
    {
        for (int col=0; col<num_cols; col++)
        {
            auto position = Position(row*2*rmax, col*2*rmax, 0.0*mm);

            Tube            capillary_solid(0.0*mm, capillary_outer_r, z_half);
            std::string     capillary_name = "capillary_"+std::to_string(row)+"_"+std::to_string(col);
            Volume          capillary_volume(capillary_name, capillary_solid, capillary_material);
            if (x_capillary.isSensitive())
            {
                capillary_volume.setSensitiveDetector(sens);
            }
            PlacedVolume    capillary_placed = module_volume.placeVolume(capillary_volume, position);
            capillary_placed.addPhysVolID("module", (row+1)*(col+1));

            Tube            scinfibre_solid(0.0*mm, scinfibre_outer_r, z_half);
            std::string     scinfibre_name = "scinfibre_"+std::to_string(row)+"_"+std::to_string(col);
            Volume          scinfibre_volume(scinfibre_name, scinfibre_solid, scinfibre_material);
            if (x_scinfibre.isSensitive())
            {
                scinfibre_volume.setSensitiveDetector(sens);
            }
            PlacedVolume    scinfibre_placed = capillary_volume.placeVolume(scinfibre_volume);
            scinfibre_placed.addPhysVolID("layer", (row+1)*(col+1));
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