#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
                                    xml_h entities,
                                    dd4hep::SensitiveDetector sens) {

    // XML Detector Element (confusingly also XML::DetElement)
    xml_det_t x_det = entities;

    // DetElement of our detector instance, attach additional information, sub-elements...
    // uses name of detector and ID number as defined in the XML detector tag
    std::string detName = x_det.nameStr();
    sens.setType("calorimeter");
    dd4hep::DetElement sdet (detName, x_det.id());

    //get the dimensions tag
    xml_dim_t dim = x_det.dimensions();
    //read its attributes
    double rmin = dim.rmin();
    double rmax = dim.rmax();
    double zmax = dim.zmax();

    //Make a Cylinder
    dd4hep::Tube envelope(rmin, rmax, zmax);
    dd4hep::Material air = theDetector.air();
    dd4hep::Volume envelopeVol(detName+"_envelope", envelope, air);
    dd4hep::PlacedVolume physvol = theDetector.pickMotherVolume(sdet).placeVolume(envelopeVol);

    // add system ID and identify as barrel (as opposed to endcap +/-1)
    physvol.addPhysVolID("system", sdet.id()).addPhysVolID(_U(side),0);
    sdet.setPlacement(physvol);

    double currentInnerRadius = rmin; // running inner radius
    dd4hep::Layering layering(x_det); // convenience class
    int layerNum = 0;

    for(xml_coll_t c(x_det,_U(layer)); c; ++c, ++layerNum) {
        xml_comp_t x_layer = c;
        const dd4hep::Layer* lay = layering.layer(layerNum); // Get the layer from the layering engine.
        const double layerThickness = lay->thickness();

        //loop over the number of repetitions
        for(int i=0, repeat=x_layer.repeat(); i<repeat; ++i, ++layerNum) {
            std::string layerName = detName + dd4hep::_toString(layerNum,"_layer%d");
            
            //make a volume for the layer
            dd4hep::Tube layerTube(currentInnerRadius, currentInnerRadius + layerThickness, zmax);
            dd4hep::Volume layerVol(layerName, layerTube, air);
            dd4hep::DetElement layerElement(sdet, layerName, layerNum);
            dd4hep::PlacedVolume layerVolPlaced = envelopeVol.placeVolume(layerVol);
            layerVolPlaced.addPhysVolID("layer",layerNum);

            //loop over slices
            int sliceNum = 0;
            for(xml_coll_t slice(x_layer,_U(slice)); slice; ++slice, ++sliceNum) {
                xml_comp_t x_slice = slice;
                double sliceThickness = x_slice.thickness();
                dd4hep::Material sliceMat = theDetector.material(x_slice.materialStr());
                std::string sliceName = layerName + dd4hep::_toString(sliceNum,"slice%d");

                dd4hep::Tube sliceTube(currentInnerRadius,
                currentInnerRadius + sliceThickness, zmax);
                dd4hep::Volume sliceVol(sliceName, sliceTube, sliceMat);
                
                if ( x_slice.isSensitive() ) {
                    sliceVol.setSensitiveDetector(sens);
                }

                //place the slice in the layer
                layerVol.placeVolume(sliceVol);
                currentInnerRadius += sliceThickness;

            } // slices

        } //repetitions

    }//layers


    return sdet;
}

DECLARE_DETELEMENT(DRCaloTubes,create_detector)