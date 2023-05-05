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
    double      z_half      = x_dim.zhalf();
    double      rotation    = x_dim.angle();
    int         num_rows    = x_dim.number();
    int         num_cols    = x_dim.count();

    DetElement    s_detElement(det_name, det_id);
    Volume        mother_volume = description.pickMotherVolume(s_detElement);
    Assembly      module_volume(det_name+"_module");

    sens.setType("calorimeter");


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

    /*
    // Workaround for sensitive tubes:
    Tube        air_solid(0.0*mm, capillary_outer_r, z_half);
    std::string air_name = "air";
    Volume      scin_air_volume("scin_"+air_name, air_solid, air);
    Volume      cher_air_volume("cher_"+air_name, air_solid, air);
    scin_air_volume.setVisAttributes(description, x_capillary.visStr());
    cher_air_volume.setVisAttributes(description, x_capillary.visStr()); 

    // Scintillation capillary
    Tube        scin_tube_solid(0.0*mm, capillary_outer_r, z_half);
    std::string scin_tube_name = "scin_tube";
    Volume      scin_tube_volume(scin_tube_name, scin_tube_solid, capillary_material);
    if (x_capillary.isSensitive())
    {
        scin_tube_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume scin_tube_placed = scin_air_volume.placeVolume(scin_tube_volume);    
    scin_tube_volume.setVisAttributes(description, x_capillary.visStr());

    // Cherenkov capillary
    Tube        cher_tube_solid(0.0*mm, capillary_outer_r, z_half);
    std::string cher_tube_name = "cher_tube";
    Volume      cher_tube_volume(cher_tube_name, cher_tube_solid, capillary_material);
    if (x_capillary.isSensitive())
    {
        cher_tube_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume cher_tube_placed = cher_air_volume.placeVolume(cher_tube_volume); 
    cher_tube_volume.setVisAttributes(description, x_capillary.visStr()); 
    */   

    /* Code for tube construction before workaround for sensitive tubes
    // Construct volumes for tubes
    Tube        capillary_solid(0.0*mm, capillary_outer_r, z_half);
    std::string capillary_name = "capillary";//_"+std::to_string(row)+"_"+std::to_string(col);
    // Need two volumes for scintillation and Cherenkov channels
    Volume      scin_tube_volume("scin_"+capillary_name, capillary_solid, capillary_material);
    Volume      cher_tube_volume("cher_"+capillary_name, capillary_solid, capillary_material);
    if (x_capillary.isSensitive())
    {
        scin_tube_volume.setSensitiveDetector(sens);
        cher_tube_volume.setSensitiveDetector(sens);
    }
    scin_tube_volume.setVisAttributes(description, x_capillary.visStr());
    cher_tube_volume.setVisAttributes(description, x_capillary.visStr());    
    */
    /*
    // Scintillation cladding
    Tube        scin_clad_solid(0.0*mm, scin_clad_outer_r, z_half);
    std::string scin_clad_name = "scin_clad";
    Volume      scin_clad_volume(scin_clad_name, scin_clad_solid, scin_clad_material);
    if (x_scin_clad.isSensitive())
    {
        scin_clad_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume scin_clad_placed = scin_tube_volume.placeVolume(scin_clad_volume);
    scin_clad_volume.setVisAttributes(description, x_scin_clad.visStr());
    //scin_clad_placed.addPhysVolID("clad", 1);

    // Scintillation fibre
    Tube        scin_fibre_solid(0.0*mm, scin_fibre_outer_r, z_half);
    std::string scin_fibre_name = "scin_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
    Volume      scin_fibre_volume(scin_fibre_name, scin_fibre_solid, scin_fibre_material);
    if (x_scin_fibre.isSensitive())
    {
        scin_fibre_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume    scin_fibre_placed = scin_clad_volume.placeVolume(scin_fibre_volume);
    scin_fibre_volume.setVisAttributes(description, x_scin_fibre.visStr());
    //scin_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 0);

    // Cherenkov cladding
    Tube        cher_clad_solid(0.0*mm, cher_clad_outer_r, z_half);
    std::string cher_clad_name = "cher_clad";
    Volume      cher_clad_volume(cher_clad_name, cher_clad_solid, cher_clad_material);
    if (x_cher_clad.isSensitive())
    {
        cher_clad_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume cher_clad_placed = cher_tube_volume.placeVolume(cher_clad_volume);
    cher_clad_volume.setVisAttributes(description, x_cher_clad.visStr());
    //cher_clad_placed.addPhysVolID("clad", 1);

    // Chrerenkov fibre
    Tube        cher_fibre_solid(0.0*mm, cher_fibre_outer_r, z_half);
    std::string cher_fibre_name = "cher_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
    Volume      cher_fibre_volume(cher_fibre_name, cher_fibre_solid, cher_fibre_material);
    if (x_cher_fibre.isSensitive())
    {
        cher_fibre_volume.setSensitiveDetector(sens);
    }
    //PlacedVolume    cher_fibre_placed = cher_clad_volume.placeVolume(cher_fibre_volume);
    cher_fibre_volume.setVisAttributes(description, x_cher_fibre.visStr());
    //cher_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 1);
    */

    int tube_id = 0;
    double x_avg = 0.0*mm;
    double y_avg = 0.0*mm;

    for (int row=0; row<num_rows; row++)
    {
        for (int col=0; col<num_cols; col++)
        {
            // Definition of Air Volume which is created newly in each loop such that daughters can be 
            // repeatedly placed with identical copy numbers
            // TODO: Check what objects can be moved outside of loop (string _name, Tube _solid, etc.)

            // Capillary tube
            Tube        capillary_solid(0.0*mm, capillary_outer_r, z_half);
            std::string capillary_name = "capillary";
            Volume      capillary_volume(capillary_name, capillary_solid, capillary_material);
            if (x_capillary.isSensitive()) capillary_volume.setSensitiveDetector(sens);
            capillary_volume.setVisAttributes(description, x_capillary.visStr()); 

            if (row & 1) // Cherenkov row
            {
                // Cherenkov cladding
                Tube        cher_clad_solid(0.0*mm, cher_clad_outer_r, z_half);
                std::string cher_clad_name = "cher_clad";
                Volume      cher_clad_volume(cher_clad_name, cher_clad_solid, cher_clad_material);
                if (x_cher_clad.isSensitive()) cher_clad_volume.setSensitiveDetector(sens);
                PlacedVolume cher_clad_placed = capillary_volume.placeVolume(cher_clad_volume, tube_id);
                cher_clad_volume.setVisAttributes(description, x_cher_clad.visStr());
                cher_clad_placed.addPhysVolID("clad", 1);

                // Chrerenkov fibre
                Tube        cher_fibre_solid(0.0*mm, cher_fibre_outer_r, z_half);
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
                Tube        scin_clad_solid(0.0*mm, scin_clad_outer_r, z_half);
                std::string scin_clad_name = "scin_clad";
                Volume      scin_clad_volume(scin_clad_name, scin_clad_solid, scin_clad_material);
                if (x_scin_clad.isSensitive()) scin_clad_volume.setSensitiveDetector(sens);
                PlacedVolume scin_clad_placed = capillary_volume.placeVolume(scin_clad_volume, tube_id);
                scin_clad_volume.setVisAttributes(description, x_scin_clad.visStr());
                scin_clad_placed.addPhysVolID("clad", 1);

                // Scintillation fibre
                Tube        scin_fibre_solid(0.0*mm, scin_fibre_outer_r, z_half);
                std::string scin_fibre_name = "scin_fibre";//_"+std::to_string(row)+"_"+std::to_string(col);
                Volume      scin_fibre_volume(scin_fibre_name, scin_fibre_solid, scin_fibre_material);
                if (x_scin_fibre.isSensitive()) scin_fibre_volume.setSensitiveDetector(sens);
                PlacedVolume    scin_fibre_placed = scin_clad_volume.placeVolume(scin_fibre_volume, tube_id);
                scin_fibre_volume.setVisAttributes(description, x_scin_fibre.visStr());
                scin_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 0);
            }
            //auto tube_to_be_placed = (row & 1) ? &cher_air_volume : &scin_air_volume;

            
            // Configuration for placing the tube
            double offset = (row & 1) ? -capillary_outer_r : 0.0*mm;
            double x = col*2*capillary_outer_r + offset;
            double D = 4.0*capillary_outer_r/sqrt(3.0);     // Long diagonal of hexagaon with capillary_outer_r as inradius
            double y = row*D*3.0/4.0;                       // Vertical spacing for hexagonal grid (pointy-top)
            auto position = Position(x, y, 0.0*mm);

            //auto sensitive_to_be_placed = (row & 1) ? Volume(cher_tube_volume) : Volume(scin_tube_volume);

            //sensitive_fibre_placed.addPhysVolID("fibre", 1).addPhysVolID("cherenkov", 1);


            PlacedVolume    tube_placed = module_volume.placeVolume(capillary_volume, tube_id, position);
            // Axial coordinate conversion following https://www.redblobgames.com/grids/hexagons/#conversions-offset
            // Slighty changed to fit my q and r directions
            int q = col + (row - (row&1))/2;
            int r = row;
            //std::cout<<"(row, col) -> (r, q) : (" <<row<<", "<<col<<") -> (" << r<<", " <<q<<")" <<std::endl;
            tube_placed.addPhysVolID("fibre", 0).addPhysVolID("q", q).addPhysVolID("r", r);

            x_avg += x;
            y_avg += y;
            tube_id++;
        }
    }

    x_avg /= tube_id;
    y_avg /= tube_id;

    Transform3D tr(RotationZYX(0,rotation,0),Position(-x_avg,-y_avg,0));
    PlacedVolume module_placed = mother_volume.placeVolume(module_volume, tr);
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