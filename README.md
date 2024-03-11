# DRCaloTubes

Implementing the IDEA dual-readout calorimeter using a capillary tube geometry (Bucatini structure) 

**Installation**:
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install/ ..
make install -j6
cd ../install/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib64
```
(assuming `$PWD=<path_to_install_directory>`)

**Running geoDisplay**
```
cd compact
geoDisplay DDDRCaloTubes.xml
```

**Running geometry overlap check**
```
cd compact
ddsim --compactFile DDDRCaloTubes.xml --runType run --macroFile overlap.mac --part.userParticleHandler=''
```

**Running example simulation**
```
cd compact
ddsim --compactFile DDDRCaloTubes.xml --runType=batch -G -N=5 --steeringFile steering.py --outputFile=test.root --gun.position "14.5*mm 16.45*mm 100.0*cm" --gun.direction "0.0 0.0 -1.0" --gun.energy "30*GeV" --part.userParticleHandler=""   --gun.particle "e-"
```