# DRCaloTubes

Implementing the IDEA dual-readout calorimeter using a capillary tube geometry (Bucatini structure) 

## Installation

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
**The export command needs to be executed in each new session!**

## Running the simulation

To run the following commands, make sure to change into the compact directory containing the DDDRCaloTubes.xml or adjust the paths in the command accordingly

### Running geoDisplay

```
geoDisplay DDDRCaloTubes.xml
```

### Running geometry overlap check

```
ddsim --compactFile DDDRCaloTubes.xml --runType run --macroFile overlap.mac --part.userParticleHandler=''
```

### Running example simulation

```
ddsim --compactFile DDDRCaloTubes.xml -G --steeringFile steering.py --outputFile=test.root --part.userParticleHandler=""
```