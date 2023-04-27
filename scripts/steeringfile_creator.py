#! /usr/bin/env python

import numpy as np
#from module_layout import create_points
from pathlib import Path
from random import randrange
from itertools import cycle
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
pdg_code_to_geant4_D = {"11":"e-",
                        "-11":"e+",
                        "12":"nu_e",
                        "13":"mu-",
                        "211":"pi+",
                        "-211":"pi-",
                        "0":"geantino",}


def create_points(xstart=-4.5, npoints=19):
        
        pointlist = []
        # shift for certain beam spots. Variable toggles between +1 and -1
        shift = cycle([1, -1])
        # stepsize for symmetric points around (0,0)
        step = np.abs(xstart*2)/(npoints-1) if npoints>1 else 0
        for i in range(npoints):
                m = 26.0/45.0
                x = i*step + xstart
                y = x*m
                # As it happens the following points need a slight shift (only works if start point is fibre centre)
                if (x.is_integer()) and (x%3!=0):
                        #print(x, " correction needed")
                        angle = np.arctan2(26.0, 45.0)
                        # 0.172542824 is roughly the shift needed in radial direction w.r.t. the closest fibre to get to the center of the tube
                        r = x/np.cos(angle) + next(shift)*0.172542824 
                        # calculate new x and y
                        x = r*np.cos(angle)
                        y = x*m
                        
                pointlist.append((x, y))
                
        return pointlist

def create_file(pdg, momentum, x, y, outdir='runcards_dump/', prefix='runcard'):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    if y==None:
        y=x/2.0
    with open(outdir+prefix+'_pdg'+pdg+'_p'+str(momentum)+'_x'+str(int(round(x*100)))+'_y'+str(int(round(y*100)))+'_50k.py', 'w') as f:
        f.write('from DDSim.DD4hepSimulation import DD4hepSimulation\n')
        f.write('from g4units import mm, GeV, MeV\n')
        f.write('SIM = DD4hepSimulation()\n\n')

        f.write('SIM.enableGun = True\n')

        f.write('SIM.runType = "batch"\n\n')

        f.write('SIM.action.calo = "DRCaloTubesSDAction"\n')
        f.write('SIM.action.calorimeterSDTypes = [u\'calorimeter\']\n')
        f.write('SIM.action.mapActions[\'DDDRCaloTubes\'] = "DRCaloTubesSDAction"\n')
        f.write('SIM.filter.calo = ""\n')
        f.write(f'SIM.gun.particle = "{pdg_code_to_geant4_D[pdg]}"\n')
        f.write(f'SIM.gun.position = (\'{x+0.87275324641}*cm\', \'{y}*cm\', \'-100.0*cm\')\n')
        f.write(f'SIM.gun.momentumMin = {momentum}*GeV\n')
        f.write(f'SIM.gun.momentumMax = {momentum}*GeV\n\n')

        f.write('def setupCerenkov(kernel):\n')
        f.write('\tfrom DDG4 import PhysicsList\n')
        f.write('\tseq = kernel.physicsList()\n')
        f.write('\tcerenkov = PhysicsList(kernel, \'Geant4CerenkovPhysics/CerenkovPhys\')\n')
        f.write('\tcerenkov.MaxNumPhotonsPerStep = 10\n')
        f.write('\tcerenkov.MaxBetaChangePerStep = 10.0\n')
        f.write('\tcerenkov.TrackSecondariesFirst = True\n')
        f.write('\tcerenkov.VerboseLevel = 0\n')
        f.write('\tcerenkov.enableUI()\n')
        f.write('\tseq.adopt(cerenkov)\n')
        f.write('\tph = PhysicsList(kernel, \'Geant4OpticalPhotonPhysics/OpticalGammaPhys\')\n')
        f.write('\tph.addParticleConstructor(\'G4OpticalPhoton\')\n')
        f.write('\tph.VerboseLevel = 0\n')
        f.write('\tph.enableUI()\n')
        f.write('\tseq.adopt(ph)\n')
        f.write('\treturn None\n\n')

        f.write('SIM.physics.setupUserPhysics(setupCerenkov)\n')

        f.write('def exampleUserPlugin(dd4hepSimulation):\n')
        f.write('\tfrom DDG4 import EventAction, Kernel\n')
        f.write('\tdd = dd4hepSimulation  # just shorter variable name\n')
        f.write('\tevt_root = EventAction(Kernel(), \'Geant4Output2ROOT/\' + dd.outputFile, True)\n')
        f.write('\tevt_root.HandleMCTruth = False\n')
        f.write('\tevt_root.Control = False\n')
        #f.write('\tevt_root.Output = None\n')
        #f.write('\tevt_root.enableUI()\n')
        #f.write('\tKernel().eventAction().add(evt_root)\n')
        f.write('\treturn None\n\n')

        f.write('SIM.outputConfig.userOutputPlugin = exampleUserPlugin\n\n')


        f.write(f'SIM.random.seed = {randrange(999999999)}\n')

        #f.write('/gps/position '+str(x+0.87275324641)+' '+str(y)+' -200 cm\n')

        #f.write(f'/gps/energy {energy} GeV\n')

        f.write('SIM.numberOfEvents = 10000\n')
    
    f.close()

#for i in np.linspace(-0.45,0.45,19):
#    create_file(i)
#    print(i)

#pointlist = [(x, -1.25) for x in np.linspace(35, 50, 31)]
#print(pointlist)
pointlist = create_points()
pdgcodelist = [-11]
momentumlist = [20]
category = 'impact_comparison'
outdir = '/its/home/al723/dualrocalo/simulation/DD4hep/dremtubes_comparison/'+category+'/'+timestr+'/' # format '/absolute/path/'
id = 0
for point in pointlist:
    for m in momentumlist:
        for pdg in pdgcodelist:
            prefix = category+f"_{id:02}"
            print(point)
            # conversion from mm to cm
            newpoint = tuple([x/10.0 for x in point])
            print(newpoint[0], newpoint[1])
            create_file(str(pdg), m, newpoint[0], newpoint[1], outdir=outdir+f'steering_{id:05d}/', prefix=prefix)
            id += 1
""" """