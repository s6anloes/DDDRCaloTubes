from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import GeV, deg
SIM = DD4hepSimulation()

SIM.enableGun = True
SIM.runType = "batch"

SIM.action.calo = "DRCaloTubesSDAction"
SIM.action.calorimeterSDTypes = [u'calorimeter']
SIM.action.mapActions['DDDRCaloTubes'] = "DRCaloTubesSDAction"
SIM.filter.calo = ""
# # configure regex SD
# SIM.geometry.regexSensitiveDetector['DDDRCaloTubes'] = {'Match':['DRBT'],
#                                        'OutputLevel':3
#                                       }

SIM.gun.particle = "e+"
SIM.gun.position = ('0.0*cm', '0.0*cm', '0.0*cm')
SIM.gun.distribution = 'uniform'
SIM.gun.energy = 10*GeV
SIM.gun.multiplicity = 1
SIM.gun.phiMin = -0.5*deg
SIM.gun.phiMax = 0.5*deg
SIM.gun.thetaMin = 109.0*deg
SIM.gun.thetaMax = 110.0*deg

def setupCerenkov(kernel):
	from DDG4 import PhysicsList
	seq = kernel.physicsList()
	cerenkov = PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
	cerenkov.MaxNumPhotonsPerStep = 1000
	# cerenkov.MaxBetaChangePerStep = 10.0
	# cerenkov.TrackSecondariesFirst = True
	cerenkov.VerboseLevel = 0
	cerenkov.enableUI()
	seq.adopt(cerenkov)
	ph = PhysicsList(kernel, 'Geant4OpticalPhotonPhysics/OpticalGammaPhys')
	ph.addParticleConstructor('G4OpticalPhoton')
	ph.VerboseLevel = 0
	ph.enableUI()
	seq.adopt(ph)
	return None

SIM.physics.setupUserPhysics(setupCerenkov)

# Suppress default output file (due to size)
def exampleUserPlugin(dd4hepSimulation):
	from DDG4 import EventAction, Kernel
	dd = dd4hepSimulation  # just shorter variable name
	evt_root = EventAction(Kernel(), 'Geant4Output2ROOT/' + dd.outputFile, True)
	evt_root.HandleMCTruth = False
	evt_root.Control = False
	return None

SIM.outputConfig.userOutputPlugin = exampleUserPlugin

SIM.random.seed = 1234567890
SIM.numberOfEvents = 1
