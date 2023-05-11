from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV
SIM = DD4hepSimulation()

SIM.enableGun = True
SIM.runType = "batch"

SIM.action.calo = "DRCaloTubesSDAction"
SIM.action.calorimeterSDTypes = [u'calorimeter']
SIM.action.mapActions['DDDRCaloTubes'] = "DRCaloTubesSDAction"
SIM.filter.calo = ""
SIM.gun.particle = "e+"
SIM.gun.position = ('0.87275324641*cm', '0.0*cm', '-100.0*cm')
SIM.gun.momentumMin = 20*GeV
SIM.gun.momentumMax = 20*GeV

def setupCerenkov(kernel):
	from DDG4 import PhysicsList
	seq = kernel.physicsList()
	cerenkov = PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
	cerenkov.MaxNumPhotonsPerStep = 10
	cerenkov.MaxBetaChangePerStep = 10.0
	cerenkov.TrackSecondariesFirst = True
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
SIM.numberOfEvents = 10
