# -*- coding: utf-8 -*-

from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers


# Simulation population and outputs
saveEvery = 100
renderEvery = 100
max_population = 100000

# cell types and properties
number_of_types = 2
# 0 (red) and 1 (green)
cell_colors = {0:[1.0, 0.0, 0.0], 1:[0.0, 1.0, 0.0]}
    
# Properties of cell type 0
rad0 = 0.75 #radius
len0 = 3 # length
initialVol0 = len0
targetVol0 = 2*initialVol0 # Dividing length
g0 = 3

# Properties of cell type 0
rad1 =  0.75
len1 = 3
initialVol1 = len1
targetVol1 = 2*initialVol1
g1 = 1



def setup(sim):
    global TimeStep

    # Set biophysics and regulation models
    biophys = CLBacterium(sim, max_cells=max_population, jitter_z=False)
    regul = ModuleRegulator(sim)

    # add biophysics, regulation, CDK and EPS objects to simulator
    sim.init(biophys, regul, None, None)
    
    # specify cell starting positions and orientations
    sim.addCell(cellType=0, rad=rad0, length=len0, pos=(-25,0,0), dir=(0,0,0))
    sim.addCell(cellType=1, rad=rad1, length=len1, pos=(25,0,0), dir=(0,0,0))

    # add some objects to draw the models
    mainRenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(mainRenderer)

    # specify how often we output data to the GUI / to file
    sim.renderEveryNSteps = renderEvery
    sim.pickleSteps = saveEvery
    TimeStep = sim.dt
    print('Timestep:',TimeStep)

def init(cell):
    if cell.cellType == 0:
        cell.targetVol = targetVol0 #+ random.uniform(0.0, 0.09*targetVol0)
        cell.growthRate = g0
        cell.length = len0
        cell.rad = rad0
    if cell.cellType == 1:
        cell.targetVol = targetVol1 #+ random.uniform(0.0, 0.9*targetVol1)
        cell.growthRate = g1
        cell.length = len1
        cell.rad = rad1
        
    cell.color = cell_colors[cell.cellType]


def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    if parent.cellType == 0:
        d1.targetVol = targetVol0 #+ random.uniform(0.0, 0.9*targetVol0)
        d1.length = len0 #+ random.uniform(0.0, 0.9*len0)
        d1.rad = rad0 #+ random.uniform(0.0, 0.9*rad0)
        d1.growthRate = g0
        d2.targetVol = targetVol0 #+ random.uniform(0.0, 0.9*targetVol0)
        d2.length = len0 #+ random.uniform(0.0, 0.9*len0)
        d2.rad = rad0 #+ random.uniform(0.0, 0.9*rad0)
        d2.growthRate = g0
    if parent.cellType == 1:
        d1.targetVol = targetVol1 #+ random.uniform(0.0, 0.9*targetVol1)
        d1.length = len1 #+ random.uniform(0.0, 0.9*len1)
        d1.rad = rad1 #+ random.uniform(0.0, 0.9*rad1)
        d1.growthRate = g1
        d2.targetVol = targetVol1 #+ random.uniform(0.0, 0.9*targetVol1)
        d2.length = len1 #+ random.uniform(0.0, 0.9*len1)
        d2.rad = rad1 #+ random.uniform(0.0, 0.9*rad1)
        d2.growthRate = g1


