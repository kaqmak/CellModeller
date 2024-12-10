# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:15:30 2024

@author: cordermi
"""
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import random

# Simulation population and outputs
saveEvery = 20
renderEvery = 20
max_population = 100000
 
muA = 10

radD = 0.34
lenD = 2.68
initialVolD= lenD
targetVolD = 2*initialVolD
growD = 1

radR = 0.34
lenR = 2.68
initialVolR = lenR
targetVolR = 2*initialVolR
growR = 1

# Do we want the transconjugant to change something? Maybe it should grow slower.
radT = radR
lenT = lenR
initialVolT = lenT
targetVolT = 2*initialVolT
growT = growR

con_prob = 0.8

# cell types and properties
number_of_types = 3
# Type 0 cells correspond to the donor, Green
# Type 1 cells correspond to the recipiemt, Red
# Type 2 cells correspond to the transconjugant, Blue
cell_colors = {0:[0.0, 1.0, 0.0], 1:[1.0, 0.0, 0.0],2:[0.0, 0.0, 1.0]}
cell_lengths = {0:lenD, 1:lenR,2:lenT}
cell_radius = {0:radD, 1:radR, 2:radT}
cell_tarVol = {0:targetVolD, 1:targetVolR, 2:targetVolT}
cell_growthRate = {0:growD, 1:growR, 2:growT}


def randomXY(L,R,ktype):
    R2 = R * R
    L2 = L * L
    if ktype == 1:    ####### inside
        while True:
            x = numpy.random.uniform(-1, 1)
            y = numpy.random.uniform(-1, 1)
            x = x * L
            y = y * L
            r2 = x * x + y * y
            if r2 < R2 :
                break
    else:
        while True:
            x = numpy.random.uniform(-1, 1)
            y = numpy.random.uniform(-1, 1)
            x = x * L
            y = y * L
            r2 = x * x + y * y
            if r2 < L2 and r2 > R2:
                break
    return x,y

def setup(sim):
    global TimeStep

    # Set biophysics and regulation models
    biophys = CLBacterium(sim, max_cells=max_population, jitter_z=False,compNeighbours=True)
    regul = ModuleRegulator(sim)

    # add biophysics, regulation, CDK and EPS objects to simulator
    sim.init(biophys, regul, None, None)

    sim.init(biophys, regul, None, None)
    R1 =35.6
    R2 =35.5
    for i in range(300):
        x1,y1 = randomXY(R1,R2,1)
        x2,y2 = randomXY(R1,R2,2)
        x = numpy.random.uniform(-1, 1)
        y = numpy.random.uniform(-1, 1)
        sim.addCell(cellType=0, rad=radD, length=lenD, pos=(x1,y1,0), dir=(x,y,0))
    for i in range(300):
        x1,y1 = randomXY(R1,R2,1)
        x2,y2 = randomXY(R1,R2,2)
        x = numpy.random.uniform(-1, 1)
        y = numpy.random.uniform(-1, 1)
        sim.addCell(cellType=1, rad=radR, length=lenR, pos=(x1,y1,0), dir=(x,y,0))
   # Add cells into the coffee ring to see the effect
    for i in range(30):
        x1,y1 = randomXY(R1,R2,1)
        x2,y2 = randomXY(R1,R2,2)
        x = numpy.random.uniform(-1, 1)
        y = numpy.random.uniform(-1, 1)
        sim.addCell(cellType=0, rad=radD, length=lenD, pos=(x1,y1,0), dir=(x,y,0))
    for i in range(30):
        x1,y1 = randomXY(R1,R2,1)
        x2,y2 = randomXY(R1,R2,2)
        x = numpy.random.uniform(-1, 1)
        y = numpy.random.uniform(-1, 1)
        sim.addCell(cellType=1, rad=radR, length=lenR, pos=(x1,y1,0), dir=(x,y,0))
    # specify cell starting positions and orientations

    # add some objects to draw the models
    mainRenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(mainRenderer)

    # specify how often we output data to the GUI / to file
    sim.renderEveryNSteps = renderEvery
    sim.pickleSteps = saveEvery
    TimeStep = sim.dt
    print('Timestep:',TimeStep)

def init(cell):
    cell.rad = cell_radius[cell.cellType]
    cell.length = cell_lengths[cell.cellType]
    cell.color = cell_colors[cell.cellType]
    cell.targetVol = cell_tarVol[cell.cellType]
    cell.growthRate = cell_growthRate[cell.cellType]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
    # Add conjugation 
    if cell.cellType == 1: #if a cell is an recipient
       for index in cell.neighbours: #loop through all contacts
           if cells[index].cellType != 1: #if donor or transconjugant is in contact
               if random.random() < con_prob: 
                   cell.cellType = 2 #become transconjugant
    cell.rad = cell_radius[cell.cellType]
    cell.length = cell_lengths[cell.cellType]
    cell.color = cell_colors[cell.cellType]
    cell.targetVol = cell_tarVol[cell.cellType]
    cell.growthRate = cell_growthRate[cell.cellType]

def divide(parent, d1, d2):
    if parent.cellType == 0:
        d1.targetVol = targetVolD 
        d1.length = lenD
        d1.rad = radD
        d1.growthRate = growD
        d2.targetVol = targetVolD
        d2.length = lenD
        d2.rad = radD
        d2.growthRate = growD
    if parent.cellType == 1:
        d1.targetVol = targetVolR
        d1.length = lenR
        d1.rad = radR 
        d1.growthRate = growR
        d2.targetVol = targetVolR 
        d2.length = lenR
        d2.rad = radR
        d2.growthRate = growR
    if parent.cellType == 2:
        d1.targetVol = targetVolT
        d1.length = lenT
        d1.rad = radT 
        d1.growthRate = growT
        d2.targetVol = targetVolT 
        d2.length = lenT
        d2.rad = radT
        d2.growthRate = growT



