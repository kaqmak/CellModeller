import random

from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium

# Import Euler integrator for solving ODE system of chemical species inside the cells
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator

max_cells = 2**15


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)

    integ = CLEulerIntegrator(sim, 2, max_cells)

    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, None, integ)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0, 0, 0))

    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers

        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 10


def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 2.5 + random.uniform(0.0, 0.5)
    # Specify growth rate of cells
    cell.growthRate = 2.0
    # Specify initial concentration of chemical species
    cell.species[:] = [10, 0]


def specRateCL():  # Add
    return """
    const float k1 = 0.f;
    const float k2 = 1.f; 
    const float k3 = 2.f;
    float x0 = species[0];
    float x1 = species[1];
    rates[0] = k1;
    rates[1] = k2*k3*k3/(k3*k3 + x0*x0);
   
    """
    # k1 = production rate of x0
    #


def update(cells):
    # Iterate through each cell and flag cells that reach target size for division
    for id, cell in cells.items():
        cell.color = [
            0.1 + cell.species[0],
            0.1 + cell.species[1] * 0.1,
            0.1,
        ]  # Add/change
        if cell.volume > cell.targetVol:
            cell.divideFlag = True


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 2.5 + random.uniform(0.0, 0.5)
    d2.targetVol = 2.5 + random.uniform(0.0, 0.5)
