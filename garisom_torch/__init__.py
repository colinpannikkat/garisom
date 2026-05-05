"""
GARISOM-Torch: PyTorch Implementation of the Carbon Gain vs Hydraulic Risk 
Stomatal Optimization Model

A fully differentiable implementation of the GARISOM plant physiology model
using PyTorch, enabling gradient-based optimization and sensitivity analysis.

Original C++ model by German Vargas G. (german.vargas@utah.edu)
PyTorch implementation for differentiable plant modeling.

Model Components:
    - Soils: Van Genuchten soil hydraulics and rhizosphere
    - Hydraulics: Plant hydraulic conductance and water transport
    - Carbon: Farquhar-type photosynthesis model
    - Xylem: Root, stem, and leaf components with Weibull vulnerability
    - Plant: Main model integration and timestep iteration
"""

from .constants import *
from .component import Component
from .soils import RhizosphereComponent, SoilLayer
from .xylem import RootComponent, StemComponent, LeafComponent, XylemComponent
from .carbon import CarbonAssimilationModel
from .hydraulics import HydraulicsModel
from .plant import Plant, Parameters
from .data_utils import (
    load_parameters, 
    load_climate_data, 
    load_configuration,
    ClimateData,
    ModelConfig,
    save_results,
    create_data_loader,
)

__version__ = "3.0.0-torch"
__author__ = "PyTorch port of GARISOM by German Vargas G."

__all__ = [
    # Core model
    'Plant',
    'Parameters',
    
    # Components
    'Component',
    'RhizosphereComponent',
    'SoilLayer',
    'RootComponent',
    'StemComponent',
    'LeafComponent',
    'XylemComponent',
    'CarbonAssimilationModel',
    'HydraulicsModel',
    
    # Data utilities
    'load_parameters',
    'load_climate_data',
    'load_configuration',
    'ClimateData',
    'ModelConfig',
    'save_results',
    'create_data_loader',
    
    # Constants
    'PI',
    'SBC',
    'GAS',
    'OA',
    'SOLAR',
    'CURVE_MAX',
    'GMAX',
    'VAN_GENUCHTEN_PARAMS',
]
