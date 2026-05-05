"""
Plant Model for GARISOM-Torch.

This module implements the main Plant class that integrates:
    - Xylem components (root, stem, leaf, soil layers)
    - Hydraulics model (supply/demand)
    - Carbon assimilation (Farquhar photosynthesis)
    - Parameters storage

The Plant class orchestrates the timestep simulation following
Sperry et al. (2016, 2017) gain-risk optimization framework.
"""

import torch
import torch.nn as nn
from typing import Optional, Dict, Tuple, List
import math

from .constants import (
    PI, SBC, GAS, OA, SOLAR, CURVE_MAX, GMAX, MIN_WIND_THRESH,
    VAN_GENUCHTEN_PARAMS, tensor, zeros, get_device, get_dtype
)
from .component import Component
from .soils import SoilLayer, RhizosphereComponent
from .xylem import XylemComponent, LeafComponent, StemComponent, RootComponent
from .carbon import CarbonAssimilationModel
from .hydraulics import HydraulicsModel


class Parameters(nn.Module):
    """
    Storage class for model parameters.
    
    Holds site, atmospheric, soil, stand, hydraulic, and
    photosynthesis parameters needed for the Plant model.
    """
    
    def __init__(self):
        super().__init__()
        
        # Site parameters
        self.lat = nn.Parameter(tensor(0.0), requires_grad=False)  # Latitude (radians)
        self.longitude = nn.Parameter(tensor(0.0), requires_grad=False)  # Longitude (degrees)
        self.slope = nn.Parameter(tensor(0.0), requires_grad=False)  # Slope inclination
        self.slope_asp = nn.Parameter(tensor(0.0), requires_grad=False)  # Slope aspect
        self.alt = nn.Parameter(tensor(0.0), requires_grad=False)  # Elevation (m)
        
        # Atmosphere
        self.tau = nn.Parameter(tensor(0.65), requires_grad=False)  # Atmospheric transmittance
        self.tsn_corr = nn.Parameter(tensor(0.0), requires_grad=False)  # Solar noon correction
        self.ca = nn.Parameter(tensor(400e-6), requires_grad=False)  # Ambient CO2 (mol/mol)
        self.emiss = nn.Parameter(tensor(0.97), requires_grad=False)  # Long-wave emissivity
        self.p_atm = nn.Parameter(tensor(101.325), requires_grad=False)  # Atmospheric pressure (kPa)
        
        # Soil
        self.layers = nn.Parameter(tensor(5.0), requires_grad=False)  # Number of soil layers
        self.rock_frac = nn.Parameter(tensor(0.0), requires_grad=False)  # Rock fraction
        self.rhizo_targ = nn.Parameter(tensor(0.5), requires_grad=False)  # Target rhizosphere resistance fraction
        self.ffc = nn.Parameter(tensor(1.0), requires_grad=False)  # Field capacity fraction
        self.field_cap_frac = nn.Parameter(tensor(0.5), requires_grad=False)  # Field capacity / saturation
        self.p_ground = nn.Parameter(tensor(-0.5), requires_grad=False)  # Ground water pressure (MPa)
        self.ground_distance = nn.Parameter(tensor(10.0), requires_grad=False)  # Distance to GW (m)
        self.soil_abs_sol = nn.Parameter(tensor(0.8), requires_grad=False)  # Soil absorptivity
        
        # Stand parameters
        self.tree_to_photo_lai = nn.Parameter(tensor(1.0), requires_grad=False)
        self.lai = nn.Parameter(tensor(3.0), requires_grad=False)  # Leaf area index
        self.leaf_angle_param = nn.Parameter(tensor(1.0), requires_grad=False)
        self.leaf_per_basal = nn.Parameter(tensor(100.0), requires_grad=False)  # LA per BA
        self.height = nn.Parameter(tensor(10.0), requires_grad=False)  # Tree height (m)
        self.aspect = nn.Parameter(tensor(2.0), requires_grad=False)  # Root radius/depth ratio
        self.root_depth = nn.Parameter(tensor(2.0), requires_grad=False)  # Max rooting depth (m)
        self.beta = nn.Parameter(tensor(0.97), requires_grad=False)  # Root biomass distribution
        self.leaf_width = nn.Parameter(tensor(0.05), requires_grad=False)  # Characteristic dimension (m)
        self.ba_per_ga = nn.Parameter(tensor(0.001), requires_grad=False)  # Basal area per ground area
        self.xh = nn.Parameter(tensor(10.0), requires_grad=False)  # Height above soil for wind (cm)
        
        # Hydraulics
        self.ksatp = nn.Parameter(tensor(10.0), requires_grad=False)  # Plant kmax (kg/hr/m2/MPa)
        self.rsatp = nn.Parameter(tensor(0.1), requires_grad=False)  # Plant resistance
        self.lsc = nn.Parameter(tensor(1.0), requires_grad=False)  # Leaf specific conductance
        self.k_sat_root = nn.Parameter(tensor(20.0), requires_grad=False)  # Root kmax
        self.p_inc = nn.Parameter(tensor(0.00075), requires_grad=False)  # Pressure increment
        self.k_min = nn.Parameter(tensor(0.005), requires_grad=False)  # Minimum conductance
        self.p_grav = nn.Parameter(tensor(0.1), requires_grad=False)  # Gravity pressure drop
        self.e_inc = nn.Parameter(tensor(0.02), requires_grad=False)  # E increment
        self.g_max = nn.Parameter(tensor(GMAX), requires_grad=False)  # Max conductance
        self.g_max_l = nn.Parameter(tensor(500.0), requires_grad=False)  # Max g per leaf area
        
        # Photosynthesis (Farquhar model)
        self.light_comp = nn.Parameter(tensor(10.0), requires_grad=False)  # Light compensation (µmol/m2/s)
        self.q_max = nn.Parameter(tensor(0.3), requires_grad=False)  # Quantum yield
        self.v_max25 = nn.Parameter(tensor(60.0), requires_grad=False)  # Vmax at 25°C
        self.j_max25 = nn.Parameter(tensor(120.0), requires_grad=False)  # Jmax at 25°C
        self.kc_25 = nn.Parameter(tensor(404.9e-6), requires_grad=False)  # Kc at 25°C
        self.ko_25 = nn.Parameter(tensor(278.4e-3), requires_grad=False)  # Ko at 25°C
        self.comp_25 = nn.Parameter(tensor(42.75e-6), requires_grad=False)  # Gamma at 25°C
        self.theta_c = nn.Parameter(tensor(0.98), requires_grad=False)  # Colimitation shape
        self.ha_vmax = nn.Parameter(tensor(65330.0), requires_grad=False)  # Ha for Vmax
        self.hd_vmax = nn.Parameter(tensor(200000.0), requires_grad=False)  # Hd for Vmax
        self.sv_vmax = nn.Parameter(tensor(485.0), requires_grad=False)  # Sv for Vmax
        self.ha_jmax = nn.Parameter(tensor(43540.0), requires_grad=False)  # Ha for Jmax
        self.hd_jmax = nn.Parameter(tensor(200000.0), requires_grad=False)  # Hd for Jmax
        self.sv_jmax = nn.Parameter(tensor(495.0), requires_grad=False)  # Sv for Jmax
        self.light_curv = nn.Parameter(tensor(0.9), requires_grad=False)  # Light response curvature
        
    def set_from_dict(self, params: Dict[str, float]):
        """Set parameters from dictionary."""
        for key, value in params.items():
            if hasattr(self, key):
                param = getattr(self, key)
                if isinstance(param, nn.Parameter):
                    param.data = tensor(value)
                    
    def get(self, name: str) -> torch.Tensor:
        """Get parameter value by name."""
        if hasattr(self, name):
            return getattr(self, name)
        raise KeyError(f"Parameter '{name}' not found")


class Plant(nn.Module):
    """
    Main Plant model integrating all components.
    
    The Plant class orchestrates:
        1. Xylem hydraulics (root, stem, leaf conductances)
        2. Soil-plant-atmosphere continuum
        3. Carbon assimilation (Farquhar model)
        4. Stomatal optimization (gain-risk framework)
        
    The model timestep iteration follows Sperry et al. (2016, 2017).
    
    Attributes:
        params: Model parameters
        xylem: Xylem component with soil layers
        carbon: Carbon assimilation model
        hydraulics: SPAC hydraulics model
        num_layers: Number of soil layers
    """
    
    def __init__(self, num_layers: int = 5, texture: str = "loam"):
        """
        Initialize Plant model.
        
        Args:
            num_layers: Number of soil layers
            texture: Soil texture for Van Genuchten parameters
        """
        super().__init__()
        
        self.num_layers = num_layers
        self.texture = texture
        
        # Initialize components
        self.params = Parameters()
        self.xylem = XylemComponent(num_layers=num_layers, texture=texture)
        self.carbon = CarbonAssimilationModel()
        self.hydraulics = HydraulicsModel()
        
        # State variables (registered as buffers for persistence)
        self.register_buffer('predawn', tensor(-0.5))  # Initial predawn (MPa)
        self.register_buffer('midday', tensor(-1.0))
        self.register_buffer('transpiration', tensor(0.0))
        self.register_buffer('assimilation', tensor(0.0))
        self.register_buffer('conductance', tensor(0.0))
        self.register_buffer('ecrit', tensor(0.0))
        self.register_buffer('pcrit', tensor(0.0))
        
        # Water content per layer - initialize at field capacity (~0.3 m³/m³)
        water_init = torch.full((num_layers + 1,), 0.30)
        self.register_buffer('water', water_init)
        
        # Field capacity
        fc_init = torch.full((num_layers + 1,), 0.35)
        self.register_buffer('fc', fc_init)
        
        # E(P) curve storage
        self.register_buffer('e_p', zeros(CURVE_MAX))
        
        # Fatigue tracking for hysteresis
        self.register_buffer('b_fatigue', zeros(3, 10))  # [component, year]
        self.register_buffer('max_plc', tensor(0.0))
        
    def set_parameters(self, params_dict: Dict[str, float]):
        """
        Set model parameters from dictionary.
        
        Args:
            params_dict: Dictionary of parameter names and values
        """
        self.params.set_from_dict(params_dict)
        
        # Update derived quantities
        self._update_derived_params()
        
    def _update_derived_params(self):
        """Update derived parameters after setting base parameters."""
        # Atmospheric pressure from elevation
        alt = self.params.alt
        p_atm = 101.325 * torch.pow(
            1 - 0.0065 * alt / (288.15 + 0.0065 * alt), 
            5.257
        )
        self.params.p_atm.data = p_atm
        
        # Gravity pressure drop from height
        p_grav = self.params.height * 0.01
        self.params.p_grav.data = p_grav
        
        # E increment from ksatp
        e_inc = self.params.ksatp / 500.0
        self.params.e_inc.data = e_inc
        
        # K minimum cutoff
        k_min = self.params.ksatp / 2000.0
        self.params.k_min.data = k_min
        
    def solar_geometry(self, 
                       julian_day: torch.Tensor,
                       hour: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Calculate solar geometry for radiation partitioning.
        
        Based on Campbell & Norman equations for solar position.
        
        Args:
            julian_day: Julian day of year (1-365)
            hour: Hour of day (0-24, local standard time)
            
        Returns:
            Tuple of (solar_altitude, direct_beam_fraction)
        """
        lat = self.params.lat
        
        # Solar declination (radians)
        delta = 0.4102 * torch.sin(2 * PI * (julian_day - 80) / 365)
        
        # Hour angle
        tsn = 12.0 + self.params.tsn_corr  # Solar noon
        hour_angle = PI * (hour - tsn) / 12.0
        
        # Solar altitude
        sin_alt = (torch.sin(lat) * torch.sin(delta) + 
                   torch.cos(lat) * torch.cos(delta) * torch.cos(hour_angle))
        solar_alt = torch.asin(torch.clamp(sin_alt, -1.0, 1.0))
        
        # Direct beam fraction (simple approximation)
        tau = self.params.tau
        air_mass = 1.0 / (torch.clamp(sin_alt, min=0.01))
        direct_frac = torch.pow(tau, air_mass)
        
        return solar_alt, direct_frac
    
    def partition_radiation(self,
                            solar: torch.Tensor,
                            solar_alt: torch.Tensor,
                            direct_frac: torch.Tensor,
                            lai: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Partition solar radiation into sunlit and shaded fractions.
        
        Based on Campbell & Norman two-stream approximation.
        
        Args:
            solar: Total incident solar radiation (W/m²)
            solar_alt: Solar altitude angle (radians)
            direct_frac: Fraction of direct beam radiation
            lai: Leaf area index
            
        Returns:
            Tuple of (PPFD_sunlit, PPFD_shaded) in µmol/m²/s
        """
        # Convert W/m² to PPFD (approximate: 4.6 µmol/J for PAR)
        ppfd = solar * 4.6 * 0.5  # 50% of solar is PAR
        
        # Beam and diffuse components
        beam = ppfd * direct_frac
        diffuse = ppfd * (1 - direct_frac)
        
        # Extinction coefficient for spherical leaf angle distribution
        sin_alt = torch.sin(torch.clamp(solar_alt, min=0.01))
        kb = 0.5 / sin_alt  # Beam extinction
        
        # Sunlit and shaded LAI fractions
        lai_sun = (1 - torch.exp(-kb * lai)) / kb
        lai_shade = lai - lai_sun
        
        # Sunlit receives both beam and diffuse
        ppfd_sun = beam + diffuse * 0.5  # Simplified
        
        # Shaded receives only diffuse (attenuated)
        ppfd_shade = diffuse * 0.2  # Simplified
        
        return ppfd_sun, ppfd_shade, lai_sun, lai_shade
    
    def vapor_pressure_deficit(self,
                                tair: torch.Tensor,
                                vpd_input: Optional[torch.Tensor] = None,
                                rh: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Calculate or validate vapor pressure deficit.
        
        Args:
            tair: Air temperature (°C)
            vpd_input: Measured VPD (kPa), optional
            rh: Relative humidity (0-1), optional
            
        Returns:
            VPD in mole fraction
        """
        # Saturation vapor pressure (Tetens equation)
        esat = 0.611 * torch.exp(17.502 * tair / (tair + 240.97))
        max_vpd = esat  # Maximum possible VPD
        
        if vpd_input is not None:
            vpd = vpd_input
        elif rh is not None:
            vpd = esat * (1 - rh)
        else:
            vpd = esat * 0.5  # Default 50% RH
            
        # Clamp to maximum
        vpd = torch.clamp(vpd, min=0.0, max=max_vpd)
        
        # Convert to mole fraction
        vpd_mol = vpd / self.params.p_atm
        
        return vpd_mol
    
    def get_predawns(self) -> torch.Tensor:
        """
        Calculate predawn water potential from soil water content.
        
        Returns weighted average of soil layer water potentials.
        
        Returns:
            Predawn water potential (MPa)
        """
        predawn = tensor(0.0)
        total_weight = tensor(0.0)
        
        for k in range(1, self.num_layers + 1):
            soil_layer = self.xylem.soils[k]
            
            # Get soil water potential from water content
            psi = soil_layer.rhizosphere.pressure_from_theta(self.water[k])
            
            # Weight by root conductance
            weight = soil_layer.root.k_max
            
            predawn = predawn + psi * weight
            total_weight = total_weight + weight
            
        predawn = predawn / (total_weight + 1e-10)
        
        return predawn
    
    def composite_curve(self,
                        predawn: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, int]:
        """
        Calculate the E(P) composite supply curve.
        
        Integrates conductance from predawn to critical pressure,
        building the supply function for the whole plant.
        
        Args:
            predawn: Predawn water potential (MPa)
            
        Returns:
            Tuple of (E_crit, P_crit, curve_length)
        """
        p_inc = self.params.p_inc
        k_min = self.params.k_min
        e_inc = self.params.e_inc
        
        # Initialize curve
        e = tensor(0.0)
        p = 0
        
        # Build E(P) curve by incrementing E
        max_iters = min(CURVE_MAX - 1, 1000)
        
        for p in range(max_iters):
            e = e + e_inc
            
            # Get whole-plant conductance at this E
            # Simplified: use Weibull vulnerability
            pressure = predawn + e / (self.params.ksatp + 1e-10)
            
            k_plant = self.xylem.whole_plant_conductance(pressure)
            
            self.e_p[p] = e
            
            # Check for critical point
            if k_plant < k_min:
                break
                
        ecrit = e
        pcrit = pressure
        
        self.ecrit = ecrit
        self.pcrit = pcrit
        
        return ecrit, pcrit, p
    
    def canopy_pressure_optimization(self,
                                      predawn: torch.Tensor,
                                      vpd: torch.Tensor,
                                      ppfd_sun: torch.Tensor,
                                      ppfd_shade: torch.Tensor,
                                      tleaf: torch.Tensor,
                                      lai_sun: torch.Tensor,
                                      lai_shade: torch.Tensor
                                      ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Find optimal canopy pressure maximizing gain minus risk.
        
        Implements the Sperry et al. (2017) gain-risk optimization:
            max(A' - θ') where A' = A/Amax, θ' = 1 - K/Kmax
            
        This version uses a differentiable soft-argmax to enable
        gradient-based optimization of plant parameters.
            
        Args:
            predawn: Predawn water potential (MPa)
            vpd: Vapor pressure deficit (mol fraction)
            ppfd_sun: PPFD for sunlit leaves (µmol/m²/s)
            ppfd_shade: PPFD for shaded leaves (µmol/m²/s)
            tleaf: Leaf temperature (°C)
            lai_sun: Sunlit LAI
            lai_shade: Shaded LAI
            
        Returns:
            Tuple of (E_optimal, P_optimal, A_optimal, G_optimal)
        """
        # Get parameters
        ca = self.params.ca
        g_max = self.params.g_max_l
        ksatp = self.params.ksatp
        e_inc = self.params.e_inc
        vmax25 = self.params.v_max25
        jmax25 = self.params.j_max25
        
        # Number of discrete steps for optimization (reduced for speed)
        n_steps = 20
        
        # Build E vector from 0 to critical E
        # Use a range based on expected maximum E
        e_max = ksatp * 0.5  # Approximate maximum E
        e_vals = torch.linspace(0.0, e_max.item(), n_steps, device=get_device(), dtype=get_dtype())
        
        # Pressure at each E
        pressure_vals = predawn + e_vals / (ksatp + 1e-10)
        
        # Get maximum assimilation (at maximum conductance)
        gcanc_max = (g_max / 1.6) * 1000.0
        a_max_sun = self.carbon.forward(tleaf, ppfd_sun, gcanc_max, ca, vmax25, jmax25)
        a_max_shade = self.carbon.forward(tleaf, ppfd_shade, gcanc_max, ca, vmax25, jmax25)
        a_max = (a_max_sun * lai_sun + a_max_shade * lai_shade) / (lai_sun + lai_shade + 1e-10)
        
        # Compute profit at each E (vectorized)
        profits = []
        a_vals = []
        g_vals = []
        
        for i in range(n_steps):
            e = e_vals[i]
            pressure = pressure_vals[i]
            
            # Conductance at this pressure
            k = self.xylem.whole_plant_conductance(pressure)
            
            # Stomatal conductance from E and VPD
            g = e / (vpd + 1e-10)
            g = torch.clamp(g, max=g_max)
            
            # Convert to CO2 conductance
            gcanc = (g / 1.6) * 1000.0
            
            # Assimilation at this g
            a_sun = self.carbon.forward(tleaf, ppfd_sun, gcanc, ca, vmax25, jmax25)
            a_shade = self.carbon.forward(tleaf, ppfd_shade, gcanc, ca, vmax25, jmax25)
            a = (a_sun * lai_sun + a_shade * lai_shade) / (lai_sun + lai_shade + 1e-10)
            
            # Normalized gain and risk
            a_prime = a / (a_max + 1e-10)
            theta_prime = 1 - k / (ksatp + 1e-10)
            
            # Profit = gain - risk
            profit = a_prime - theta_prime
            
            profits.append(profit)
            a_vals.append(a)
            g_vals.append(g)
        
        # Stack into tensors
        profits = torch.stack(profits)
        a_vals = torch.stack(a_vals)
        g_vals = torch.stack(g_vals)
        
        # Differentiable soft-argmax: use softmax weights to select optimal point
        # Temperature parameter controls sharpness (higher = closer to hard argmax)
        temperature = 10.0
        weights = torch.softmax(profits * temperature, dim=0)
        
        # Weighted sum gives differentiable approximation to argmax
        best_e = torch.sum(weights * e_vals)
        best_a = torch.sum(weights * a_vals)
        best_g = torch.sum(weights * g_vals)
        best_p = predawn + best_e / (ksatp + 1e-10)
        
        return best_e, best_p, best_a, best_g
    
    def forward(self,
                solar: torch.Tensor,
                vpd: torch.Tensor,
                tair: torch.Tensor,
                wind: torch.Tensor,
                julian_day: torch.Tensor,
                hour: torch.Tensor,
                soil_water: Optional[torch.Tensor] = None
                ) -> Dict[str, torch.Tensor]:
        """
        Forward pass: compute plant gas exchange for one timestep.
        
        This is the main entry point for model predictions. Given
        environmental forcing, returns transpiration, assimilation,
        and canopy pressure.
        
        Args:
            solar: Incident solar radiation (W/m²)
            vpd: Vapor pressure deficit (kPa)
            tair: Air temperature (°C)
            wind: Wind speed (m/s)
            julian_day: Julian day of year
            hour: Hour of day (0-24)
            soil_water: Soil water content per layer (optional)
            
        Returns:
            Dictionary with model outputs:
                - transpiration: E (kg/hr/m² BA)
                - assimilation: A (µmol/m²/s leaf area)
                - conductance: G (mmol/m²/s)
                - predawn: Predawn pressure (MPa)
                - midday: Midday pressure (MPa)
                - leaf_temp: Leaf temperature (°C)
        """
        # Minimum wind threshold
        wind = torch.clamp(wind, min=MIN_WIND_THRESH)
        
        # Update soil water if provided
        if soil_water is not None:
            self.water = soil_water
            
        # Solar geometry
        solar_alt, direct_frac = self.solar_geometry(julian_day, hour)
        
        # Check if daytime
        lai = self.params.lai
        is_day = solar_alt > 0.0
        
        # Partition radiation
        ppfd_sun, ppfd_shade, lai_sun, lai_shade = self.partition_radiation(
            solar, solar_alt, direct_frac, lai
        )
        
        # Check light compensation
        light_comp = self.params.light_comp
        is_active = ppfd_sun > light_comp
        
        # VPD in mole fraction
        vpd_mol = self.vapor_pressure_deficit(tair, vpd)
        
        # Get predawn
        predawn = self.get_predawns()
        
        # Leaf temperature (simplified: assume close to air temp initially)
        tleaf = tair + 2.0  # Slight warming due to radiation
        
        if is_active:
            # Optimize canopy pressure
            e_opt, p_opt, a_opt, g_opt = self.canopy_pressure_optimization(
                predawn, vpd_mol, ppfd_sun, ppfd_shade, tleaf, lai_sun, lai_shade
            )
            
            # Update leaf temperature with energy balance
            tleaf = self.xylem.leaf.get_leaf_temp_simple(
                tair, e_opt, vpd, wind, self.params.leaf_width, self.params.p_atm
            )
            
        else:
            # Night or too dark - stomata closed
            e_opt = tensor(0.0)
            p_opt = predawn
            a_opt = tensor(0.0)  # Respiration would be negative
            g_opt = tensor(0.0)
            
        # Update state
        self.predawn = predawn
        self.midday = p_opt
        self.transpiration = e_opt
        self.assimilation = a_opt
        self.conductance = g_opt
        
        return {
            'transpiration': e_opt,
            'assimilation': a_opt,
            'conductance': g_opt,
            'predawn': predawn,
            'midday': p_opt,
            'leaf_temp': tleaf,
            'ppfd_sun': ppfd_sun,
            'ppfd_shade': ppfd_shade,
            'lai_sun': lai_sun,
            'lai_shade': lai_shade,
            'vpd': vpd_mol * self.params.p_atm,  # Back to kPa
            'is_active': is_active,
        }
    
    def update_soil_water(self,
                          transpiration: torch.Tensor,
                          timestep_hours: float = 1.0):
        """
        Update soil water content after transpiration.
        
        Simple bucket model for soil water balance.
        
        Args:
            transpiration: Total transpiration (kg/hr/m² BA)
            timestep_hours: Timestep duration in hours
        """
        # Convert transpiration to volume per ground area
        ba_per_ga = self.params.ba_per_ga
        vol_water = transpiration * timestep_hours * ba_per_ga / 1000.0  # m³/m²
        
        # Distribute extraction according to root distribution
        total_k = tensor(0.0)
        for k in range(1, self.num_layers + 1):
            total_k = total_k + self.xylem.soils[k].root.k_max
            
        for k in range(1, self.num_layers + 1):
            layer = self.xylem.soils[k]
            frac = layer.root.k_max / (total_k + 1e-10)
            
            # Extract water from layer
            extract = vol_water * frac / layer.depth
            self.water[k] = torch.clamp(self.water[k] - extract, min=0.01)
            
    def __repr__(self) -> str:
        return (f"Plant(layers={self.num_layers}, texture='{self.texture}', "
                f"lai={self.params.lai.item():.1f}, height={self.params.height.item():.1f}m)")
