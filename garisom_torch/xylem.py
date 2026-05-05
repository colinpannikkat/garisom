"""
Xylem components for GARISOM-Torch.

This module implements the xylem transport system including root, stem,
and leaf components. Each component uses Weibull vulnerability curves
to model the relationship between xylem pressure and hydraulic conductance.

The xylem system forms the water transport pathway:
    Soil -> Rhizosphere -> Roots -> Stem -> Leaves -> Atmosphere
"""

import torch
import torch.nn as nn
from typing import List, Optional, Tuple
from .constants import (
    CURVE_MAX, MAX_LAYERS, PI, SBC, SHA, ABS_SOLAR, ABS_PAR, ABS_NIR,
    MIN_WIND_THRESH, tensor, zeros, get_device, get_dtype
)
from .component import Component


class RootComponent(Component):
    """
    Root hydraulic component.
    
    Roots connect the rhizosphere to the stem, transporting water
    from each soil layer. Root vulnerability is modeled using
    Weibull curves specific to the species.
    """
    
    def __init__(self):
        super().__init__()
        # Default root Weibull parameters (species-specific)
        self.b_wb.data = tensor(1.0)
        self.c_wb.data = tensor(2.0)


class StemComponent(Component):
    """
    Stem hydraulic component.
    
    The stem transports water from the root crown to the leaves.
    Stem vulnerability is typically less than roots due to 
    wider vessels and lower cavitation risk.
    """
    
    def __init__(self):
        super().__init__()
        # Default stem Weibull parameters
        self.b_wb.data = tensor(2.0)
        self.c_wb.data = tensor(3.0)


class LeafComponent(Component):
    """
    Leaf hydraulic component with energy balance calculations.
    
    Leaves are the terminus of the water transport pathway and
    the site of transpiration and photosynthesis. This component
    includes leaf energy balance for temperature calculations.
    
    Additional attributes for leaf-specific processes:
        - Leaf temperature calculations
        - Vapor pressure deficit at leaf surface
        - Radiation absorption
    """
    
    def __init__(self):
        super().__init__()
        
        # Default leaf Weibull parameters (typically most vulnerable)
        self.b_wb.data = tensor(1.5)
        self.c_wb.data = tensor(3.0)
        
        # Leaf area parameters
        self.lai = nn.Parameter(tensor(2.0), requires_grad=False)  # Total LAI
        self.laisl = nn.Parameter(tensor(1.0), requires_grad=False)  # Sunlit LAI
        self.laish = nn.Parameter(tensor(1.0), requires_grad=False)  # Shaded LAI
        
        # Radiation parameters
        self.emiss = nn.Parameter(tensor(0.97), requires_grad=False)  # Emissivity
        
        # Solar calculations
        self.ssun = nn.Parameter(tensor(0.0), requires_grad=False)
        self.sshade = nn.Parameter(tensor(0.0), requires_grad=False)
        self.sref = nn.Parameter(tensor(0.0), requires_grad=False)
        self.sbottom = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # Boundary layer parameters
        self.la = nn.Parameter(tensor(0.0), requires_grad=False)  # Longwave from atmosphere
        self.lg = nn.Parameter(tensor(0.0), requires_grad=False)  # Longwave from ground
        self.lambda_val = nn.Parameter(tensor(44000.0), requires_grad=False)  # Latent heat
        self.grad = nn.Parameter(tensor(0.0), requires_grad=False)  # Radiation gradient
        self.gha = nn.Parameter(tensor(0.0), requires_grad=False)  # Heat transfer coeff
        
        # Transpiration and temperature curves
        self.register_buffer('lavpd', zeros(CURVE_MAX))  # Leaf-air VPD
        self.register_buffer('lavpdsh', zeros(CURVE_MAX))  # Shade leaf VPD
        self.register_buffer('leaftemp', zeros(CURVE_MAX))  # Leaf temperature
        self.register_buffer('leaftempsh', zeros(CURVE_MAX))  # Shade leaf temp
        self.register_buffer('eplantl', zeros(CURVE_MAX))  # Transpiration per leaf area
        
        # Midday values
        self.emd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.lavpdmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.lavpdshmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.leaftmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.leaftshmd = nn.Parameter(tensor(0.0), requires_grad=False)
    
    def calculate_leaf_temp(self, p: int, airtemp: torch.Tensor,
                            eplant: torch.Tensor, vpd: torch.Tensor,
                            wind: torch.Tensor, laperba: torch.Tensor,
                            leafwidth: torch.Tensor, patm: torch.Tensor) -> None:
        """
        Calculate leaf temperature from energy balance.
        
        Uses iterative energy balance to find leaf temperature where
        energy inputs (radiation) equal outputs (convection + transpiration).
        
        The energy balance equation:
            Rabs = H + λE + εσT⁴
        
        where:
            Rabs = Absorbed radiation
            H = Sensible heat flux
            λE = Latent heat flux (transpiration)
            ε = Emissivity
            σ = Stefan-Boltzmann constant
            T = Leaf temperature
        
        Args:
            p: Pressure index
            airtemp: Air temperature (°C)
            eplant: Plant transpiration rate (mmol m⁻² s⁻¹)
            vpd: Vapor pressure deficit (kPa)
            wind: Wind speed (m s⁻¹)
            laperba: Leaf area per basal area
            leafwidth: Characteristic leaf width (m)
            patm: Atmospheric pressure (kPa)
        """
        # Ensure minimum wind speed
        wind = torch.clamp(wind, min=MIN_WIND_THRESH)
        
        # Convert transpiration to energy flux
        # eplant is in mmol m⁻² s⁻¹ (leaf area)
        # Convert to kg m⁻² s⁻¹: mmol * 0.018 g/mmol * 0.001 kg/g = kg m⁻² s⁻¹
        e_kg = eplant * 0.018 * 0.001  # kg m⁻² s⁻¹
        
        # Boundary layer conductance for heat (mol m⁻² s⁻¹)
        # gHa = 0.135 * sqrt(wind/d) where d is characteristic dimension
        gha = 0.135 * torch.sqrt(wind / (leafwidth + 1e-6))
        
        # Initial guess: leaf temp = air temp
        tleaf = airtemp.clone()
        
        # Absorbed radiation (stored from solar calculations)
        rabs = self.ssun  # W m⁻²
        
        # Iterate to find equilibrium leaf temperature
        for _ in range(10):
            # Saturation vapor pressure at leaf temperature (kPa)
            esat_leaf = 0.611 * torch.exp(17.502 * tleaf / (tleaf + 240.97))
            
            # Actual vapor pressure in air
            eair = esat_leaf - vpd
            eair = torch.clamp(eair, min=0.01)
            
            # Leaf-to-air VPD
            lavpd_val = esat_leaf - eair
            lavpd_val = torch.clamp(lavpd_val, min=0.0)
            
            # Sensible heat flux: H = gHa * Cp * (Tleaf - Tair)
            # Cp for air ≈ 29.3 J mol⁻¹ K⁻¹
            H = gha * SHA * (tleaf - airtemp)  # W m⁻²
            
            # Latent heat flux: λE
            lambda_E = e_kg * self.lambda_val  # W m⁻²
            
            # Longwave emission: εσT⁴
            tk = tleaf + 273.15
            lw_emit = self.emiss * SBC * torch.pow(tk, 4)
            
            # Longwave absorbed (from sky and ground)
            lw_abs = self.la + self.lg
            
            # Energy balance residual
            # Rabs + LWabs = H + λE + LWemit
            residual = rabs + lw_abs - H - lambda_E - lw_emit
            
            # Derivative of residual with respect to Tleaf
            # d(H)/dT = gHa * Cp
            # d(LWemit)/dT = 4εσT³
            dH_dT = gha * SHA
            dLW_dT = 4.0 * self.emiss * SBC * torch.pow(tk, 3)
            
            dR_dT = -dH_dT - dLW_dT
            
            # Newton-Raphson update
            if torch.abs(dR_dT) > 1e-10:
                delta_T = residual / (-dR_dT)
                delta_T = torch.clamp(delta_T, -5.0, 5.0)  # Limit step size
                tleaf = tleaf + delta_T
            
            # Check convergence
            if torch.abs(residual) < 0.1:
                break
        
        # Store results
        self.leaftemp[p] = tleaf
        
        # Calculate final VPD at leaf surface
        esat_leaf = 0.611 * torch.exp(17.502 * tleaf / (tleaf + 240.97))
        eair = esat_leaf - vpd
        self.lavpd[p] = (esat_leaf - eair) / patm  # Normalized by atmospheric pressure
    
    def get_leaf_temp_simple(self, airtemp: torch.Tensor,
                              eplant: torch.Tensor, vpd: torch.Tensor,
                              wind: torch.Tensor, leafwidth: torch.Tensor,
                              patm: torch.Tensor) -> torch.Tensor:
        """
        Calculate leaf temperature (simplified, returns value directly).
        
        A simplified energy balance for quick computation.
        
        Args:
            airtemp: Air temperature (°C)
            eplant: Plant transpiration rate
            vpd: Vapor pressure deficit (kPa)
            wind: Wind speed (m s⁻¹)
            leafwidth: Characteristic leaf width (m)
            patm: Atmospheric pressure (kPa)
            
        Returns:
            Leaf temperature (°C)
        """
        # Ensure minimum wind speed
        wind = torch.clamp(wind, min=MIN_WIND_THRESH)
        
        # Boundary layer conductance for heat
        gha = 0.135 * torch.sqrt(wind / (leafwidth + 1e-6))
        
        # Convert transpiration to latent heat flux
        # Assume eplant is in reasonable units
        e_kg = eplant * 0.018 * 0.001 / 3600.0  # Rough conversion
        lambda_E = e_kg * 44000.0  # Latent heat of vaporization
        
        # Simplified: leaf temp = air temp + radiation heating - evaporative cooling
        # Delta_T ≈ (Rabs - λE) / (gha * Cp + 4εσT³)
        
        tk = airtemp + 273.15
        denominator = gha * SHA + 4.0 * self.emiss * SBC * torch.pow(tk, 3)
        
        # Assume absorbed radiation proportional to incoming solar (stored in ssun)
        rabs = self.ssun if hasattr(self, 'ssun') and self.ssun > 0 else tensor(200.0)
        
        delta_t = (rabs - lambda_E) / (denominator + 1e-6)
        delta_t = torch.clamp(delta_t, -10.0, 15.0)
        
        tleaf = airtemp + delta_t
        
        return tleaf
    
    def calculate_shade_leaf_temp(self, p: int, airtemp: torch.Tensor,
                                   patm: torch.Tensor, vpd: torch.Tensor) -> None:
        """
        Calculate shade leaf temperature (simplified version).
        
        Shade leaves have lower radiation load and are typically
        closer to air temperature.
        """
        # Shade leaves are closer to air temperature
        # Use absorbed shade radiation
        rabs = self.sshade
        
        # Simpler calculation for shade
        # Assume small deviation from air temp
        delta_t = rabs / (self.gha * SHA + 4.0 * self.emiss * SBC * 
                         torch.pow(airtemp + 273.15, 3) + 1e-6)
        delta_t = torch.clamp(delta_t, -10.0, 10.0)
        
        tleaf = airtemp + delta_t * 0.3  # Damped response
        
        self.leaftempsh[p] = tleaf
        
        # VPD at shade leaf
        esat_leaf = 0.611 * torch.exp(17.502 * tleaf / (tleaf + 240.97))
        eair = esat_leaf - vpd
        self.lavpdsh[p] = (esat_leaf - eair) / patm
    
    def clean_parameters(self) -> None:
        """Reset all parameters including leaf-specific ones."""
        super().clean_parameters()
        
        self.lai.data.fill_(2.0)
        self.laisl.data.fill_(1.0)
        self.laish.data.fill_(1.0)
        
        self.ssun.data.fill_(0.0)
        self.sshade.data.fill_(0.0)
        self.sref.data.fill_(0.0)
        self.sbottom.data.fill_(0.0)
        
        self.la.data.fill_(0.0)
        self.lg.data.fill_(0.0)
        self.grad.data.fill_(0.0)
        self.gha.data.fill_(0.0)
        
        self.emd.data.fill_(0.0)
        self.lavpdmd.data.fill_(0.0)
        self.lavpdshmd.data.fill_(0.0)
        self.leaftmd.data.fill_(0.0)
        self.leaftshmd.data.fill_(0.0)
        
        self.lavpd.zero_()
        self.lavpdsh.zero_()
        self.leaftemp.zero_()
        self.leaftempsh.zero_()
        self.eplantl.zero_()


class XylemComponent(nn.Module):
    """
    Complete xylem system combining all plant hydraulic elements.
    
    The xylem component integrates:
        - Multiple soil layers with roots and rhizosphere
        - Single stem element
        - Single leaf element
        
    It computes the composite E(P) supply curve by solving for
    pressures throughout the hydraulic pathway.
    
    The model follows Sperry et al. (2016) where:
        - Root and rhizosphere in each layer are in series
        - All layers connect in parallel to root crown
        - Stem connects root crown to leaves
        
    Attributes:
        leaf: Leaf hydraulic component
        stem: Stem hydraulic component
        soils: List of soil layers with root+rhizosphere pairs
        num_layers: Number of active soil layers
    """
    
    def __init__(self, num_layers: int = 5, texture: str = "loam"):
        super().__init__()
        
        self.num_layers = num_layers
        
        # Main xylem components
        self.leaf = LeafComponent()
        self.stem = StemComponent()
        
        # Soil layers (0 is surface layer without roots)
        # Import here to avoid circular dependency
        from .soils import SoilLayer
        self.soils = nn.ModuleList([
            SoilLayer(texture) for _ in range(num_layers + 1)
        ])
        
        # Set surface layer (0) to have no roots
        self.soils[0].root.k_max.data.fill_(0.0)
        
        # Composite curves
        self.register_buffer('e_p', zeros(CURVE_MAX))  # Whole plant E(P)
        self.register_buffer('k_curve', zeros(CURVE_MAX))  # Whole plant K(P)
        self.register_buffer('root_pressure', zeros(CURVE_MAX))  # Root crown pressure
        
        # Boundary layer parameters
        self.rough = nn.Parameter(tensor(0.01), requires_grad=False)  # Soil roughness
        self.zdispl = nn.Parameter(tensor(0.065), requires_grad=False)  # Displacement height
        self.zh = nn.Parameter(tensor(0.002), requires_grad=False)  # Temperature roughness
    
    def calc_flow(self, p_inc: torch.Tensor, k_min: torch.Tensor) -> None:
        """
        Calculate E(P) curves for all components.
        
        This initializes the vulnerability curves for each hydraulic
        element in the plant.
        
        Args:
            p_inc: Pressure increment for curve calculation (MPa)
            k_min: Minimum conductance threshold
        """
        # Calculate curves for each component
        self.leaf.calc_flow_rate(p_inc, k_min)
        self.stem.calc_flow_rate(p_inc, k_min)
        
        for i in range(1, self.num_layers + 1):
            self.soils[i].root.calc_flow_rate(p_inc, k_min)
            self.soils[i].rhizosphere.calc_flow_rate(p_inc, k_min)
    
    def fatigue(self, b_wb: torch.Tensor, sapwood_t: torch.Tensor,
                conduit_d: torch.Tensor, max_plc: torch.Tensor) -> torch.Tensor:
        """
        Calculate cavitation fatigue effect on Weibull b parameter.
        
        Cavitation fatigue represents the accumulation of damage to
        xylem conduits over multiple stress events. This modifies
        the vulnerability curve for subsequent years.
        
        Per Hacke et al. (2001), repeated cavitation cycles reduce
        the pressure at which a given level of conductivity loss occurs.
        
        Args:
            b_wb: Current Weibull b parameter
            sapwood_t: Sapwood thickness (m)
            conduit_d: Conduit diameter (μm)
            max_plc: Maximum percent loss of conductivity from previous year
            
        Returns:
            Modified Weibull b parameter accounting for fatigue
        """
        # Fatigue coefficient (empirical)
        # Based on relationship between conduit size and fatigue susceptibility
        fatigue_coef = 0.1 * conduit_d / 50.0  # Normalized to 50 μm diameter
        
        # Fatigue effect proportional to previous damage
        fatigue_effect = fatigue_coef * max_plc / 100.0
        
        # Reduce b parameter (makes vulnerability curve shift left)
        b_new = b_wb * (1.0 - fatigue_effect * 0.1)
        
        # Don't let b go below 50% of original
        b_new = torch.clamp(b_new, min=0.5 * b_wb)
        
        return b_new
    
    def whole_plant_conductance(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate whole-plant hydraulic conductance at given pressure.
        
        The whole plant is modeled as:
            - Parallel root+rhizosphere elements (one per soil layer)
            - In series with single stem and leaf
            
        Total conductance: 1/K_plant = 1/K_roots + 1/K_stem + 1/K_leaf
        
        where K_roots is the parallel sum of all layer conductances.
        
        Args:
            pressure: Xylem pressure (MPa)
            
        Returns:
            Whole-plant hydraulic conductance (kg hr⁻¹ MPa⁻¹ m⁻² BA)
        """
        # Calculate root system conductance (parallel layers)
        k_root_total = tensor(0.0)
        
        for i in range(1, self.num_layers + 1):
            if not self.soils[i].cavitated:
                # Each layer has rhizosphere and root in series
                k_rhizo = self.soils[i].rhizosphere.vg(pressure)
                k_root = self.soils[i].root.weibull(pressure)
                
                # Series combination for this layer
                k_layer = 1.0 / (1.0 / (k_rhizo + 1e-10) + 1.0 / (k_root + 1e-10))
                
                # Add to parallel sum
                k_root_total = k_root_total + k_layer
        
        # Stem and leaf conductances at this pressure
        k_stem = self.stem.weibull(pressure)
        k_leaf = self.leaf.weibull(pressure)
        
        # Combine all in series
        # 1/K_plant = 1/K_roots + 1/K_stem + 1/K_leaf
        r_plant = (1.0 / (k_root_total + 1e-10) + 
                   1.0 / (k_stem + 1e-10) + 
                   1.0 / (k_leaf + 1e-10))
        
        k_plant = 1.0 / (r_plant + 1e-10)
        
        return k_plant
    
    def clean_parameters(self) -> None:
        """Reset all xylem parameters."""
        self.leaf.clean_parameters()
        self.stem.clean_parameters()
        
        for soil in self.soils:
            soil.clean_parameters()
        
        self.e_p.zero_()
        self.k_curve.zero_()
        self.root_pressure.zero_()
    
    def forward(self, predawn_pressures: torch.Tensor,
                transpiration: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass computes canopy pressure for given transpiration.
        
        This is a simplified forward model that computes the steady-state
        canopy pressure given soil water potentials and transpiration rate.
        
        Args:
            predawn_pressures: Soil water potential in each layer (MPa)
            transpiration: Transpiration rate (kg hr⁻¹ m⁻² basal area)
            
        Returns:
            Tuple of (canopy_pressure, whole_plant_conductance)
        """
        # This is a placeholder for the full Newton-Raphson solution
        # In practice, this requires iterating to find consistent pressures
        
        # Simple estimate: assume linear pressure drop
        # Total resistance = sum of series resistances
        
        # Get conductances at predawn
        k_total = tensor(0.0)
        for i in range(1, self.num_layers + 1):
            if not self.soils[i].cavitated:
                k_rhizo = self.soils[i].rhizosphere(predawn_pressures[i-1])
                k_root = self.soils[i].root(predawn_pressures[i-1])
                
                # Series: 1/k_total = 1/k_rhizo + 1/k_root
                k_layer = 1.0 / (1.0/k_rhizo + 1.0/k_root + 1e-10)
                k_total = k_total + k_layer
        
        # Add stem and leaf in series
        k_stem = self.stem(tensor(0.0))
        k_leaf = self.leaf(tensor(0.0))
        
        k_plant = 1.0 / (1.0/k_total + 1.0/k_stem + 1.0/k_leaf + 1e-10)
        
        # Pressure drop from transpiration
        # E = K * ΔP  =>  ΔP = E / K
        if k_plant > 1e-10:
            delta_p = transpiration / k_plant
        else:
            delta_p = tensor(10.0)  # Maximum pressure drop
        
        # Canopy pressure
        mean_predawn = predawn_pressures.mean()
        p_canopy = mean_predawn + delta_p
        
        return p_canopy, k_plant
    
    def __repr__(self) -> str:
        return (f"XylemComponent(layers={self.num_layers}, "
                f"leaf={self.leaf.p_crit.item():.2f} MPa, "
                f"stem={self.stem.p_crit.item():.2f} MPa)")
