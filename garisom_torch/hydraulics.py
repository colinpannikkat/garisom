"""
Hydraulics Model for GARISOM-Torch.

This module implements the soil-plant-atmosphere continuum (SPAC) hydraulics,
including the calculation of transpiration supply and demand functions.

The hydraulics model follows Sperry et al. (2016, 2017) where:
    - E(P) supply function describes water availability as canopy pressure drops
    - Demand function locates plant on supply curve via stomatal regulation
    - Optimization finds canopy pressure maximizing gain minus risk
"""

import torch
import torch.nn as nn
from typing import Optional, Tuple
from .constants import (
    CURVE_MAX, GMAX, PI, tensor, zeros, get_device, get_dtype
)


class HydraulicsModel(nn.Module):
    """
    Soil-Plant-Atmosphere Continuum (SPAC) hydraulics model.
    
    Implements the hydraulic transport calculations following Sperry et al.:
        - Transpiration supply: E(P_canopy) from soil to leaves
        - Transpiration demand: E from atmospheric demand (VPD)
        - Hydraulic conductance: K(P) vulnerability curves
        
    The model integrates conductance curves to determine water supply
    and uses stomatal regulation to balance carbon gain with hydraulic risk.
    
    Key equations:
        E = ∫K(P)dP  (supply from predawn to canopy pressure)
        E = G × D    (demand from stomata and VPD)
        
    where:
        E = Transpiration rate
        K = Hydraulic conductance
        P = Xylem pressure
        G = Canopy diffusive conductance
        D = Vapor pressure deficit
    """
    
    def __init__(self):
        super().__init__()
        
        # Maximum canopy conductance (kg hr⁻¹ MPa⁻¹ m⁻² basal area)
        self.gmax = nn.Parameter(tensor(GMAX), requires_grad=False)
        
        # Current states
        self.e_supply = nn.Parameter(tensor(0.0), requires_grad=False)
        self.e_demand = nn.Parameter(tensor(0.0), requires_grad=False)
        self.p_canopy = nn.Parameter(tensor(0.0), requires_grad=False)
        self.k_plant = nn.Parameter(tensor(0.0), requires_grad=False)
    
    def transpiration_supply(self, p_canopy: torch.Tensor,
                              predawn: torch.Tensor,
                              k_plant: torch.Tensor) -> torch.Tensor:
        """
        Calculate transpiration supply from hydraulic pathway.
        
        The supply function describes the potential rate of water supply
        as a function of canopy xylem pressure. It is derived by integrating
        the K(P) vulnerability curve from predawn to canopy pressure:
        
            E(P_c) = ∫K(P)dP from P_pd to P_c
        
        For a simple linear approximation:
            E ≈ K × (P_c - P_pd)
        
        Args:
            p_canopy: Canopy xylem pressure (MPa)
            predawn: Predawn pressure / soil water potential (MPa)
            k_plant: Whole-plant hydraulic conductance (kg hr⁻¹ MPa⁻¹ m⁻²)
            
        Returns:
            Transpiration rate (kg hr⁻¹ m⁻² basal area)
        """
        # Pressure drop from predawn to canopy
        delta_p = p_canopy - predawn
        delta_p = torch.clamp(delta_p, min=0.0)
        
        # Supply = K × ΔP
        e_supply = k_plant * delta_p
        
        return e_supply
    
    def transpiration_demand(self, gcanopy: torch.Tensor,
                              vpd: torch.Tensor) -> torch.Tensor:
        """
        Calculate transpiration demand from atmospheric conditions.
        
        The demand function specifies the transpiration rate given
        stomatal opening and vapor pressure deficit:
        
            E = G × D
        
        where G is canopy diffusive conductance and D is VPD.
        
        Unit conversions:
            - G typically in mmol m⁻² s⁻¹ (leaf area)
            - D in kPa
            - E in kg hr⁻¹ m⁻² (basal area)
        
        Args:
            gcanopy: Canopy diffusive conductance (mmol m⁻² s⁻¹)
            vpd: Vapor pressure deficit (kPa)
            
        Returns:
            Transpiration rate (mmol m⁻² s⁻¹ leaf area)
        """
        # E = G × D (mmol m⁻² s⁻¹)
        e_demand = gcanopy * vpd
        
        return e_demand
    
    def optimal_stomatal_conductance(self, 
                                      dedp: torch.Tensor,
                                      dadp: torch.Tensor,
                                      gmax: torch.Tensor) -> torch.Tensor:
        """
        Calculate optimal stomatal conductance from marginal gain/risk.
        
        Stomatal optimization theory (Sperry et al. 2017) predicts that
        stomata regulate to the point where marginal carbon gain equals
        marginal hydraulic risk:
        
            dA'/dP = dθ'/dP
        
        This determines the optimal operating point on the supply curve.
        
        Args:
            dedp: Marginal transpiration response (∂E/∂P)
            dadp: Marginal assimilation response (∂A/∂P)
            gmax: Maximum stomatal conductance
            
        Returns:
            Optimal stomatal conductance (mmol m⁻² s⁻¹)
        """
        # The marginal cost of water = dedp / k_max
        # The marginal benefit = dadp / amax
        
        # At optimum, d(A - θ)/dP = 0
        # Simplified: g_opt ∝ sqrt(A × k)
        
        g_opt = torch.sqrt(dedp * dadp + 1e-10)
        g_opt = torch.clamp(g_opt, max=gmax)
        
        return g_opt
    
    def profit_maximization(self,
                            e_curve: torch.Tensor,
                            a_curve: torch.Tensor,
                            p_curve: torch.Tensor,
                            e_crit: torch.Tensor,
                            a_max: torch.Tensor) -> Tuple[int, torch.Tensor, torch.Tensor]:
        """
        Find operating point that maximizes gain minus risk.
        
        The gain-risk framework (Sperry et al. 2017) defines:
            - Gain: A' = A / A_max (normalized assimilation)
            - Risk: θ' = 1 - K/K_max (normalized conductance loss)
            
        The optimal point maximizes (A' - θ').
        
        Args:
            e_curve: Transpiration at each pressure increment
            a_curve: Assimilation at each pressure increment
            p_curve: Canopy pressure at each increment
            e_crit: Critical transpiration (at hydraulic failure)
            a_max: Maximum assimilation
            
        Returns:
            Tuple of (optimal_index, optimal_E, optimal_A)
        """
        # Normalize curves
        a_prime = a_curve / (a_max + 1e-10)  # Normalized gain
        e_prime = e_curve / (e_crit + 1e-10)  # Normalized transpiration
        
        # Risk function: θ' = 1 - K/K_max = E/E_crit (approximately)
        # This simplification assumes linear relationship
        theta_prime = e_prime
        
        # Profit = A' - θ'
        profit = a_prime - theta_prime
        
        # Find maximum profit
        max_idx = torch.argmax(profit)
        
        optimal_e = e_curve[max_idx]
        optimal_a = a_curve[max_idx]
        
        return int(max_idx.item()), optimal_e, optimal_a
    
    def calculate_k_plant(self,
                          e: torch.Tensor,
                          p_canopy: torch.Tensor,
                          p_predawn: torch.Tensor,
                          p_grav: torch.Tensor) -> torch.Tensor:
        """
        Calculate whole-plant hydraulic conductance.
        
        K_plant = E / (P_canopy - P_predawn - P_grav)
        
        Args:
            e: Transpiration rate
            p_canopy: Canopy pressure (MPa)
            p_predawn: Predawn pressure (MPa)
            p_grav: Gravitational pressure drop (MPa)
            
        Returns:
            Whole-plant conductance (kg hr⁻¹ MPa⁻¹ m⁻² basal area)
        """
        delta_p = p_canopy - p_predawn - p_grav
        delta_p = torch.clamp(delta_p, min=1e-6)
        
        k_plant = e / delta_p
        
        return k_plant
    
    def percent_loss_conductance(self,
                                  k_current: torch.Tensor,
                                  k_max: torch.Tensor) -> torch.Tensor:
        """
        Calculate percent loss of conductance (PLC).
        
        PLC = 100 × (1 - K/K_max)
        
        Args:
            k_current: Current conductance
            k_max: Maximum (initial) conductance
            
        Returns:
            Percent loss of conductance (0-100)
        """
        plc = 100.0 * (1.0 - k_current / (k_max + 1e-10))
        plc = torch.clamp(plc, 0.0, 100.0)
        
        return plc
    
    def vapor_pressure_deficit(self,
                                tair: torch.Tensor,
                                rh: Optional[torch.Tensor] = None,
                                eair: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Calculate vapor pressure deficit.
        
        VPD = e_sat(T) - e_air
        
        Args:
            tair: Air temperature (°C)
            rh: Relative humidity (0-1), optional
            eair: Actual vapor pressure (kPa), optional
            
        Returns:
            Vapor pressure deficit (kPa)
        """
        # Saturation vapor pressure (Tetens equation)
        esat = 0.611 * torch.exp(17.502 * tair / (tair + 240.97))
        
        if eair is not None:
            vpd = esat - eair
        elif rh is not None:
            vpd = esat * (1.0 - rh)
        else:
            # Assume 50% RH as default
            vpd = esat * 0.5
        
        vpd = torch.clamp(vpd, min=0.0)
        
        return vpd
    
    def forward(self, 
                predawn: torch.Tensor,
                vpd: torch.Tensor,
                k_max: torch.Tensor,
                gmax: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward pass computes steady-state transpiration and canopy pressure.
        
        Finds the operating point where supply meets demand for given
        soil water potential and atmospheric conditions.
        
        Args:
            predawn: Soil water potential / predawn pressure (MPa)
            vpd: Vapor pressure deficit (kPa)
            k_max: Maximum plant hydraulic conductance
            gmax: Maximum stomatal conductance
            
        Returns:
            Tuple of (transpiration, canopy_pressure, conductance)
        """
        # Maximum possible transpiration (at critical point)
        # For now, use simplified estimate
        e_max = k_max * 5.0  # Assume 5 MPa pressure drop at critical
        
        # Demand at maximum stomatal opening
        e_demand = gmax * vpd
        
        # Actual transpiration is minimum of supply and demand
        e_actual = torch.min(e_max, e_demand)
        
        # Canopy pressure from supply equation
        # E = K × ΔP  =>  ΔP = E / K
        delta_p = e_actual / (k_max + 1e-10)
        p_canopy = predawn + delta_p
        
        # Actual conductance (may be reduced from max)
        if vpd > 1e-6:
            g_actual = e_actual / vpd
        else:
            g_actual = gmax
        
        return e_actual, p_canopy, g_actual
    
    def __repr__(self) -> str:
        return f"HydraulicsModel(gmax={self.gmax.item():.2e})"
