"""
Soil and Rhizosphere components for GARISOM-Torch.

This module implements Van Genuchten soil hydraulics for the rhizosphere,
which models the soil-to-root water pathway. The rhizosphere is treated
as a distinct hydraulic element with its own vulnerability characteristics.
"""

import torch
import torch.nn as nn
from typing import Optional, Tuple, Dict
from dataclasses import dataclass

from .constants import (
    VAN_GENUCHTEN_PARAMS, CURVE_MAX, RHIZO_EPSX, TRAP_ITER_MAX,
    tensor, zeros, get_device, get_dtype
)
from .component import Component


class RhizosphereComponent(nn.Module):
    """
    Rhizosphere hydraulic component using Van Genuchten functions.
    
    The rhizosphere represents the soil-to-root interface and uses
    Van Genuchten equations to describe the relationship between
    soil water potential and hydraulic conductivity.
    
    Van Genuchten equations:
        θ/θs = [1 + (α|ψ|)^n]^(-m)  where m = 1 - 1/n
        K/Ks = Se^0.5 * [1 - (1 - Se^(1/m))^m]^2
    
    where:
        θ = volumetric water content
        θs = saturated water content
        ψ = soil water potential (MPa)
        α = Van Genuchten alpha parameter (MPa^-1)
        n = Van Genuchten n parameter
        K = hydraulic conductivity
        Ks = saturated hydraulic conductivity
        Se = effective saturation
        
    Attributes:
        van_gen_alpha: α parameter (MPa^-1)
        van_gen_n: n parameter (dimensionless)
        theta_sat: Saturated water content (m^3/m^3)
        k_max: Maximum hydraulic conductance
    """
    
    def __init__(self, texture: str = "loam"):
        super().__init__()
        
        # Get Van Genuchten parameters from texture
        if texture in VAN_GENUCHTEN_PARAMS:
            alpha, n, ksat, theta_sat = VAN_GENUCHTEN_PARAMS[texture]
        else:
            # Default to loam
            alpha, n, ksat, theta_sat = VAN_GENUCHTEN_PARAMS["loam"]
        
        # Van Genuchten parameters
        self.van_gen_alpha = nn.Parameter(tensor(alpha), requires_grad=False)
        self.van_gen_n = nn.Parameter(tensor(n), requires_grad=False)
        self.theta_sat = nn.Parameter(tensor(theta_sat), requires_grad=False)
        
        # Conductance parameters
        self.k_max = nn.Parameter(tensor(ksat), requires_grad=False)
        self.k_sat = nn.Parameter(tensor(ksat), requires_grad=False)
        self.p_crit = nn.Parameter(tensor(10.0), requires_grad=False)
        self.k_min = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # Current state
        self.pressure = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # Curves for composite calculations
        self.register_buffer('e_p', zeros(CURVE_MAX))
        self.register_buffer('k_curve', zeros(CURVE_MAX))
        self.register_buffer('pressure_comp', zeros(CURVE_MAX))
        self.register_buffer('e_comp', zeros(CURVE_MAX))
        self.register_buffer('k_comp', zeros(CURVE_MAX))
    
    def van_genuchten_se(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate effective saturation using Van Genuchten equation.
        
        Se = [1 + (α|ψ|)^n]^(-m)  where m = 1 - 1/n
        
        Args:
            pressure: Soil water potential (MPa, absolute value)
            
        Returns:
            Effective saturation (0-1)
        """
        p = torch.abs(pressure)
        
        # Van Genuchten m parameter
        m = 1.0 - 1.0 / self.van_gen_n
        
        # Effective saturation
        term = torch.pow(self.van_gen_alpha * p, self.van_gen_n)
        se = torch.pow(1.0 + term, -m)
        
        return se
    
    def vg(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate hydraulic conductance using Van Genuchten function.
        
        K = Ks * Se^0.5 * [1 - (1 - Se^(1/m))^m]^2
        
        Args:
            pressure: Soil water potential (MPa)
            
        Returns:
            Hydraulic conductance (kg hr^-1 MPa^-1 m^-2)
        """
        p = torch.abs(pressure) + 1e-10  # Small offset to avoid zero
        
        # Van Genuchten m parameter
        m = 1.0 - 1.0 / self.van_gen_n
        
        # Effective saturation
        se = self.van_genuchten_se(p)
        
        # Hydraulic conductivity (Mualem-Van Genuchten)
        # K/Ks = Se^0.5 * [1 - (1 - Se^(1/m))^m]^2
        se_inv_m = torch.pow(se + 1e-10, 1.0 / m)
        term = 1.0 - torch.pow(1.0 - se_inv_m + 1e-10, m)
        k_rel = torch.pow(se, 0.5) * torch.pow(term, 2.0)
        
        # Absolute conductance
        k = self.k_max * k_rel
        
        return k
    
    def rvg(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Reverse Van Genuchten: calculate pressure from water content.
        
        Args:
            pressure: Actually theta/theta_sat ratio here
            
        Returns:
            Corresponding pressure (MPa)
        """
        # This inverts the VG equation
        theta_frac = torch.clamp(pressure, 0.001, 0.999)
        m = 1.0 - 1.0 / self.van_gen_n
        
        # Invert Se = [1 + (αψ)^n]^(-m)
        # (αψ)^n = Se^(-1/m) - 1
        se_inv_m = torch.pow(theta_frac, -1.0 / m)
        term = se_inv_m - 1.0
        term = torch.clamp(term, min=1e-10)
        
        alpha_psi_n = term
        psi = torch.pow(alpha_psi_n, 1.0 / self.van_gen_n) / self.van_gen_alpha
        
        return psi
    
    def swc(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate soil water content fraction from pressure.
        
        Args:
            pressure: Soil water potential (MPa)
            
        Returns:
            θ/θs ratio (0-1)
        """
        return self.van_genuchten_se(pressure)
    
    def theta_from_pressure(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate volumetric water content from pressure.
        
        Args:
            pressure: Soil water potential (MPa)
            
        Returns:
            Volumetric water content (m^3/m^3)
        """
        se = self.van_genuchten_se(pressure)
        return se * self.theta_sat
    
    def pressure_from_theta(self, theta: torch.Tensor) -> torch.Tensor:
        """
        Calculate soil water potential from volumetric water content.
        
        Inverts the Van Genuchten equation:
            Se = θ/θs
            ψ = (1/α) * (Se^(-1/m) - 1)^(1/n)
        
        Args:
            theta: Volumetric water content (m^3/m^3)
            
        Returns:
            Soil water potential (MPa, negative)
        """
        # Effective saturation
        se = theta / (self.theta_sat + 1e-10)
        se = torch.clamp(se, 0.001, 0.999)  # Avoid numerical issues
        
        # Van Genuchten m parameter
        m = 1.0 - 1.0 / self.van_gen_n
        
        # Invert: Se = [1 + (αψ)^n]^(-m)
        # => (αψ)^n = Se^(-1/m) - 1
        # => ψ = (1/α) * (Se^(-1/m) - 1)^(1/n)
        se_inv_m = torch.pow(se, -1.0 / m)
        term = se_inv_m - 1.0
        term = torch.clamp(term, min=1e-10)
        
        psi = torch.pow(term, 1.0 / self.van_gen_n) / self.van_gen_alpha
        
        # Return as negative (suction)
        return -psi
    
    def trapzd(self, p1: torch.Tensor, p2: torch.Tensor,
               s: torch.Tensor, n: int) -> Tuple[torch.Tensor, int]:
        """Trapezoidal integration step for Van Genuchten function."""
        if n == 1:
            k1 = self.vg(p1)
            k2 = self.vg(p2)
            s = 0.5 * (p2 - p1) * (k1 + k2)
            it = 1
        else:
            it = 2 ** (n - 2)
            delta = (p2 - p1) / it
            x = p1 + 0.5 * delta
            
            sum_k = torch.zeros_like(s)
            for j in range(it):
                sum_k = sum_k + self.vg(x)
                x = x + delta
            
            s = 0.5 * (s + (p2 - p1) * sum_k / it)
            
        return s, it
    
    def qtrap(self, p1: torch.Tensor, p2: torch.Tensor,
              eps: float = RHIZO_EPSX) -> torch.Tensor:
        """Adaptive trapezoidal integration for rhizosphere flow."""
        s = tensor(0.0)
        olds = tensor(-1e30)
        
        for n in range(1, TRAP_ITER_MAX + 1):
            s, _ = self.trapzd(p1, p2, s, n)
            
            if n > 5:
                if torch.abs(s - olds) < eps * torch.abs(olds) + 1e-10:
                    break
                if torch.abs(s) < 1e-10 and torch.abs(olds) < 1e-10:
                    break
                    
            olds = s.clone()
            
        return s
    
    def calc_flow_rate(self, p_inc: torch.Tensor, k_min: torch.Tensor) -> None:
        """Calculate E(P) curve for rhizosphere."""
        e = tensor(0.0)
        p = tensor(0.0)
        i = 0
        
        self.e_p.zero_()
        self.k_curve.zero_()
        
        while i < CURVE_MAX - 1:
            k = self.vg(p)
            
            if k < k_min:
                self.p_crit.data = p.clone()
                break
            
            self.k_curve[i] = k
            
            if i > 0:
                p_prev = tensor((i - 1) * p_inc.item())
                e = e + self.qtrap(p_prev, p)
            
            self.e_p[i] = e
            
            i += 1
            p = p + p_inc
        
        if i > 0:
            self.p_crit.data = tensor(i * p_inc.item())
    
    def forward(self, pressure: torch.Tensor) -> torch.Tensor:
        """Forward pass returns hydraulic conductance."""
        return self.vg(pressure)
    
    def clean_parameters(self) -> None:
        """Reset parameters to default values."""
        self.pressure.data.fill_(0.0)
        self.k_min.data.fill_(0.0)
        self.e_p.zero_()
        self.k_curve.zero_()
        self.pressure_comp.zero_()
        self.e_comp.zero_()
        self.k_comp.zero_()
    
    def __repr__(self) -> str:
        return (f"RhizosphereComponent("
                f"α={self.van_gen_alpha.item():.2f}, "
                f"n={self.van_gen_n.item():.2f}, "
                f"θs={self.theta_sat.item():.3f})")


@dataclass
class SoilLayerState:
    """Container for soil layer state variables."""
    cavitated: bool = False
    failure: str = "ok"


class SoilLayer(nn.Module):
    """
    A single soil layer containing root and rhizosphere components.
    
    Each soil layer represents a horizontal slice of the rooting zone,
    with its own root biomass fraction and hydraulic properties.
    
    Attributes:
        root: Root component for this layer
        rhizosphere: Rhizosphere (soil) component
        biomass_fraction: Fraction of total root biomass in this layer
        layer_depth: Lower depth of this layer (m)
        depth: Thickness of this layer (m)
        vert_distance: Vertical distance to biomass center (m)
        radius: Radial extent of roots in this layer (m)
        length: Total transport distance (m)
    """
    
    def __init__(self, texture: str = "loam"):
        super().__init__()
        
        # Import here to avoid circular import
        from .xylem import RootComponent
        
        self.root = RootComponent()
        self.rhizosphere = RhizosphereComponent(texture)
        
        # Layer geometry
        self.biomass_fraction = nn.Parameter(tensor(0.0), requires_grad=False)
        self.layer_depth = nn.Parameter(tensor(0.0), requires_grad=False)
        self.depth = nn.Parameter(tensor(0.0), requires_grad=False)
        self.vert_distance = nn.Parameter(tensor(0.0), requires_grad=False)
        self.radius = nn.Parameter(tensor(0.0), requires_grad=False)
        self.length = nn.Parameter(tensor(1.0), requires_grad=False)
        
        # Hydraulic parameters
        self.kkmax = nn.Parameter(tensor(0.0), requires_grad=False)
        self.swclimit = nn.Parameter(tensor(0.0), requires_grad=False)
        self.soilredist = nn.Parameter(tensor(0.0), requires_grad=False)
        self.flow = nn.Parameter(tensor(0.0), requires_grad=False)
        self.predawn_pressure = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # State tracking
        self.cavitated = False
        self.cavitated_t = False
        self.failure = "ok"
        self.failure_t = "ok"
    
    def set_van_genuchten_params(self, alpha: float, n: float,
                                  ksat: float, theta_sat: float) -> None:
        """Set Van Genuchten parameters for this layer."""
        self.rhizosphere.van_gen_alpha.data = tensor(alpha)
        self.rhizosphere.van_gen_n.data = tensor(n)
        self.rhizosphere.k_max.data = tensor(ksat)
        self.rhizosphere.k_sat.data = tensor(ksat)
        self.rhizosphere.theta_sat.data = tensor(theta_sat)
        self.kkmax.data = tensor(ksat)
    
    def reset_state(self) -> None:
        """Reset layer failure state."""
        self.cavitated = False
        self.failure = "ok"
    
    def forward(self, pressure: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass computes root and rhizosphere conductances.
        
        Args:
            pressure: Input pressure (MPa)
            
        Returns:
            Tuple of (root_conductance, rhizosphere_conductance)
        """
        k_root = self.root(pressure)
        k_rhizo = self.rhizosphere(pressure)
        return k_root, k_rhizo
    
    def clean_parameters(self) -> None:
        """Reset all parameters."""
        self.root.clean_parameters()
        self.rhizosphere.clean_parameters()
        self.biomass_fraction.data.fill_(0.0)
        self.layer_depth.data.fill_(0.0)
        self.depth.data.fill_(0.0)
        self.vert_distance.data.fill_(0.0)
        self.radius.data.fill_(0.0)
        self.length.data.fill_(1.0)
        self.kkmax.data.fill_(0.0)
        self.swclimit.data.fill_(0.0)
        self.soilredist.data.fill_(0.0)
        self.flow.data.fill_(0.0)
        self.predawn_pressure.data.fill_(0.0)
        self.cavitated = False
        self.failure = "ok"


def get_vg_params_from_texture(texture: str) -> Tuple[float, float, float, float]:
    """
    Get Van Genuchten parameters for a soil texture class.
    
    Args:
        texture: Soil texture name (e.g., "loam", "sand", "clay")
        
    Returns:
        Tuple of (alpha, n, ksat, theta_sat)
    """
    if texture.lower() in VAN_GENUCHTEN_PARAMS:
        return VAN_GENUCHTEN_PARAMS[texture.lower()]
    else:
        print(f"Warning: Unknown texture '{texture}', using 'loam' defaults")
        return VAN_GENUCHTEN_PARAMS["loam"]
