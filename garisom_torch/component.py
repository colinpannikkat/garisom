"""
Base Component class for plant hydraulic elements.

This module implements the fundamental hydraulic component that serves as 
the base class for roots, stems, and leaves. It includes Weibull vulnerability
curves and flow rate calculations using differentiable PyTorch operations.
"""

import torch
import torch.nn as nn
from typing import Optional, Tuple
from .constants import (
    CURVE_MAX, XYLEM_EPSX, TRAP_ITER_MAX, FLOW_ITER_LIMIT,
    tensor, zeros, get_device, get_dtype
)


class Component(nn.Module):
    """
    Base class for plant hydraulic components (roots, stems, leaves).
    
    Implements Weibull vulnerability curves for hydraulic conductance and
    integral transforms for computing flow rates as a function of pressure.
    
    All operations are fully differentiable to enable gradient-based
    optimization and sensitivity analysis.
    
    Attributes:
        b_wb: Weibull b parameter (scale parameter related to P50)
        c_wb: Weibull c parameter (shape parameter related to curve steepness)
        k_max: Maximum hydraulic conductance (kg hr^-1 MPa^-1 m^-2)
        p_crit: Critical pressure at which conductance approaches zero
        k_min: Minimum recorded conductance (for tracking cavitation)
        res_percent: Percent resistance at k_sat
    """
    
    def __init__(self):
        super().__init__()
        
        # Weibull parameters (can be made learnable)
        self.b_wb = nn.Parameter(tensor(1.0), requires_grad=False)
        self.c_wb = nn.Parameter(tensor(1.0), requires_grad=False)
        
        # Conductance parameters
        self.k_max = nn.Parameter(tensor(1.0), requires_grad=False)
        self.p_crit = nn.Parameter(tensor(10.0), requires_grad=False)
        self.k_min = nn.Parameter(tensor(0.0), requires_grad=False)
        self.res_percent = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # Current pressure state
        self.pressure = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # E(P) and K(P) curves - stored for composite calculations
        self.register_buffer('e_p', zeros(CURVE_MAX))
        self.register_buffer('e_pv', zeros(CURVE_MAX))  # Virgin curve
        self.register_buffer('e_comp', zeros(CURVE_MAX))  # Composite curve
        self.register_buffer('e_pt', zeros(CURVE_MAX))  # Historical curve
        
        self.register_buffer('k_curve', zeros(CURVE_MAX))
        self.register_buffer('k_v', zeros(CURVE_MAX))  # Virgin K curve
        self.register_buffer('k_comp', zeros(CURVE_MAX))  # Composite K curve
        self.register_buffer('k_t', zeros(CURVE_MAX))  # Historical K curve
        
        self.register_buffer('pressure_comp', zeros(CURVE_MAX))
        self.register_buffer('pressure_v', zeros(CURVE_MAX))  # Virgin pressure
        
        # Fatigue tracking for multi-year simulations
        self.register_buffer('wb_fatigue', zeros(MAX_YEARS := 90))
        
    def weibull(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Calculate hydraulic conductance using Weibull vulnerability curve.
        
        The Weibull function describes how conductance declines with increasingly
        negative xylem pressure (water potential):
        
            K(P) = K_max * exp(-(P/b)^c)
        
        where:
            - K_max is maximum conductance
            - b is the scale parameter (related to P50)
            - c is the shape parameter (curve steepness)
        
        Args:
            pressure: Xylem pressure (MPa, typically negative or absolute value)
            
        Returns:
            Hydraulic conductance (kg hr^-1 MPa^-1 m^-2)
        """
        # Ensure pressure is non-negative for Weibull calculation
        p = torch.abs(pressure)
        
        # Weibull function: k = k_max * exp(-(p/b)^c)
        # Use safe operations to avoid numerical issues
        ratio = p / (self.b_wb + 1e-10)
        exponent = -torch.pow(ratio, self.c_wb)
        
        # Clamp to avoid extreme values
        exponent = torch.clamp(exponent, min=-50.0)
        
        k = self.k_max * torch.exp(exponent)
        
        return k
    
    def wb(self, pressure: torch.Tensor) -> torch.Tensor:
        """Alias for weibull() to match C++ API."""
        return self.weibull(pressure)
    
    def trapzd(self, p1: torch.Tensor, p2: torch.Tensor, 
               s: torch.Tensor, n: int) -> Tuple[torch.Tensor, int]:
        """
        Trapezoidal integration step for computing flow rate.
        
        Computes the n-th refinement stage of the extended trapezoidal rule
        for integrating the conductance function K(P) from p1 to p2.
        
        Args:
            p1: Lower pressure bound (MPa)
            p2: Upper pressure bound (MPa)
            s: Previous integral estimate
            n: Refinement stage (1 for first call)
            
        Returns:
            Tuple of (refined integral estimate, iteration count)
        """
        if n == 1:
            # First call - simple trapezoid
            k1 = self.weibull(p1)
            k2 = self.weibull(p2)
            s = 0.5 * (p2 - p1) * (k1 + k2)
            it = 1
        else:
            # Refined estimate with 2^(n-2) additional interior points
            it = 2 ** (n - 2)
            delta = (p2 - p1) / it
            x = p1 + 0.5 * delta
            
            # Sum contributions from interior points
            sum_k = torch.zeros_like(s)
            for j in range(it):
                sum_k = sum_k + self.weibull(x)
                x = x + delta
            
            # Update integral estimate
            s = 0.5 * (s + (p2 - p1) * sum_k / it)
            
        return s, it
    
    def qtrap(self, p1: torch.Tensor, p2: torch.Tensor, 
              eps: float = XYLEM_EPSX) -> torch.Tensor:
        """
        Compute integral of K(P) from p1 to p2 using adaptive trapezoidal rule.
        
        Returns the flow rate E for a given pressure drop:
            E = integral(K(P), P=p1 to p2)
        
        Args:
            p1: Lower pressure bound (MPa)
            p2: Upper pressure bound (MPa)
            eps: Desired relative accuracy
            
        Returns:
            Flow rate (kg hr^-1 m^-2)
        """
        s = tensor(0.0)
        olds = tensor(-1e30)
        
        for n in range(1, TRAP_ITER_MAX + 1):
            s, _ = self.trapzd(p1, p2, s, n)
            
            if n > 5:
                # Check for convergence
                if torch.abs(s - olds) < eps * torch.abs(olds) + 1e-10:
                    break
                if torch.abs(s) < 1e-10 and torch.abs(olds) < 1e-10:
                    break
                    
            olds = s.clone()
            
        return s
    
    def calc_flow_rate(self, p_inc: torch.Tensor, k_min: torch.Tensor,
                       virgin: bool = False) -> None:
        """
        Calculate the E(P) supply curve by integrating K(P).
        
        Increments pressure from 0 to P_crit and computes cumulative flow
        rate at each pressure step.
        
        Args:
            p_inc: Pressure increment (MPa)
            k_min: Minimum conductance threshold
            virgin: If True, store results in virgin curves
        """
        e = tensor(0.0)
        p = tensor(0.0)
        i = 0
        
        # Reset the arrays
        curve_e = self.e_pv if virgin else self.e_p
        curve_k = self.k_v if virgin else self.k_curve
        
        curve_e.zero_()
        curve_k.zero_()
        
        while i < CURVE_MAX - 1:
            # Get conductance at current pressure
            k = self.weibull(p)
            
            # Check termination
            if k < k_min:
                self.p_crit.data = p.clone()
                break
            
            # Store values
            curve_k[i] = k
            
            # Integrate to get flow
            if i > 0:
                p_prev = (i - 1) * p_inc
                e = e + self.qtrap(p_prev, p)
            
            curve_e[i] = e
            
            # Increment
            i += 1
            p = p + p_inc
            
            if i >= FLOW_ITER_LIMIT:
                break
        
        # Update p_crit
        if i > 0:
            self.p_crit.data = tensor(i * p_inc.item())
    
    def calc_through_flow(self, p1: torch.Tensor, p2: torch.Tensor,
                          p_inc: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Calculate flow through component between two pressures.
        
        Args:
            p1: Upstream pressure (MPa)
            p2: Downstream pressure (MPa)
            p_inc: Pressure increment for lookup
            
        Returns:
            Tuple of (flow, k_lower, k_upper)
        """
        if torch.abs(p2 - p1) < 1e-10:
            flow = tensor(0.0)
            k = self.weibull(p1)
            return flow, k, k
        
        # Compute flow via integration
        flow = self.qtrap(p1, p2)
        
        # Get conductances at boundaries
        k_lower = self.weibull(p1)
        k_upper = self.weibull(p2)
        
        return flow, k_lower, k_upper
    
    def calc_pressure(self, e: torch.Tensor, bottom_pressure: torch.Tensor,
                      p_grav: torch.Tensor, p_inc: torch.Tensor) -> torch.Tensor:
        """
        Calculate downstream pressure given flow rate and upstream pressure.
        
        Uses Newton-Raphson iteration to find pressure P such that:
            E = integral(K(P'), P'=bottom_pressure to P)
        
        Args:
            e: Target flow rate (kg hr^-1 m^-2)
            bottom_pressure: Upstream pressure (MPa)
            p_grav: Gravitational pressure drop (MPa)
            p_inc: Pressure increment for initial guess
            
        Returns:
            Downstream pressure (MPa)
        """
        if e <= 0:
            return bottom_pressure + p_grav
        
        # Initial guess from E(P) curve
        p = bottom_pressure + p_grav
        
        # Newton-Raphson iteration
        for _ in range(50):
            # Current flow estimate
            e_current = self.qtrap(bottom_pressure, p)
            
            # Residual
            residual = e_current - e
            
            if torch.abs(residual) < 1e-8:
                break
            
            # Derivative (K at current pressure)
            k = self.weibull(p)
            
            if k < 1e-10:
                break
            
            # Update pressure
            p = p - residual / k
            
            # Ensure pressure doesn't go below upstream
            p = torch.clamp(p, min=bottom_pressure + p_grav)
        
        return p
    
    def update_kmin(self) -> None:
        """Update minimum conductance based on current pressure exposure."""
        k_current = self.weibull(self.pressure)
        self.k_min.data = torch.min(self.k_min, k_current)
    
    def store_curves(self) -> None:
        """Store current curves to historical buffers."""
        self.e_pt.copy_(self.e_p)
        self.k_t.copy_(self.k_curve)
    
    def restore_curves(self) -> None:
        """Restore curves from historical buffers."""
        self.e_p.copy_(self.e_pt)
        self.k_curve.copy_(self.k_t)
    
    def use_virgin_curves(self) -> None:
        """Reset to virgin (uncavitated) curves."""
        self.e_p.copy_(self.e_pv)
        self.k_curve.copy_(self.k_v)
    
    def clean_parameters(self) -> None:
        """Reset all parameters to default values."""
        self.b_wb.data.fill_(1.0)
        self.c_wb.data.fill_(1.0)
        self.k_max.data.fill_(1.0)
        self.p_crit.data.fill_(10.0)
        self.k_min.data.fill_(0.0)
        self.res_percent.data.fill_(0.0)
        self.pressure.data.fill_(0.0)
        
        self.e_p.zero_()
        self.e_pv.zero_()
        self.e_comp.zero_()
        self.e_pt.zero_()
        self.k_curve.zero_()
        self.k_v.zero_()
        self.k_comp.zero_()
        self.k_t.zero_()
        self.pressure_comp.zero_()
        self.pressure_v.zero_()
        self.wb_fatigue.zero_()
    
    def forward(self, pressure: torch.Tensor) -> torch.Tensor:
        """
        Forward pass computes conductance for given pressure.
        
        Args:
            pressure: Input pressure(s) (MPa)
            
        Returns:
            Hydraulic conductance (kg hr^-1 MPa^-1 m^-2)
        """
        return self.weibull(pressure)
    
    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"b_wb={self.b_wb.item():.4f}, "
                f"c_wb={self.c_wb.item():.4f}, "
                f"k_max={self.k_max.item():.4f}, "
                f"p_crit={self.p_crit.item():.4f})")
