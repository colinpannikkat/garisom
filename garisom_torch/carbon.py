"""
Carbon Assimilation Model for GARISOM-Torch.

This module implements the Farquhar-type photosynthesis model following
Medlyn et al. (2002). It calculates net carbon assimilation based on
stomatal conductance, light, temperature, and CO2 concentration.

The model computes:
    - Stomatal conductance from transpiration and VPD
    - Temperature-corrected Rubisco and electron transport parameters
    - Light-limited and Rubisco-limited photosynthesis
    - Net carbon assimilation accounting for respiration
"""

import torch
import torch.nn as nn
from typing import Optional, Tuple
from .constants import (
    CURVE_MAX, GAS, OA, tensor, zeros, get_device, get_dtype
)


class CarbonAssimilationModel(nn.Module):
    """
    Farquhar-von Caemmerer-Berry photosynthesis model.
    
    Implements the carbon assimilation equations following:
        - Farquhar et al. (1980) for C3 photosynthesis
        - Medlyn et al. (2002) for temperature corrections
        - Ball-Berry type stomatal conductance
    
    The model calculates net photosynthesis A_net as the minimum of
    Rubisco-limited (Ac) and electron transport-limited (Aj) rates:
    
        A_net = min(Ac, Aj) - Rd
    
    where:
        Ac = Vcmax * (Ci - Γ*) / (Ci + Kc*(1 + O/Ko))
        Aj = J * (Ci - Γ*) / (4*(Ci + 2*Γ*))
        Rd = Day respiration
    
    Attributes:
        ca: Atmospheric CO2 concentration (mol/mol)
        vmax25: Maximum carboxylation rate at 25°C (μmol m⁻² s⁻¹)
        jmax25: Maximum electron transport rate at 25°C (μmol m⁻² s⁻¹)
        comp25: CO2 compensation point at 25°C (μmol/mol)
        kc25: Michaelis constant for CO2 at 25°C (μmol/mol)
        ko25: Michaelis constant for O2 at 25°C (mmol/mol)
    """
    
    def __init__(self):
        super().__init__()
        
        # Atmospheric CO2 (mol/mol, converted from ppm)
        self.ca = nn.Parameter(tensor(390e-6), requires_grad=False)
        
        # Kinetic parameters at 25°C
        self.vmax25 = nn.Parameter(tensor(50.0), requires_grad=False)  # μmol m⁻² s⁻¹
        self.jmax25 = nn.Parameter(tensor(90.0), requires_grad=False)  # μmol m⁻² s⁻¹
        self.comp25 = nn.Parameter(tensor(42.75e-6), requires_grad=False)  # mol/mol (Γ*)
        self.kc25 = nn.Parameter(tensor(404.9e-6), requires_grad=False)  # mol/mol
        self.ko25 = nn.Parameter(tensor(278.4e-3), requires_grad=False)  # mol/mol
        
        # Temperature response parameters
        self.havmax = nn.Parameter(tensor(73637.0), requires_grad=False)  # J mol⁻¹
        self.hdvmax = nn.Parameter(tensor(149252.0), requires_grad=False)  # J mol⁻¹
        self.svvmax = nn.Parameter(tensor(486.0), requires_grad=False)  # J mol⁻¹ K⁻¹
        
        self.hajmax = nn.Parameter(tensor(50300.0), requires_grad=False)  # J mol⁻¹
        self.hdjmax = nn.Parameter(tensor(152044.0), requires_grad=False)  # J mol⁻¹
        self.svjmax = nn.Parameter(tensor(495.0), requires_grad=False)  # J mol⁻¹ K⁻¹
        
        # Light response parameters
        self.qmax = nn.Parameter(tensor(0.3), requires_grad=False)  # Quantum yield
        self.thetac = nn.Parameter(tensor(0.98), requires_grad=False)  # Curvature
        self.lightcurv = nn.Parameter(tensor(0.9), requires_grad=False)  # Light curve
        self.lightcomp = nn.Parameter(tensor(30.0), requires_grad=False)  # Light comp point
        
        # Current light levels (from solar calculations)
        self.qsl = nn.Parameter(tensor(0.0), requires_grad=False)  # PPFD sunlit
        self.qsh = nn.Parameter(tensor(0.0), requires_grad=False)  # PPFD shaded
        
        # Calculated assimilation values
        self.psynmax = nn.Parameter(tensor(0.0), requires_grad=False)  # Max A sunlit
        self.psynmaxsh = nn.Parameter(tensor(0.0), requires_grad=False)  # Max A shaded
        self.psynact = nn.Parameter(tensor(0.0), requires_grad=False)  # Actual A sunlit
        self.psynactsh = nn.Parameter(tensor(0.0), requires_grad=False)  # Actual A shaded
        self.atree = nn.Parameter(tensor(0.0), requires_grad=False)  # Tree-level A
        
        # Intercellular CO2 and conductances
        self.cinc = nn.Parameter(tensor(0.0), requires_grad=False)  # Ci sunlit
        self.cincsh = nn.Parameter(tensor(0.0), requires_grad=False)  # Ci shaded
        
        # Conductances at midday
        self.gcanwmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.gcanwshmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.gcancmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.gcancshmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.gcmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.gcmdsh = nn.Parameter(tensor(0.0), requires_grad=False)
        
        self.cinmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.cinshmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.psynmaxmd = nn.Parameter(tensor(0.0), requires_grad=False)
        self.psynmaxshmd = nn.Parameter(tensor(0.0), requires_grad=False)
        
        # Curves for all pressure increments
        self.register_buffer('cin', zeros(CURVE_MAX))  # Ci at each E increment
        self.register_buffer('cinsh', zeros(CURVE_MAX))
        self.register_buffer('psyn', zeros(CURVE_MAX))  # A at each E increment
        self.register_buffer('psynsh', zeros(CURVE_MAX))
        self.register_buffer('psynmd', zeros(CURVE_MAX))
        self.register_buffer('psynshmd', zeros(CURVE_MAX))
        self.register_buffer('psync', zeros(CURVE_MAX))
        self.register_buffer('gcanw', zeros(CURVE_MAX))  # Water conductance
        self.register_buffer('gcanc', zeros(CURVE_MAX))  # CO2 conductance
        self.register_buffer('gcanwsh', zeros(CURVE_MAX))
        self.register_buffer('gcancsh', zeros(CURVE_MAX))
    
    def temperature_correction_vmax(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Temperature correction for Vmax using peaked Arrhenius function.
        
        Vmax(T) = Vmax25 * f(T) where f is the peaked function
        from Medlyn et al. (2002).
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Vmax at leaf temperature (μmol m⁻² s⁻¹)
        """
        tk = tleaf + 273.15
        tk25 = 298.15
        
        # Peaked Arrhenius function
        numerator = (1.0 + torch.exp((self.svvmax * tk25 - self.hdvmax) / (GAS * tk25))) * \
                    torch.exp((self.havmax / (GAS * tk25)) * (1.0 - tk25 / tk))
        denominator = 1.0 + torch.exp((self.svvmax * tk - self.hdvmax) / (GAS * tk))
        
        vmax = self.vmax25 * numerator / (denominator + 1e-10)
        
        return vmax
    
    def temperature_correction_jmax(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Temperature correction for Jmax using peaked Arrhenius function.
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Jmax at leaf temperature (μmol m⁻² s⁻¹)
        """
        tk = tleaf + 273.15
        tk25 = 298.15
        
        numerator = (1.0 + torch.exp((self.svjmax * tk25 - self.hdjmax) / (GAS * tk25))) * \
                    torch.exp((self.hajmax / (GAS * tk25)) * (1.0 - tk25 / tk))
        denominator = 1.0 + torch.exp((self.svjmax * tk - self.hdjmax) / (GAS * tk))
        
        jmax = self.jmax25 * numerator / (denominator + 1e-10)
        
        return jmax
    
    def temperature_correction_kc(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Temperature correction for Kc (Michaelis constant for CO2).
        
        Uses Arrhenius equation per Bernacchi et al. (2001).
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Kc at leaf temperature (mol/mol)
        """
        tk = tleaf + 273.15
        tk25 = 298.15
        
        # Activation energy for Kc = 79430 J/mol
        kc = self.kc25 * torch.exp((79430.0 * (tk - tk25)) / (tk25 * GAS * tk))
        
        return kc
    
    def temperature_correction_ko(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Temperature correction for Ko (Michaelis constant for O2).
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Ko at leaf temperature (mol/mol)
        """
        tk = tleaf + 273.15
        tk25 = 298.15
        
        # Activation energy for Ko = 36380 J/mol
        ko = self.ko25 * torch.exp((36380.0 * (tk - tk25)) / (tk25 * GAS * tk))
        
        return ko
    
    def temperature_correction_gamma(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Temperature correction for Γ* (CO2 compensation point).
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Γ* at leaf temperature (mol/mol)
        """
        tk = tleaf + 273.15
        tk25 = 298.15
        
        # Activation energy = 37830 J/mol per Bernacchi
        gamma = self.comp25 * torch.exp((37830.0 * (tk - tk25)) / (tk25 * GAS * tk))
        
        return gamma
    
    def day_respiration(self, tleaf: torch.Tensor) -> torch.Tensor:
        """
        Calculate day respiration rate.
        
        Rd = 0.01 * Vmax25 * 2^((T-25)/10) with high-temp inhibition.
        
        Args:
            tleaf: Leaf temperature (°C)
            
        Returns:
            Day respiration rate (μmol m⁻² s⁻¹)
        """
        rday25 = self.vmax25 * 0.01
        rday = rday25 * torch.pow(tensor(2.0), (tleaf - 25.0) / 10.0)
        
        # High temperature inhibition (Collatz et al.)
        inhibition = torch.pow(1.0 + torch.exp(1.3 * (tleaf - 55.0)), -1.0)
        rday = rday * inhibition
        
        return rday
    
    def electron_transport_rate(self, ppfd: torch.Tensor, 
                                 jmax: torch.Tensor) -> torch.Tensor:
        """
        Calculate actual electron transport rate from light and Jmax.
        
        Uses the non-rectangular hyperbola:
            θ*J² - (I + Jmax)*J + I*Jmax = 0
        
        where I = quantum yield * PPFD
        
        Args:
            ppfd: Photosynthetic photon flux density (μmol m⁻² s⁻¹)
            jmax: Maximum electron transport rate (μmol m⁻² s⁻¹)
            
        Returns:
            Actual electron transport rate J (μmol m⁻² s⁻¹)
        """
        alpha_q = self.qmax * ppfd  # Light-limited electron transport
        
        # Quadratic solution
        a = self.lightcurv
        b = -(alpha_q + jmax)
        c = alpha_q * jmax
        
        discriminant = b * b - 4.0 * a * c
        discriminant = torch.clamp(discriminant, min=0.0)
        
        # Take smaller root (actual J)
        j = (-b - torch.sqrt(discriminant)) / (2.0 * a + 1e-10)
        
        return j
    
    def rubisco_limited_rate(self, ci: torch.Tensor, vmax: torch.Tensor,
                              kc: torch.Tensor, ko: torch.Tensor,
                              gamma: torch.Tensor) -> torch.Tensor:
        """
        Calculate Rubisco-limited assimilation rate (Ac).
        
        Ac = Vmax * (Ci - Γ*) / (Ci + Kc*(1 + O/Ko))
        
        Args:
            ci: Intercellular CO2 (mol/mol)
            vmax: Maximum carboxylation rate (μmol m⁻² s⁻¹)
            kc: Michaelis constant for CO2 (mol/mol)
            ko: Michaelis constant for O2 (mol/mol)
            gamma: CO2 compensation point (mol/mol)
            
        Returns:
            Rubisco-limited rate Ac (μmol m⁻² s⁻¹)
        """
        numerator = vmax * (ci - gamma)
        denominator = ci + kc * (1.0 + OA / ko)
        
        ac = numerator / (denominator + 1e-10)
        
        return ac
    
    def electron_transport_limited_rate(self, ci: torch.Tensor, j: torch.Tensor,
                                         gamma: torch.Tensor) -> torch.Tensor:
        """
        Calculate electron transport-limited assimilation rate (Aj).
        
        Aj = J * (Ci - Γ*) / (4 * (Ci + 2*Γ*))
        
        Args:
            ci: Intercellular CO2 (mol/mol)
            j: Actual electron transport rate (μmol m⁻² s⁻¹)
            gamma: CO2 compensation point (mol/mol)
            
        Returns:
            RuBP regeneration-limited rate Aj (μmol m⁻² s⁻¹)
        """
        numerator = j * (ci - gamma)
        denominator = 4.0 * (ci + 2.0 * gamma)
        
        aj = numerator / (denominator + 1e-10)
        
        return aj
    
    def co_limited_rate(self, ac: torch.Tensor, aj: torch.Tensor) -> torch.Tensor:
        """
        Calculate co-limited assimilation using smooth minimum.
        
        Uses quadratic smoothing to avoid discontinuities:
            θ*A² - (Ac + Aj)*A + Ac*Aj = 0
        
        Args:
            ac: Rubisco-limited rate (μmol m⁻² s⁻¹)
            aj: Electron transport-limited rate (μmol m⁻² s⁻¹)
            
        Returns:
            Co-limited gross assimilation (μmol m⁻² s⁻¹)
        """
        a = self.thetac
        b = -(ac + aj)
        c = ac * aj
        
        discriminant = b * b - 4.0 * a * c
        discriminant = torch.clamp(discriminant, min=0.0)
        
        # Take smaller root
        agross = (-b - torch.sqrt(discriminant)) / (2.0 * a + 1e-10)
        
        return agross
    
    def assimilation(self, p: int, gmax: torch.Tensor, tleaf: torch.Tensor,
                     lavpd: torch.Tensor, eplant: torch.Tensor,
                     is_night: bool = False) -> None:
        """
        Calculate net photosynthesis at a given pressure index.
        
        This is the main assimilation calculation following Medlyn et al. (2002).
        It solves for the intercellular CO2 concentration that balances
        supply (through stomata) with demand (from photosynthesis).
        
        Args:
            p: Pressure index (for storing in curves)
            gmax: Maximum stomatal conductance (mmol m⁻² s⁻¹)
            tleaf: Leaf temperature (°C)
            lavpd: Leaf-to-air VPD (kPa/kPa normalized)
            eplant: Transpiration rate (mmol m⁻² s⁻¹)
            is_night: Whether it is nighttime
        """
        # Get stomatal conductance from transpiration and VPD
        if lavpd > 1e-10:
            gcanw = eplant / lavpd  # mmol m⁻² s⁻¹
        else:
            gcanw = gmax
        
        # Convert to CO2 conductance (mol m⁻² s⁻¹)
        # Factor 1.6 is the ratio of diffusivities (H2O/CO2)
        gcanc = (gcanw / 1.6) * 1000.0  # μmol m⁻² s⁻¹
        
        # Store conductances
        self.gcanw[p] = gcanw
        self.gcanc[p] = gcanc
        
        if is_night or gcanc <= 0:
            # Nighttime: only respiration
            rday = self.day_respiration(tleaf)
            self.psyn[p] = -rday
            self.cin[p] = self.ca
            return
        
        # Temperature corrections
        vmax = self.temperature_correction_vmax(tleaf)
        jmax = self.temperature_correction_jmax(tleaf)
        kc = self.temperature_correction_kc(tleaf)
        ko = self.temperature_correction_ko(tleaf)
        gamma = self.temperature_correction_gamma(tleaf)
        rday = self.day_respiration(tleaf)
        
        # Electron transport rate from light
        j = self.electron_transport_rate(self.qsl, jmax)
        
        # Solve for Ci iteratively
        # At steady state: A_biochem(Ci) = A_diffusion(Ci) = gcanc * (Ca - Ci)
        
        # Start from compensation point
        ci = gamma.clone()
        ca = self.ca
        
        # Iterate to find equilibrium Ci
        for _ in range(100):
            # Rubisco-limited rate
            ac = self.rubisco_limited_rate(ci, vmax, kc, ko, gamma)
            
            # Electron transport-limited rate
            aj = self.electron_transport_limited_rate(ci, j, gamma)
            
            # Co-limited gross assimilation
            agross = self.co_limited_rate(ac, aj)
            
            # Net assimilation
            anet = agross - rday
            
            # Diffusion-limited rate (supply through stomata)
            a_supply = gcanc * (ca - ci)
            
            # Check convergence
            if torch.abs(anet - a_supply) < 0.01:
                break
            
            # Adjust Ci
            # If A_biochem > A_supply, Ci is too high, decrease
            # If A_biochem < A_supply, Ci is too low, increase
            ci = ci + (a_supply - anet) * 1e-9
            ci = torch.clamp(ci, min=gamma, max=ca)
        
        # Store results
        self.psyn[p] = anet
        self.cin[p] = ci
        
        # Update maximum
        if anet > self.psynmax:
            self.psynmax.data = anet.clone()
    
    def assimilation_shade(self, p: int, gmax: torch.Tensor, tleaf: torch.Tensor,
                           lavpd: torch.Tensor, eplant: torch.Tensor,
                           is_night: bool = False) -> None:
        """
        Calculate net photosynthesis for shade leaves.
        
        Same as assimilation() but uses shade light level (qsh).
        """
        if lavpd > 1e-10:
            gcanw = eplant / lavpd
        else:
            gcanw = gmax
        
        gcanc = (gcanw / 1.6) * 1000.0
        
        self.gcanwsh[p] = gcanw
        self.gcancsh[p] = gcanc
        
        if is_night or gcanc <= 0:
            rday = self.day_respiration(tleaf)
            self.psynsh[p] = -rday
            self.cinsh[p] = self.ca
            return
        
        vmax = self.temperature_correction_vmax(tleaf)
        jmax = self.temperature_correction_jmax(tleaf)
        kc = self.temperature_correction_kc(tleaf)
        ko = self.temperature_correction_ko(tleaf)
        gamma = self.temperature_correction_gamma(tleaf)
        rday = self.day_respiration(tleaf)
        
        # Use shade light level
        j = self.electron_transport_rate(self.qsh, jmax)
        
        ci = gamma.clone()
        ca = self.ca
        
        for _ in range(100):
            ac = self.rubisco_limited_rate(ci, vmax, kc, ko, gamma)
            aj = self.electron_transport_limited_rate(ci, j, gamma)
            agross = self.co_limited_rate(ac, aj)
            anet = agross - rday
            a_supply = gcanc * (ca - ci)
            
            if torch.abs(anet - a_supply) < 0.01:
                break
            
            ci = ci + (a_supply - anet) * 1e-9
            ci = torch.clamp(ci, min=gamma, max=ca)
        
        self.psynsh[p] = anet
        self.cinsh[p] = ci
        
        if anet > self.psynmaxsh:
            self.psynmaxsh.data = anet.clone()
    
    def forward(self, tleaf: torch.Tensor, ppfd: torch.Tensor,
                gcanc: torch.Tensor, ca: Optional[torch.Tensor] = None,
                vmax25: Optional[torch.Tensor] = None,
                jmax25: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass computes net photosynthesis.
        
        This is a simplified interface for batch computation of
        net assimilation.
        
        Args:
            tleaf: Leaf temperature (°C)
            ppfd: Photosynthetic photon flux density (μmol m⁻² s⁻¹)
            gcanc: Canopy CO2 conductance (μmol m⁻² s⁻¹)
            ca: Atmospheric CO2 (mol/mol), uses self.ca if None
            vmax25: Maximum carboxylation rate at 25°C (μmol m⁻² s⁻¹), 
                   uses self.vmax25 if None
            jmax25: Maximum electron transport rate at 25°C (μmol m⁻² s⁻¹),
                   uses self.jmax25 if None
            
        Returns:
            Net assimilation rate (μmol m⁻² s⁻¹)
        """
        if ca is None:
            ca = self.ca
        if vmax25 is None:
            vmax25 = self.vmax25
        if jmax25 is None:
            jmax25 = self.jmax25
        
        # Temperature corrections using provided vmax25/jmax25
        vmax = self._temperature_correction_vmax_with_param(tleaf, vmax25)
        jmax = self._temperature_correction_jmax_with_param(tleaf, jmax25)
        kc = self.temperature_correction_kc(tleaf)
        ko = self.temperature_correction_ko(tleaf)
        gamma = self.temperature_correction_gamma(tleaf)
        rday = self.day_respiration(tleaf)
        
        # Electron transport
        j = self.electron_transport_rate(ppfd, jmax)
        
        # Solve for Ci (simplified) - reduced iterations for speed
        # Initial guess: 70% of Ca
        ci = 0.7 * ca
        
        for _ in range(10):
            ac = self.rubisco_limited_rate(ci, vmax, kc, ko, gamma)
            aj = self.electron_transport_limited_rate(ci, j, gamma)
            agross = self.co_limited_rate(ac, aj)
            anet = agross - rday
            a_supply = gcanc * (ca - ci)
            
            if torch.abs(anet - a_supply).max() < 0.1:
                break
            
            ci = ci + (a_supply - anet) * 1e-8
            ci = torch.clamp(ci, min=gamma, max=ca)
        
        return anet
    
    def _temperature_correction_vmax_with_param(self, tleaf: torch.Tensor, 
                                                 vmax25: torch.Tensor) -> torch.Tensor:
        """Temperature correction for Vmax using provided vmax25 parameter."""
        tk = tleaf + 273.15
        tk25 = 298.15
        
        numerator = (1.0 + torch.exp((self.svvmax * tk25 - self.hdvmax) / (GAS * tk25))) * \
                    torch.exp((self.havmax / (GAS * tk25)) * (1.0 - tk25 / tk))
        denominator = 1.0 + torch.exp((self.svvmax * tk - self.hdvmax) / (GAS * tk))
        
        vmax = vmax25 * numerator / (denominator + 1e-10)
        return vmax
    
    def _temperature_correction_jmax_with_param(self, tleaf: torch.Tensor,
                                                 jmax25: torch.Tensor) -> torch.Tensor:
        """Temperature correction for Jmax using provided jmax25 parameter."""
        tk = tleaf + 273.15
        tk25 = 298.15
        
        numerator = (1.0 + torch.exp((self.svjmax * tk25 - self.hdjmax) / (GAS * tk25))) * \
                    torch.exp((self.hajmax / (GAS * tk25)) * (1.0 - tk25 / tk))
        denominator = 1.0 + torch.exp((self.svjmax * tk - self.hdjmax) / (GAS * tk))
        
        jmax = jmax25 * numerator / (denominator + 1e-10)
        return jmax
    
    def clear_parameters(self) -> None:
        """Reset assimilation curves and states."""
        self.psynmax.data.fill_(0.0)
        self.psynmaxsh.data.fill_(0.0)
        self.psynact.data.fill_(0.0)
        self.psynactsh.data.fill_(0.0)
        self.atree.data.fill_(0.0)
        
        self.cinc.data.fill_(0.0)
        self.cincsh.data.fill_(0.0)
        
        self.qsl.data.fill_(0.0)
        self.qsh.data.fill_(0.0)
        
        self.cin.zero_()
        self.cinsh.zero_()
        self.psyn.zero_()
        self.psynsh.zero_()
        self.psynmd.zero_()
        self.psynshmd.zero_()
        self.psync.zero_()
        self.gcanw.zero_()
        self.gcanc.zero_()
        self.gcanwsh.zero_()
        self.gcancsh.zero_()
    
    def __repr__(self) -> str:
        return (f"CarbonAssimilationModel("
                f"Vmax25={self.vmax25.item():.1f}, "
                f"Jmax25={self.jmax25.item():.1f}, "
                f"Ca={self.ca.item()*1e6:.0f} ppm)")
