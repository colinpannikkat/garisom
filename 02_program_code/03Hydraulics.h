#ifndef HYDRAULICS_H
#define HYDRAULICS_H

class HydraulicsModel {
    
    public:

        /*
            Per Sperry and Love 2015.

            Calculates supply transpiration rate, E.

            Potential rate of water supply for transpiration (E), is a function
            of canopy xylem pressure (potential).
        
            The E(P_canopy) supply function depends on how the hydraulic 
            conductance (K) of soil and plant declines with more negative 
            sap pressure (P).

            We must derive, soil-canopy conductance (k) from the hydraulic 
            conductivity (K). Conductance can be expressed relative to a
            reference cross-sectional area. All k, E, and G are expressed per
            basal ground area (BAGA).

            If there is no k_max provided, we use allometric equations to assume
            that an individual leaf specific conductance can extrapolate to the 
            whole conductance of the tree.

            k(P_canopy) = \dd{E}{P_canopy}
            E_i = \int_{p_0}^{p_i}{k(P)}dP

            This function takes in canopy pressure and soil conductance and
            performs an integration over the full range of P_i, holding P_0
            constant.

            E will monotonically increase from zero to E_crit, as k approaches
            zero.

            k(P) decline is represented by a two-parameter Weibull function,
            since approximation may not reach zero the minimum relevant 
            K_cont = 0.01 kg h^{-1} MPa^{-1} m^{-2}. It falls from a K_max as
            P_canopy falls.

            P_0 holds the predawn canopy sap pressure (assume no nocturnal 
            transpiration). Intercept also incorporates rooted soil mositure
            profile (-0.5MPa), and gravitational gradient (-0.01 MPA m^{-1}).

            E/(\delta{P}) slope, where \delta{P} = (predawn - P_canopy) provides 
            the hydraulic conductance of the soil-canopy path during transpiration 
            when there is a \delta{P} soil-to-canopy drop. This is greater than 
            k(P_canopy) as sap pressure is less negative than the extreme across 
            continuum. The \dd{E}{P_canopy} provides local loss at the canopy 
            end of flow path.

            Parameters:
                - Canopy xylem sap pressure (-MPa): P_canopy
                - Soil-canopy hydraulic conductance (kg h^{-1} MPa^{-1} m^{-2}): K_cont

            Outputs:
                - Transpiration rate (kg h^{-1} m^{-2}): E
        */
        double& transpiration_supply(double &P_canopy, double &K_cont);

        /*
            Per Sperry and Love 2015.

            Calculates canopy diffusive conductance, G and transpiration rate
            demand, E.

            Stomatal and canopy diffusive conductances decrease in response
            to water stress.

            Transpiration loss function specifies where plant regulates the
            transpiration rate along E(P_canopy) supply function. Regulation
            is achieved via control of canopy diffusive conductance (G_canopy).

            Runaway cavitation (hydraulic failure) where E > E_crit only occurs
            under extreme soil drought conditions. Stable cavitation is where
            E <= E_crit.

            The transpiration loss function is dependent on supply. We use the
            derivative, \dd{E}{P_canopy}, to drive loss. This is also the 
            soil-canopy conductance. Maxmium canopy diffusive conductance G_max 
            represents maximal stomatal opening.

            Canopy diffusive conductance (kg h^{-1} m^{-2}), G, is a function of 
            stomatal aperture, boundary layer resistance, and other properties 
            of the leaf such as surface area.

            Vapor pressure deficit (kPa), D, is the difference between moisture in the
            air and the moisture air can hold when saturated. Can be defined as
            e_leaf (saturated vapor pressure inside of leaf) - e_air (actual
            vapor pressure of surrounding air).

            E' = G_max * D, this is unregulated transpiration
            E = G * D, this is regulated transpiration

            Drop in regulated pressure \delta{P} is calculated via the fact that
            stomatal regulation is assumed to reduce unregulated drop in pressure 
            \delta{P'}. Using this \delta{P} we can find E, and from E / D, find
            G.

            \delta{P} = \delta{P'} * (\dd{E'}{P_canopy'})/(\dd{E}{P_max}),
            this \delta{P} reaches a maxmium before  E' = E_crit.

            Since plant wants to keep E and P_canopy constant, stomata are
            sensitive towards D, therein changing stomatal opening to maintain
            balance.

            G_max obtained via yielding observed E under k_max conditions.

            Inputs:
                - Maxmium canopy diffusive conductance (kg h^{-1} m^{-2}): G_max
                - Atmospheric vapor pressure deficit (kPA): D

            Outputs:
                - Transpiration rate (kg h^{-1} m^{-2}): E
                - Canopy diffusive conductance (kg h^{-1} m^{-2}): G
        
        */
        double& transpiration_demand(double &G_max);

        /*
            Per Sperry et. al. 2016.

            Models relationship between transpiration rate (E) and canopy xylem
            pressure (P_canopy), given a constant soil-water potential.

            Stomatal demand function locates plant on the supply function.

            Demand calculations can be broken down into five steps:
                1. Find unregulated E' = DG_max on supply function
                2. Find ratio of derivatives for conductance as unregulated P'
                   and P_max
                3. Find regulated pressure drop given unregulated pressure and 
                   ratio of derivatives
                4. Regulated E corresponding to regulated \delta{P} is found on 
                   supply function
                5. G is then solved via E / D

            Determines transpiration gain and loss functions, and then finds the
            intersection.

            Is able to run in both reversible and irreversible cavitation modes.
                - Reversible: use original vulnerability curves
                - Irreversible: vulnerability curves change after drop in conductance

            For irreversible, k = k(P_min) for P = 0 to P_min. For P more negative
            than P_min, k(P) is recalculated. Even if set to irreversible, the
            demand function uses original uncavitated supply for calculating
            different in regulated pressure, \delta{P}.

            Regulated E, E, is calculated from current supply function.

            When in irreversible mode, model is initialized to prior exposure,
            and Rhizosphere vulnerability curves are still fully reversible.

            Since hydraulic conductance is expressed per trunk basal area,
            transpiration rate and canopy diffusive conductance is also expressed
            per trunk area.

            Roots and rhizosphere components are partioned into N paths draining
            horizontal soil layers. Layer depths are set as follows:

            \Beta = 1 - \beta^d, where \Beta is fraction of biomass above depth
            d in cm, and 0 < \beta < 1. Max root depth set at \Beta = 0.995.

            "The rhizosphere k_max for the whole root system was
            partitioned equally among the N layers. Total root system k_max
            was divided among layers in proportion to the inverse of the
            transport distance to each layer. Transport distance was depth to
            the center of layer biomass plus the radial spread of roots within
            each layer. The radial spread for the top layer was calculated by
            multiplying maximum root depth by an inputed aspect ratio of maximum 
            radial spread divided by maximum root depth. The rooted soil volume 
            in the top layer was calculated as a cylinder from spread and layer 
            thickness. By assuming this volume was constant for each layer, the 
            radial root spread in deeper soil layers was calculated." - Sperry
            et al. 2016

            Inputs:
                - Predawn soil water potential: P_soil
                - Vapor pressure deficit: D
                - Maxmium canopy diffusive conductance (kg h^{-1} m^{-2}): G_max

            Outputs:
                - Canopy transpiration (kg h^{-1} m^{-2}): E
                - Canopy diffusive conductance (kg h^{-1} m^{-2}): G
                - Canopy xylem pressure (-MPa): P_canopy

        */
       double& soil_plant_atmosphere();
};

#endif