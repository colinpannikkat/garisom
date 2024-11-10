#ifndef CARBON_H
#define CARBON_H

class CarbonAssimilationModel {

    private:

    public:

        double lai,
               laish,
               laisl,
               ca,
               cinc,
               atree;

        /*
            Per Venturas et al. 2018.

            This function calculates the net photosynthesis gain, A_net, at every E 
            increment.
            
            This represents the gain of stomatal opening, by the amount of carbon
            ingested.

            The function calculates in the following sequence:

                1. Leaf temperature, T_l, is obtained from E, W, and leaf width
                L_w via the leaf energy budget.
                2. Leaf-air vapor pressure deficit, D_l, is calculated from T_l and
                D_air.
                3. Diffusive conductance to water vapor (G_w) is given by G_w = E/D,
                this is obtained from the E(P_c) curve.
                4. Diffusive conductance to CO2, G_c, is estimated as G_w / 1.6
                5. A_net is then calculated from a Farquhar-type model implemented by
                Medlyn et al. 2002. More details for this can be found in Sperry
                et al. 2017.

            Leaf respiration at 25c is set to V_max25 * 0.01

            We then derive A' as follows

                A' = \frac{A_net}{A_max}

            Such that A_max is the maximum A_net at that time step, usually at E_crit.

            There are separate gain functions computed for sun and shade layers
            following the light model of Campbell & Normal 1998, this requires L_ai.

            Inputs:
                - Maximum rates of carboxylation (at 25c): V_max25
                - Maxmimum rates of electron transport (at 25c): J_max25
                - Canopy leaf area index: L_ai
            
        */
    double &net_photosynthesis();


};

#endif