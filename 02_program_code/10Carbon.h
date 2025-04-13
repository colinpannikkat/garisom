#ifndef CARBON_H
#define CARBON_H

#include "01Utilities.h"

class CarbonAssimilationModel {

    private:

    public:

        double ca,
               cinc,
               cincsh,
               atree,
               qsl,
               qsh,
               psynmax,
               psynmaxsh,
               psynact,
               psynactsh;

        double gcanwmd,
               gcanwshmd,
               gcancmd,
               gcancshmd,
               cinmd,
               cinshmd,
               psynmaxmd,
               psynmaxshmd,
               gcmd,
               gcmdsh;

        double cin[CURVE_MAX],
               cinsh[CURVE_MAX],
               psyn[CURVE_MAX],
               psynsh[CURVE_MAX],
               psynmd[CURVE_MAX],
               psynshmd[CURVE_MAX],
               psync[CURVE_MAX],
               gcanw[CURVE_MAX],
               gcanc[CURVE_MAX],
               gcanwsh[CURVE_MAX],
               gcancsh[CURVE_MAX];

        void clearParameters();

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
        double &net_photosynthesis(); // this function is empty for now, the actual calculations are done in canopypressure()

        void assimilation(const int &p,
                          const double &gmax,
                          const double &qmax,
                          const double &comp25,
                          const double &thetac,
                          const double &vmax25,
                          const double &jmax25,
                          const double &kc25,
                          const double &ko25,
                          const double &svvmax,
                          const double &svjmax,
                          const double &havmax,
                          const double &hdvmax,
                          const double &hdjmax,
                          const double &hajmax,
                          const double &lightcurv,
                          const std::string &night,
                          const double eplantl[],
                          const double lavpd[],
                          const double leaftemp[]);
                          
        void assimilationShade(const int &p,
                               const double &gmax,
                               const double &qmax,
                               const double &comp25,
                               const double &thetac,
                               const double &vmax25,
                               const double &jmax25,
                               const double &kc25,
                               const double &ko25,
                               const double &svvmax,
                               const double &svjmax,
                               const double &havmax,
                               const double &hdvmax,
                               const double &hdjmax,
                               const double &hajmax,
                               const double &lightcurv,
                               const std::string &night,
                               const double eplantl[],
                               const double lavpd[],
                               const double leaftemp[]);

        void assimilationMd(const int &p,
                                           const double &gmax,
                                           const double &qmax,
                                           const double &comp25,
                                           const double &thetac,
                                           const double &vmax25,
                                           const double &jmax25,
                                           const double &kc25,
                                           const double &ko25,
                                           const double &svvmax,
                                           const double &svjmax,
                                           const double &havmax,
                                           const double &hdvmax,
                                           const double &hdjmax,
                                           const double &hajmax,
                                           const double &lightcurv,
                                           const std::string &night,
                                           const double emd,
                                           const double lavpdmd,
                                           const double leaftmd);

        void assimilationShadeMd(const int &p,
                                const double &gmax,
                                const double &qmax,
                                const double &comp25,
                                const double &thetac,
                                const double &vmax25,
                                const double &jmax25,
                                const double &kc25,
                                const double &ko25,
                                const double &svvmax,
                                const double &svjmax,
                                const double &havmax,
                                const double &hdvmax,
                                const double &hdjmax,
                                const double &hajmax,
                                const double &lightcurv,
                                const std::string &night,
                                const double emd,
                                const double lavpdshmd,
                                const double leaftshmd);
};

#endif