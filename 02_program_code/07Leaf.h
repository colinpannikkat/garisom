#ifndef LEAF_H
#define LEAF_H

#include "04Component.h"

class LeafComponent : public Component {
    
    public:

        double lavpd[CURVE_MAX] = {0.0},        // Leaf-to-air vapor pressure deficit
               lavpdsh[CURVE_MAX] = {0.0},      // Above but shade
               leaftemp[CURVE_MAX] = {0.0},     // Leaf temp in C
               leaftempsh[CURVE_MAX] = {0.0},
               eplantl[CURVE_MAX] = {0.0};      // Transpiration of plant in leaf area

        // Midday attributes, same as above
        double emd = 0.0,
               lavpdmd = 0.0,
               lavpdshmd = 0.0,
               leaftmd = 0.0,
               leaftshmd = 0.0;
        
        // Parameters used in solarcalc()
        double ssun = 0.0,
               sshade = 0.0,
               sref = 0.0,
               sbottom = 0.0,
               emiss = 0.0,
               la = 0.0,
               lg = 0.0,
               lai = 0.0,
               laisl = 0.0,
               laish = 0.0,
               lambda = 0.0,
               grad = 0.0,
               gha = 0.0;

        void cleanParameters() override;

        double calc_diffusion_coefficient(const double &leaftemp,
                                          const double &airtemp,
                                          const double &p,
                                          const double &D_0);

        double calc_nu(const double &leaftemp,
                       const double &airtemp,
                       const double &windspeed,
                       const double &leafwidth);
        
        void temp(const int p,
                  const double &r_sw,
                  const double &al_s,
                  const double &a_s,
                  const double &a_l,
                  const double &airtemp,
                  const double eplant[],
                  const double &vpd,
                  const double &wind,
                  const double &laperba,
                  const double &leafwidth,
                  const double &patm);
        void tempShade(const int &p,
                       const double &airtemp,
                       const double &patm,
                       const double &vpd);
        void tempMd(const int p,
                  const double &r_sw,
                  const double &al_s,
                  const double &a_s,
                  const double &a_l,
                  const double &e,
                  const double &airtemp,
                  const double &vpd,
                  const double &wind,
                  const double &laperba,
                  const double &leafwidth,
                  const double &patm);
        void tempShadeMd(const int &p,
                           const double &e,
                           const double &airtemp,
                           const double &vpd,
                           const double &patm); 
};

#endif