#ifndef LEAF_H
#define LEAF_H

#include "04Component.h"

class LeafComponent : public Component {
    
    public:

        double lavpd[CURVE_MAX],        // Leaf-to-air vapor pressure deficit
                                        // 'saturated mole fraction'
               lavpdsh[CURVE_MAX],      // Above but shade
               leaftemp[CURVE_MAX],     // Leaf temp in C
               leaftempsh[CURVE_MAX],
               eplantl[CURVE_MAX];      // Transpiration of plant in leaf area
                                        // mol m-2s-1 (per leaf area)

        // Midday attributes, same as above
        double emd,
               lavpdmd,
               lavpdshmd,
               leaftmd,
               leaftshmd;
        
        // Parameters used in solarcalc()
        double ssun,
               sshade,
               sref,
               sbottom,
               emiss,
               la,
               lg,
               lai,
               laisl,
               laish,
               lambda,
               grad,
               gha;

        void cleanParameters() override;
        
        void temp(const int p,
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
        void tempMd(const int &p,
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