#ifndef LEAF_H
#define LEAF_H

#include "04Component.h"

class LeafComponent : public Component {
    
    public:

        double lavpd[CURVE_MAX],
               lavpdsh[CURVE_MAX],
               leaftemp[CURVE_MAX],
               leaftempsh[CURVE_MAX],
               eplantl[CURVE_MAX];

        double emd,
               lavpdmd,
               lavpdshmd,
               leaftmd,
               leaftshmd;
        
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