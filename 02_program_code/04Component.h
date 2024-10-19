#ifndef COMPONENT_H
#define COMPONENT_H

#include <stdio.h>
#include <cmath>    // math utility functions
#include <string>   // the C++ String Class, easier to deal with than char arrays for this application
#include <ctime>    // timers for performance testing
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <assert.h>

/*
    This implementation follows Sperry et al. 2016

    Implemented by Colin Pannikkat
    Last updated: 10/05/24
*/

class Component {
    
    private:

        double  b_wb,           // b for weibull
                c_wb,           // c for weibull
                res_percent,    // percent resistances at ksat
                res_wb,         // weibull resistance, 1/cond_wb
                cond_wb,        // weibull conductance, 1/res_wb           
                k_max;          // kmax, conductance per basal area
                                // ksat used interchangably
        std::vector<double> wb_fatigue;

    public:

        /* Getters */
        double getBwb();
        double getCwb();
        double getResPercent();
        double getResWb();
        double getCondWb();
        double getKmax();
        double& getFatigue(int index);

        /* Setters */
        void setBwb(double value);
        void setCwb(double value);
        void setResPercent(double value);
        void setResWb(double value);
        void setCondWb(double value);
        void setKmax(double value);
        void setFatigue(int index, double value);

        /* 
            Weibull function.
        */
       double wb(double pressure);

        /* 

            Calculates hydraulic conductance of component, via a weibull function.

            Inputs:
                - Maximum flow rate per pressure drop: K_max
                - Negative sap pressure (-MPa): P
                - Curve shift: b_wb
                - Shape: c_wb
                    Exponential corresponds to c <= 1, sigmoidal has c > 1

            Outputs:
                - Hydraulic conductance of component: k

        */
        double& calc_hydraulic_conductance(double *K_max, double &P, double &b_wb, double &c_wb);

        /*
        
            Calculates steady-state flow rate, E_i, for a component

            E_i = \int_{P_up}^{P_down}{k(P)_i}dP

            Inputs:
                - Downstream pressure: P_down
                - Upstream pressure: P_up

            Outputs:
                - Steady-state flow rate: E_i
        */ 
       double &calc_flow_rate(double &P_down, double &P_up);
};

#endif