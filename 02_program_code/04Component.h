#ifndef COMPONENT_H
#define COMPONENT_H

#include <stdio.h>
#include <cmath>    // math utility functions
#include <string>   // the C++ String Class, easier to deal with than char arrays for this application
#include <cstring>  // for memcpy
#include <cstdint>
#include <ctime>    // timers for performance testing
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <assert.h>
#include "01Utilities.h"

/*
    This implementation follows Sperry et al. 2016

    Implemented by Colin Pannikkat
    Last updated: 10/05/24
*/

class Component {
    
    protected:

        double  b_wb = 0,           // b for weibull
                c_wb = 0,           // c for weibull
                res_percent = 0,    // percent resistances at ksat
                res_wb = 0,         // weibull resistance, 1/cond_wb
                cond_wb = 0,        // weibull conductance, 1/res_wb           
                k_max = 0,          // kmax, conductance per basal area
                                // ksat used interchangably
                p_crit = 0,         // critical pressure, of VC
                k_min = 0,
                pressure = 0;       // downstream pressure used in composite calculations
        std::vector<double> wb_fatigue;
        double e_p[CURVE_MAX] = {0};    // E(P) curve for component
        double e_pv[CURVE_MAX] = {0};   // Virgin E(P) curve
        double e_comp[CURVE_MAX] = {0};
        double e_pt[CURVE_MAX] = {0};   // Historical curve
        double k[CURVE_MAX] = {0};      // K (conductivity) curve for component
        double k_v[CURVE_MAX] = {0};    // Virgin K (conductivity) curve for component
        double k_comp[CURVE_MAX] = {0}; // New(conductivity) curve for component based on other components
        double k_t[CURVE_MAX] = {0};    // Historical curve
        double pressure_comp[CURVE_MAX] = {0};
        double pressure_v[CURVE_MAX] = {0}; // virgin pressure curve

    public:

        /* Getters */
        double getBwb();
        double getCwb();
        double getResPercent();
        double getResWb();
        double getCondWb();
        double getKmax();
        double getPcrit();
        double getKmin();
        double getPressure();
        double& getFatigue(int index);
        double& getPressureComp(int index);
        double& getPressureVirgin(int index);
        double& getEp(int index);
        double& getEpVirgin(int index);
        double& getEComp(int index);
        double& getK(int index);
        double& getKVirgin(int index);
        double& getKComp(int index);

        /* Setters */
        void setBwb(double value);
        void setCwb(double value);
        void setResPercent(double value);
        void setResWb(double value);
        void setCondWb(double value);
        void setKmax(double value);
        void setPcrit(double value);
        void setKmin(double value);
        void setPressure(double value);
        void setFatigue(int index, double value);
        void setPressureComp(int index, double value);
        void setPressureVirgin(int index, double value);
        void setEp(int index, double value);
        void setEComp(int index, double value);
        void setEpVirgin(int index, double value);
        void setK(int index, double value);
        void setKVirgin(int index, double value);
        void setKComp(int index, double value);

        /* Storage */
        void storeTranspirationCurve();
        void storeConductivityCurve();
        void storeCurves();
        void storeTranspirationCurveAndUseVirgin();
        void storeConductivityCurveAndUseVirgin();
        void storeCurvesAndUseVirgin();

        void restoreTranspirationCurve();
        void restoreConductivityCurve();
        void restoreCurves();
        void clearHistoricalCurves();

        /* Clean parameters */
        virtual void cleanParameters();

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
        double wb(const double &pressure);

        /*
        
            Calculates steady-state flow rate, E_i(P), for a component

            E_i = \int_{P_up}^{P_down}{k(P)_i}dP

            Inputs:
                - p_inc: Increment for pressure in Reimann integration
                - k_min: Minimum conductivity amount

            Outputs:
                - Steady-state flow rate: e_p
        */ 
        void calc_flow_rate(const double &p_inc, const double &k_min, bool virgin=false);
        void trapzd(const double &p1, const double &p2, double &s, const int &t, int &it);
        void virtual qtrap(double &p1, double &p2, double &s);
        void virtual calc_through_flow(const double &p1, 
                                       const double &p2, 
                                       const double &p_inc,
                                       double &flow,
                                       double &klower,
                                       double &kupper);
        int virtual calc_pressure(const double &e, const double &bottom_pressure, const double &p_grav, const double &p_inc);

        void printCurveToFile(const double &p_inc, const std::string &filename) const {
            std::ofstream outFile(filename);
            if (!outFile) {
                std::cerr << "Error opening file: " << filename << std::endl;
                return;
            }

            outFile << "p_inc,E(P)" << std::endl;
            for (int i = 0; p_inc * i <= this->p_crit; ++i) {
                outFile << p_inc * i << "," << e_p[i] << std::endl;
            }

            outFile.close();
        }
};

#endif