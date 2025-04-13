#ifndef PARAMETERS_H
#define PARAMETERS_H

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
#include "01Utilities.h"

#define MAX_YEARS 90

class Parameters {
    private:
        /* General model parameters, copied from Plant::param_data */
        std::unordered_map<std::string, double> model_parameters;

        /* Growing season data */
        std::vector<int> gs_ar_year;
        std::vector<int> gs_ar_start;
        std::vector<int> gs_ar_end;
        std::vector<double> gs_ar_ppm;

        /* Saved weibull params */
        double stem_b_wb[MAX_YEARS];
        double root_b_wb[MAX_YEARS];
    
    public:

        Parameters();

        double& getModelParam(const std::string &param_name);
        double& operator()(const std::string &param_name);

        void setModelParam(double value, const std::string &param_name);

        int& getGsArYear(int index);
        void setGsArYear(int index, int value);

        int& getGsArStart(int index);
        void setGsArStart(int index, int value);

        int& getGsArEnd(int index);
        void setGsArEnd(int index, int value);

        double& getGsArPpm(int index);
        void setGsArPpm(int index, double value);

        double &getStemBWb(int index);
        void setStemBWb(int index, double value);

        double &getRootBWb(int index);
        void setRootBWb(int index, double value);
};

#endif