#ifndef UTILITIES_H
#define UTILITIES_H

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
#include "11Parameters.h"

// input file names
#define CONFIG_FILE_PATH "../03_test_data/configuration.csv"
#define PARAMETER_FILE_PATH "../03_test_data/parameters.csv"

// output precision
#define FIO_PRECISION 12

// staging
#define STAGE_ID_NONE 0 // use this mode to skip all of the stage code and just run based on param sheet settings (a "normal" run, the default config)
#define STAGE_ID_HIST_OPT 1
#define STAGE_ID_HIST_STRESS 2
#define STAGE_ID_FUT_OPT 3
#define STAGE_ID_FUT_STRESS 4
#define STAGE_ID_FUT_STRESS_NOACCLIM 5

// data arrays
#define MAX_SUMMARY_COLS 121
#define MAX_SUMMARY_ROWS 2001
#define PARAMFILE_MAXROWS 2001
#define PARAMFILE_MAXCOLS 101
#define DATAFILE_MAXROWS 2000001
#define DATAFILE_MAXCOLS 101
#define CONFIGFILE_MAXROWS 10
#define CONFIGFILE_MAXCOLS 101

// constants
#define PROFT_MAX_RUN_MEAN 1        // running mean for profit maximization
#define DPA_MAX_CUTOFF 1            // cutoff for stopping dpamax search
#define MIN_WIND_THRESH 0.4515      // m s-1'minimum wind threshold
#define TRAP_ITER_MAX 70            // itmax for trapzd routine used for xylem only
#define EPSX 0.0001                 // acceptable error for e integral for xylem
#define GMAX 1000000                // maximum G, wet soil, vpd = 0, in kg m - 2 hr - 1 basal area
#define PI 3.14159
#define SBC 0.0000000567            // Stefan-Boltzmann constant in W m-2 K-4
#define SHA 29.3                    // Specific heat of air in J mol-1C-1
#define GAS 8.3144598               // Universal gas constant J mol-1K-1
#define OA 0.21                     // Mole fraction of O2
#define SOLAR 1362                  // Solar constant W m-2
#define ABS_OLAR 0.5                // Absorptivity of solar for leaves
#define ABS_PAR 0.8                 // Absorptivity of PAR for leaves
#define ABS_NIR 0.2                 // Absorptivity of near infrared for leaves

class CSVData {
    private:
        std::unordered_map<std::string, int> header;
        std::vector<std::vector<std::string>> data;

    public:
        CSVData() {};
        CSVData(const std::string &data_file_name);
        CSVData(const std::string &data_file_name, const std::string &header_file_name, bool skip_first=true);

        void readHeader(std::ifstream &dataFile);
        void readData(std::ifstream &dataFile);

        bool empty();
        void print(size_t page_size);
        int row_size();
        int col_size();
        
        std::string getColumnValue(const std::string& column_name, int row=0);
        bool getColumnValue(std::string &store, const std::string &column_name, int row=0);
        bool getColumnValue(double &store, const std::string &column_name, int row=0);
        bool getColumnValue(bool &store, const std::string& column_name, int row=0);
        bool getColumnValue(int &store, const std::string &column_name, int row=0);

        bool setColumnValue(std::string value, int row, const std::string &column_name);

};

bool locateRanges(CSVData &config_data, CSVData &param_data);
void readGSSheet(CSVData &gs_data, std::string &gs_file_name);
void readGrowSeasonData(Parameters &param, CSVData &gs_data);
void readDataSheet(CSVData &data, CSVData &sum_data, std::string &data_file_name, std::string &header_file_name, std::string &sum_head_file_name);

#endif