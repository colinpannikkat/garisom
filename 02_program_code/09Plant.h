#ifndef PLANT_H
#define PLANT_H

#include "01Utilities.h"
#include "08Xylem.h"
#include "03Hydraulics.h"
#include "10Carbon.h"
#include "11Parameters.h"

/*

    Per Sperry et al. 2017.

    "The plant is represented by canopy (all leaves in parallel) and stem (all 
    stems in parallel) elements in via in-series rhizosphere elements to bulk 
    soil in five horizontal series, connected to five root elements in parallel, 
    each connected layers."

    The in-series components are contained within XylemComponent.

*/

class Plant {

    public:

        Plant(int stage_id);

        // Model Components
        XylemComponent xylem;
        HydraulicsModel hydraulics;
        CarbonAssimilationModel carbon;
        Parameters param;       // all running parameters
        CSVData<std::string> config_data;
        CSVData<std::string> param_data;
        CSVData<double> gs_data;
        CSVData<double> data;           // up to 2000k rows and 100 columns of input data

        // Model Variables
        std::string stageNames[STAGE_ID_FUT_STRESS_NOACCLIM + 1] = {"standard"},
                    growing_season_limits_data_path, 
                    climate_forcing_data_path, 
                    data_header_file_path, 
                    night;      // know whether is night

        int         species_no = 0, 
                    stage_id = 0,
                    unknowns = 0,
                    halt = 0,
                    haltsh = 0,         // track when hydraulics fail
                    layers = 0,
                    tod = 0;

        bool        ground = false,
                    oldground = false, 
                    soilred = false, 
                    sevap = false, 
                    raining = false,
                    oldraining = false, 
                    useGSData = false, 
                    mode_predawns = false, 
                    refilling = false, 
                    hysteresis = false,
                    stem_only = false,
                    iter_gwEnable = false,
                    iter_ffcEnable = false,
                    iter_bagaEnable = false,
                    iter_useAreaTable = false,
                    iter_yearsAsCount = false,
                    iter_runSupplyCurve = false,
                    isNewYear = false,
                    gs_inGrowSeason = false,
                    gs_doneFirstDay = false, //done all first day of grow season calculations?
                    rainEnabled = false,
                    oldrainEnabled = false;

        double      iter_gwInc = 0.0,
                    iter_gwStart = 0.0,
                    iter_gwEnd = 0.0,
                    iter_ffcInc = 0.0,
                    iter_ffcStart = 0.0,
                    iter_ffcEnd = 0.0,
                    iter_bagaInc = 0.0,
                    iter_bagaStart = 0.0,
                    iter_bagaEnd = 0.0,
                    iter_bagaRef = 0.0,
                    iter_bagaCutoff = 0.0,
                    max_plc_x = 0.0,
                    gwflow = 0.0,             // Used like twice, needed for soilwetness and deepflow
                    drainage = 0.0,           // used in deepflow
                    runoff = 0.0,
                    soilevap = 0.0,           // set in soilevaporation(), used in soilwetness and soilevaporation, maybe just declare local in modelTimeStepIter and pass?
                    iter_refK = 0.0,          // used in getsoilwetness, never initialized though?
                    kpday1 = 0.0,
                    kxday1 = 0.0,             // this and above used in getsoilwetness
                    transpiration_tree = 0.0, // used in multiple functions, need to preserve old value
                    transpiration = 0.0, 
                    transpirationsh = 0.0, 
                    md = 0.0, 
                    mdsh = 0.0,                // need to preserve this and three above so during night we still retain values
                    rmean = 0.0,
                    // waterold,        // read from previous dd
                    b_fatigue[3][10] = {{0.0}},   // [0, :] is roots, [1, :] is stem, [2, :] is leaves
                    ecritsystem = 0.0,
                    pcritsystem = 0.0,
                    kmin = 0.0,               // k-min for the plant
                    e_p[CURVE_MAX] = {0.0},    // whole plant transpiration curve, likely not needed now since only one xylem is typically used
                    albedo_soil = 0.2,
                    absorption_short = 0.5;

        long        gs_yearIndex = 0,       // this is a counter from 0 (for the first year) indicating how many years have passed
                        // get the actual year from gs_ar_years(gs_yearIndex)
                        // this avoids having to make year an input column in the model -- will just count how many years have passed when running
                    gs_prevDay = 0,
                    year_cur = 0, 
                    year_start = 0,         // not to be confused with the year array index, which is year_cur - year_start
                    yearVal = 0;            // the temporary variable where we hold the year read from the sheet

        std::vector<double> water = {},
                            fc = {};
        std::vector<std::vector<double>> jmatrix = {};

        /* 
            Per Venturas et al. 2018.

            At every time step, optimization algorithm finds the canopy pressure
            P_c that maximizes difference between carbon gain (A') and hydraulic
            risk (\theta').

            This model assumes that plants are profit-maximizers, attempting to
            maximize their net gain out of carbon gain and hydraulic risk via
            transpiration and corresponding canopy pressure P_c. Therein, the
            stomata maintain a P_c where the marginal normalized carbon gain
            \dd{A'}{P_c} is equal to the marginal normalized risk of hydraulic
            failure \dd{\theta'}{\P_c}.

            The model is iterated in timesteps designated in the provided
            model parameters.

            For each timestep, the following boundary conditions were used.
                - Bulk soil water potential: P_s
                - Vapor pressure deficit: D_air
                - Solar radiation: W
                - Windspeed: U_wind
                - Air temperature: T_air
                - Soil surface temperature: T_soil

            At each timestep, the transpiration rate (E) is increment from zero
            to E_crit (hydraulic failure). At each of these increments the
            transpiration rate/pressure drop is calculated via integral transform
            based on the conductivity function of rhizosphere or xylem element
            vulnerability curves. By incrementing E, it yields the steady-state 
            relatioship with canopy pressure, P_c.

            All fluxes and conductances were scaled from leaf-to-basal-to ground
            area using leaf area per basal area, La_Ba, and basal area per
            ground area, Ba_Ga.

            The derivative of the E(P_c) supply function is proportional to
            hydraulic conductance of the downstream end of flow path. We then
            calculate the risk function as:

                \theta' = 1 - \frac{K_c}{K_cmax}

            Within each timestep, the model also calculates net photosynthesis,
            A_net, at every E increment. For more information see 
            CarbonAssimilationModel::net_photosynthesis() in '10Carbon.h'. The
            function outputs a normalized gain function A'.

            Then via the provided \theta' and A' the model calculates:

                max(A' - \theta') w.r.t to P_c

            This maximizes the gain-risk difference w.r.t to P_c.

            With xylem refilling on, K(P) functions have no hysteresis. With it 
            off fluxes and conductances are reduced to a permanent loss of 
            conductance.

            After each time step, bulk soil water content and P_s boundary
            conditions were recalculated for each soil layer. See paper for
            further explanation.

            Each time step outputs:
                - Canopy pressure: P_c
                - Transpiration rate: E
                - Net photosynthesis: A_net
                - Leaf diffusive conductance: G_w
                - Leaf temperature: T_l
                - Tree sapflow per basal area, S_f
                - Soil-canopy hydraylic conductance per basal area: K_p

            Inputs:
                - Solar radiation: W
                - Air temperature: T_air
                - Air vapor pressure deficit: D_air
                - Wind speed: U_wind
                - Soil temperature: T_soil
                - Predawn pressure: P_pd
                - Leaf area index: L_ai
                - Max carboxylation rate: V_max25
                - Max electron transport rate: J_max25
                - Percent resistance in leaf: R_l
                - Leaf vulnerability curve: VC_l
                - Stem vulnerability curve: VC_s
                - Leaf to basal area ratio: La_Ba
                - Maximum plant conductance: K_max
                - Basal area to ground area ratio: Ba_Ga
                - Avg percent resistance in rhizosphere: R_rh
                - Root vulenrability curve: VC_r

            Outputs:
                - Midday pressure: P_md
                - Transpiration: E
                - Net assimilation: A_net
                - Leaf diffusive conductance: G_w
                - Leaf temperature: T_l
                - Sap flow: S_f
                - Whole-plant hydraulic conductance: K_p
        
        */
        double &model();

        void setConfig(int config_setting);
        void initModelVars();
        void cleanModelVars();
        void readin();
        void resetLayerStatus();
        void componentPCrits();
        int modelTimestepIter(int &dd);
        double getCarbonByYear(int yearVal);
        void modelProgramNewYear();
        bool isInGrowSeasonSimple(const int &jd);
        void getsoilwetness(const int &dd, 
                            const int &timestep,
                            const double &lai,
                            const double &laish,
                            const double &laisl,
                            const double &laperba,
                            const double &atree,
                            const double &cinc,
                            const double &ca);
        void solarcalc(const int &dd,
                       const int &jd,
                       const double &obssolar,
                       const double &maxvpd,
                       const double &airtemp,
                       const double &vpd,
                       double &lai,
                       double &laisl,
                       double &laish,
                       double &qsl,
                       double &qsh,
                       double &ssun,
                       double &sshade,
                       double &sref,
                       double &la,
                       double &lg);
        int getpredawns(const int &dd);
        int newtonrhapson(const int &dd, const double &p_inc, const double &e);
        int compositeCurve(const double &e, const int &p, int &total);
        void canopypressure(const int &dd,
                           const int &total,
                           double &transpiration,
                           double &transpirationsh,
                           double &md,
                           double &mdsh);
        void updatecurves(const int &halt);
        void soilflow();
        void deepflow(const double &timestep);
        void soilevaporation(const double &soiltemp,
                             const double &maxvpd,
                             const double &vpd,
                             const double &airtemp,
                             const double &us);
        void storehistory();
        void gethistory();
};

#endif