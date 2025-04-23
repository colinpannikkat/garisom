#include "09Plant.h"

Plant::Plant(int stage) {
    stage_id = stage;
    std::cout.precision(FIO_PRECISION);
    srand(42);
}

/**
 * @brief Configures the Plant model by setting up various parameters and modes
 *        based on the input configuration data. This includes settings for 
 *        plant community, soil, climate, hydraulics, BA:GA optimization, and 
 *        file locations.
 * 
 * @details
 * The function performs the following tasks:
 * 
 * - **Plant Community**:
 *   - Configures multi-species mode and sets the number of species.
 * 
 * - **Soil**:
 *   - Configures groundwater flow, soil water redistribution, and soil water evaporation.
 * 
 * - **Climate**:
 *   - Configures rainfall inputs, multiple growing seasons, and growing season data usage.
 * 
 * - **Hydraulics**:
 *   - Configures pre-dawn water potential mode, xylem refilling mode, xylem hysteresis, 
 *     and cavitation fatigue in roots.
 * 
 * - **File Locations**:
 *   - Sets paths for climate forcing data, growing season data, and time-step header file.
 * 
 * @note The function outputs the configuration status to the console for debugging purposes.
 * 
 * @param None
 * @return void
 */
void Plant::setConfig(int config_setting) { // sets up model configuration
    // // Model Configurations
    // Plant Community
    std::cout << "Plant Community ----------------------------------" << std::endl;
    std::cout << "            Multi-Species Mode: ";
    // Are we working with multiple species? Values: y; n
    if (config_data.getColumnValue("i_multipleSP", config_setting) == "y")
    {
        std::cout << "On" << std::endl;
        config_data.getColumnValue(species_no, "i_speciesN", config_setting);
        species_no -= 1;
        std::cout << "MODE: Setting species number to: " << species_no + 1 << std::endl; 
    } else {
        std::cout << "Off" << std::endl;
        species_no = 0; // default
        std::cout << "MODE: Setting species number to: " << species_no + 1 << std::endl;
    }
    
    // Soil
    std::cout << std::endl;
    std::cout << "Soil ---------------------------------------------" << std::endl;
    std::cout << "             Ground water flow: "; //i_gWaterEnable
    // turns on/off groundwater flow. Values: on: y; off: n
    config_data.getColumnValue(ground, "i_gWaterEnable", config_setting);
    if (ground){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    std::cout << "     Soil water redistribution: "; //i_soilRedEnable
    // turns on/off soil redistribution routine. Values: on: y; off: n
    config_data.getColumnValue(soilred, "i_soilRedEnable", config_setting);
    if(soilred){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    std::cout << "        Soil water evaporation: "; //i_soilEvapEnable
    // turns on/off soil evaporation routine. Values: on: y; off: n
    config_data.getColumnValue(sevap, "i_soilEvapEnable", config_setting);
    if(sevap){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    // // Climate
    std::cout << std::endl;
    std::cout << "Climate ------------------------------------------" << std::endl;
    std::cout << "               Rainfall inputs: "; //i_rainEnable
    config_data.getColumnValue(raining, "i_rainEnable", config_setting);
    if (!raining) {
        rainEnabled = false;
        std::cout << "Off" << std::endl;
    }
    else {
        //enabled is made the "else" case becasue the "default" state is to process rain, so any input OTHER than //n// enables
        rainEnabled = true;
        std::cout << "On" << std::endl;
    } // endif

    std::cout << "               Multiple years: "; 
    // turns on/off if we are working with multiple growing seasons. Values: n(off); y(on)
    if (stage_id == STAGE_ID_NONE || stage_id == STAGE_ID_HIST_STRESS || stage_id == STAGE_ID_FUT_STRESS || stage_id == STAGE_ID_FUT_STRESS_NOACCLIM)
    {
        if (config_data.getColumnValue("i_useGSDataStress") == "y") {
            useGSData = true;
            std::cout << "On" << std::endl;
        }
        else {
            useGSData = false;
            std::cout << "Off" << std::endl;
        } // endif
    }
    else if (stage_id == STAGE_ID_HIST_OPT || stage_id == STAGE_ID_FUT_OPT)
    {
        if (config_data.getColumnValue("i_useGSDataOpt", config_setting) == "y") { // different parameter name
            useGSData = true;
            std::cout << "On" << std::endl;
        }
        else {
            useGSData = false;
            std::cout << "Off" << std::endl;
        } // endif
    }
    else { // impossible unknown stage failsafe
        useGSData = false;
        std::cout << "Off" << std::endl;
    }

    // Hydraulics
    std::cout << std::endl;
    std::cout << "Hydraulics ---------------------------------------" << std::endl;
    std::cout << "                Pre-dawns mode: ";
    // turns on/off if model should consider measured pre-dawn water potential values. Values: n (off); y (on)
    if (config_data.getColumnValue("i_predawnsMode", config_setting) == "y")
    {
        mode_predawns = true;
        std::cout << "On" << std::endl;
        std::cout << "MODE: Running from predawn inputs (in rain column, MPa). Soil sim disabled." << std::endl;
    }
    else
    {
        mode_predawns = false;// default
        std::cout << "Off" << std::endl;
        std::cout << "MODE: Running from soil simulation (calculated predawns)." << std::endl;
    }
    
    std::cout << "          Xylem refilling mode: ";
    config_data.getColumnValue(refilling, "i_refilling", config_setting);
    if(refilling){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }
    
    std::cout << "              Xylem hysteresis: "; 
    // turns on/off xylem hysteresis from previous growing season. Values: n(off); y(on) 
    if(config_data.getColumnValue("i_cavitFatigue", config_setting) == "y"){// "i_cavitFatigue"
        hysteresis = true;
        std::cout << "On" << std::endl;
    } else {
        hysteresis = false;// default
        std::cout << "Off" << std::endl;
    }

    std::cout << "   Cavitation fatigue in roots: ";
    if (config_data.getColumnValue("i_stemOnly", config_setting) == "n")
    {
        stem_only = false;
        std::cout << "On" << std::endl;
    } else {
        stem_only = true; // Default
        std::cout << "Off" << std::endl;
    }

    // BA:GA optimization routine
    std::cout << std::endl;
    std::cout << "BAGA Optimization --------------------------------" << std::endl;
    std::cout << "          Iterate Ground Water: ";
    // Iterating ground water in for each stand. Values: y; n
    if (config_data.getColumnValue("i_iter_gwEnable", config_setting) == "y")
        iter_gwEnable = true;
    else
        iter_gwEnable = false;
    if (iter_gwEnable)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl;//default
    }
    
    param_data.getColumnValue(iter_gwInc, "i_iter_gwInc", species_no);
    param_data.getColumnValue(iter_gwStart, "i_iter_gwStart", species_no);
    param_data.getColumnValue(iter_gwEnd, "i_iter_gwEnd", species_no);

    std::cout << "  Iterate Field Capacity (FFC): ";
    // Iterating field capacity for each stand. Values: y; n
    if (config_data.getColumnValue("i_iter_ffcEnable") == "y")
        iter_ffcEnable = true;
    else
        iter_ffcEnable = false;
    if (iter_ffcEnable)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl;//default
    }

    param_data.getColumnValue(iter_ffcInc, "i_iter_ffcInc",species_no);
    param_data.getColumnValue(iter_ffcStart, "i_iter_ffcStart",species_no);
    param_data.getColumnValue(iter_ffcEnd, "i_iter_ffcEnd",species_no);
    
    std::cout << "           Iterate Field BA:GA: ";
    // Iterating BA:GA for each stand. Values: off; on
    if (config_data.getColumnValue("i_iter_bagaEnable", config_setting) == "y")
        iter_bagaEnable = true;
    else
        iter_bagaEnable = false;
    if (iter_bagaEnable)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl;//default
    }

    param_data.getColumnValue(iter_bagaInc, "i_iter_bagaInc",species_no);
    param_data.getColumnValue(iter_bagaStart, "i_iter_bagaStart",species_no);
    param_data.getColumnValue(iter_bagaEnd, "i_iter_bagaEnd",species_no);
    param_data.getColumnValue(iter_bagaRef, "i_iter_bagaRef",species_no);
    param_data.getColumnValue(iter_bagaCutoff, "i_iter_bagaCutoff",species_no);
    
    iter_bagaRef = 1.0; //156.269; // 1.0; // TODO TEMP can remove after param sheets updated
    iter_bagaEnd = 500.0; // TODO TEMP " -> for bisection method, allow extreme range

    std::cout << "                Use Area Table: ";
    // Are we using BA:GA values from another csv file. Values: y; n
    if (config_data.getColumnValue("i_iter_useAreaTable") == "y")
        iter_useAreaTable = true;
    else
        iter_useAreaTable = false;
    if (iter_useAreaTable)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl;// default
    }
    
    std::cout << "       Iterate Years as Counts: ";
    // Iterating yearly ground water for each stand. Values: off; on
    if (config_data.getColumnValue("i_iter_yearsAsCount") == "y")
        iter_yearsAsCount = true;
    else
        iter_yearsAsCount = false;
    if (iter_yearsAsCount)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl; // default
    }

    std::cout << "        Iteration Supply Curve: ";
    // Turns on the iterations in the BAGA optimization routine. Values: y; n
    if (config_data.getColumnValue("i_iter_runSupplyCurve") == "y")
        iter_runSupplyCurve = true;
    else
        iter_runSupplyCurve = false;
    if (iter_runSupplyCurve)
    {
        std::cout << "On" << std::endl;
    } else
    {
        std::cout << "Off" << std::endl;// default
    }

    // File locations
    std::cout << std::endl;
    std::cout << "File locations -----------------------------------" << std::endl;
    std::cout << "  Path to climate forcing data: " << std::endl;
    climate_forcing_data_path = config_data.getColumnValue("i_ClimateData", config_setting);
    std::cout << climate_forcing_data_path << std::endl;
    
    std::cout << "   Path to growing season data: " << std::endl;
    growing_season_limits_data_path = config_data.getColumnValue("i_GSData", config_setting);
    std::cout << growing_season_limits_data_path << std::endl;

    std::cout << "  Path to time-step header file: " << std::endl;
    data_header_file_path = DATA_HEADER_FILE_PATH;
    std::cout << data_header_file_path << std::endl;
}

/**
 * @brief Initializes the model variables for the Plant class.
 * 
 * This method is responsible for resetting and initializing all the variables
 * used in the Plant model. It is typically called at the start of every iteration,
 * but it should not be called on new years. The method ensures that all variables
 * are set to their default or initial values, preparing the model for a new simulation
 * or iteration cycle.
 * 
 * Key functionalities:
 * - Resets general state variables such as `isNewYear`, `gs_yearIndex`, and `gs_inGrowSeason`.
 * - Initializes water-related variables like `runoff`, `drainage`, and `transpiration`.
 * - Resets arrays such as `b_fatigue` and `e_p` to zero.
 * - Sets various flags and parameters to their default states.
 * - Prepares iteration-specific variables for groundwater, ffc, and baga simulations.
 * 
 * Note:
 * - This method should not be called on new years to avoid overwriting year-specific data.
 */
void Plant::initModelVars() {
    // this can be called on the start of every iteration BUT NOT ON NEW YEARS
    isNewYear = true;
    // this is a good place for general initialization too
    gs_yearIndex = 0;
    gs_prevDay = 0;
    gs_inGrowSeason = false;
    gs_doneFirstDay = false;
    runoff = 0;
    drainage = 0;
    year_cur = 0;
    year_start = 0;
    yearVal = 0;

    transpiration = 0.0;
    transpirationsh = 0.0;
    md = 0.0;
    mdsh = 0.0;
    rmean = 0.0;
    ecritsystem = 0.0;
    pcritsystem = 0.0;
    kmin = 0.0;
    iter_refK = 0.0;
    kpday1 = 0.0;
    kxday1 = 0.0;
    transpiration_tree = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 10; ++j) {
            b_fatigue[i][j] = 0.0;
        }
    }
    for (int i = 0; i < CURVE_MAX; ++i) {
        e_p[i] = 0.0;
    }
    stage_id = 0;
    unknowns = 0;
    halt = 0;
    haltsh = 0;
    layers = 0;
    tod = 0;
    mode_predawns = false;
    hysteresis = false;
    stem_only = false;
    iter_gwEnable = false;
    iter_ffcEnable = false;
    iter_bagaEnable = false;
    iter_useAreaTable = false;
    iter_yearsAsCount = false;
    iter_runSupplyCurve = false;
    iter_gwInc = 0.0;
    iter_gwStart = 0.0;
    iter_gwEnd = 0.0;
    iter_ffcInc = 0.0;
    iter_ffcStart = 0.0;
    iter_ffcEnd = 0.0;
    iter_bagaInc = 0.0;
    iter_bagaStart = 0.0;
    iter_bagaEnd = 0.0;
    iter_bagaRef = 0.0;
    iter_bagaCutoff = 0.0;
    max_plc_x = 0.0;
    gwflow = 0.0;
    drainage = 0.0;
    soilevap = 0.0;
}

/**
 * @brief Resets all member variables of the Plant class to their default values.
 * 
 * This function is responsible for initializing or resetting the state of the Plant object.
 * It sets numerical variables to zero, boolean variables to false, clears containers, 
 * and resets parameters in associated objects. Additionally, it ensures that certain 
 * conditions are met before clearing specific data structures.
 * 
 * Key operations:
 * - Resets numerical variables such as species_no, stage_id, and various iteration parameters.
 * - Sets boolean flags like mode_predawns, hysteresis, and others to false.
 * - Clears containers such as water, fc, and jmatrix.
 * - Resets parameters in associated objects like carbon and xylem.
 * - Initializes arrays like e_p to zero.
 * 
 * Note:
 * - The function includes a conditional check to avoid clearing water and fc containers 
 *   if useGSData is true and gs_yearIndex is greater than 0.
 */
void Plant::cleanModelVars() {
    species_no = 0;
    stage_id = 0;
    unknowns = 0;
    halt = 0;
    haltsh = 0;
    layers = 0;
    tod = 0;

    mode_predawns = false;
    hysteresis = false;
    stem_only = false;
    iter_gwEnable = false;
    iter_ffcEnable = false;
    iter_bagaEnable = false;
    iter_useAreaTable = false;
    iter_yearsAsCount = false;
    iter_runSupplyCurve = false;

    iter_gwInc = 0.0;
    iter_gwStart = 0.0;
    iter_gwEnd = 0.0;
    iter_ffcInc = 0.0;
    iter_ffcStart = 0.0;
    iter_ffcEnd = 0.0;
    iter_bagaInc = 0.0;
    iter_bagaStart = 0.0;
    iter_bagaEnd = 0.0;
    iter_bagaRef = 0.0;
    iter_bagaCutoff = 0.0;
    max_plc_x = 0.0;
    gwflow = 0.0;
    drainage = 0.0;
    soilevap = 0.0;
    iter_refK = 0.0;
    transpiration_tree = 0.0;
    transpiration = 0.0;
    transpirationsh = 0.0;
    md = 0.0;
    mdsh = 0.0;
    ecritsystem = 0.0;
    pcritsystem = 0.0;

    for (int i = 0; i < CURVE_MAX; ++i) {
        e_p[i] = 0.0;
    }

    if (!(useGSData && gs_yearIndex > 0)) {
        water.clear();
        fc.clear();
    }
    
    // reset jacobian matrix
    jmatrix.clear();

    carbon.clearParameters();
    xylem.cleanParameters();

}

/**
 * @brief Sets Van Genuchten parameters for soil layers based on soil texture.
 *
 * This function assigns Van Genuchten parameters (alpha, n, saturated hydraulic 
 * conductivity, and saturated water content) to a collection of soil layers 
 * based on the specified soil texture. If an invalid texture is provided, 
 * the function will terminate the program with an error message.
 *
 * @param texture A string representing the soil texture category. 
 *                Valid options include:
 *                "sand", "loamy sand", "sandy loam", "loam", "silt", 
 *                "silt loam", "sandy clay loam", "clay loam", 
 *                "silty clay loam", "sandy clay", "silty clay", "clay".
 * @param layers The number of soil layers to process.
 * @param soils A vector of pointers to SoilLayer objects, representing the soil layers.
 *
 * @throws std::runtime_error If an invalid soil texture is provided, the program 
 *         will terminate with an error message.
 *
 * @note The function modifies the rhizosphere properties of each soil layer 
 *       in the provided vector.
 * @note Alpha is MPa-1 instead of normal cm-1, and n is unitless
 * @note Thetasat is unitless as well, representing volume/volume
 *       Can be estimated via porosity at saturation (volume of water at saturation/total volume of soil)
 */
void get_vgparams(std::string texture, const double &layers, std::vector<SoilLayer*> &soils) {
    double a, n, soilkmax, thetasat;
    if (texture == "sand"){
        a = 1479.5945; 
        n = 2.68;
        soilkmax = 30305.88; 
        thetasat = 0.43;
    } else if (texture == "loamy sand"){
        a = 1265.3084;
        n = 2.28;
        soilkmax = 14897.84;
        thetasat = 0.41;
    } else if (texture == "sandy loam") {
        a = 765.3075;
        n = 1.89;
        soilkmax = 4510.168;
        thetasat = 0.41;
    } else if (texture == "loam") {
        a = 367.3476;
        n = 1.56;
        soilkmax = 1061.216;
        thetasat = 0.43;
    } else if (texture == "silt") {
        a = 163.2656;
        n = 1.37;
        soilkmax = 255.1;
        thetasat = 0.46;
    } else if (texture == "silt loam") {
        a = 204.082;
        n = 1.41;
        soilkmax = 459.18;
        thetasat = 0.45;
    } else if (texture == "sandy clay loam") {
        a = 602.0419;
        n = 1.48;
        soilkmax = 1336.724;
        thetasat = 0.39;
    } else if (texture == "clay loam") {
        a = 193.8779;
        n = 1.31;
        soilkmax = 265.304;
        thetasat = 0.41;
    } else if (texture == "silty clay loam") {
        a = 102.041;
        n = 1.23;
        soilkmax = 71.428;
        thetasat = 0.43;
    } else if (texture == "sandy clay") {
        a = 275.5107;
        n = 1.23;
        soilkmax = 122.448;
        thetasat = 0.38;
    } else if (texture == "silty clay") {
        a = 51.0205;
        n = 1.09;
        soilkmax = 20.408;
        thetasat = 0.36;
    } else if (texture == "clay") {
        a = 81.6328;
        n = 1.09;
        soilkmax = 204.08;
        thetasat = 0.38;
    } else if (texture == "pofr-dbg") { // using estimate of all populations (CCR/JLA/NRV/TSZ)
        a = 775.4012357;
        n = 1.471392609;
        soilkmax = 1336.724; // taken from sandy clay loam, since that is closest soil profile
        thetasat = 0.47; // from soil tin estimates, 0.47 +/- 0.03
    } else {
        std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
        std::cout << "SOURCE: Incorrect soil texture category" << std::endl;
        std::cout << "ACTION: Model stops " << std::endl;
        std::cout << std::endl;
        abort();
    }
    for (int k = 0; k <= layers; k++){
        soils[k]->rhizosphere.setVanGenAlpha(a);
        soils[k]->rhizosphere.setVanGenN(n);
        soils[k]->kkmax = soilkmax;
        soils[k]->rhizosphere.setThetasat(thetasat);
    }
}

/**
 * @brief Initializes and calculates all parameters for the Plant model at the 
 * start of the simulation.
 * 
 * This function reads input parameters from a data source, performs necessary 
 * calculations, and sets up the initial conditions for the Plant model. It 
 * handles parameters related to site, atmosphere, soil, stand, hydraulics, and 
 * photosynthesis. Additionally, it calculates initial conditions and adjusts 
 * parameters for soil and plant hydraulics.
 * 
 * @details
 * - **Site Parameters**: Latitude, longitude, slope inclination, slope aspect, 
 *   elevation.
 * - **Atmosphere Parameters**: Atmospheric transmittance, solar noon correction, 
 *   ambient CO2, long wave emissivity.
 * - **Soil Parameters**: Number of soil layers, rock fraction, field capacity, 
 *   ground water properties, soil absorptivity.
 * - **Stand Parameters**: Leaf area index, tree height, root depth, root biomass 
 *   distribution.
 * - **Hydraulics**: Weibull parameters for roots, stems, and leaves, accounting 
 *   for cavitation fatigue.
 * - **Photosynthesis Parameters**: Light compensation point, quantum yield, 
 *   temperature dependencies, and other photosynthetic constants.
 * - **Initial Conditions**: Calculates conductance, resistance, and transport 
 *   distances for soil and plant components.
 * 
 * @note The function also adjusts parameters for cavitation fatigue in roots 
 * and stems based on previous years' drought stress.
 * 
 * @warning This function assumes that all required input parameters are 
 * available in the data source.
 * 
 * @param None
 * 
 * @return void
 */
void Plant::readin() { //'inputs and calculates all parameters at the start

    std::cout << std::endl;
    std::cout << "INIT: Setting up model parameters. (year count = " << gs_yearIndex << ")" << std::endl;
    std::cout << std::endl;

    double temp;
    // site identifiers & parameters
    param_data.getColumnValue(temp, "i_latitude", species_no); // latitude in degree fraction north
    temp = PI / 180 * temp; // converted degrees to radians
    param.setModelParam(temp, "lat");

    param_data.getColumnValue(temp, "i_longitude", species_no); // longitude in degree fraction west
    param.setModelParam(temp, "longitude");

    param_data.getColumnValue(temp, "i_slopeI", species_no); // slope inclination, degrees from horizontal
    param.setModelParam(temp, "slope");
    param_data.getColumnValue(temp, "i_slopeA", species_no); // slope aspect, counterclockwise degrees from south
    param.setModelParam(temp, "slope_asp");
    param_data.getColumnValue(temp, "i_elevation", species_no); // elevation in m
    param.setModelParam(temp, "alt");

    // atmosphere
    param_data.getColumnValue(temp, "i_atmTrans", species_no); // atmospheric transmittance from weather data
    param.setModelParam(temp, "tau");
    param_data.getColumnValue(temp, "i_solarNoon", species_no); // solar noon correction from weather data
    param.setModelParam(temp, "tsn_corr");
    param_data.getColumnValue(temp, "i_co2AmbPPM", species_no); // ambient co2 in ppm, first year only
    param.setModelParam(temp, "param_ca");
    temp = temp * 0.000001; // ambient co2 in moles per mole
    carbon.ca = temp;
    param_data.getColumnValue(temp, "i_emiss", species_no); // long wave emissivity
    param.setModelParam(temp, "emiss");
    
    // soil
    param_data.getColumnValue(temp, "i_layers", species_no); // number of soil layers
    if (mode_predawns)// predawns mode, only 1 layer as we don't need to calculate pre-dawns
    {
        temp = 1;
    }
    param.setModelParam(temp, "layers");
    this->layers = temp; // easier access

    water.resize(layers+1, 0);
    fc.resize(layers+1, 0);

    param_data.getColumnValue(temp, "i_rockFrac", species_no); // fraction of soil volume as rocks
    param.setModelParam(temp, "rock_frac");
    // TO-DO: update here to be able to input VG parameters from soil water retention curves
    param_data.getColumnValue(temp, "i_rhizoPer", species_no);
    temp /= 100.0; // average fraction of whole plant resistance in rhizosphere(maximum soil limitation)
    param.setModelParam(temp, "rhizo_targ");
    param_data.getColumnValue(temp, "i_fieldCapPercInit", species_no);
    temp /= 100.0; // fraction of field capacity for starting the season
    param.setModelParam(temp, "ffc");
    param_data.getColumnValue(temp, "i_fieldCapFrac", species_no); // fraction that field capacity is of saturation(minus residual)
    param.setModelParam(temp, "field_cap_frac");
    param_data.getColumnValue(temp, "i_gWaterP", species_no); // ground water pressure
    param.setModelParam(temp, "p_ground");
    param_data.getColumnValue(temp, "i_gWaterDist", species_no); // distance to ground water source in m
    param.setModelParam(temp, "ground_distance");
    param_data.getColumnValue(temp, "i_soilAbsSol", species_no); // absorptivity of soil surface for solar
    param.setModelParam(temp, "soil_abs_sol");
    
    // stand
    param_data.getColumnValue(temp, "i_treeToPhotoLAI", species_no); // tree LAI / Photo LAI
    param.setModelParam(temp, "tree_to_photo_lai");
    param_data.getColumnValue(temp, "i_leafAreaIndex", species_no); // canopy lai
    param.setModelParam(temp, "lai");
    param_data.getColumnValue(temp, "i_leafAngleParam", species_no); // leaf angle parameter, CN 15.4
    param.setModelParam(temp, "leaf_angle_param");
    param_data.getColumnValue(temp, "i_leafPerBasal", species_no); // initial leaf area per basal area, m2 m - 2
    param.setModelParam(temp, "leaf_per_basal");
    param_data.getColumnValue(temp, "i_height", species_no); // average tree height in m
    param.setModelParam(temp, "height");
    param_data.getColumnValue(temp, "i_aspect", species_no); // max radius of root system per max depth
    param.setModelParam(temp, "aspect");
    param_data.getColumnValue(temp, "i_rootDepth", species_no); // maximum rooting depth
    param.setModelParam(temp, "root_depth");
    
    /* Root biomass distribution is allocated based on the equation reported in Love et al (2019):
    M = 1 - Beta^d, where M is the fraction of biomass above depth d expressed in cm. We find the
    Beta that provides an M of 0.995 for the maximum rooting depth. */
    //rootdepth = 1 / (rootdepth*100);
    //beta = pow(0.005,rootdepth); // Calculates Beta as: Beta = (1-M)^(1/rootdepth); where M = 0.995 and rootdepth is cm
    param_data.getColumnValue(temp, "i_rootBeta", species_no); // root beta for Y = 1 - B ^ d
    param.setModelParam(temp, "beta");
    param_data.getColumnValue(temp, "i_leafWidth", species_no);
    temp *= 0.72; // leaf width x factor = characteristic dimension(campbell and norman)
    param.setModelParam(temp, "leaf_width");
    param_data.getColumnValue(temp, "i_baperga", species_no);
    temp *= 0.0001; // basal area per ground area converted from m2 / Ha to m2 / m2
    param.setModelParam(temp, "ba_per_ga");
    param_data.getColumnValue(temp, "i_soilXHeight", species_no);
    temp *= 100.0; // height above soil surface for understory wind and gh in cm
    param.setModelParam(temp, "xh");

    /* Hydraulics */
    
    // ROOTS
    double c_temp, b_temp;
    param_data.getColumnValue(c_temp, "i_cr", species_no);
    param_data.getColumnValue(b_temp, "i_br", species_no);
    std::cout << "Weibull parameter before accounting for previous year PLC: " << std::endl;       
    std::cout << "root b: " << b_temp << " root c: " << c_temp << std::endl;
    if (hysteresis == true && stem_only == false)
    {
        param_data.getColumnValue(temp, "i_sapwoodT",species_no);
        param.setModelParam(temp, "sapwood_t");
        param_data.getColumnValue(temp, "i_conduitDiam",species_no);
        param.setModelParam(temp, "conduit_d");

        if (gs_yearIndex == 1)
        {
            b_fatigue[0][0] = b_temp; // current year "fresh xylem"
            // accounting for drought stress in previous year (1 year old)
            b_fatigue[0][1] = xylem.fatigue(b_temp, 
                                            param.getModelParam("sapwood_t"), 
                                            param.getModelParam("conduit_d"),
                                            max_plc_x);
            b_temp = ((b_fatigue[0][0] * 1.0) + (b_fatigue[0][1] * 0.75)) / (1.0 + 0.75); // updated parameter b for roots
        } else if (gs_yearIndex == 2) 
        {
            b_fatigue[0][0] = b_temp; // current year "fresh xylem"
            b_fatigue[0][2] = b_fatigue[0][1]; // previous year ring now 1 year older (2 years old)
            b_fatigue[0][1] = (xylem.fatigue(b_temp, 
                                             param.getModelParam("sapwood_t"), 
                                             param.getModelParam("conduit_d"),
                                             max_plc_x)); // accounting for drought stress in previous year (1 year old)
            b_temp = ((b_fatigue[0][0] * 1.0) + (b_fatigue[0][1] * 0.75) + (b_fatigue[0][2] * 0.50)) / (1.0 + 0.75 + 0.50);// updated parameter b for roots
        } else if (gs_yearIndex >= 3)
        {   
            b_fatigue[0][0] = b_temp;// current year "fresh xylem"
            b_fatigue[0][3] = b_fatigue[0][2]; // 2 years old ring now 3 years old
            b_fatigue[0][2] = b_fatigue[0][1]; // previous year ring now 1 year older (2 years old)
            b_fatigue[0][1] = xylem.fatigue(b_temp, 
                                            param.getModelParam("sapwood_t"), 
                                            param.getModelParam("conduit_d"),
                                            max_plc_x); // accounting for drought stress in previous year (1 year old)
            b_temp = ((b_fatigue[0][0] * 1) + (b_fatigue[0][1] * 0.75) + (b_fatigue[0][2] * 0.5) + (b_fatigue[0][3] * 0.25))/(1.0 + 0.75 + 0.5 + 0.25);// updated parameter b for roots
        } else 
        {
            // we don't account for cavitation fatigue in the first year unless we now about it
        }
        std::cout << "Weibull parameter after accounting for previous year PLC: " << std::endl;
        std::cout << "root b: " << b_temp << " root c: " << c_temp << std::endl;
    }

    for (int k = 0; k <= layers; k++) {
        SoilLayer *t = new SoilLayer; // Xylem destructor deletes the allocated SoilLayers
        t->root.setCwb(c_temp);
        t->root.setBwb(b_temp);
        xylem.soils.push_back(t);

        for (int j = 0; j < sizeof(b_fatigue[0]) / sizeof(b_fatigue[0][0]); j++) {
            xylem.soils[k]->root.setFatigue(b_fatigue[0][j], j);
        }
    }
    get_vgparams(param_data.getColumnValue("i_texture", species_no), layers, xylem.soils); // extract VG parameters for a given texture

    // STEMS
    param_data.getColumnValue(c_temp, "i_cs",species_no);
    param_data.getColumnValue(b_temp, "i_bs",species_no);
    std::cout << "Weibull parameter before accounting for previous year PLC: " << std::endl;
    std::cout << "stem b: " << b_temp << " stem c: " << c_temp << std::endl;
    if (hysteresis == true && stem_only == false)
    {
        if(gs_yearIndex == 1)
        {
            b_fatigue[1][0] = b_temp; // current year "fresh xylem"
            b_fatigue[1][1] = xylem.fatigue(b_temp, 
                                                    param.getModelParam("sapwood_t"), 
                                                    param.getModelParam("conduit_d"),
                                                    max_plc_x);// accounting for drought stress in previous year (1 year old)
            b_temp = ((b_fatigue[1][0] * 1.0) + (b_fatigue[1][1] * 0.75)) / (1.0+0.75); // updated parameter b for stems
        } else if(gs_yearIndex == 2) 
        {
            b_fatigue[1][0] = b_temp; // current year "fresh xylem"
            b_fatigue[1][2] = b_fatigue[1][1]; // previous year ring now 1 year older (2 years old)
            b_fatigue[1][1] = xylem.fatigue(b_temp, 
                                                    param.getModelParam("sapwood_t"), 
                                                    param.getModelParam("conduit_d"),
                                                    max_plc_x); // accounting for drought stress in previous year (1 year old)
            b_temp = ((b_fatigue[1][0] * 1.0) + (b_fatigue[1][1] * 0.75) + (b_fatigue[1][2] * 0.50)) / (1.0+0.75+0.50);// updated parameter b for stems
        } else if (gs_yearIndex >= 3)
        {   
            b_fatigue[1][0] = b_temp;// current year "fresh xylem"
            b_fatigue[1][3] = b_fatigue[1][2]; // 2 years old ring now 3 years old
            b_fatigue[1][2] = b_fatigue[1][1]; // previous year ring now 1 year older (2 years old)
            b_fatigue[1][1] = xylem.fatigue(b_temp, 
                                                    param.getModelParam("sapwood_t"), 
                                                    param.getModelParam("conduit_d"),
                                                    max_plc_x); // accounting for drought stress in previous year (1 year old)
            b_temp = ((b_fatigue[1][0] * 1) + (b_fatigue[1][1] * 0.75) + (b_fatigue[1][2] * 0.5) + (b_fatigue[1][3] * 0.25))/(1.0+0.75+0.5+0.25);// updated parameter b for stems
        } else 
        {
            // we don't account for cavitation fatigue in the first year unless we now about it
        }
        std::cout << "Weibull parameter after accounting for previous year PLC: " << std::endl;
        std::cout << "stem b: " << b_temp << " stem c: " << c_temp << std::endl;
    }
    xylem.stem.setCwb(c_temp);
    xylem.stem.setBwb(b_temp);
    
    // LEAVES, no histeresis in leaves. Leaf xylem can be replaced//refilled more actively.
    //p12_l = param_data.getColumnValue(temp, "i_leafP12",species_no);
    //p50_l = param_data.getColumnValue(temp, "i_leafP50",species_no);
    //cl = get_cweibull(p12_l,p50_l);// weibull c for each root element
    //bl = get_bweibull(p12_l,cl);// weibull b for each root element
    param_data.getColumnValue(temp, "i_cl",species_no);
    xylem.leaf.setCwb(temp);
    c_temp = temp;
    param_data.getColumnValue(temp, "i_bl",species_no);
    xylem.leaf.setBwb(temp);
    b_temp = temp;
    std::cout << "leaf b: " << b_temp << " leaf c: " << c_temp << std::endl;
    param_data.getColumnValue(temp, "i_leafPercRes",species_no); // saturated % of tree R in leaves
    xylem.leaf.setResPercent(temp);
    param_data.getColumnValue(temp, "i_kmaxTree",species_no); // kmax of tree in kg hr - 1m - 2MPa - 1 per basal area
    param.setModelParam(temp, "ksatp");
    std::string ksatpcalc = "old";
    if (param.getModelParam("ksatp") <= 0)
    {
        param_data.getColumnValue(temp, "i_LSC", species_no);// LSC
        param.setModelParam(temp, "lsc_input");

        param_data.getColumnValue(temp, "i_LSCpref", species_no);// LSCPref
        param.setModelParam(temp, "lsc_pref");

        temp = param.getModelParam("lsc_input") / (std::exp(-(std::pow((param.getModelParam("lsc_pref") / b_temp), c_temp))));//--//
        param.setModelParam(temp, "lsc_max");
        ksatpcalc = "new";
    }
    
    param_data.getColumnValue(temp, "i_pinc",species_no); // MPa increment for global K(P) curves...has to be small enough for NR convergence.
    if (temp <= 0)
    {
        temp = 0.00075; // default, override if it's zero
    }
    param.setModelParam(temp, "p_inc");
    
    // photosynthesis
    param_data.getColumnValue(temp, "i_lightComp",species_no); // light compensation point in ppfd
    param.setModelParam(temp, "light_comp");
    param_data.getColumnValue(temp, "i_qMax",species_no); // quantum yield of electron transport, moles e per mols photons
    param.setModelParam(temp, "q_max");
    param_data.getColumnValue(temp, "i_vmax25",species_no); // vmax at 25C
    param.setModelParam(temp, "v_max25");
    param_data.getColumnValue(temp, "i_jmax25",species_no); // jmax at 25C
    param.setModelParam(temp, "j_max25");
    param_data.getColumnValue(temp, "i_kc25",species_no); // m - m constant for CO2 in mole fraction at 25C
    param.setModelParam(temp, "kc_25");
    param_data.getColumnValue(temp, "i_ko25",species_no); // m - m constant for O2 in mole fraction at 25C
    param.setModelParam(temp, "ko_25");
    param_data.getColumnValue(temp, "i_comp25",species_no); // photorespiratory compensation point in mole fraction at 25C
    param.setModelParam(temp, "comp_25");
    param_data.getColumnValue(temp, "i_thetaC",species_no); // shape factor for A - ci colimitation, 0.98
    param.setModelParam(temp, "theta_c");
    param_data.getColumnValue(temp, "i_havmax",species_no); // these are all temp - dependency parameters from Leunig 2002
    param.setModelParam(temp, "ha_vmax");
    param_data.getColumnValue(temp, "i_hdvmax",species_no); //
    param.setModelParam(temp, "hd_vmax");
    param_data.getColumnValue(temp, "i_svvmax",species_no); //
    param.setModelParam(temp, "sv_vmax");
    param_data.getColumnValue(temp, "i_hajmax",species_no); //
    param.setModelParam(temp, "ha_jmax");
    param_data.getColumnValue(temp, "i_hdjmax",species_no); //
    param.setModelParam(temp, "hd_jmax");
    param_data.getColumnValue(temp, "i_svjmax",species_no); //
    param.setModelParam(temp, "sv_jmax");
    param_data.getColumnValue(temp, "i_lightCurv",species_no); //
    param.setModelParam(temp, "light_curv");

    // initial conditions calculations ---------------------------------------------------------------------------------------------------------
    temp = 1000000;
    param.setModelParam(temp, "g_max");
    temp = temp * (1.0 / param.getModelParam("leaf_per_basal")) * (1 / 3600.0) * 55.56 * 1000.0; //convert to param.getModelParam("g_max") per leaf area in mmol m-2s-1
    param.setModelParam(temp, "g_max_l");
    temp = 101.325 * pow((1 - 0.0065 * param.getModelParam("alt") / (288.15 + 0.0065 * param.getModelParam("alt"))), 5.257); //atmospheric pressure, T = 15 C, average sealevel patm; approximation
    param.setModelParam(temp, "p_atm");
    if (ksatpcalc == "old")
    {   
        //[HNT] calculate and output the stem and root percent resistances @ ksat
        xylem.stem.setResPercent(1.0 / 3.0 * (100.0 - xylem.leaf.getResPercent())); // saturated % of tree R in stem
        double rootpercent = 2.0 / 3.0 * (100.0 - xylem.leaf.getResPercent()); // saturated % of tree R in roots

        temp = param.getModelParam("ksatp") * (100.0 / xylem.leaf.getResPercent());
        xylem.leaf.setKmax(temp); // leaf conductance per basal area

        temp = temp * 1.0 / param.getModelParam("leaf_per_basal"); // lsc in kg hr-1m-2MPa-1//lsc per leaf area in kg hr - 1
        param.setModelParam(temp, "lsc");

        temp = 1.0 / param.getModelParam("ksatp"); // convert to resistance
        param.setModelParam(temp, "rsatp");
        temp = 1.0 / ((rootpercent / 100.0) * temp); // kmax of root system; assumes zero % rhizosphere resistance in WET soil
        param.setModelParam(temp, "k_sat_root");
        for (int k = 1; k <= layers; k++) {
            xylem.soils[k]->root.setResPercent(rootpercent);
        }

        temp = 1.0 / ((xylem.stem.getResPercent() / 100.0) * param.getModelParam("rsatp")); // kmax of stem system
        xylem.stem.setKmax(temp);
        std::cout << "Conductance calculations (old mode): kMaxTree = " << param.getModelParam("ksatp") << " lsc = " << param.getModelParam("lsc") << " kMaxStem = " << xylem.stem.getKmax() << " kMaxRoot = " << param.getModelParam("k_sat_root") << std::endl;
    } else {
        temp = param.getModelParam("lsc_max") * 3600.0 * 0.000018;
        param.setModelParam(temp, "lsc_max");
        temp = temp * param.getModelParam("leaf_per_basal");
        xylem.leaf.setKmax(temp);

        double rSatLeaf = 1.0 / temp; // 25%
        double rSatStem = 1.0 * rSatLeaf; // 25%
        double rSatRoot = 2.0 * rSatLeaf; // 50%

        double rSatTree = rSatLeaf + rSatStem + rSatRoot;
        temp = 1.0 / rSatTree;
        param.setModelParam(temp, "ksatp");
        temp = rSatTree;
        param.setModelParam(temp, "rsatp");
        param.setModelParam(param.getModelParam("lsc_max"), "lsc");
        temp = 1.0 / rSatRoot;
        param.setModelParam(temp, "k_sat_root");
        xylem.stem.setKmax(1.0 / rSatStem);

        // any chance this worked? (later me: It did, once I remembered the unit conversion for lscMax)
        std::cout << "Conductance calculations: kMaxTree = " << param.getModelParam("ksatp") << " lsc = " << param.getModelParam("lsc") << " kMaxStem = " << xylem.stem.getKmax() << " kMaxRoot = " << param.getModelParam("k_sat_root") << std::endl;
    }
    temp = param.getModelParam("height") * 0.01; //pressure drop from gravity in MPa
    param.setModelParam(temp, "p_grav");
    temp = param.getModelParam("ksatp") / 500.0; // e increment in kg hr-1 m-2 basal area for composite curve
    param.setModelParam(temp, "e_inc");
    temp = param.getModelParam("ksatp") / 2000.0; //"instantaneous K" DPA_MAX_CUTOFF for global K(P) curves for each element
    param.setModelParam(temp, "k_min");

    xylem.num_layers = layers;
    //for this soil data we want to use the original anchor-offset system
    for (int k = 1; k <= layers; k++) { // set layer depths and % root ksat
        xylem.soils[k]->layer_depth = 0.01 * log(1.0 - k * 0.995 / layers) / log(param.getModelParam("beta")); // lower depth of each layer converted to m
    }
    double depthmax = xylem.soils[layers]->layer_depth;
    // calculate transport distances
    std::vector<double> temp_vert_distance((layers * 2) + 1, 0);
    // first get vertical distance to biomass center of each layer
    for (int k = 1; k <= layers * 2.0; k++)
    {
        temp_vert_distance[k] = 0.01 * log(1 - k * 0.995 / (layers * 2.0)) / log(param.getModelParam("beta")); // get half depths
    }
    int i = 0;
    for (int k = 1; k <= layers * 2; k += 2) // To layers * 2 Step 2
    {
        i = i + 1;
        xylem.soils[i]->vert_distance = temp_vert_distance[k]; //take every other vertdistance
    }
    // now get radial distances
    for (int k = 1; k <= layers; k++) //To layers //get thicknesses
    {
        if (k == 1)
        {
            xylem.soils[k]->depth = xylem.soils[k]->layer_depth; //Cells(8 + k, 11)
        }
        else {
            xylem.soils[k]->depth = xylem.soils[k]->layer_depth - xylem.soils[k - 1]->layer_depth; // depth is layer thickness in meters
        } // endif
    }
    temp = xylem.soils[1]->depth * PI * pow((depthmax * param.getModelParam("aspect")), 2.0); //volume of first, and hence all, layers
    param.setModelParam(temp, "vol");
    double shallow = 0;
    // get radial widths of each layer and transport length
    for (int k = 1; k <= layers; k++) {
        xylem.soils[k]->radius = pow(temp / (xylem.soils[k]->depth * PI), 0.5); // width in m
        xylem.soils[k]->length = xylem.soils[k]->radius + xylem.soils[k]->vert_distance; // transport distance
        if (k == 1)
            shallow = xylem.soils[k]->length;
        xylem.soils[k]->length = xylem.soils[k]->length / shallow;    
    }

    unknowns = layers + 1; // number of unknowns to be solved for; also dimensions of matrix
    // Instantiate jmatrix with zeros
    for (int k = 0; k <= unknowns; k++) {
        std::vector<double> temp;
        for (int j = 0; j <= unknowns; j++)
            temp.push_back(0);
        jmatrix.push_back(temp);
    }

    temp = 1 - param.getModelParam("rock_frac"); // fraction of volume with no rocks
    for (int k = 1; k <= layers; k++) { // read in soil properties
        xylem.soils[k]->rhizosphere.setThetasat(xylem.soils[k]->rhizosphere.getThetaSat() * temp); // reduce for actual rock-free fraction of soil
    }
    param.setModelParam(temp, "rock_frac");
    /* Now add toplayer (layer 0) of rootless soil 2 cm thick w. same properties as layer 1 */
    xylem.soils[0]->depth = 0.02; //sets top layer to 2 cm
    xylem.soils[0]->root.setKmax(0); //set to 0
    xylem.soils[0]->rhizosphere.setKmax(0);
    xylem.soils[0]->kkmax = xylem.soils[1]->kkmax;
    xylem.soils[0]->rhizosphere.setThetasat(xylem.soils[1]->rhizosphere.getThetaSat());

    /* Now solve for kmax rhizosphere that gives the desired ave % rhizosphere resistance */
    int z = 1; //use layer 1 as stand in for whole root system
    xylem.soils[1]->root.setKmax(param.getModelParam("k_sat_root")); //set to whole root system
    double x = 0.5; //start by finding kmaxrh at 0.5 MPa...a deliberate under-shoot
    double rootr, rstem, rleaf, rplant, rhizor, vp, vgterm, kinc, sum, rrhizofrac;
    rootr = 1.0 / xylem.soils[1]->root.wb(x);
    rstem = 1.0 / xylem.stem.wb(x);
    rleaf = 1.0 / xylem.leaf.wb(x);
    rplant = rootr + rstem + rleaf; //rplant here is just the xylem part
    rhizor = rplant * (param.getModelParam("rhizo_targ") / (1.0 - param.getModelParam("rhizo_targ"))); //solve for what rhizor has to be at the target
    vp = 1.0 / (pow((xylem.soils[z]->rhizosphere.getVanGenAlpha() * x), xylem.soils[z]->rhizosphere.getVanGenN()) + 1); //van genuchten terms // vp = 1 / ((a(z) * x) ^ n(z) + 1) 
    vgterm = pow(vp, ((xylem.soils[z]->rhizosphere.getVanGenN() - 1) / (2.0 * xylem.soils[z]->rhizosphere.getVanGenN()))) * pow((pow((1 - vp), ((xylem.soils[1]->rhizosphere.getVanGenN() - 1) / xylem.soils[z]->rhizosphere.getVanGenN())) - 1), 2.0); //van genuchten terms // vgterm = vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
    xylem.soils[1]->rhizosphere.setKmax((1.0 / rhizor) / vgterm); //solve for kmaxrh[1]
    kinc = xylem.soils[1]->rhizosphere.getKmax() * 0.1;
    do //loop through rhizosphere kmax
    {
        xylem.soils[1]->rhizosphere.setKmax(xylem.soils[1]->rhizosphere.getKmax() + kinc); //increase from deliberate undershoot
        x = 0;
        sum = 0;
        do //loop through pressures
        {
            x = x + 0.1;
            rootr = 1.0 / xylem.soils[1]->root.wb(x);
            rstem = 1.0 / xylem.stem.wb(x);
            rleaf = 1.0 / xylem.leaf.wb(x);
            rhizor = 1.0 / xylem.soils[z]->rhizosphere.vg(x);
            rplant = rootr + rstem + rleaf + rhizor;
            rrhizofrac = rhizor / rplant; //fraction of resistance in rhizosphere
            sum = sum + rrhizofrac; //add up fractions
        } while (!((1.0 / rplant) < param.getModelParam("k_min"))); //Loop Until 1 / rplant < kmin //average over full range
        sum = sum / (x / 0.1); //average fraction
    } while (!(sum < param.getModelParam("rhizo_targ"))); // Until sum < rhizotarg //loop until desired soil limitation is reached
    
    xylem.soils[1]->rhizosphere.setKmax(xylem.soils[1]->rhizosphere.getKmax() / layers);  //divide whole root rhizokmax into equal portions for each layer
    // end of soil limitation adjustment
    // now set soil layer parameters based on aroot of entire root system
    for (int k = 1; k <= layers; k++)
    {
        xylem.soils[k]->rhizosphere.setKmax(xylem.soils[1]->rhizosphere.getKmax()); // soil to root MAXIMUM conductance in kg hr-1 MPa-1//re - set for individual layers
                                                                                    // note: this is kMAX...at P=0; not at saturated PD
    }

    // t = 0;
    //loop to find ksatroot for each layer
    double coef = 0.0;
    do
    {
        coef = coef + 0.01;
        sum = 0.0;
        for (int k = 1; k <= layers; k++)//k = 0 To layers //soil layers from top to bottom
        {
            xylem.soils[k]->root.setKmax(coef / xylem.soils[k]->length); //assumes ksatr proportional to biomass/length
            sum = sum + xylem.soils[k]->root.getKmax();
        }
    } while (!(sum > param.getModelParam("k_sat_root"))); //Loop Until sum > ksatroot //loop until each layer adds to total
    // for (k = 1; k <= layers; k++) {
    //         soillayersTable[rowLR + k][colLR + 2] = std::to_string(100.0 * (ksatr[k] / ksatroot)); // % root saturated k
    //         soillayersTable[rowLR + k][colLR + 9] = std::to_string(ksatr[k]); // root kmax
    // }

    xylem.rough = 0.01; //soil Zm, eqn 14.9, using table 5.1 for smooth surface, cm
    xylem.zdispl = 6.5 * xylem.rough; // soil d, eqn 14.9, using d = 6.5 Zm, eq 5.2,5.3
    xylem.zh = 0.2 * xylem.rough; // roughness for temperature
}

/**
 * @brief Resets the status of all soil layers in the xylem.
 * 
 * This function iterates through all the soil layers in the xylem and resets
 * their status by setting the `cavitated` property to `false` and the `failure`
 * property to `"ok"`.
 * 
 * @note Assumes that the `xylem.soils` array is properly initialized and has
 *       at least `layers + 1` elements.
 */
void Plant::resetLayerStatus() {
    for (int i = 0; i <= layers; i++) {
        xylem.soils[i]->cavitated = false;
        xylem.soils[i]->failure = "ok";
    }
}

/**
 * @brief Calculates the flow and pcrit for the plant's xylem component.
 *
 * This function invokes the `calc_flow` method of the `xylem` object to compute
 * flow characteristics. It retrieves the required.
 *
 * @details
 * - "p_inc" represents the pressure increment parameter.
 * - "k_min" represents the minimum conductivity parameter.
 *
 * The calculated flow parameters are used to model the behavior of the plant's
 * xylem under varying conditions.
 */
void Plant::componentPCrits() {
    std::cout << "Calculating critical points for components" << std::endl;
    xylem.calc_flow(param.getModelParam("p_inc"), param.getModelParam("k_min"));
}

/**
 * @brief
 */
/**
 * @brief Simulates a single timestep iteration for the plant model.
 *
 * This function performs a comprehensive simulation of plant processes for a given timestep.
 * It handles initialization for new years, updates environmental parameters, calculates 
 * photosynthesis, transpiration, and hydraulic properties, and manages soil and canopy processes.
 *
 * @param dd Reference to the current timestep index.
 * @return int Returns 0 on successful completion, a positive integer for a year index reset, 
 *         or -1 in case of failure.
 *
 * @details
 * - Initializes model parameters and variables at the start of a new year or the first timestep.
 * - Updates atmospheric CO2 concentration based on the year and input data.
 * - Retrieves environmental parameters such as solar radiation, vapor pressure deficit (VPD),
 *   air temperature, wind speed, and soil temperature.
 * - Simulates soil wetness and photosynthesis processes for sun and shade layers.
 * - Handles hydraulic processes, including soil pressure updates, canopy pressure, and water flow.
 * - Manages soil evaporation, groundwater flow, and redistribution of soil water.
 * - Outputs various plant and environmental metrics, including transpiration, photosynthesis,
 *   hydraulic conductance, and failure status for each soil layer.
 * - Handles convergence failures and retries calculations with adjusted parameters.
 * - Resets layer failure statuses and updates critical soil water content at the end of the timestep.
 *
 * @note This function is computationally intensive and involves multiple nested loops and conditionals.
 *       It also includes mechanisms to handle failures and ensure convergence of calculations.
 */
int Plant::modelTimestepIter(int &dd) {
    bool failure = false;
    int check = 0,
        reset_guess = 0,
        total = 0,
        test_count = 0;
    std::string failspot;
    int    timestep;
    double obssolar,
           vpd,
           airtemp,
           maxvpd,
           wind,
           us,
           soiltemp,
           chalk,
           lightcomp,
           kplantold,
           patm = param.getModelParam("p_atm");

    double laperba = param.getModelParam("leaf_per_basal");
    xylem.leaf.emiss = param.getModelParam("emiss");

    if (dd == 1 || isNewYear)
    {
        failure = 0;
        failspot = "no failure";
        componentPCrits();//'gets pcrits for each component
        failspot = "no failure";

        for (int i = 1; i <= layers; i++)// k = 1 To layers ;//'exclude the top layer
        {
            xylem.soils[i]->root.setKmin(xylem.soils[i]->root.getKmax());
        }

        xylem.stem.setKmin(xylem.stem.getKmax());
        xylem.leaf.setKmin(xylem.leaf.getKmax());
        kmin = param.getModelParam("ksatp");

        gwflow = 0; //'inflow to bottom of root zone
        drainage = 0; //'drainage from bottom of root zone

        gs_prevDay = 0;
        gs_inGrowSeason = false;
        gs_doneFirstDay = false; // prevents PLC misses on first year
                                //[/HNT]
    }

    yearVal = std::lround(data.getColumnValue("year", dd));
    if (yearVal != year_cur)
    {
        // all of these year values get zeroed out during initModelVars, so we can assume they will be zero at the start of a new iteration and test against 0 meaning "not set"
        if (yearVal > year_cur && yearVal > year_start)
        {
            // if the start year is zero, this is the first year we've processed
            if (year_start == 0) // model runs starting in the year zero are not supported, I guess?
                year_start = yearVal;
            year_cur = yearVal;

            gs_yearIndex = year_cur - year_start;
            if (gs_yearIndex >= 100)
                gs_yearIndex = 99; //safety check

            param.setGsArYear(gs_yearIndex, year_cur); // correct the year listing in the "growing seasons" array 
                                                // TODO make that input data able to handle the new system?? Otherwise just eliminate
                                                // on VBA side it might be sufficient to pull the start year from the first line of data
                                                // though this still assumes that our data is in linear time order ...

            if (gs_yearIndex > 0)
            {
                // if we've set the year and it was anythign other than the first timestep of this iteration (in which case we set the year index from zero TO zero)
                // then return to VBA and let it handle the model reset
                // it will re-run this timestep afterwards
                gs_prevDay = 0;
                gs_inGrowSeason = false;
                gs_doneFirstDay = false;

                isNewYear = true; // this can be used to override the dd==/>1 cases -- it will be set to false upon successful completion of a timestep
                return gs_yearIndex;
            }
            // if we detected a year change, but the yearIndex comes out to zero, then this must be the first timestep
            // of the first year in the data -- so we can continue without returning and re-setting the model
        }
        else // going back in time means something went wrong -- may be able to handle out of order years later though
        {
            std::cout << "FAILURE: Went back in time!" << std::endl;
            return -1; // failure
        }
    }
    
    if ( dd == 1 || isNewYear){ // Get CO2 for current year
        if (useGSData)
            carbon.ca = getCarbonByYear(yearVal); // get current year atmospheric CO2
        else
            carbon.ca = param.getModelParam("param_ca"); // if no GS data use ca from param file
        std::cout << "Atmospheric CO2 concentration for " << yearVal << ": " << carbon.ca << std::endl;
        carbon.ca = carbon.ca * 0.000001;
        param.setStemBWb(gs_yearIndex, xylem.stem.getBwb());
        param.setRootBWb(gs_yearIndex, xylem.soils[1]->root.getBwb());
    }

    int jd = data.getColumnValue("julian-day", dd); //'julian day

    if (dd > 1 && !isNewYear) { //if// //'set timestep
        if (tod < data.getColumnValue("standard-time", dd)) { //if// //'same old day getting older
            timestep = data.getColumnValue("standard-time", dd) - tod;
        }
        else {
            timestep = (24 - tod) + data.getColumnValue("standard-time", dd); //'a new day has started
                                                                            //[HNT] multi-year support
                                                                            // following method of new year detection has been replaced with the method above
                                                                            // now that we are including the year as a data input
            gs_prevDay = jd;
            gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's a new day, so let's see if this is in the growing season or not
                                    //[/HNT]
        }// //'tod if
    }// //'dd>1 if
        //[HNT] multi-year support
    else
    {
        gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's the first data point, test if we're starting in growing season
    }

    gs_inGrowSeason = isInGrowSeasonSimple(jd); // just always call this!
                                                //[/HNT]
    lightcomp = param.getModelParam("light_comp");
    tod = data.getColumnValue("standard-time", dd); //'time of day, standard local time in hour fraction
    obssolar = data.getColumnValue("solar", dd); //'observed total solar, horizontal, wm-2
    vpd = data.getColumnValue("D-MD", dd); //'midday vpd in kPa
    vpd = vpd / param.getModelParam("p_atm"); //'vpd in mole fraction
    airtemp = data.getColumnValue("T-air", dd); //'in C
    maxvpd = (101.3 / param.getModelParam("p_atm")) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    wind = data.getColumnValue("wind", dd); //'wind speed
    if (wind < MIN_WIND_THRESH) { //if//
        data.setColumnValue(MIN_WIND_THRESH, dd, "wind"); //'set to minimum wind
        wind = MIN_WIND_THRESH;
    }//
    us = wind * 0.1; //'understory windspeed in m s-1
    soiltemp = data.getColumnValue("T-soil", dd); //'surface temp of soil
    if (vpd > maxvpd) { //if//
        vpd = maxvpd;
        data.setColumnValue(maxvpd * param.getModelParam("p_atm"), dd, "D-MD"); //'print out maximum vpd
    }//

    xylem.leaf.lai = param.getModelParam("lai");

    getsoilwetness(dd, timestep, xylem.leaf.lai, xylem.leaf.laish, xylem.leaf.laisl, laperba, carbon.atree, carbon.cinc, carbon.ca); //'after initializing, start updating water contents of soil layers
    solarcalc(dd, jd, obssolar, maxvpd, airtemp, vpd, xylem.leaf.lai, xylem.leaf.laisl, xylem.leaf.laish, carbon.qsl, carbon.qsh, xylem.leaf.ssun, xylem.leaf.sshade, xylem.leaf.sref, xylem.leaf.la, xylem.leaf.lg); //'get radiation for timestep
    if (carbon.qsl > lightcomp /*[HNT] multiyear*/ && gs_inGrowSeason /*[/HNT]*/) { //if//
        night = "n"; //'it//'s light enough for sun layer to do business
    }
    else {
        night = "y"; //'too dark
    }// //'end night if

    gwflow = 0; //'re-set inflow to bottom of root zone
    drainage = 0; //'re-set drainage from bottom of root zone
                    //'Call getpredawns //'update soil pressure of each layer
                    //'if failure = 1 { //if// Exit do
    chalk = 0; //'einc counter

twentyMarker:

    failure = getpredawns(dd); //'passed initializing...update soil pressure of each layer

    if (failure == 1) {
        return -1;
    }
    int test = 1; //'=0 if stem or leaf fails
    int p = -1; //'E(Pleaf) point counter
    double einc = param.getModelParam("e_inc");
    double e = -einc; //'total e: einc and e are still in kg hr-1 m-2 basal area, as are conductances
    carbon.psynmax = -100;
    carbon.psynmaxsh = -100;
    int skip = 0; //'this turns off psynthesis
    double sum = 0;

    do // This loop obtains and stores the entire composite curve by incrementing e from zero
    {
        p = p + 1;
        e = e + einc;

        check = newtonrhapson(dd, param.getModelParam("p_inc"), e); // solves root and rhizosphere pressures

        if (check > 500) {
            std::cout << "500 restarts for newtonrhapson" << std::endl;
            break;
        }

        // test for total failure
        sum = 0;
        for (int k = 1; k <= layers; k++)
        {
            sum = sum + xylem.soils[k]->cavitated;
        }
        if (sum == layers) {
            failspot = "below ground";
            break;
        }
        test = xylem.stem.calc_pressure(e, xylem.soils[1]->root.getPressure(), param.getModelParam("p_grav"), param.getModelParam("p_inc")); //'gets stem and leaf pressures
        test = xylem.leaf.calc_pressure(e, xylem.stem.getPressure(), 0, param.getModelParam("p_inc")); // apparently gravity doesn't affect leaves downstream pressure??
        if (test == 1) {
            break;
        }
        test = compositeCurve(e, p, total); //'stores the entire composite curve for every P_c = p
        xylem.leaf.temp(p, airtemp, xylem.e_p, vpd, wind, laperba, param.getModelParam("leaf_width"), param.getModelParam("p_atm")); //'gets sun layer leaf temperature from energy balance
        xylem.leaf.tempShade(p, airtemp, param.getModelParam("p_atm"), vpd); //'gets shade layer leaf temperature
        carbon.assimilation(p, 
                            param.getModelParam("g_max"), 
                            param.getModelParam("q_max"), 
                            param.getModelParam("comp_25"), 
                            param.getModelParam("theta_c"), 
                            param.getModelParam("v_max25"), 
                            param.getModelParam("j_max25"), 
                            param.getModelParam("kc_25"), 
                            param.getModelParam("ko_25"), 
                            param.getModelParam("sv_vmax"), 
                            param.getModelParam("sv_jmax"), 
                            param.getModelParam("ha_vmax"), 
                            param.getModelParam("hd_vmax"), 
                            param.getModelParam("hd_jmax"), 
                            param.getModelParam("ha_jmax"), 
                            param.getModelParam("light_curv"), 
                            night, 
                            xylem.leaf.eplantl, 
                            xylem.leaf.lavpd, 
                            xylem.leaf.leaftemp); //'gets sun layer photosynthesis
        carbon.assimilationShade(p, 
                                 param.getModelParam("g_max"), 
                                 param.getModelParam("q_max"), 
                                 param.getModelParam("comp_25"), 
                                 param.getModelParam("theta_c"), 
                                 param.getModelParam("v_max25"), 
                                 param.getModelParam("j_max25"), 
                                 param.getModelParam("kc_25"), 
                                 param.getModelParam("ko_25"), 
                                 param.getModelParam("sv_vmax"), 
                                 param.getModelParam("sv_jmax"), 
                                 param.getModelParam("ha_vmax"), 
                                 param.getModelParam("hd_vmax"), 
                                 param.getModelParam("hd_jmax"), 
                                 param.getModelParam("ha_jmax"), 
                                 param.getModelParam("light_curv"), 
                                 night, 
                                 xylem.leaf.eplantl, 
                                 xylem.leaf.lavpdsh, 
                                 xylem.leaf.leaftempsh); //'gets shade layer photosynthesis

    } while (!(sum == layers || test == 1 || night == "y" && (dd > 1 && !isNewYear) || check >= 500)); //'loop to complete failure unless it//'s night

    if (chalk > 0) { //if//
        reset_guess = 0; //'done our best
        failspot = "convergence";
        for (int z = 1; z <= layers; z++)//z = 0 To layers //'restore layers to functioning if they//'ve been turned off by convergence failure
        {
            if (xylem.soils[z]->root.getKmin() != 0)
                xylem.soils[z]->cavitated = false;
        } //end for z

        goto fortyMarker; //'got as much of the composite curve as is going to happen
    }//
    if (dd == 1 || isNewYear || night == "n") { //if//

        if (check >= 500) { //if// //'try once more
                            //'Stop
            chalk = chalk + 1;

            if (ecritsystem == 0)
            {
                einc = param.getModelParam("ksatp") / 500.0;
                param.setModelParam(einc, "e_inc");
                std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }
            goto twentyMarker;
        }//

        if (total > 500 || total < 400) { //if//
            // std::cout << "total > 500 or < 400" << std::endl;
            einc = ecritsystem / 450.0; //'re-set Einc
            param.setModelParam(einc, "e_inc");
            if (ecritsystem == 0)
            {
                einc = param.getModelParam("ksatp") / 500.0;
                param.setModelParam(einc, "e_inc");
                std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }
            test_count++; // [DEBUG]
            if (test_count > 10)
            {
                    std::cout << "test count > 10" << std::endl;
                test_count = 0;
                goto fortyMarker;
            }
            goto twentyMarker; //'recalculate the composite curve
        }// //'total ok

    }// //'night <>"n" or it//'s not the first round

fortyMarker:

    bool isNight = true;
    if (night == "n")
        isNight = false;

    if (night == "n" && carbon.psynmax > 0 && carbon.psynmaxsh > 0 && reset_guess == 0) { //if//
                                                                        //DoEvents //'[HNT] this was required to prevent a hard lock -- this portion of the loop is the most intensive, so let Excel take a "breath" by processing system events to prevent lockup
        canopypressure(dd, total, transpiration, transpirationsh, md, mdsh); //'returns canopy P and associated output
                          // 'if check >= 2000 { //if// GoTo 60:
        if (refilling == false)
            updatecurves(halt); //'updates element E(P) curves as required for midday exposure for no refilling
    }// //'night <> "n", psynmax IF

    if (soilred == true) { //if//
        soilflow(); //'gets vertical soil flow between layers in m3/m2
    }
    else {
        for (int z = 0; z <= layers; z++)//z = 0 To layers
        {
            xylem.soils[z]->soilredist = 0;
        } //end for z
    }// //'soil red <> y
    if (ground == true)
        deepflow(timestep); //'gets groundwater flow into bottom layer in m3/m2
                    // '}// //'pet <> y or n
    if (gs_inGrowSeason && sevap == true) { //if//
        soilevaporation(soiltemp, maxvpd, vpd, airtemp, us); //'gets soil evaporation rate
    }
    else {
        soilevap = 0;
    }//

    if (failure == 0 || reset_guess == 1) { //if//
                                        //Debug.Print "DOING A LOOP-8 " & dd

        if (night == "y" || carbon.psynmax == 0 || carbon.psynmaxsh == 0) { //if// //'set everything to starting point
            int k = 0;
            transpiration = xylem.leaf.eplantl[k]; //'all gas exchange values are for closed stomata
            md = xylem.leaf.getPressureComp(k);
            carbon.psynact = carbon.psyn[k];
            carbon.gcmd = carbon.gcanw[k]; //'g for water in mmol
            xylem.leaf.lavpdmd = xylem.leaf.lavpd[k] * patm;
            carbon.cinc = carbon.cin[k];
            halt = k;
            transpirationsh = xylem.leaf.eplantl[k]; //'all gas exchange values are from most recent historical values
            mdsh = xylem.leaf.getPressureComp(k);
            carbon.psynactsh = carbon.psynsh[k];
            carbon.gcmdsh = carbon.gcanwsh[k]; //'g for water in mmol
            xylem.leaf.lavpdshmd = xylem.leaf.lavpdsh[k] * patm;
            carbon.cincsh = carbon.cinsh[k];
            haltsh = k; //'halt is index of midday datum
        }// //'night<>y

        data.setColumnValue(xylem.leaf.getPressureComp(0), dd, "P-PD"); //'the predawn
                                                                    //'SUN LAYER OUTPUT
        data.setColumnValue(md, dd, "P-MD"); //'the midday
        data.setColumnValue(transpiration, dd, "E-MD"); //'midday transpiration, mmol s-1 m-2 leaf area
        data.setColumnValue(carbon.gcmd, dd, "GW"); //'midday canopy diffusive conductance to water, mmol s-1m-2
        data.setColumnValue(xylem.leaf.lavpdmd, dd, "leaf-air-vpd"); //'leaf-to-air vpd
        data.setColumnValue(xylem.leaf.leaftemp[halt], dd, "leaftemp"); //'leaf temp
        data.setColumnValue(carbon.psynact, dd, "Anet-la"); //'net A in umol s-1m-2 leaf area
        data.setColumnValue(carbon.cinc * patm * 1000, dd, "ci"); //'partial pressure of CO2 in Pa
        data.setColumnValue(carbon.qsl, dd, "PPFD"); //'umol s-1m-2 photon flux density
        data.setColumnValue(mdsh, dd, "S-P-MD"); //'the midday
        data.setColumnValue(transpirationsh, dd, "S-E-MD"); //'midday transpiration, mmol s-1 m-2 leaf areas
        data.setColumnValue(carbon.gcmdsh, dd, "S-GW"); //'midday canopy diffusive conductance to water, mmol s-1m-2
        data.setColumnValue(xylem.leaf.lavpdshmd, dd, "S-leaf-air-vpd"); //'leaf-to-air vpd
        data.setColumnValue(xylem.leaf.leaftempsh[haltsh], dd, "S-leaftempt"); //'leaf temp
        data.setColumnValue(carbon.psynactsh, dd, "S-Anet-la"); //'A in umol s-1m-2 leaf area
        data.setColumnValue(carbon.cincsh * patm * 1000, dd, "S-ci"); //'partial pressure of CO2 in Pa
        data.setColumnValue(carbon.qsh, dd, "S-PPFD"); //'umol s-1m-2 photon flux density

        if (night == "n")
            transpiration_tree = xylem.leaf.laisl / xylem.leaf.lai * transpiration + xylem.leaf.laish / xylem.leaf.lai * transpirationsh; //'weighted mean
        if (night == "y")
            transpiration_tree = transpiration;
        if (night == "n")
            carbon.atree = xylem.leaf.laisl / xylem.leaf.lai * carbon.psynact + xylem.leaf.laish / xylem.leaf.lai * carbon.psynactsh; //'weighted mean
        if (night == "y")
            carbon.atree = (carbon.psynact + carbon.psynactsh) / 2.0; //'simple average at night when there//'s no sun or shade leaves

        data.setColumnValue(transpiration_tree, dd, "S-E-tree"); //'weighted mean
        data.setColumnValue(carbon.atree, dd, "Anet-tree"); // Anet Tree per Leaf Area (umol s-1m-2)
        //'Cells(16 + dd, o + 35) = dpamax //'shade leaf dpa
        //'HYDRAULIC OUTPUT (BASED ON SUN MD)
        data.setColumnValue(pcritsystem, dd, "Pcrit");
        data.setColumnValue(ecritsystem * (1 / laperba) * (1.0 / 3600.0) * 55.4 * 1000, dd, "Ecrit");
        data.setColumnValue(xylem.leaf.getPressureComp(halt), dd, "P-leaf");
        data.setColumnValue(xylem.stem.getPressureComp(halt), dd, "P-stem");
        data.setColumnValue(xylem.root_pressure[halt], dd, "P-root");
        data.setColumnValue(xylem.stem.getKComp(halt), dd, "K-stem");
        data.setColumnValue(xylem.leaf.getKComp(halt), dd, "K-leaf");

        kplantold = 0;
        if (!isNewYear)
            kplantold = data.getColumnValue("K-plant", dd - 1); // set to last just in case

        if (transpiration > 0) { //if//
            kplantold = xylem.e_p[halt] / (xylem.leaf.getPressureComp(halt) - xylem.leaf.getPressureComp(0));  //'whole plant k at midday in kg hr-1 m-2 basal area...sun value
            data.setColumnValue(kplantold, dd, "K-plant");
        }//
        if (transpiration == 0)
            data.setColumnValue(kplantold, dd, "K-plant"); //'use most recent kplant

        if (kplantold < gs_data.getColumnValue("K-plant", gs_yearIndex) || gs_data.getColumnValue("K-plant", gs_yearIndex) == 0)
            gs_data.setColumnValue(kplantold, gs_yearIndex, "K-plant");

        // k = o + dColF_CP_kroot1 - 1;//43;
        sum = 0;
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            std::ostringstream oss;
            oss << "K-root-" << z;
            data.setColumnValue(xylem.soils[z]->root.getKComp(halt), dd, oss.str()); //'root k at midday, sun fluxes
            sum = sum + xylem.soils[z]->root.getKComp(halt);
        } //end for z
        data.setColumnValue(sum, dd, "K-root-all"); //'total root k at midday
        if (failure == 0) { //if//
            double tempDouble = 0.0;
            tempDouble = 1 / (1 / xylem.leaf.getKmin() + 1 / xylem.stem.getKmin() + 1 / data.getColumnValue("K-root-all", dd)); //'total xylem k
            data.setColumnValue(tempDouble, dd, "K-xylem"); //1 / (1 / kminleaf + 1 / kminstem + 1 / dSheet.Cells(rowD + dd, colD + k + 1 + layers)); //'total xylem k
            if (tempDouble < gs_data.getColumnValue("K-xylem", gs_yearIndex) || gs_data.getColumnValue("K-xylem", gs_yearIndex) == 0)
                gs_data.setColumnValue(tempDouble, gs_yearIndex, "K-xylem");
        }//
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            std::ostringstream oss;
            oss << "E-root-" << z;
            data.setColumnValue(xylem.soils[z]->root.getEComp(halt) * (1 / laperba) * (1.0 / 3600.0) * 55.56 * 1000, dd, oss.str()); //'uptake in mmol s-1m-2 leaf area...sun rate
        } //end for z
        //'Cells(16 + dd, k + 2 + 2 * layers) = failspot //'position of failure at critical point

        // TODO since all IO is being handled as double currently, cannot add these failure notes
        // temporarily putting in an obvious number to flag failure
        for (int z = 1; z <= layers; z++)//z = 1 To 1
        {
            std::ostringstream oss;
            oss << "Layer-" << z << "-failure";
            if (xylem.soils[z]->cavitated == 1)
                data.setColumnValue(-1137, dd, oss.str()); // xylem.soils[z]->failure; //'layers failed at critical point
        } //end for z

        //Debug.Print "DOING A LOOP-9 " & dd
    }// //'failure IF (basically...failure can//'t happen!)


    if (dd == 1 || isNewYear) { //if// //'NOTE: must be sure that pcritsystem is computed for dd=1!!! (i.e., it//'s not computed at night)
        int x = pcritsystem; //'estimate of "permanent wilting point"
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            xylem.soils[z]->swclimit = xylem.soils[z]->rhizosphere.swc(x); //'theta/thetasat at critical point
            xylem.soils[z]->swclimit = xylem.soils[z]->swclimit * xylem.soils[z]->rhizosphere.getThetaSat(); //'convert to water content
            xylem.soils[z]->swclimit = xylem.soils[z]->swclimit * xylem.soils[z]->depth; //'water content left over in m3/m2 ground area
                                                //'sumsoil = sumsoil + (fc[z] - xylem.soils[z]->swclimit) //'sum is total m3 water per m2 ground withdrawn
        } //end for z

        //Debug.Print "DOING A LOOP-11 " & dd
    }//

    //'now...need to reset layer failure status to midday values (not critical point)
    for (int z = 1; z <= layers; z++)//z = 0 To layers
    {
        if (xylem.soils[z]->cavitated == 1) { //if// //'check to see if kminroot[z]=0; otherwise, restore it to function
            if (xylem.soils[z]->failure == "root" && xylem.soils[z]->root.getKmin() > 0) { //if//
                xylem.soils[z]->cavitated = 0;
                xylem.soils[z]->failure = "ok";
            }//
            if (xylem.soils[z]->failure == "rhizosphere" && xylem.soils[z]->root.getKmin() > 0) { //if//
                xylem.soils[z]->cavitated = 0;
                xylem.soils[z]->failure = "ok";
            }//
        }//
    } //end for z

    if (isNewYear)
        isNewYear = false; // always set this

    return 0;
}

/**
 * @brief Retrieves the atmospheric CO2 concentration (in ppm) for a specified year.
 * 
 * This function iterates through a range of years up to MAX_YEARS and checks if the 
 * given year matches the year stored in the parameter object. If a match is found 
 * the function returns the CO2 concentration for that year.
 * 
 * @param yearVal The year for which the atmospheric CO2 concentration is requested.
 * @return double The atmospheric CO2 concentration (in ppm) for the specified year.
 *                Returns 0.0 if no matching year is found or if the CO2 concentration
 *                is not greater than 0.
 */
double Plant::getCarbonByYear(int yearVal) {
    double ca_year = 0.0;
    for (long gsC = 0; gsC <= MAX_YEARS; gsC++)
    {
        if (yearVal == param.getGsArYear(gsC) && param.getGsArPpm(gsC) > 0) {
            ca_year = param.getGsArPpm(gsC); // get current year atmospheric CO2 in ppm
        }
    }
    return ca_year;
}

/**
 * @brief Prepares the plant model for a new year by resetting variables, 
 *        saving previous states, and reinitializing parameters.
 * 
 * This function performs the following steps:
 * 1. Outputs a message indicating the start of a new year.
 * 2. Saves the current states of the water system (ground, raining, rainEnabled).
 * 3. Clears all model variables and deletes soil layers to avoid memory effects 
 *    from previous runs.
 * 4. Retrieves the previous year's mean Xylem PLC (Percent Loss of Conductivity).
 * 5. Reads in all necessary parameters for the model again (need to optimize)
 * 6. Resets growth season variables (e.g., day and growth season status).
 * 7. Configures iteration details for supply curve runs, if applicable.
 * 8. Initializes soil layer properties, such as source pressures and failure states.
 * 
 * @note This function is critical for ensuring the model starts with a clean 
 *       state and accurate parameters for the new simulation year.
 */
void Plant::modelProgramNewYear()
{   
    std::cout << std::endl;
    std::cout << "Starting new year" << std::endl;
    // save the iterate-able water system states
    oldground = ground;
    oldraining = raining;
    oldrainEnabled = rainEnabled;
    
    // clear all the model variables
    /* Deletes all soil layers as well */
    cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
    max_plc_x = gs_data.getColumnValue("PLC-x", gs_yearIndex - 1); // previous year mean PLC?
    std::cout << "Previous year Xylem PLC: " << max_plc_x << std::endl;

    readin(); // get all the parameters again

    gs_prevDay = 0;
    gs_inGrowSeason = false;

    // set up the iteration details based on what we saved
    if (iter_runSupplyCurve) {
        ground = false; //disable ground water for the first two iterations
        raining = false; //disable rain for the first iteration
        rainEnabled = false; //disable rain for the first iteration
        ground = oldground;
        raining = oldraining;
        rainEnabled = oldrainEnabled;

    }

    for (int k = 0; k <= layers; k++) // To layers //assign source pressures, set layer participation
    {
        xylem.soils[k]->failure = "ok";
        xylem.soils[k]->cavitated = false; //1 if out of function
    }
}

/**
 * @brief Determines if the given Julian day (jd) falls within the growing season.
 * 
 * This function checks whether the specified Julian day is within the growing season
 * based on the growing season data (GSData) and parameters. If the growing season
 * limits are disabled, the function always returns true.
 * 
 * @param jd The Julian day to check.
 * @return true If the Julian day is within the growing season or if growing season
 *         limits are disabled.
 * @return false If the Julian day is outside the growing season.
 * 
 * @note The function uses the following conditions:
 *       - If `useGSData` is true, it checks the growing season data for the current year index.
 *       - The year index (`gs_yearIndex`) must be valid (0 <= gs_yearIndex < 100).
 *       - The growing season start and end days are retrieved using `param.getGsArStart` 
 *         and `param.getGsArEnd` respectively.
 *       - If the Julian day falls within the start and end days, the function returns true.
 *       - If `useGSData` is false, the function assumes the growing season is always active.
 */
bool Plant::isInGrowSeasonSimple(const int &jd) {
    if (useGSData)
    {
        if (gs_yearIndex >= 0 && gs_yearIndex < 100)
        {
            if (param.getGsArYear(gs_yearIndex) > 0)
            {
                if (jd >= param.getGsArStart(gs_yearIndex) && jd <= param.getGsArEnd(gs_yearIndex))
                {
                    return true;
                }
            }
        }
        return false;
    }
    else
    {
        return true; // if the GS limits are disabled, we're always in the growing season
    }
}

//'gets predawns accounting for rain, groundwater flow, transpiration, and redistribution via soil and roots
/**
 * @brief Calculates and updates soil wetness and water content dynamics for a plant model.
 * 
 * This function simulates the water content in soil layers, including processes such as 
 * transpiration, soil redistribution, rain infiltration, groundwater input, and runoff. 
 * It also tracks water content changes over time and records various metrics related to 
 * plant and soil water dynamics.
 * 
 * @param dd The current day of the simulation.
 * @param timestep The time step of the simulation in hours.
 * @param lai Leaf area index (total).
 * @param laish Leaf area index for shaded leaves.
 * @param laisl Leaf area index for sunlit leaves.
 * @param laperba Leaf area per basal area.
 * @param atree Net assimilation rate of the tree.
 * @param cinc Internal CO2 concentration.
 * @param ca Atmospheric CO2 concentration.
 * 
 * @details
 * - Initializes soil water content at the start of the year or simulation.
 * - Handles water redistribution between soil layers during day and night.
 * - Simulates rain infiltration and adjusts soil water content accordingly.
 * - Tracks groundwater input and runoff.
 * - Updates water content metrics such as field capacity, saturation, and extraction limits.
 * - Records transpiration, soil evaporation, and total evapotranspiration (ET).
 * - Tracks plant hydraulic conductance (K-plant) and xylem conductance (K-xylem).
 * - Calculates and records plant and xylem percent loss of conductivity (PLC).
 * - Handles growing season and off-season water input/output tracking.
 * 
 * @note This function interacts with multiple external data structures and parameters, 
 * including soil properties, plant parameters, and simulation data tables.
 * 
 * @warning Ensure that all required external data and parameters are properly initialized 
 * before calling this function. Misconfigured inputs may lead to undefined behavior.
 */
void Plant::getsoilwetness(const int &dd, 
                           const int &timestep,
                           const double &lai,
                           const double &laish,
                           const double &laisl,
                           const double &laperba,
                           const double &atree,
                           const double &cinc,
                           const double &ca) 
{
    double tempDouble = 0.0;
    double layerflow = 0;

    std::vector<double> thetafracres(layers+1),
                        thetafracfc(layers+1),
                        thetafc(layers+1);
    if (dd == 1 || isNewYear) { //if// //'every layer starts at initial % of field capacity

        double waterold = 0; //'total root zone water content (mmol m-2 ground area)

        if (dd != 1) // if not first iteration and not new year, then get the previous year water content
            waterold = data.getColumnValue("water-content", dd - 1) / 1000;

        if (!(useGSData && gs_yearIndex > 0))
        {
            waterold = 0;
            
            for (int z = 0; z <= layers; z++)//z = 0 To layers //'
            {
                // int x = 10; //'MPa water potential for getting residual thetafrac
                thetafracres[z] = xylem.soils[z]->rhizosphere.swc(10); //'residual thetafrac
                thetafracfc[z] = (1 - thetafracres[z]) * param.getModelParam("field_cap_frac") + thetafracres[z]; //'thetafrac at field capacity
                thetafc[z] = xylem.soils[z]->rhizosphere.getThetaSat() * thetafracfc[z]; //'water content at field capacity
            } //for//z //'

            for (int z = 0; z <= layers; z++)//z = 1 To layers
            {
                water[z] = thetafc[z] * xylem.soils[z]->depth; //'field capacity estimated as 1/2 saturated capacity, water content of layer in m3 water per m2 ground area
                fc[z] = water[z]; //'records field capacity in m3 water volume per m2 ground area.
                water[z] = param.getModelParam("ffc") * water[z]; //'start off with initial fraction of field capacity
                                            //'Cells(16 + dd, 42 + z) = water[z]
                waterold = waterold + water[z]; //'in m3/m2 ground area
            } //for//z
                                                           
            data.setColumnValue(waterold * 1000, dd, "water-content");  //'root zone water content in mm m-2 ground area
                                                                        //'waterold = waterold * 55555556# //'convert m3 water per ground area to mmol water per ground area
                                                                        // [HNT] starting water now counts as an input
            gs_data.setColumnValue(gs_data.getColumnValue("input") + waterold * 1000, gs_yearIndex, "input");
        }

        if (gs_yearIndex == 0) // if it's the first year, store this as the off-season starting water because we don't have a real value
            gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "water-initial-off");

        // store the initial water content to check how much we consume at the end
        gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "initial");
        // [/HNT]
    }// //'dd=1 if
        //'if pet = "y" Or pet = "n" { //if// //'do the original routine...doesn//'t run for PET scenario
    if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'get flows that happened during previous timestep
        for (int z = 0; z < layers; z++)//z = 0 To layers - 1 //'transpiration, root and soil redistribution
        {
            if (night == "n") { // if it's day, must adjust elayer for sun vs. shade weighting
                layerflow = xylem.soils[z]->root.getEComp(halt) * laisl / lai + xylem.soils[z]->root.getEComp(haltsh) * laish / lai; //'weighted flow; NOTE: elayer = 0 for layer 0
            }
            else {
                layerflow = xylem.soils[z]->root.getEComp(halt); //'no adjustment necessary at night; NOTE: elayer = 0 for layer 0
            }// //'night if
            layerflow = layerflow * param.getModelParam("ba_per_ga") * 1.0 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
            layerflow = layerflow + xylem.soils[z]->soilredist * 1.0 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow). NOTE: xylem.soils(0)->soilredist includes soil evaporation for layer 0
            water[z] = water[z] - layerflow; //'subtracts rootflow from layer on per ground area basis
        } //for//z
        //'now do the bottom layer and potential groundwater input
        if (night == "n") { // if it's day, must adjust layerflow for sun vs. shade weighting
            layerflow = xylem.soils[layers]->root.getEComp(halt) * laisl / lai + xylem.soils[layers]->root.getEComp(haltsh) * laish / lai; //'weighted flow
        }
        else {
            layerflow = xylem.soils[layers]->root.getEComp(halt); //'no adjustment necessary at night
        }// //'night if
        layerflow = layerflow * param.getModelParam("ba_per_ga") * 1 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
        layerflow = layerflow + xylem.soils[layers]->soilredist * 1 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow)
                                                                            //'water(layers) = water(layers) - layerflow //'subtracts rootflow from layer on per ground area basis
        if (layerflow < 0) { //if// //'water is added
            for (int z = layers; z >= 0; z--)//z = layers To 0 Step -1 //'start at bottom, go up
            {
                double deficit = xylem.soils[z]->rhizosphere.getThetaSat() * xylem.soils[z]->depth - water[z]; //'m of water required to wet up layer to SATURATION
                if (-1 * layerflow - deficit >= 0) { //if// //'there//'s enough to wet the layer...remember, negative flow is flow into the layer
                    water[z] = xylem.soils[z]->rhizosphere.getThetaSat() * xylem.soils[z]->depth; //'m water at saturation in layer
                    layerflow = layerflow + deficit; //'reduce what//'s left over for next layers
                }
                else { //'just soak up all the groundwater
                    water[z] = water[z] - layerflow; //'add to bottom layer
                    layerflow = 0; //' all gone
                }// //'wetting if
            } //for//z
            
            runoff = runoff - layerflow; //'add what//'s left to runoff...
        }
        else { //'groundwater is positive...bottom layer is losing water
            water[layers] = water[layers] - layerflow; //'subtract from bottom layer
        }// //'layerflow if

        //'now reset any exhausted layers to extraction limit
        if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            if (water[z] < xylem.soils[z]->swclimit)
                water[z] = xylem.soils[z]->swclimit; //'water at limit
        } //for//z

        bool rainOverride = false;
        if (tod == 23 && (stage_id == STAGE_ID_HIST_OPT || stage_id == STAGE_ID_FUT_OPT))
            rainOverride = true;

        //predawns mode
        if (mode_predawns)
            rainOverride = true;
        // end predawns mode

        //'now check for rain during PREVIOUS TIME STEP
        if (data.getColumnValue("rain", dd - 1) >= 0.0 || rainOverride == true) { //if// //'there//'s been rain!
            double rain = data.getColumnValue("rain", dd - 1) * 0.001; //'convert mm of rain to m depth...equivalent to m3/m2 water volume per ground area
                                                                        //'if raining(xx) = "n" { //if// rain = 0 //'rain turned off

            if (rainOverride) // just force everything to FC on every timestep
            {
                for (int z = 0; z <= layers; z++)//z = 0 To layers //'add in rain
                {
                    water[z] = fc[z];
                }
            } // do this before the rain addition, so that we see all rain as runoff

            if (rainEnabled == false || mode_predawns) // predawns mode also turns off rain
                rain = 0;
            // [HNT] temp
            //if (rain > 0.1) // more than 1/10th meter per hour is a data anomoly, however this is a poor assumption and these should be checked in the weather curator
            //   rain = 0.0;
            // [/HNT]
            //'sumrain = rain * 1000 + sumrain //'total rain input in mm m-2
            for (int z = 0; z <= layers; z++)//z = 0 To layers //'add in rain
            {
                if (rain <= 0)
                    break; //'rain//'s used up
                double deficit = fc[z] - water[z]; //'m3 of water per m2 to wet up layer to fc
                                            //if (deficit >= -0.001) { //if// //'layer is at or below field capacity
                if (rain - deficit >= 0) { //if// //'there//'s enough to wet the layer
                    water[z] = fc[z];
                    rain = rain - deficit;
                    //'sumrain = sumrain + deficit //'absorbed rain
                    drainage = rain * 1000.0; //'any left over will drain out the bottom unless layer is rising above field capacity
                }
                else {
                    water[z] = water[z] + rain; //'it all goes into the first layer
                                                //'sumrain = sumrain + rain
                    rain = 0; //'rain used up
                    drainage = 0;
                }// //'wetting up to field capacity "if"
                //}
            } //for//z

            // If there's drainage, and the ground water is on, that should be used to fill up to saturation
            if (rain > 0.0 && ground) // the remaining "drainage" is also still stored in the rain variable
            {
                // this is kind of inefficient, but the rain routine actually drained all the layers to FC even if GW was on and we should have been filling to sat
                // now we start at the bottom and fill the layers to saturation using the drainage
                for (int j = layers; j >= 0; j--) //j = z To 0 Step -1 //'go back up to fill profile to saturation
                {
                    if (rain <= 0) // if rain is no more
                        break;
                    double deficit = xylem.soils[j]->rhizosphere.getThetaSat() * xylem.soils[j]->depth - water[j];
                    if (deficit >= 0) { //if// //'got capacity
                        if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                            water[j] = xylem.soils[j]->rhizosphere.getThetaSat() * xylem.soils[j]->depth; //'saturate the layer
                            rain = rain - deficit; //'reduce rain
                        }
                        else { //'rain absorbed by layer
                            water[j] = water[j] + rain;
                            rain = 0; //'use up rain
                        }// //'deficit=>0 "if"
                    }
                    else { //'deficit<0...layer//'s saturated
                        rain = rain - deficit; //'increase rain by super-saturated amount (deficit is negative)
                        water[j] = xylem.soils[j]->rhizosphere.getThetaSat() * xylem.soils[j]->depth; //'reset to saturation
                    }// //'deficit <>0 if
                } //for//j

                runoff = runoff + rain; //'whatever is left over will run off
                drainage = 0; //'no drainage if any layer is rising above field capacity
                rain = 0; //'reset rain to zero
            }
            //'sumdrain = sumdrain + drainage //'total drainage
        }// //'rain if

        //'now check for exhausted layers
        if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            if (water[z] < xylem.soils[z]->swclimit)
                xylem.soils[z]->cavitated = true; //'water exhausted
        } //for//z

        //'now get water content change over PREVIOUS time step
        if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'now get updated new water content

            double waternew = 0;
            double waterold = data.getColumnValue("water-content", dd - 1) / 1000;

            for (int z = 0; z <= layers; z++)//z = 0 To layers //'check for exhausted layers
            {
                waternew += water[z];
            } //for//z

            double waterchange = waternew - waterold; //'total water change in mm m-2 ground area
                                            //'waterchange = waterchange * 1 / param.getModelParam("ba_per_ga") * 1 / laperba //'total water in mmol per m2 leaf area
            data.setColumnValue(waternew * 1000, dd, "water-content"); //'root zone water content in mm

                                                                                    // always store the water content as the "final" -- not worth testing if it's really the last day of the year
            gs_data.setColumnValue(waternew * 1000, gs_yearIndex, "water-final");
            data.setColumnValue(waterchange * 1000, dd, "water-content-delta");             //'change in water content over PREVIOUS timestep
                                                                                    //'if raining(xx) = "y" { //if// Cells(16 + dd, 59) = Cells(16 + dd - 1, 4) //'rain input per previous timestep
            if (rainEnabled == true)// && data.getColumnValue("rain", dd - 1) < 100.0)
                data.setColumnValue(data.getColumnValue("rain", dd - 1), dd, "end-rain"); //'rain input per previous timestep
                                                                    //'Cells(16 + dd, 61) = transpiration_tree * 3600 * timestep * laperba * param.getModelParam("ba_per_ga") * 0.000000018 * 1000 //'transpiration per ground area in mm m-2 per timestop
            data.setColumnValue(gwflow, dd, "end-ground-water"); //'groundwater input in mm per timestep
            data.setColumnValue(drainage, dd, "end-drainage"); //'total drainage in mm per timestep
                                                                                                                                                                //'Cells(16 + dd, 64) = water(0) * 1000 //'water in top layer in mm
            data.setColumnValue(data.getColumnValue("end-rain", dd) + data.getColumnValue("end-ground-water", dd), dd, "end-total-water-input"); //'total input per timestep in mm
            
            gs_data.setColumnValue(gs_data.getColumnValue("input", gs_yearIndex) + data.getColumnValue("end-total-water-input", dd), gs_yearIndex, "input");

            data.setColumnValue(runoff * 1000, dd, "end-runoff"); //'excess root zone water per timestep in mm
        }// //'dd>1 if
    }// //'dd>1 if
        //'}////'pet if
    if (dd > 1 && !isNewYear) { //if//
        tempDouble = transpiration_tree * 3600 * timestep * laperba * param.getModelParam("ba_per_ga") * 0.000000018 * 1000;
        data.setColumnValue(tempDouble, dd, "end-E"); //transpiration_tree * 3600 * timestep * laperba * param.getModelParam("ba_per_ga") * 0.000000018 * 1000; //'transpiration per ground area in mm m-2 per timestop
        if (gs_inGrowSeason) // only record growing season E (should not be any non-GS E, but just for safety)
            gs_data.setColumnValue(gs_data.getColumnValue("E", gs_yearIndex) + tempDouble, gs_yearIndex, "E");

        data.setColumnValue(soilevap * 1 / 998.2 * timestep * 1000, dd, "end-soil-evap");  //'evaporative water loss in mm per timestep
        
        tempDouble = transpiration_tree * 3600 * timestep * laperba * param.getModelParam("ba_per_ga") * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
        data.setColumnValue(tempDouble, dd, "end-ET"); //transpiration_tree * 3600 * timestep * laperba * param.getModelParam("ba_per_ga") * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
        gs_data.setColumnValue(gs_data.getColumnValue("ET", gs_yearIndex) + tempDouble, gs_yearIndex, "ET");

        tempDouble = atree * timestep * 3600 * 0.001;
        data.setColumnValue(tempDouble, dd, "end-Anet-la"); //atree * timestep * 3600 * 0.001; //'Anet per timestep in mmoles per leaf area
        if (!std::isnan(tempDouble) && gs_inGrowSeason) // only record growing season A
        {
            gs_data.setColumnValue(gs_data.getColumnValue("Anet", gs_yearIndex) + tempDouble, gs_yearIndex, "Anet");
            //anytime we record A, also record the ci-related outputs
            // TODO CRIT units???
            // Anet in mmoles per leaf area as in final output columns? (calc above)
            if (night == "n") // daytime only
            {
                gs_data.setColumnValue(gs_data.getColumnValue("cica", gs_yearIndex) + cinc / ca, gs_yearIndex, "cica");
                gs_data.setColumnValue(gs_data.getColumnValue("cica-N", gs_yearIndex) + 1, gs_yearIndex, "cica-N");

                gs_data.setColumnValue(gs_data.getColumnValue("Aci", gs_yearIndex) + tempDouble * cinc, gs_yearIndex, "Aci");
                gs_data.setColumnValue(gs_data.getColumnValue("AnetDay", gs_yearIndex) + tempDouble, gs_yearIndex, "AnetDay");
            }
        }
    }// //'dd>1 if

    if (tod == 16 && !gs_doneFirstDay && gs_inGrowSeason && data.getColumnValue("K-plant", dd - 3) > 0.000000001) { //if// //'get midday k//'s for day 1
                                                                                                                      // VPD zero case -- if the stomata did not open on the first day of the GS, kplant won't have been set and will be zero... in which case, try again tomorrow
        gs_doneFirstDay = true;

        double sum = 0.;
        for (int z = 1; z <= 3; z++)//z = 1 To 3
        {
            sum = sum + data.getColumnValue("K-plant", dd - z);
        } //for//z
        kpday1 = sum / 3.0; //'average midday kplant on day one
        sum = 0;
        for (int z = 1; z <= 3; z++)//z = 1 To 3
        {
            sum = sum + data.getColumnValue("K-xylem", dd - z);
        } //for//z
        kxday1 = sum / 3.0; //'average midday kxylem on day one

                            // [HNT] first day of the new growing season! Record the starting water
        gs_data.setColumnValue(data.getColumnValue("water-content", dd), gs_yearIndex, "initial"); // no matter what happened above, this has been set by this point
                                                                                                    // and record this as the FINAL water for THIS YEAR's off season -- remember that this year's off season is the PRECEDING winter
        gs_data.setColumnValue(data.getColumnValue("water-content", dd), gs_yearIndex, "final-off");
        // [/HNT]
    }// //'dd=16 if
    
    if (gs_doneFirstDay) { //if// //'calculate plc relative to midday of day 1
        // if (iter_refK < 0.000000001) // || iter_Counter == 0 // no longer appropriate to test for iter_Counter == 0 ... may be doing a seperate stress profile that refers to a saved refK, which will have been loaded at start of modelProgramMain
        tempDouble = 100 * (1 - data.getColumnValue("K-plant", dd - 1) / kpday1); //' done to avoid repeating this calculation
        // else // if we haven't loaded a ref K it will be set to zero, and we need to fall back to the old method.. otherwise use refK to calculate PLC
        // {
        //     tempDouble = 100 * (1 - data.getColumnValue("K-plant", dd - 1) / iter_refK); // PLCp calculated from refK if it exists
        //     if (tempDouble < 0.0) // if we're using the refK, it's possible for this to be negative briefly -- should be considered zero
        //         tempDouble = 0.0;
        // }
        data.setColumnValue(tempDouble, dd, "end-PLC-plant"); //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1) 'plc plant...prior timestep

                                                                        // no matter what, add it to the tally for calculating the mean over-season PLCp
        gs_data.setColumnValue(gs_data.getColumnValue("PLC-sum", gs_yearIndex) + tempDouble, gs_yearIndex, "PLC-sum"); // yearly running total of PLC values
        gs_data.setColumnValue(gs_data.getColumnValue("PLC-sum-N", gs_yearIndex) + 1, gs_yearIndex, "PLC-sum-N"); // total hours in GS
                                                                        // now test for highest PLC and hours > 85
        if (tempDouble > gs_data.getColumnValue("PLC-p", gs_yearIndex))
            gs_data.setColumnValue(tempDouble, gs_yearIndex, "PLC-p");

        if (tempDouble > 85.0)
            gs_data.setColumnValue(gs_data.getColumnValue("PLC-85", gs_yearIndex) + 1, gs_yearIndex, "PLC-85");

        tempDouble = 100 * (1 - data.getColumnValue("K-xylem", dd - 1) / kxday1);
        data.setColumnValue(tempDouble, dd, "end-PLC-xylem"); //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1) 'plc xylem...prior timestep
        if (tempDouble > gs_data.getColumnValue("PLC-x", gs_yearIndex))
            gs_data.setColumnValue(tempDouble, gs_yearIndex, "PLC-x");

        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1)); //'plc plant...prior timestep
        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1)); //'plc xylem...prior timestep

        // [HNT] keep track of in-season input
        if (gs_inGrowSeason)
        {
            // if we've done the first GS day and we're in the growing season then it's "this year" during GS
            gs_data.setColumnValue(gs_data.getColumnValue("water-input", gs_yearIndex) + data.getColumnValue("end-total-water-input", dd), 
                                   gs_yearIndex, 
                                   "water-input");
        }
        else
        {
            // Make sure that if we are going to next year that this exists in gs_data sheet
            if ((gs_yearIndex + 1) >= gs_data.row_size()) {
                    gs_data.setColumnValue(0.0, gs_yearIndex + 1, "water-initial-off");
            }
            // if we've done the first GS day but we're NOT in the growing season, it's the winter following GS.
            // this is considered NEXT YEAR's off-season input! First check if we already have a value for the FINAL water for THIS YEAR, because if it's zero then
            // this is the first timestep of the winter and we need to store it
            if (gs_data.getColumnValue("water-final", gs_yearIndex) <= 0.0 && gs_data.getColumnValue("water-initial-off", gs_yearIndex + 1) <= 0.0) // ok to use == or <= with double here because this will be memset to 0 if it hasn't been set
            {
                gs_data.setColumnValue(data.getColumnValue("water-content", dd), 
                                       gs_yearIndex, 
                                       "water-final"); // ok to overshoot by 1 hour, I think (instead of using dd - 1)
                gs_data.setColumnValue(data.getColumnValue("water-content", dd), 
                                       gs_yearIndex + 1, 
                                       "water-initial-off"); // also the off-season initial for next year
            }
            else // otherwise, we're in the middle of NEXT YEAR's off season ... note this +1 on an array index is super lazy and bad. Make sure to never run exactly this array size # of years
            {
                gs_data.setColumnValue(gs_data.getColumnValue("water-initial-off", gs_yearIndex + 1) + data.getColumnValue("end-total-water-input", dd), 
                                       gs_yearIndex + 1, 
                                       "water-initial-off"); // add the stored input to the input tally, for NEXT YEAR
            }
        }
        // [/HNT]
    }// //'dd>16 if
    else if (!gs_inGrowSeason)// have NOT done first day and are NOT in growing season
    {
        // we must be in the pre-GS winter of what we called the NEXT year above... so gs_yearIndex is now that year, and we add to the off-season input for THIS year
        gs_data.setColumnValue(gs_data.getColumnValue("water-input-off", gs_yearIndex) + data.getColumnValue("end-total-water-input", dd), 
                               gs_yearIndex, 
                               "water-input-off");
    }
}

/**
 * @brief Calculates the soil predawn water potential for each layer and updates the state of the soil layers.
 * 
 * This function determines the predawn water potential for each soil layer, taking into account the 
 * participation of rooted layers and their cavitation status. It also updates the pressures in the 
 * rhizosphere and root systems, and handles cases where layers are disconnected due to cavitation.
 * 
 * @param dd The index of the current day or time step for which the predawn water potential is calculated.
 * @return int Returns 0 if the calculation is successful, or 1 if one or more layers are disconnected (failure).
 * 
 * @details
 * - The function first checks the participation of each soil layer and updates their cavitation status.
 * - For each layer, it calculates the predawn pressure based on rain or estimated via the Van Genuchten 
 *   function if mode_predawns is disabled.
 * - If the predawn pressure exceeds critical thresholds for the rhizosphere or root, the layer is marked as 
 *   cavitated and disconnected.
 * - The function computes an average predawn pressure for connected layers and assigns it to the root system.
 * - Finally, it stores the calculated predawn pressures for each layer in the output data structure.
 * 
 * @note
 * - Layer 0 is not considered for root participation.
 * - The function uses external parameters and data columns for calculations, such as "rain" and "p_grav".
 * - The function outputs a message to the console for each disconnected layer.
 * 
 */
int Plant::getpredawns(const int &dd) // gets soil predawn water potential for each layer
{
    double theta = 0,
           sum = 0,
           pr = 0;
    int    t = 0,
           failure = 0;

    //'first check for layer participation...only rooted layers, not layer 0
    for (int k = 1; k <= layers; k++)//k = 1 To layers //'assign source pressures, set layer participation
    {
        if (xylem.soils[k]->failure == "root") { //if//
            if (xylem.soils[k]->root.getKmin() == 0) {
                xylem.soils[k]->cavitated = true; //'gone from scene if roots cavitated at midday
            }
            if (xylem.soils[k]->root.getKmin() != 0) { //if// //'roots still around
                xylem.soils[k]->failure = "ok";
                xylem.soils[k]->cavitated = false;
            }//
        }//
        if (xylem.soils[k]->failure == "rhizosphere")
            xylem.soils[k]->cavitated = false; //'layer can come back to life
    } //for//k
        //'after getting water[z] and layer participation, get predawns
        
    for (int z = 0; z <= layers; z++)//z = 0 To layers
    {
        if (xylem.soils[z]->cavitated == false) { //if//
            theta = water[z] / xylem.soils[z]->depth; //'convert m3 water per m2 ground back to m3 water / m3 soil
            double x = theta / xylem.soils[z]->rhizosphere.getThetaSat(); //'remember, VG function takes theta/thetasat as input
            //predawns mode
            if (mode_predawns)
                xylem.soils[z]->predawn_pressure = data.getColumnValue("rain", dd) - param.getModelParam("p_grav"); // read the predawns+pgrav from the "rain" column
            else
                xylem.soils[z]->predawn_pressure = xylem.soils[z]->rhizosphere.rvg(x); //'soil pressure of layer

            //end predawns mode changes
            xylem.soils[z]->rhizosphere.setPressure(xylem.soils[z]->predawn_pressure);
            if (xylem.soils[z]->predawn_pressure >= xylem.soils[z]->rhizosphere.getPcrit() && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
                xylem.soils[z]->cavitated = true;
                xylem.soils[z]->failure = "rhizosphere";
            }//
            if (xylem.soils[z]->predawn_pressure >= xylem.soils[z]->root.getPcrit() && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
                xylem.soils[z]->cavitated = true;
                xylem.soils[z]->failure = "root";
                xylem.soils[z]->root.setKmin(0);
            }//
        }
        else { //'layer//'s disconnected
            std::cout << z << " is disconnected." << std::endl;
            xylem.soils[z]->predawn_pressure = xylem.soils[z]->root.getPcrit();
            xylem.soils[z]->rhizosphere.setPressure(xylem.soils[z]->root.getPcrit());
        }//
    } //for//z
        // exit(1);
        //'now get guess of proot
    sum = 0;
    t = 0;
    for (int k = 1; k <= layers; k++)//k = 1 To layers
    {  
        if (xylem.soils[k]->cavitated == false) { //if//
            sum = sum + xylem.soils[k]->predawn_pressure;
        }
        else { //'predawn is not seen by the roots
            t = t + 1;
        }//
    } //for//k
    // failspot = "no failure";
    if (t < layers) { //if//
        pr = sum / (layers - t); //'set unknown proot to average pd
        xylem.soils[1]->root.setPressure(pr); //'set unknown proot to average pd
    }
    else
        failure = 1;

    for (int z = 0; z <= layers; z++)//z = 1 To layers
    {
        std::ostringstream oss;
        oss << "P" << z;
        data.setColumnValue(xylem.soils[z]->predawn_pressure, dd, oss.str()); //'soil pressures by layer (only for rooted layers)
    }

    return failure;
}

/**
 * @brief Performs LU decomposition on a given square matrix (Jacobian matrix) 
 *        to prepare it for solving linear equations using forward and backward 
 *        substitution.
 * 
 * @param unknowns The number of unknowns (size of the square matrix).
 * @param jmatrix A reference to the square matrix (Jacobian matrix) to be 
 *                decomposed. On output, it contains the LU decomposition of 
 *                the original matrix.
 * @param indx A reference to a vector that stores the permutation indices 
 *             resulting from partial pivoting during the decomposition.
 * 
 * @details
 * The function decomposes the input matrix `jmatrix` into a lower triangular 
 * matrix (L) and an upper triangular matrix (U) such that the product of L and 
 * U equals the original matrix. Partial pivoting is used to improve numerical 
 * stability, and the permutation indices are stored in the `indx` vector. The 
 * diagonal elements of L are assumed to be 1 and are not stored explicitly.
 * 
 * If the matrix is singular (i.e., a row is entirely zero), the function exits 
 * early. To avoid division by zero, a small value (1E-25) is assigned to 
 * diagonal elements that are zero.
 * 
 * @note The matrix indices are assumed to be 1-based as the top soil layer is 
 *       not included, the input matrix `jmatrix` should be resized 
 *       accordingly. The `indx` vector should also be resized to accommodate 
 *       the number of unknowns.
 */
void ludcmp(const int &unknowns, std::vector<std::vector<double>> &jmatrix, std::vector<double> &indx) { //'does LU decomposition on the jacobian prior to solution by lubksb
    int imax = 0;
    double aamax = 0,
           sum = 0,
           dum = 0;
    std::vector<double> vv(unknowns + 1, 0.0);

    for (int i = 1; i <= unknowns; i++)//i = 1 To unknowns
    {
        aamax = 0;
        for (int j = 1; j <= unknowns; j++)//j = 1 To unknowns
        {
            if (std::abs(jmatrix[i][j]) > aamax)
            aamax = std::abs(jmatrix[i][j]);
        }
        if (aamax == 0)
            return;
        vv[i] = 1 / aamax;
    }
    for (int j = 1; j <= unknowns; j++)//j = 1 To unknowns
    {
        for (int i = 1; i < j; i++) //i = 1 To j - 1
        {
            sum = jmatrix[i][j];
            for (int k = 0; k < i; k++) //k = 1 To i - 1
            {
                sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
        }
        aamax = 0;
        for (int i = j; i <= unknowns; i++)//i = j To unknowns
        {
            sum = jmatrix[i][j];
            for (int k = 0; k < j; k++) //k = 1 To j - 1
            {
                sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
            dum = vv[i] * std::abs(sum);
            if (dum > aamax)
            {
                imax = i;
                aamax = dum;
            }
        }
        if (j != imax) // j <> imax
        {
            for (int k = 1; k <= unknowns; k++) //k = 1 To unknowns
            {
                dum = jmatrix[imax][k];
                jmatrix[imax][k] = jmatrix[j][k];
                jmatrix[j][k] = dum;
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (jmatrix[j][j] == 0)
            jmatrix[j][j] = 1E-25;
        if (j != unknowns)
        {
            dum = 1 / jmatrix[j][j];
            for (int i = j + 1; i <= unknowns; i++) //i = j + 1 To unknowns - 1
            {
                jmatrix[i][j] = jmatrix[i][j] * dum;
            }
        }
    }
}

/**
 * @brief Solves the decomposed Jacobian matrix for delta values using back substitution.
 * 
 * This function performs forward and backward substitution to solve a system of linear equations
 * represented by a decomposed Jacobian matrix. It modifies the input vector `func` to contain
 * the solution of the system.
 * 
 * @param unknowns The number of unknowns in the system (size of the matrix).
 * @param jmatrix A 2D vector representing the decomposed Jacobian matrix (LU decomposition).
 * @param func A vector representing the right-hand side of the equation. It is modified in-place
 *             to store the solution of the system.
 * @param indx A vector containing the permutation indices from the LU decomposition.
 *             It is used to rearrange the rows of the matrix during substitution.
 * 
 * @note The input matrix `jmatrix` is assumed to be in LU-decomposed form, and the `indx` vector
 *       must be generated from the decomposition process (e.g., using a function like `ludcmp`).
 *       The function modifies the `func` vector in-place to store the solution.
 */
void lubksb(const int &unknowns, std::vector<std::vector<double>> &jmatrix, std::vector<double> &func, std::vector<double> const &indx) { //'solves the decomposed jacobian for delta p's
    int ii = 0;
    int ll = 0;
    double sum = 0;
    for (int i = 1; i <= unknowns; i++)//i = 1 To unknowns
    {
        ll = int(indx[i]); //'indx array comes from ludcmp
        sum = func[ll]; //'the func array input is the right-hand vector
        func[ll] = func[i];
        if (ii != 0)
        {
            for (int j = ii; j < i; j++)//j = ii To i - 1
            {
                sum = sum - jmatrix[i][j] * func[j];
            }
        }
        else
        {
            if (sum != 0)
                ii = i;
        }
        func[i] = sum;
    }
    for (int i = unknowns; i >= 1; i--) //i = unknowns To 1 Step -1
    {
        sum = func[i];
        for (int j = i + 1; j <= unknowns; j++) //j = i + 1 To unknowns
        {
            sum = sum - jmatrix[i][j] * func[j];
        }
        func[i] = sum / jmatrix[i][i];
    }
}

/**
 * @brief Implements the Newton-Raphson method to calculate rhizosphere pressures 
 *        and root pressure as a function of predawn pressures and transpiration rate.
 * 
 * If root and rhizosphere components are partitioned into N parallel 
 * paths that drain given a known P_soil, there are then N + 1 unknown 
 * pressures:
 * 
 * E_{i(rhizosphere)} - E_{i(root)} = 0, root surface pressure
 * 
 * \sum{E_{i(root)}} - E = 0, root crown pressure
 * 
 * This is solved via multidimensional Newton-Raphson.
 * 
 * @param dd The day index or identifier for the current simulation step.
 * @param p_inc Incremental pressure step used in flow calculations.
 * @param e Transpiration rate or water demand.
 * @return int The number of restarts required for convergence or failure.
 * 
 * This function iteratively solves for the pressures in the rhizosphere and root 
 * using the Newton-Raphson method. It adjusts guesses for pressures and recalculates 
 * flows and derivatives until convergence is achieved or a maximum number of retries 
 * is reached. The function also handles cavitation events and tracks failures.
 * 
 * Key steps:
 * - Initializes or resets pressure guesses.
 * - Iteratively updates pressures using the Jacobian matrix and flow conservation equations.
 * - Detects cavitation in soil layers and adjusts calculations accordingly.
 * - Tracks convergence thresholds and restarts the process if necessary.
 * - Updates failure statistics in case of non-convergence.
 * 
 * Convergence criteria:
 * - The threshold of the flow conservation equations must be below 0.02.
 * - The maximum number of retries is limited to 500.
 * 
 * Failure handling:
 * - If convergence fails, the function logs failure statistics and resets the system state.
 * - Cavitated layers are marked as non-functional and excluded from further calculations.
 * 
 * Preconditions:
 * - The xylem and soil layers must be properly initialized with valid pressure and flow properties.
 * 
 * Postconditions:
 * - The pressures in the rhizosphere and root are updated to steady-state values.
 * - Cavitated layers are identified and marked.
 */
int Plant::newtonrhapson(const int &dd, const double &p_inc, const double &e) { //returns rhizosphere pressures and root pressure, pr, as function of pd's and e
    bool reset_guess = false; //tracks pr estimate
    int heck = 0; //restart loop counter
    int check = 0;
    int ticks = 0;

    std::string failspot;
    
    double pr = xylem.soils[1]->root.getPressure(), // root pressure from top layer
           frt = 0,
           dfrdpr = 0,
           p1 = 0,
           p2 = 0,
           flow = 0,
           klower = 0,
           kupper = 0,
           initialthreshold = 0,
           threshold = 0;

    std::vector<double> func(unknowns + 1, 0.0),
                        dfrhdprh(unknowns + 1, 0.0),
                        dfrhdpr(unknowns + 1, 0.0),
                        dfrdprh(unknowns + 1, 0.0);

    do //loop to reset guesses
    {
        //'restore layer functioning
        for (int z = 1; z <= layers; z++)
        {  
            if (xylem.soils[z]->root.getKmin() != 0)
            {
                xylem.soils[z]->cavitated = false; //make sure to start with all layers functioning
                xylem.soils[z]->failure = "no failure"; //reset
            }
        }

        // failspot = "no failure";
        if (reset_guess == 1) //reset guesses
        {
            // int k = (rand() % (layers)) + 1; // new rand used
            double rFloat = (double)rand() / (double)RAND_MAX;
            int k = int((layers - 1 + 1) * rFloat + 1); // old rand in line with garisom 2.0.5
            pr = xylem.soils[k]->predawn_pressure; //random choice of pd
            xylem.soils[1]->root.setPressure(pr);

            if (false) // can enable this is running into erroneous solutions -- but allowing these to vary instead of resetting results in more frequent solutions in my experience
            { // alternatively could randomize them properly
                for (k = 1; k <= layers; k++) //reset prhz(k)
                {
                    xylem.soils[k]->rhizosphere.setPressure(xylem.soils[k]->predawn_pressure);
                }
            }

            reset_guess = 0; //reset DPA_MAX_CUTOFF
        } //end reset guesses loop
        check = check + 1; //number of restarts
        ticks = 0; //convergence counter

        do //loop to seek convergence
        {
            ticks = ticks + 1;

            if (ticks > 1000)
            {
                reset_guess = 1;
                std::cout << "NR ticks exceeded 1000 -- setting reset_guess for retry. pinc too high? dd = " << dd << std::endl;
            }
            //'get top row of matrix and right-hand func vector
            //'zero out the jacobian first
            for (int k = 1; k <= unknowns; k++)
            {
                for (int j = 1; j <= unknowns; j++)
                    jmatrix[k][j] = 0;
            }

            //'fill up four arrays:
            //'func(i) is zero flow function...the right-hand flow vector
            //'dfrhdprh(i) is partial derivative of frh(i) for prh(i)(rhizo pressure)...the diagonal of the jacobian
            //'dfrhdpr(i) is partial derivative of frh(i)for pr (root pressure)...the last column of the jacobian
            //'dfrdprh(i) is partial derivative of fr for prh(i)...the last row of the jacobian
            frt = 0; //this is the last row of right-hand flow vector
            dfrdpr = 0; //this is lower right-hand partial for jacobian
            for (int z = 1; z <= layers; z++)
            {
                /* If predawn pressure or soil pressure is past PCrit, then cavitation has occured */
                if (xylem.soils[z]->predawn_pressure >= xylem.soils[z]->rhizosphere.getPcrit() || xylem.soils[z]->rhizosphere.getPressure() >= xylem.soils[z]->rhizosphere.getPcrit())
                {
                    xylem.soils[z]->cavitated = 1; //layer's just gone out of function
                    xylem.soils[z]->failure = "rhizosphere";
                    //reset_guess = 1; //could be result of non-convergence
                }
                if (xylem.soils[z]->cavitated == 0)  //it's functional
                {
                    p1 = xylem.soils[z]->predawn_pressure; //p1 is not a guess
                    p2 = xylem.soils[z]->rhizosphere.getPressure(); //prh[z] IS a guess and it is initially set prior to the RN routine
                    xylem.soils[z]->rhizosphere.calc_through_flow(p1, p2, p_inc, flow, klower, kupper); //gets flows through rhizosphere element from p2-p1 and E(p)curve; gets K's at p1 and p2 as well
                    func[z] = flow;
                    dfrhdprh[z] = kupper;
                }
                if (xylem.soils[z]->rhizosphere.getPressure() >= xylem.soils[z]->root.getPcrit() || pr >= xylem.soils[z]->root.getPcrit())
                {
                    xylem.soils[z]->cavitated = 1; //layer's just gone out of function
                    xylem.soils[z]->failure = "root";
                    //reset_guess = 1; //could be result of non-convergence
                }
                if (xylem.soils[z]->cavitated == 0) //it's functional
                {
                    p1 = xylem.soils[z]->rhizosphere.getPressure(); //now re-set p1 to prh[z]...the guess
                    p2 = pr; //guess comes from average predawn pressure across layers, set at soils[0]->root
                    xylem.soils[z]->root.calc_through_flow(p1, p2, p_inc, flow, klower, kupper); //gets flows through root element, and K's at upper and lower bounds
                    func[z] = func[z] - flow;
                    dfrhdprh[z] = dfrhdprh[z] + klower;
                    dfrhdpr[z] = -kupper;
                    dfrdprh[z] = -klower;
                }
                if (xylem.soils[z]->cavitated == 1)  //layer's out of function...zero out values
                {
                    func[z] = 0;
                    dfrhdprh[z] = 1E-25;
                    dfrdprh[z] = 1E-25;
                    dfrhdpr[z] = 1E-25;
                    kupper = 1E-25;
                    flow = 0;
                }

                dfrdpr = dfrdpr + kupper;
                frt = frt + flow;
            }
            frt = frt - e;
                //'now load jacobian
            for (int k = 1; k <= layers; k++)
            {
                jmatrix[k][unknowns] = dfrhdpr[k]; //last column with dFrh/dPr partials
            }
            for (int k = 1; k <= layers; k++)
            {
                jmatrix[unknowns][k] = dfrdprh[k]; //last row with dFr/dPrh partials
            }
            for (int k = 1; k <= layers; k++)
            {
                jmatrix[k][k] = dfrhdprh[k]; //diagonal of dFrh/dPrh partials
            }
            jmatrix[unknowns][unknowns] = dfrdpr; //lower right corner with dFr/dPr partial
            func[unknowns] = frt; //last position in right-hand flow vector

                                //'ok, jacobian and righthand vector are loaded
                                //'test for total failure
            double sum = 0;
            for (int k = 1; k <= layers; k++)
            {
                sum = sum + xylem.soils[k]->cavitated;
            }
            if (sum == layers)  //total failure
            {
                failspot = "belowground";
                reset_guess = 1; //trigger a restart
            }
            //'test for flow conservation (steady-state)
            threshold = 0;
            for (int k = 1; k <= unknowns; k++) //k = 1 To unknowns
            {
                threshold = threshold + std::abs(func[k]);
            }
            if (ticks == 1)
                initialthreshold = threshold;
            //'remember to replace "n" with "unknowns" in ludcmp and lubksb
            std::vector<double> indx(unknowns + 1, 0.0);
            ludcmp(unknowns, jmatrix, indx); //numerical recipe for doing LU decomposition of jacobian prior to solving
            lubksb(unknowns, jmatrix, func, indx); //solves the decomposed jacobian for delta p's
                    //'print out solution vector of pressures
                    //'revise unknown pressures
            for (int k = 1; k <= layers; k++)//k = 1 To layers
            {
                xylem.soils[k]->rhizosphere.setPressure(xylem.soils[k]->rhizosphere.getPressure() - func[k]); //NOTE lubksb replaces original right-side func()vector with the solution vector
            }
            pr = pr - func[unknowns];
            //'check for jumping lower bound
            for (int k = 1; k <= layers; k++)//k = 1 To layers
            {
                if (xylem.soils[k]->rhizosphere.getPressure() < 0)
                    xylem.soils[k]->rhizosphere.setPressure(0);
            }
            if (pr < 0)
                pr = 0;
            xylem.soils[1]->root.setPressure(pr);
            //'if pr > pcritr Then
            //'pr = prinitial
            //'reset_guess = 1 //trigger a re-start
            //'}
            if (ticks > 1)  //check for non convergence
            {
                if (threshold > initialthreshold) {
                    reset_guess = 1;//pr is spiraling, restart NR with new guesses
                }
                if (pr >= xylem.stem.getPcrit()) {
                    reset_guess = 1; //
                }
            }
        } while (!(threshold < 0.02 || reset_guess == 1));
        //Loop Until threshold < 0.01 Or reset_guess = 1 //reset_guess = 1 restarts NR
    } while (!((threshold < 0.02 && reset_guess == 0) || check > 500));

    if (check > 500)
    {
        double waterold = data.getColumnValue("water-content", dd) / 1000;
        // disable this output if it's causing too much spam
        //std::cout << "NR Failure " << threshold << " check = " << check << " dd = " << dd << " watercontent = " << waterold * 1000.0 << " reset_guess = " << reset_guess << std::endl;
        // keep track of the frequency of NR failures
        gs_data.setColumnValue(gs_data.getColumnValue("fail-converge", gs_yearIndex) + 1,
                               gs_yearIndex,
                               "fail-converge"); // non-convergent failure
        gs_data.setColumnValue(gs_data.getColumnValue("fail-converge-water", gs_yearIndex) + waterold * 1000.0,
                               gs_yearIndex,
                               "fail-converge-water"); // non-convergent failure
        if (waterold * 1000.0 > gs_data.getColumnValue("fail-converge-water-max", gs_yearIndex)) 
            gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "fail-converge-water-max");

    }

    // //Loop Until threshold < 0.01 And reset_guess = 0 Or check > 500 //give up after 2000 restarts
    // //'if check >= 500 Then Stop

    // //final step -- recheck the layers
    for (int z = 1; z <= layers; z++)
    {
        if (xylem.soils[z]->root.getKmin() != 0)
        {
            xylem.soils[z]->cavitated = false; //make sure to start with all layers functioning
            xylem.soils[z]->failure = "no failure"; //reset
        }

        if (xylem.soils[z]->predawn_pressure >= xylem.soils[z]->rhizosphere.getPcrit() || xylem.soils[z]->rhizosphere.getPressure() >= xylem.soils[z]->rhizosphere.getPcrit())
        {
            xylem.soils[z]->cavitated = true; //layer's just gone out of function
            xylem.soils[z]->failure = "rhizosphere";
        }
        if (xylem.soils[z]->rhizosphere.getPressure() >= xylem.soils[z]->root.getPcrit() || pr >= xylem.soils[z]->root.getPcrit())
        {
            xylem.soils[z]->cavitated = true; //layer's just gone out of function
            xylem.soils[z]->failure = "root";
        }
    }

    // xylem.soils[1]->root.setPressure(pr); // set initial root pressure

    return check;
}

/**
 * @brief Computes the composite E(P) curve and element conductances for the plant.
 * 
 * This function calculates the flow and conductance for each soil layer, root, stem, 
 * and leaf components of the plant based on the given parameters. It also handles 
 * scenarios such as cavitation, refilling, and failure of root or rhizosphere elements.
 * The results are stored in the respective components of the plant model.
 * 
 * @param e The total flow rate (E) for the plant.
 * @param p The current pressure index.
 * @param total Reference to an integer where the total pressure index will be stored.
 * @return int Returns 1 if the leaf pressure compensation does not change between 
 *         consecutive pressure indices, otherwise returns 0.
 * 
 * @details
 * - For each soil layer, the function calculates the flow through the root and 
 *   updates the root and rhizosphere conductances.
 * - Handles cavitation and refilling scenarios for soil layers.
 * - Updates the pressure and conductance for the stem and leaf components.
 * - Computes the whole plant conductance and stores the results in the plant model.
 * - If the flow rate (e) is zero, the function handles conductance based on whether 
 *   refilling is enabled or not.
 * - Updates the critical pressure and flow rate for the system.
 * 
 * @note The function assumes that the plant model and its components are properly 
 *       initialized before calling this function.
 */
int Plant::compositeCurve(const double &e, const int &p, int &total) //'stores composite E(P)curve and the element conductances
{
    double p1 = 0,
           p2 = 0,
           flow = 0,
           klower = 0,
           kupper = 0,
           x = 0;

    int test = 0;
    xylem.soils[0]->root.setEComp(p, 0); //'no root mediated flow in topmost layer
    for (int z = 1; z <= layers; z++)//z = 1 To layers
    {
        if (xylem.soils[z]->cavitated == false) { //if//
            xylem.soils[z]->rhizosphere.setPressureComp(p, xylem.soils[z]->rhizosphere.getPressure());
            p1 = xylem.soils[z]->rhizosphere.getPressure();
            p2 = xylem.soils[1]->root.getPressure();
            xylem.soils[z]->root.calc_through_flow(p1, p2, param.getModelParam("p_inc"), flow, klower, kupper);
            xylem.soils[z]->root.setEComp(p, flow); // flow through layer
            if (flow != 0) {
                xylem.soils[z]->root.setKComp(p, std::abs(xylem.soils[z]->root.getEComp(p) / (xylem.soils[1]->root.getPressure() - xylem.soils[z]->rhizosphere.getPressureComp(p))));
            }
            if (flow == 0) { //if//
                if (refilling == true) { //if// //'for refilling, starting point is always weibull
                    x = xylem.soils[z]->predawn_pressure;
                    xylem.soils[z]->root.setKComp(p, xylem.soils[z]->root.wb(x));
                }//
                if (refilling == false)
                    xylem.soils[z]->root.setKComp(p, xylem.soils[z]->root.getKmin());
            }//
        }//
        if (xylem.soils[z]->cavitated == true) { //if//
            xylem.soils[z]->root.setEComp(p, 0); //'no flow
            if (xylem.soils[z]->failure == "root") { //if// //'root element has failed
                xylem.soils[z]->root.setKComp(p, 0); //'total cavitation in root
                xylem.soils[z]->rhizosphere.setPressureComp(p, xylem.soils[z]->predawn_pressure); //'rhizosphere pressure returns to the predawn value
            }//
            if (xylem.soils[z]->failure == "rhizosphere") { //if// //'rhizosphere element has failed
                x = xylem.soils[1]->root.getPressure();
                xylem.soils[z]->root.setKComp(p, xylem.soils[z]->root.wb(x)); //'root element conductance = instantaneous conductance from weibull curve at pr
                xylem.soils[z]->rhizosphere.setPressureComp(p, xylem.soils[z]->rhizosphere.getPcrit());
            }//
        }//
    } //for//z
    xylem.root_pressure[p] = xylem.soils[1]->root.getPressure();
    xylem.stem.setPressureComp(p, xylem.stem.getPressure());
    xylem.leaf.setPressureComp(p, xylem.leaf.getPressure());
    
    if (e > 0) { //if// 
        xylem.leaf.setKComp(p, e / (xylem.leaf.getPressure() - xylem.stem.getPressure())); //'leaf element conductance
        xylem.stem.setKComp(p, e / (xylem.stem.getPressure() - xylem.soils[1]->root.getPressure() - param.getModelParam("p_grav"))); //'stem element conductance subtracting extra gravity drop
        xylem.k[p] = e / (xylem.leaf.getPressure() - xylem.leaf.getPressureComp(0) - param.getModelParam("p_grav")); //'whole plant k, subtracting extra gravity drop
    }
    else {
        if (refilling == false) { //if//
            xylem.leaf.setKComp(p, xylem.leaf.getKmin()); //'leaf element conductance
            xylem.stem.setKComp(p, xylem.stem.getKmin()); //'stem element conductance subtracting extra gravity drop
            xylem.k[p] = this->kmin; //'whole plant k, subtracting extra gravity drop
        }//
        if (refilling == true) { //if//
            xylem.leaf.setKComp(p, xylem.leaf.wb(xylem.leaf.getPressure()));
            xylem.stem.setKComp(p, xylem.stem.wb(xylem.stem.getPressure())); //'stem element conductance subtracting extra gravity drop
            //'kplant[p]=??? ignore kplant setting for e=0...probably not important
        }//
        //'dedp[p] = kminplant
    }//
    xylem.e_p[p] = e; //'total flow
    if (p > 0) { //if//
        if (xylem.leaf.getPressureComp(p) - xylem.leaf.getPressureComp(p - 1) == 0) { //if//
            test = 1;
        }
        // else { // These are never used...
        //     dedp[p] = einc / (pleaf[p] - pleaf[p - 1]); //'dedp=instantaneous K of system
        //     dedpf[p] = dedp[p] / dedp[1]; //'fractional canopy conductance
        // }//
    }//
    pcritsystem = xylem.leaf.getPressure();
    ecritsystem = e;
    total = p;

    return test;
}

/**
 * @brief Calculates radiative terms for energy balance and assimilation in a 
 * plant canopy.
 * 
 * This function computes various solar radiation and energy balance parameters, 
 * including solar declination, zenith angle, azimuth angle, beam and diffuse 
 * radiation, and their effects on leaf area index (LAI), photosynthetically 
 * active radiation (PAR), and near-infrared (NIR) radiation. It also calculates 
 * longwave irradiance from the sky and ground.
 * 
 * @param dd         Day of the year (integer).
 * @param jd         Julian day (integer).
 * @param obssolar   Observed solar radiation on the horizontal surface (W/m²).
 * @param maxvpd     Maximum vapor pressure deficit (kPa).
 * @param airtemp    Air temperature (°C).
 * @param vpd        Vapor pressure deficit (kPa).
 * @param lai        Canopy leaf area index (input).
 * @param laisl      Sunlit leaf area index (output).
 * @param laish      Shaded leaf area index (output).
 * @param qsl        Average PPFD incident on sunlit leaves (µmol m⁻² s⁻¹) 
 *                   (output).
 * @param qsh        Average PPFD incident on shaded leaves (µmol m⁻² s⁻¹) 
 *                   (output).
 * @param ssun       Total solar radiation incident on sunlit leaves (W/m²) 
 *                   (output).
 * @param sshade     Total solar radiation incident on shaded leaves (W/m²) 
 *                   (output).
 * @param sref       Reflected solar radiation (W/m²) (output).
 * @param la         Longwave irradiance from the clear sky (W/m²) (output).
 * @param lg         Longwave irradiance from the ground (W/m²) (output).
 * 
 * @details 
 * The function uses various atmospheric and geometric parameters to calculate 
 * solar radiation components and their interactions with the plant canopy. It 
 * accounts for factors such as solar declination, zenith angle, azimuth angle, 
 * atmospheric transmittance, and cloud cover. The calculations include:
 * - Direct beam and diffuse radiation.
 * - Extinction coefficients for beam and diffuse radiation.
 * - Partitioning of radiation into sunlit and shaded leaf areas.
 * - Photosynthetically active radiation (PAR) and near-infrared (NIR) 
 *   radiation.
 * - Longwave radiation from the sky and ground.
 * 
 * The results are used for energy balance and carbon assimilation modeling in 
 * plant canopies.
 */
void Plant::solarcalc(const int &dd,
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
                      double &lg) //'gets radiative terms for energy balance and assimilation
{

    /* Variables, once everything is ported will move into carbon assimilation or leaf class for use in all solar/carbon calculations */
    double fet = 0.0,
           et = 0.0,
           sm = 0.0,      // not used
           longitude = 0.0,
           tsncorr = 0.0,
           tsn = 0.0,
           sindec = 0.0,
           dec = 0.0,
           cosdec = 0.0,
           tim = 0.0,
           lat = 0.0,
           coszen = 0.0,
           zen = 0.0,
           cosaz = 0.0,
           az = 0.0,
           m = 0.0,
           patm = 0.0,
           sp = 0.0,
           tau = 0.0,
           sb = 0.0,
           sd = 0.0,
           st = 0.0,
           cloud = 0.0,
           fcd = 0.0,
           xang = 0.0,
           kbe = 0.0,
           kbezero = 0.0,
           mleafang = 0.0,
           rad = 0.0,
           sum = 0.0,
           k1 = 0.0,
           t1 = 0.0,
           t2 = 0.0,
           told = 0.0,
           kd = 0.0,
           qd = 0.0,
           qds = 0.0,
           qdt = 0.0,
           qb = 0.0,
           qbt = 0.0,
           qsc = 0.0,
           parsh = 0.0,
           parsl = 0.0,
           parbottom = 0.0,
           nirsh = 0.0,
           nirsl = 0.0,
           ssunb = 0.0, // not used
           ssund = 0.0, // not used
           par = 0.0,
           ppfd = 0.0,  // not used
           ea = 0.0,
           eac = 0.0;

    longitude = param.getModelParam("longitude");
    lat = param.getModelParam("lat");
    tau = param.getModelParam("tau");
    tsncorr = param.getModelParam("tsn_corr");
    patm = param.getModelParam("p_atm");
    xang = param.getModelParam("leaf_angle_param");

    //'j = Cells(14 + i, 3) //'julian day
    fet = 279.575 + 0.9856 * jd; //'jd is julian day, fet is factor for CN eqn 11.4
    fet = fet * PI / 180.0; //'convert to radians
    et = (-104.7 * sin(fet) + 596.2 * sin(2 * fet) + 4.3 * sin(3 * fet) - 12.7 * sin(4 * fet) - 429.3 * cos(fet) - 2 * cos(2 * fet) + 19.3 * cos(3 * fet)) / 3600.0; //'"equation of time" in fraction of HOURS C&N 11.4
                                                                                                                                                                    //'long = Cells(11, 4) //'longitude in degree fraction W
    sm = 15 * int(longitude / 15.0); //'standard meridian east of longitude
                                    //'lc = 0.0666667 * (sm - longitude) //'CN 11.3 for longitude correction in fraction of hours
    tsn = 12 - tsncorr - et; //'time of solar noon in hour fraction from midnight, CN 11.3
                            //'Cells(14 + i, 5) = tsn //'output solar noon
    sindec = PI / 180.0 * (356.6 + 0.9856 * jd);
    sindec = sin(sindec);
    sindec = PI / 180.0 * (278.97 + 0.9856 * jd + 1.9165 * sindec);
    sindec = 0.39785 * sin(sindec); //'sine of the angle of solar declination
                                    //'lat = Cells(10, 4) //'latitude, fraction of degrees N
    dec = atan(sindec / pow((-sindec * sindec + 1), 0.5)); //'arcsin of sindec in radians, dec is solar declination
                                                            //'Cells(14 + i, 6) = dec * 180 / PI //'output solar declination
    cosdec = cos(dec);
    //'timst = Cells(14 + i, 4) //'local standard time, hour fraction from midnight
    tim = 15 * (tod - tsn);
    tim = PI / 180.0 * tim; //'convert to radians
    coszen = sin(lat) * sindec + cos(lat) * cosdec * cos(tim); //'cos of zenith angle of sun from overhead, CN11.1
    zen = atan(-coszen / pow((-coszen * coszen + 1), 0.5)) + 2 * atan(1); //'zenith in radians
                                                                            //'if zen < 1.57 { Cells(14 + i, 7) = zen * 180 / PI //'output when sun//'s up (zen<90)
    cosaz = -(sindec - cos(zen) * sin(lat)) / (cos(lat) * sin(zen)); //'cos of azimuth angle measured counterclockwize from due south CN 11.5
    if (cosaz < -1)
        cosaz = -1; //'keep it in limits
    if (cosaz > 1)
        cosaz = 1; //'ditto
    if (cosaz == 1 || cosaz == -1) { //'keeps stupid acos eqn from crashing
        if (cosaz == 1)
            az = 0;
        if (cosaz == -1)
            az = 3.14159; //'180 in radians
    }
    else {
        az = atan(-cosaz / pow((-cosaz * cosaz + 1), 0.5)) + 2 * atan(1); //'solar az in radians
    } //Endif//
    if (tod > tsn)
        az = 6.28319 - az; //'correct for afternoon hours!
                            //'if zen < 1.57 { Cells(14 + i, 8) = az * 180 / PI//'output azimuth during day
                            //'dayl = (-Sin(lat) * sindec) / (Cos(lat) * Cos(dec))
                            //'dayl = Atn(-dayl / Sqr(-dayl * dayl + 1)) + 2 * Atn(1)
                            //'dayl = 2 * dayl / 15
                            //'Cells(14 + i, 9) = dayl * 180 / PI
    if (zen * 180 / PI < 90) { //'sun//'s up: calculate solar radiation
                                //'pa = Cells(12, 4) //'atmospheric p in kpa
                                //'Cells(16 + dd, 8) = zen * 180 / PI //'zenith angle in degrees
                                //'Cells(16 + dd, 9) = lat
                                //'night = "n" //'its officially day
        m = patm / (101.3 * cos(zen)); //'CN 11.12
                                        //'spo = Cells(9, 4) //'solar constant
                                        //'tau = Cells(8, 4) //'transmittance, CN 11.11
        sp = SOLAR * pow(tau, m); //'direct beam irradiance Wm-2
        sb = sp * cos(zen); //'direct beam irradiance on horizontal surface
                            //'Cells(14 + i, 11) = sb
        sd = 0.3 * (1 - pow(tau, m)) * SOLAR * cos(zen); //'clear sky diffuse radiation
                                                        //'Cells(14 + i, 12) = sd
        st = sd + sb; //'total horizontal irradiance from sun (w/o reflected radiation)
        cloud = SOLAR * pow(0.4, m) * cos(zen); //'overcast threshold
                                                //'stobs = Cells(13, 4) //'observed solar radiation on the horizontal Wm-2
        if (obssolar > 0) { //'we//'ve got solar data
            if (obssolar < st) { //'we//'ve got clouds
                if (obssolar > cloud) { //'we//'ve got partial clouds
                    fcd = 1 - (obssolar - cloud) / (st - cloud); //'fraction for converting beam to diffuse
                    sd = sd / st + fcd * sb / st; //'diffuse/total rato
                    sd = sd * obssolar; //'multiply ratio by total observed to get total diffuse
                    sb = obssolar - sd; //'leftover beam
                    st = obssolar; //'reset to stobs
                }
                else { //'its all clouds
                    sd = obssolar;
                    sb = 0;
                    st = obssolar;
                } //Endif//
            } //Endif// //'if no clouds, everything//'s already set
        } //Endif// //'if no solar data, we assume no clouds
        //'calculate reflected light as if it is equal to light at bottom of canopy
        //'xang = Cells(12, 11) //'leaf angle parameter, CN 15.4
        //'if zen < 1.57 { //'sun//'s up:
        kbe = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction coefficient CN 15.4
                                                                                                                //kbe = Sqr(xang ^ 2 + (Tan(zen)) ^ 2) / (xang + 1.774 * (xang + 1.182) ^ -0.733) 'beam extinction coefficient CN 15.4
        kbezero = xang / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction for zen=0(overhead
        mleafang = atan(-kbezero / pow((-kbezero * kbezero + 1), 0.5)) + 2 * atan(1); //'mean leaf angle in radians
                                                                                    //'Cells(14 + i, 20) = kbe
                                                                                    //'Cells(14 + i, 21) = mleafang * 180 / PI //'mean leaf angle in degrees
                                                                                    //'lai = Cells(13, 11) //'canopy leaf area index
                                                                                    //'abspar = Cells(7, 17) //'absorptivity for PAR of leaves
                                                                                    //'abssol = Cells(8, 17) //'absorptivity for total solar of leaves
                                                                                    //'gdbeamf = Exp(-Sqr(abssol) * kbe * lai) //'CN 15.6, fraction of solar beam radiation reaching ground
                                                                                    //'Cells(14 + i, 22) = gdbeamf
                                                                                    //'now solve for kd (diffuse extinction) by integrating beam over all possible zenith angles from 0 to 90, CN 15.5
        rad = 0;
        sum = 0;
        k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
        t1 = exp(-k1 * lai); //'transmittance CN 15.1
        told = t1 * sin(rad) * cos(rad); //'integral function
        do
        {
            rad = rad + 0.015708; //'0.9 degree intervals, zenith angle in radians
            k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
            t1 = exp(-k1 * lai); //'transmittance CN 15.1
            t2 = t1 * sin(rad) * cos(rad); //'integral function
            sum = sum + (t2 + told) / 2.0 * 0.015708; //'integral sum
            told = t2; //'reset
        } while (!(rad > 1.5708));
        //Loop Until rad > 1.5708 //'loop until 90 degrees
        sum = sum * 2; //'complete summing
        kd = -log(sum) / lai; //'extinction coefficient for diffuse radiation
                            //'Cells(14 + i, 23) = kd //'output
                            //'now...compute shaded leaf ppfd, q denotes a ppfd
        qd = 0.45 * sd * 4.6; //'converts total solar diffuse to PPFD diffuse in umol m-2 s-1
        qds = qd * (1 - exp(-(pow(ABS_PAR, 0.5) * kd * lai))) / (pow(ABS_PAR, 0.5) * kd * lai); //'mean diffuse irradiance for shaded leaves CN p. 261
        qdt = qd * exp(-(pow(ABS_PAR, 0.5) * kd * lai)); //'diffuse irradiance at bottom of canopy, CN p. 255, Eqn 15.6
        qb = 0.45 * sb * 4.6; //'converts total solar beam to PPFD
        qbt = qb * exp(-(pow(ABS_PAR, 0.5) * kbe * lai)); //'direct AND downscattered PPFD at bottom of canopy, CN 15.6
        qb = qb * exp(-(kbe * lai)); //'direct PPFD at bottom of canopy,CN 15.1
        qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
        qsh = qds + qsc; //'average PPFD incident (not absorbed!) on shaded leaves
        qb = 0.45 * sb * 4.6; //'re-set qb to top of canopy
                            //'now get sunlit leaf ppfd
        qsl = kbe * qb + qsh; //'average PPFD incident (not absorbed) on sunlit leaves
                            //'now get sunlit vs. shaded lai
        laisl = (1 - exp(-(kbe * lai))) / kbe; //'sunlit lai
        laish = lai - laisl; //'shaded lai
                            //'now get sun and shade PAR, NIR, & longwave in preparation for energy balance
                            //'in par range
        parsh = qsh / 4.6; //'incident PAR, Wm-2, shaded leaves, 100% diffuse
        parsl = qsl / 4.6;//'incident PAR, sunlit leaves
        parbottom = (qbt + qdt) / 4.6; //'PAR making it through the canopy
                                        //'in near-infrared range (assume same equations as for PAR, but different absorptances and incoming fluxes in Wm-2
                                        //'absnir = Cells(9, 17) //'absorptivity of leaves to NIR
        qd = 0.55 * sd; //'converts total solar diffuse to NIR diffuse in Wm-2
        qds = qd * (1 - exp(-(pow(ABS_NIR, 0.5) * kd * lai))) / (pow(ABS_NIR, 0.5) * kd * lai); //'mean diffuse nir irradiance for shaded leaves CN p. 261
        qdt = qd * exp(-(pow(ABS_NIR, 0.5) * kd * lai)); //'diffuse NIR irradiance at bottom of canopy, CN p. 255, Eqn 15.6
        qb = 0.55 * sb; //'converts total solar beam to NIR
                        //'qb = 1600
        qbt = qb * exp(-(pow(ABS_NIR, 0.5) * kbe * lai)); //'direct AND downscattered NIR at bottom of canopy, CN 15.6
        qb = qb * exp(-(kbe * lai)); //'direct NIR at bottom of canopy,CN 15.1, using same extinction coefficient as for PAR
        qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
        nirsh = (qds + qsc); //'incident NIR on shaded leaves, 100% diffuse
        qb = 0.55 * sb; //'re-set nir qb to top of canopy to get sunlit leaves
        nirsl = kbe * qb + nirsh; //'average incident NIR on sunlit leaves
        sshade = parsh + nirsh; //'total solar incident on shaded leaves, 100% diffuse
        ssun = parsl + nirsl; //'total solar incident on sunlit leaves
        xylem.leaf.sbottom = parbottom + qdt + qbt; //'total solar at bottom of canopy
        ssunb = sb / st * ssun; //'beam solar on sunlit (an approximation)
        ssund = sd / st * ssun;//'diffuse solar on sunlit (approximation)
                                //'abssolar = Cells(8, 17) //'absorptivity of leaves for total solar (0.5)
        sref = (1 - ABS_SOLAR) * sshade; //'reflected light...this used for sun/shade dichotomy
        sref = (1 - ABS_SOLAR) * xylem.leaf.sbottom; //'reflected light for monolayer version
                                        //'these below are used for monolayer version:
        par = 0.45 * st; //'wm-2 in par wavelength...45% of total solar
        ppfd = par * 4.6; //'assumes 4.6 moles photons per Joule conversion factor
    }
    else { //'sun//'s down
        sp = 0; sb = 0; sd = 0; st = 0; sref = 0; par = 0; ppfd = 0; sref = 0; ssun = 0; sshade = 0;
        qsh = 0; qsl = 0; ssunb = 0; ssund = 0; laisl = 0; laish = 0; xylem.leaf.sbottom = 0; //'sun//'s down
                                                                                    //'night = "y" //'it//'s officially night
    } //Endif//
        //'now compute long wave irradiance
    ea = patm * (maxvpd - vpd); //'vapor pressure in kPa
                                //'ta = Cells(3, 19) + 273.15 //'air temp in K
    eac = 1.72 * pow((ea / (airtemp + 273.15)), (1.0 / 7.0)); //'emissivity of clear sky CN  10.10
                                                                //'boltz = Cells(9, 11) //'boltzman constant
    la = eac * SBC * pow((airtemp + 273.15), 4); //'long wave irradiance from clear sky
    lg = 0.97 * SBC * pow((airtemp + 273.15), 4); //'long wave irradiance from ground...assumes equilibrium with air temp
                                                    //'Cells(14 + i, 16) = la
                                                    //'Cells(14 + i, 17) = lg
}

/**
 * @brief Calculates the canopy pressure and related parameters for a plant.
 *        This is the main gain-risk profit-maxmization algorithm.
 *
 * This function computes the midday pressures, transpiration rates, and other
 * physiological parameters for both sunlit and shaded layers of a plant canopy.
 * It uses historical and virgin curves to determine the optimal midday pressures
 * and associated gas exchange values. The function also handles critical system
 * states and updates historical data accordingly. Additionally, it computes the
 * optimal plant pressure using the gain-risk algorithm, optimizing profit with
 * respect to cavitation risk and assimilation gain.
 *
 * @param dd The current day index for data retrieval.
 * @param total The total number of iterations or layers to process.
 * @param transpiration Reference to store the transpiration rate for the sunlit layer.
 * @param transpirationsh Reference to store the transpiration rate for the shaded layer.
 * @param md Reference to store the midday pressure for the sunlit layer.
 * @param mdsh Reference to store the midday pressure for the shaded layer.
 *
 * @details
 * - The function retrieves environmental parameters such as vapor pressure deficit (VPD),
 *   air temperature, and wind speed from the data object.
 * - It calculates the virgin risk and gain curves for both sunlit and shaded layers,
 *   ensuring monotonicity and avoiding negative values.
 * - The function determines the midday pressures by maximizing the profit curve
 *   (difference between revenue and cost) using the gain-risk algorithm.
 * - If the system is critical or fails to find a peak in the profit curve, it resets
 *   midday pressures to predawn values.
 * - Historical data is restored at the end of the function, and if refilling is enabled,
 *   minimum conductance values are recorded.
 *
 * @note
 * - The function uses a Newton-Raphson solver for pressure calculations.
 * - It includes safeguards to handle system failures and ensure monotonicity in
 *   the calculated curves.
 * - Debugging messages are included to identify termination conditions.
 *
 * @throws std::runtime_error If critical system failure occurs.
 */
void Plant::canopypressure(const int &dd,
                           const int &total,
                           double &transpiration,
                           double &transpirationsh,
                           double &md,
                           double &mdsh) 
{
    double pr = xylem.soils[1]->root.getPressure(), // should be able to zero out, but just in case
           dedplzero = 0.0,
           dedpl = 0.0,
           maxkloss = 0.0,
           dpmax = 0.0,
           dpamax = 0.0,
           amaxmax = 0.0,
           dpamin = 0.0;
        //    rmean = 0.0;

    double klossv[CURVE_MAX] = {0},
           amaxfrac[CURVE_MAX] = {0},
           amaxfracsh[CURVE_MAX] = {0},
           dpa[CURVE_MAX] = {0};

    int check = 0,
        totalv = 0,
        failure = 0,
        k = 0,
        p = 0,
        test = 0,
        t = 0;

    double e = 0.0,
           einc = 0.0,
           sum = 0.0,
           dedplmin = 0.0;

    // Need to declare in function because sometimes definition is skipped by markers
    double vpd = data.getColumnValue("D-MD", dd); //'midday vpd in kPa
    vpd = vpd / param.getModelParam("p_atm"); //'vpd in mole fraction
    double airtemp = data.getColumnValue("T-air", dd); //'in C
    double maxvpd = (101.3 / param.getModelParam("p_atm")) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    double wind = data.getColumnValue("wind", dd); //'wind speed
    double laperba = param.getModelParam("leaf_per_basal");

    if (wind < MIN_WIND_THRESH) { //if//
        data.setColumnValue(MIN_WIND_THRESH, dd, "wind"); //'set to minimum wind
        wind = MIN_WIND_THRESH;
    }//
    if (vpd > maxvpd) { //if//
        vpd = maxvpd;
        data.setColumnValue(maxvpd * param.getModelParam("p_atm"), dd, "D-MD"); //'print out maximum vpd
    }//

    //computes carbon-based middays; pressures and cost curve from virgin supply function, gas exchange from historical values
    //store history, get MD and cost function from virgin curve
    if (ecritsystem != 0)
    {

        storehistory(); //stores xylem element curves and failure status, resets layers to full functioning

        sum = 0;
        t = 0;
        for (int k = 1; k <= layers; k++) //assign source pressures, set layer participation
        {
            if (xylem.soils[k]->cavitated == false)
            {
                xylem.soils[k]->rhizosphere.setPressure(xylem.soils[k]->predawn_pressure); //initial guess of unknown rhizosphere pressures
                sum = sum + xylem.soils[k]->predawn_pressure;
            }
            else
            {
                xylem.soils[k]->rhizosphere.setPressure(xylem.soils[k]->root.getPcrit());
                t = t + 1;
            }
        }
        if (t < layers)
        {
            pr = sum / (layers - t); //set unknown proot to average pd
            xylem.soils[1]->root.setPressure(sum / (layers - t)); //'set unknown proot to average pd
        }
        else
        {
            failure = 1; //system is critical
            return;
        }

        test = 0; //=1 if stem or leaf fails
                    //virgin gain function params
        carbon.psynmaxmd = 0;
        carbon.psynmaxshmd = 0;
        //now loop through virgin risk curve
        einc = param.getModelParam("e_inc");
        e = -einc;
        p = -1;
        dedplmin = param.getModelParam("ksatp"); //insures the kloss function is monotonic
        double predawn = 0, plold = 0;
        do
        {
            e = e + einc;
            p = p + 1;
            newtonrhapson(dd, param.getModelParam("p_inc"), e); //pd's already assigned...this solves for p's and e's in fingers of chicken foot

            if (check >= 500)
            {
                //p = p - 1
                //std::cout << "NR failed on virgin curve at e = " << e << " on timestep dd = " << dd << std::endl;
                break; //gone as far as can go
            }

            test = xylem.stem.calc_pressure(e, xylem.soils[1]->root.getPressure(), param.getModelParam("p_grav"), param.getModelParam("p_inc")); //gets stem and leaf pressures
            test = xylem.leaf.calc_pressure(e, xylem.stem.getPressure(), 0, param.getModelParam("p_inc"));
            xylem.leaf.setPressureVirgin(p, xylem.leaf.getPressure()); //pleaf from virgin curve
            if (p == 0)
            {
                predawn = xylem.leaf.getPressure(); //set the predawn...pl returned by "leaf" routine
                plold = xylem.leaf.getPressure();
            }
            if (p > 0)
            {
                if ((xylem.leaf.getPressure() - plold) == 0)
                    break; //gone to failure
                
                if (p == 1) {
                    dedplzero = einc / (xylem.leaf.getPressure() - plold); //note: pl is returned by "leaf" routine
                }
                dedpl = einc / (xylem.leaf.getPressure() - plold);
                
                if (dedpl < dedplmin)
                    dedplmin = dedpl; //insure that kloss only goes up
                klossv[p] = dedplzero - dedpl; //non-normalized kloss from virgin curve
                                            //dedpf(p) = dedpl / dedplzero //fractional k/kmax canopy from virgin curve (units don't matter!)
            }
            if (xylem.leaf.getPressure() >= pcritsystem)
                break; //gone to failure
                    //now get virgin A curve

            xylem.leaf.tempMd(p, e, airtemp, vpd, wind, laperba, param.getModelParam("leaf_width"), param.getModelParam("p_atm")); //gets virgin sun layer leaf temperature from energy balance
            xylem.leaf.tempShadeMd(p, e, airtemp, vpd, param.getModelParam("p_atm")); //gets virgin shade layer leaf temperature
            carbon.assimilationMd(p, 
                                param.getModelParam("g_max"), 
                                param.getModelParam("q_max"), 
                                param.getModelParam("comp_25"), 
                                param.getModelParam("theta_c"), 
                                param.getModelParam("v_max25"), 
                                param.getModelParam("j_max25"), 
                                param.getModelParam("kc_25"), 
                                param.getModelParam("ko_25"), 
                                param.getModelParam("sv_vmax"), 
                                param.getModelParam("sv_jmax"), 
                                param.getModelParam("ha_vmax"), 
                                param.getModelParam("hd_vmax"), 
                                param.getModelParam("hd_jmax"), 
                                param.getModelParam("ha_jmax"), 
                                param.getModelParam("light_curv"), 
                                night, 
                                xylem.leaf.emd, 
                                xylem.leaf.lavpdmd, 
                                xylem.leaf.leaftmd); //gets virgin sun layer photosynthesis
            carbon.assimilationShadeMd(p, 
                                    param.getModelParam("g_max"), 
                                    param.getModelParam("q_max"), 
                                    param.getModelParam("comp_25"), 
                                    param.getModelParam("theta_c"), 
                                    param.getModelParam("v_max25"), 
                                    param.getModelParam("j_max25"), 
                                    param.getModelParam("kc_25"), 
                                    param.getModelParam("ko_25"), 
                                    param.getModelParam("sv_vmax"), 
                                    param.getModelParam("sv_jmax"), 
                                    param.getModelParam("ha_vmax"), 
                                    param.getModelParam("hd_vmax"), 
                                    param.getModelParam("hd_jmax"), 
                                    param.getModelParam("ha_jmax"), 
                                    param.getModelParam("light_curv"), 
                                    night, 
                                    xylem.leaf.emd, 
                                    xylem.leaf.lavpdshmd, 
                                    xylem.leaf.leaftshmd); //gets virgin shade layer photosynthesis
                                    //by now we have assigned psynmd(p) and psynshmd(p) and reset psynmaxes

            plold = xylem.leaf.getPressure();
            //e = e + einc
            //p = p + 1
        } while (!(test == 1 || p > 99900)); //loop to failure or out of "p"s
                                            //Loop Until test = 1 Or p > 99900 'loop to failure or out of "p"s
                                            //If check >= 2000 Then Exit Sub

        if (p <= 2)
        {
            goto tenMarker;
        }

        klossv[0] = 0;
        maxkloss = dedplzero - dedpl; //maximum kloss...may not be kmax
        totalv = p - 1;
        klossv[totalv] = maxkloss;
        //now normalize klossv for virgin pleafv
        for (p = 0; p <= totalv; p++)
        {
            klossv[p] = klossv[p] / maxkloss;
        }

        //now, find the middday from virgin risk and gain
        //First do for sun layer
        p = -1; //p is still index for virgin curve
        dpmax = 0;
        dpamax = -100;
        amaxmax = 0; //ensures the gain function is monotonic
        // rmean = 0.0;
        dpamin = 0.0; // keep track of low values to avoid extreme negative profit curves producing a result
        do
        {
            p = p + 1;
            amaxfrac[p] = carbon.psynmd[p] / carbon.psynmaxmd; //this is the normalized revenue function from the virgin curve
            if (amaxfrac[p] > amaxmax)
            {
                amaxmax = amaxfrac[p];
            }
            else
            {
                amaxfrac[p] = amaxmax;
            } //this insures that amaxfrac monotonically increases
            if (amaxfrac[p] < 0)
                amaxfrac[p] = 0; //no negative gains
            if (klossv[p] < 0)
                klossv[p] = 0; //no negative risks
            dpa[p] = amaxfrac[p] - klossv[p]; //profit, with revenue from virgin curve and cost from virgen one
            if (p < PROFT_MAX_RUN_MEAN - 1)
                rmean = 0;
            if (p >= PROFT_MAX_RUN_MEAN - 1)  //get running mean
            {
                sum = 0;
                for (int i = p - PROFT_MAX_RUN_MEAN + 1; i <= p; i++)
                {
                    sum = sum + dpa[i];
                }
                rmean = sum / PROFT_MAX_RUN_MEAN; //the running mean
            }
            if (rmean < dpamin)
            {
                dpamin = rmean;
            }
            if (rmean < 0)
                rmean = 0; //avoid negative rmean
            if (rmean > dpamax)
            {
                dpamax = rmean;
                md = xylem.leaf.getPressureVirgin(p); //midday pressure for sun layer from virgin curves
            }

        } while (!(einc * p >= param.getModelParam("g_max") * xylem.leaf.lavpd[p] || total == 0 || (rmean < dpamax / DPA_MAX_CUTOFF && p > PROFT_MAX_RUN_MEAN && p > 15) || klossv[p] > 0.9 || p >= totalv));

        //std::cout << "DPA MIN = " << dpamin << " DPA MAX = " << dpamax << std::endl;
        if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax)
        {
            // the profit went more negative than positive, so reset mid-day to predawn
            md = xylem.leaf.getPressureVirgin(0);
        }

        // [HNT] debug
        if (!(rmean < dpamax / DPA_MAX_CUTOFF))
        {
            //std::cout << "Terminated sun layer opt without finding peak! At timestep dd = " << dd << std::endl;
            if (einc * p >= param.getModelParam("g_max") * xylem.leaf.lavpd[p])
            {
                ;// std::cout << "Terminated sun layer opt without finding peak: end case 1 einc * p >= param.getModelParam("g_max") * lavpd[p]" << std::endl;
            }
            if (total == 0)
            {
                std::cout << "Terminated sun layer opt without finding peak: end case 2 total == 0" << std::endl;
            }
            if (p >= totalv)
            {
                std::cout << "Terminated sun layer opt without finding peak: end case 3 exceeded end of virgin curve" << std::endl;
            }
            if (klossv[p] > 0.9)
            {
                std::cout << "Terminated sun layer opt without finding peak: end case 4 klossv[p] > 0.9" << std::endl;
            }
        }
        // [/HNT]

        //while (!(einc * p >= param.getModelParam("g_max") * lavpd[p] || total == 0 || rmean < dpamax / DPA_MAX_CUTOFF && p > PROFT_MAX_RUN_MEAN || klossv[p] > 0.9 || p >= totalv));
        //Loop Until einc * p >= param.getModelParam("g_max") * lavpd(p) Or total = 0 Or rmean < dpamax / DPA_MAX_CUTOFF And p > PROFT_MAX_RUN_MEAN Or klossv(p) > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1

        // dpasun = dpamax; // unused?

                        //now do for shade layer
        if (carbon.psynmaxshmd == 0) //shade layer's below light compensation point, don't open
        {
            mdsh = xylem.leaf.getPressureVirgin(0); //set midday = predawn
        }
        else
        {
            p = -1; //p is still index for virgin curve
            dpmax = 0;
            dpamax = -100;
            amaxmax = 0;

            dpamin = 0.0;
            do
            {
                p = p + 1;
                amaxfracsh[p] = carbon.psynshmd[p] / carbon.psynmaxshmd; //this is the normalized shade revenue function from the virgin curve
                if (amaxfracsh[p] > amaxmax)
                {
                    amaxmax = amaxfracsh[p];
                }
                else
                {
                    amaxfracsh[p] = amaxmax;
                }  //this insures that amaxfrac monotonically increases
                //now loop to find kloss from virgin curve that matches historical pleaf
                if (amaxfracsh[p] < 0)
                    amaxfracsh[p] = 0; //no negative gains
                if (klossv[p] < 0)
                    klossv[p] = 0; //no negative risks
                dpa[p] = amaxfracsh[p] - klossv[p]; //profit, with revenue from historical curve and cost from virgen one
                if (p < PROFT_MAX_RUN_MEAN - 1)
                    rmean = 0;
                if (p >= PROFT_MAX_RUN_MEAN - 1)  //get running mean
                {
                    sum = 0;
                    for (int i = p - PROFT_MAX_RUN_MEAN + 1; i <= p; i++)
                    {
                        sum = sum + dpa[i];
                    }
                    rmean = sum / PROFT_MAX_RUN_MEAN; //the running mean
                }
                if (rmean < dpamin)
                {
                    dpamin = rmean;
                }
                if (rmean < 0)
                    rmean = 0; //avoid negative dpa
                if (rmean > dpamax)
                {
                    dpamax = rmean;
                    mdsh = xylem.leaf.getPressureVirgin(p);  //midday pressure for shade layer from virgin curve
                }

            } while (!(einc * p >= param.getModelParam("g_max") * xylem.leaf.lavpdsh[p] || total == 0 || (rmean < dpamax / DPA_MAX_CUTOFF && p > PROFT_MAX_RUN_MEAN && p > 15) || klossv[p] > 0.9 || p >= totalv));
            //Loop Until einc * p >= param.getModelParam("g_max") * lavpdsh[p] Or total = 0 Or rmean < dpamax / DPA_MAX_CUTOFF And p > PROFT_MAX_RUN_MEAN Or klossv[p] > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1
            if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax)
            {
                // the profit went more negative than positive, so reset mid-day to predawn
                mdsh = xylem.leaf.getPressureVirgin(0);
            }
        } //psynmaxsh if

        k = -1;
        //Range("c17:f10000").ClearContents
        do
        {
            k = k + 1;
        } while (!(xylem.leaf.getPressureComp(k) >= md || xylem.leaf.getPressureComp(k) == 0));
        //Loop Until pleaf(k) >= md Or pleaf(k) = Empty //pleaf from historical curve must match pleaf from virgin curve
        transpiration = xylem.leaf.eplantl[k]; //all gas exchange values are from most recent historical values
        carbon.psynact = carbon.psyn[k];
        carbon.gcmd = carbon.gcanw[k]; //g for water in mmol
        xylem.leaf.lavpdmd = xylem.leaf.lavpd[k] * param.getModelParam("p_atm");
        carbon.cinc = carbon.cin[k];
        //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
        halt = k; //halt is index of midday datum
                    //now do shade layer
        k = -1;
        do
        {
            k = k + 1;
        } while (!(xylem.leaf.getPressureComp(k) >= mdsh || xylem.leaf.getPressureComp(k) == 0));
        //Loop Until pleaf(k) >= mdsh Or pleaf(k) = Empty //pleaf for historical curve must match pleaf from virgin curve
        transpirationsh = xylem.leaf.eplantl[k]; //all gas exchange values are from most recent historical values
        carbon.psynactsh = carbon.psynsh[k];
        carbon.gcmdsh = carbon.gcanwsh[k]; //g for water in mmol
        xylem.leaf.lavpdshmd = xylem.leaf.lavpdsh[k] * param.getModelParam("p_atm");
        carbon.cincsh = carbon.cinsh[k];
        //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
        haltsh = k; //halt is index of midday datum
    } else { // ecritsysem == 0
        goto tenMarker; // no midday
    }
    if (ecritsystem == 0)
    {
    tenMarker:     //no midday
        k = 0;
        transpiration = xylem.leaf.eplantl[k]; //all gas exchange values are from most recent historical values
        carbon.psynact = carbon.psyn[k];
        carbon.gcmd = carbon.gcanw[k]; //g for water in mmol
        xylem.leaf.lavpdmd = xylem.leaf.lavpd[k] * param.getModelParam("p_atm");
        carbon.cinc = carbon.cin[k];
        halt = k;
        transpirationsh = xylem.leaf.eplantl[k]; //all gas exchange values are from most recent historical values
        carbon.psynactsh = carbon.psynsh[k];
        carbon.gcmdsh = carbon.gcanwsh[k]; //g for water in mmol
        xylem.leaf.lavpdshmd = xylem.leaf.lavpdsh[k] * param.getModelParam("p_atm");
        carbon.cincsh = carbon.cinsh[k];
        haltsh = k; //halt is index of midday datum
    }
    gethistory(); //reinstates historical element curves and failure status prior to updating
    if (refilling == true) //need to record midday kmins, uses sunlit pressures
    {
        for (int z = 1; z <= layers; z++)
        {
            xylem.soils[z]->root.setKmin(xylem.soils[z]->root.getKComp(halt));
        }
        xylem.stem.setKmin(xylem.stem.getKComp(halt));
        xylem.leaf.setKmin(xylem.leaf.getKComp(halt));
    }
}

/**
 * @brief Updates the hydraulic conductivity curves for the plant's root, stem, 
 *        and leaf components.
 * 
 * This function recalculates the element E(P) curves for the plant's hydraulic 
 * system, ensuring that the hydraulic conductivity (K) values do not fall below 
 * their minimum thresholds (Kmin). If the conductivity is below the threshold, 
 * the function back-calculates the corresponding pressure (P) and updates the 
 * E(P) and K(P) curves accordingly.
 * 
 * @param halt An integer reference representing the current time step or index 
 *             for which the hydraulic conductivity and pressure values are 
 *             being updated.
 * 
 * @details
 * - The function iterates through all soil layers and checks the root 
 *   component's hydraulic conductivity. If the conductivity is below the 
 *   minimum threshold, it updates the Kmin value and recalculates the E(P) and 
 *   K(P) curves for the root.
 * - Similar checks and updates are performed for the stem and leaf components 
 *   of the plant.
 * - Numerical instabilities are accounted for by ensuring that the difference 
 *   between Kmin and the current conductivity is greater than a small threshold 
 *   (1e-9) before making updates.
 * - The pressure increment (pinc) is used to calculate the pressure datum and 
 *   back-calculate the E(P) values.
 * 
 */
void Plant::updatecurves(const int &halt) //'resets element E(P) curves
{

    double pinc = param.getModelParam("p_inc");
    int phigh = 0;

    //'if k<kmin, re-assign e//'s on element curve by back-calculating
    for (int z = 1; z <= layers; z++)//z = 1 To layers
    {
        RootComponent *root = &xylem.soils[z]->root;
        if (true)
        {
            // if (root->getKComp(halt) < root->getKmin())
            if ((root->getKmin() - root->getKComp(halt)) > 1e-9) // numerical instabilities means we should check if greater than some threshold value
            {
                root->setKmin(root->getKComp(halt));
                phigh = int(xylem.root_pressure[halt] / pinc) + 1; //'pressure datum just above the target
                for (int k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
                {
                    root->setEp(k, root->getEp(phigh) - root->getKmin() * pinc * (phigh - k));
                    root->setK(k, root->getKmin()); //'back-calculate KR(Z,K) too for roots (not stem or leaves)
                } //EndFor  k
            } //EndIf//
        }
    } //EndFor  z
    if (xylem.stem.getKComp(halt) < xylem.stem.getKmin())
    {
        xylem.stem.setKmin(xylem.stem.getKComp(halt));
        phigh = int(xylem.stem.getPressureComp(halt) / pinc) + 1;
        for (int k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
        {
            xylem.stem.setEp(k, xylem.stem.getEp(phigh) - xylem.stem.getKmin() * pinc * (phigh - k));
        } //EndFor  k
    } //EndIf//
    if (xylem.leaf.getKComp(halt) < xylem.leaf.getKmin())
    {
        xylem.leaf.setKmin(xylem.leaf.getKComp(halt));
        phigh = int(xylem.leaf.getPressureComp(halt) / pinc) + 1;
        for (int k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
        {
            xylem.leaf.setEp(k, xylem.leaf.getEp(phigh) - xylem.leaf.getKmin() * pinc * (phigh - k));
        } //EndFor  k
    } //EndIf//
        //'if kplant[halt] < kminplant Then kminplant = kplant[halt]NOTE: kplant CAN go up because of rhizosphere recovery!
}

void Plant::soilflow() //'gets flow out of each layer via soil, not including the groundwater basement
{
    double p1, p2, pend, e, s = 0;
    double pinc = param.getModelParam("p_inc");

    for (int z = 0; z < layers; z++)//z = 0 To layers - 1
    {
        SoilLayer *sl = xylem.soils[z];
        SoilLayer *next_sl = xylem.soils[z + 1];
        
        double store = sl->rhizosphere.getKmax(); //'store the rhizosphere kmax for now
        sl->rhizosphere.setKmax(next_sl->kkmax * 1 / (sl->depth / 2.0 + next_sl->depth / 2.0)); //'reset it to the vertical soil kmax per m2 ground using distance between midpoints of each layer. Use properties of fatter layer.
        if (sl->predawn_pressure == next_sl->predawn_pressure)
            sl->soilredist = 0; //'no redistribution (gravity ignored!)
        if (sl->predawn_pressure != next_sl->predawn_pressure) { //if//
            if (sl->predawn_pressure < next_sl->predawn_pressure) { //if// //'flow is out of layer z (positive)
                p1 = sl->predawn_pressure;
                pend = next_sl->predawn_pressure;
            }
            else { //'flow is into layer z(negative)
                p1 = next_sl->predawn_pressure;
                pend = sl->predawn_pressure;
            }//
            
            e = 0; //'flow integral
            do
            {
                p2 = p1 + pinc;
                sl->rhizosphere.qtrap(p1, p2, s);
                e = e + s;
                p1 = p2; //'reset p1 for next increment
            } while (!(p2 >= pend || p2 > 5));
            //Loop Until p2 >= pend Or p2 > 5 //'cutoff for really dry soil
            if (sl->predawn_pressure < next_sl->predawn_pressure) { //if// //'flow is out of layer z (positive)
                sl->flow = e; //'flow in kg hr-1 m-2 ground area
            }
            else { //'flow is into layer z(negative)
                sl->flow = -e;
            }//
        }//
        sl->rhizosphere.setKmax(store);
    } //for//z
    
    xylem.soils[0]->soilredist = xylem.soils[0]->flow; //'set upper layer which has only one flux
    xylem.soils[layers]->soilredist = -1 * xylem.soils[layers - 1]->flow; //'set water flowing into/out of top of the bottom layer...soilredist is net soil flow of layer
                                                //'now calculate net flows of internal layers
    for (int z = 1; z < layers; z++)//z = 1 To layers - 1
    {
        xylem.soils[z]->soilredist = (-1 * xylem.soils[z - 1]->flow) + xylem.soils[z]->flow; //'add up to get net flow...remember negative flow is flow into the layer
    } //for//z
}

/**
 * @brief Calculates the flow of water into or out of the lowermost soil layer.
 *
 * This function determines the net flow of water into or out of the bottom soil layer
 * (groundwater flow or drainage) based on the predawn pressure and groundwater pressure.
 * It integrates the flow over a pressure range and updates the respective flow variables.
 *
 * @param timestep The time step for the simulation, used to scale the flow values.
 *
 * The function performs the following steps:
 * - Stores the current rhizosphere maximum conductivity (Kmax).
 * - Resets the rhizosphere Kmax based on the vertical soil conductivity and distance to groundwater.
 * - Determines whether groundwater flows into the bottom layer (negative flow) or out of it (positive flow).
 * - Integrates the flow over the pressure range using small increments.
 * - Updates the groundwater flow or drainage totals in mm per square meter of ground area per time step.
 * - Restores the original rhizosphere Kmax.
 * - Updates the net flow into the lowermost soil layer.
 *
 * Variables:
 * - `gwflow`: Tracks the total groundwater flow into the bottom layer.
 * - `drainage`: Tracks the total drainage out of the bottom layer.
 * - `soilredist`: Tracks the net flow into the lowermost soil layer.
 */
void Plant::deepflow(const double &timestep) //'gets flow into/out of the lowermost layer
{
    int z = layers; //'just the bottom layer
    SoilLayer *sl = xylem.soils[z];
    double store = sl->rhizosphere.getKmax(); //'store the rhizosphere kmax for now

    double p1, p2, pend, e, groundflow, s = 0;
    double pinc = param.getModelParam("p_inc");

    sl->rhizosphere.setKmax(sl->kkmax / param.getModelParam("ground_distance")); //'reset it to the vertical soil kmax using distance to groundwater
    if (sl->predawn_pressure >= param.getModelParam("p_ground")) { //if// //'groundwater will flow in (negative flow)
        p1 = param.getModelParam("p_ground");
        pend = sl->predawn_pressure;
        e = 0; //'flow integral
        do
        {
            p2 = p1 + pinc;
            sl->rhizosphere.qtrap(p1, p2, s);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
        } while (!(p2 >= pend));
        //Loop Until p2 >= pend
        groundflow = -1 * e; //'negative flow (into bottom layer)
        gwflow = -1 * groundflow * 1 / 998.2 * timestep * 1000; //'groundwater flow totals in mm m-2 ground area per time step
    }
    else { //'flow is out of bottom layer(positive)
        p1 = sl->predawn_pressure;
        pend = param.getModelParam("p_ground");
        e = 0;//'flow integral
        do
        {
            p2 = p1 + pinc;
            sl->rhizosphere.qtrap(p1, p2, s);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
        } while (!(p2 >= pend));
        //Loop Until p2 >= pend
        groundflow = e; //'positive flow (out of bottom layer)
        drainage = drainage + groundflow * 1 / 998.2 * timestep * 1000; //'drainage totals in mm m-2 ground area per time step
    }//
    //'reset rhizosphere kmaxrh (though i don//'t THINK it//'s used again?)
    sl->rhizosphere.setKmax(store);
    //'now calculate net flow into lowermost layer
    sl->soilredist = sl->soilredist + groundflow; //'add up to get net flow...remember negative flow is flow into the layer
}

/* Should move this into the soils.h file */
/**
 * @brief Calculates soil evaporation based on environmental and soil parameters.
 *
 * This function computes the potential and actual soil evaporation rates using
 * various environmental inputs such as soil temperature, vapor pressure deficit,
 * air temperature, and wind speed. It also adjusts the evaporation rate based on
 * relative humidity, soil surface water potential, and basal area reduction.
 *
 * @param soiltemp The soil temperature in degrees Celsius.
 * @param maxvpd The maximum vapor pressure deficit (VPD) in kPa.
 * @param vpd The current vapor pressure deficit (VPD) in kPa.
 * @param airtemp The air temperature in degrees Celsius.
 * @param us The wind speed in m/s.
 *
 * @details
 * - The function calculates the long-wave emissivity of the soil and the thermal
 *   emissivity from the bottom of the canopy.
 * - It computes the aerodynamic conductance for heat and vapor transfer.
 * - The potential soil evaporation rate is determined based on energy balance.
 * - The actual soil evaporation rate is adjusted using relative humidity and
 *   soil surface water potential.
 * - Negative evaporation rates are set to zero.
 * - The evaporation rate is reduced based on basal area and converted to
 *   kilograms per square meter per hour.
 * - The evaporative loss is added to the soil redistribution for the top soil layer.
 */
void Plant::soilevaporation(const double &soiltemp,
                            const double &maxvpd,
                            const double &vpd,
                            const double &airtemp,
                            const double &us) //'get soil evaporation
{
    double emission = param.getModelParam("emiss") * SBC * pow((soiltemp + 273.15), 4); //'long wave emissivity of soil
    double lc = 0.97 * SBC * pow((xylem.leaf.leaftempsh[haltsh] + 273.15), 4); //'thermal emissivity from bottom of canopy
    double rabssoil = param.getModelParam("soil_abs_sol") * xylem.leaf.sbottom + 0.97 * lc; //'rabs of soil surface; does not account for slope (beam only)
    double mdensair = 44.6 * (param.getModelParam("p_atm") / 101.3) * 273.15 / (airtemp + 273.15); //'molar density of air, CN 3.3, moles/m3
    double gha = 0.16 * mdensair * us / (log((param.getModelParam("xh") - xylem.zdispl) / xylem.rough) * log((param.getModelParam("xh") - xylem.zdispl) / xylem.zh)); //'CN 7.28/14.9,heat and vapor aerodynamic conductance for soil, mol m-2 s-1, stability factors = 0
    double soilep = (rabssoil - emission - SHA * (soiltemp - airtemp) * (xylem.leaf.grad + gha)) / xylem.leaf.lambda; //'potential evaporation in moles m-2 s-1; 12.8

    if (soilep < 0)
        soilep = 0; //'set to zero if negative
    double rha = 1 - vpd / maxvpd; //'relative humidity of air
    double rhs = exp((-xylem.soils[0]->predawn_pressure * 0.00001805) / (GAS * 0.000001 * (soiltemp + 273.15))); //'relative humidity in equilibrium with soil surface water potential

    if (rha == 1) { //if//
        soilevap = 0;
    }
    else {
        soilevap = soilep * (rhs - rha) / (1 - rha); //'campbell 1985, eqn 9.14; actual soil evaporation rate in moles m-2s-1
    }//
    if (soilevap < 0)
        soilevap = 0; //'don//'t let it go negative
    soilevap = soilevap * (1 - param.getModelParam("ba_per_ga")); //'reduce for basal area
    soilevap = soilevap * 0.0180153 * 3600; //'converts from moles m-2s-1 to kg m-2 hr-1
    xylem.soils[0]->soilredist = xylem.soils[0]->soilredist + soilevap; //'add evaporative loss to redistribution for top layer
}

/**
 * @brief Stores historical element e(P) curves and updates soil states.
 * 
 * This function iterates through all soil layers and performs the following:
 * - Stores historical curves for the root system in each soil layer and switches to virgin curves.
 * - Saves the current cavitation and failure states of each soil layer.
 * - Resets the cavitation state for each soil layer.
 * - Checks if the failure state is due to the root and if the minimum hydraulic conductivity (Kmin) is zero,
 *   in which case the cavitation state is set to true.
 * 
 * Additionally, it stores the transpiration curves for the stem and leaf components and switches them to virgin curves.
 */
void Plant::storehistory() //'stores historical element e(P) curves
{
    for (int z = 1; z <= layers; z++)
    {
        xylem.soils[z]->root.storeCurvesAndUseVirgin();
        xylem.soils[z]->cavitated_t = xylem.soils[z]->cavitated;
        xylem.soils[z]->failure_t = xylem.soils[z]->failure;

        xylem.soils[z]->cavitated = false;
        if (xylem.soils[z]->failure == "root" && xylem.soils[z]->root.getKmin() == 0)
            xylem.soils[z]->cavitated = true;
}

    xylem.stem.storeTranspirationCurveAndUseVirgin();
    xylem.leaf.storeTranspirationCurveAndUseVirgin();
}

/**
 * @brief Restores the historical element e(P) curves for the plant.
 * 
 * This function iterates through all soil layers associated with the plant's xylem
 * and restores their historical curves. It also updates the cavitation and failure
 * states of the soil layers to their historical values. Additionally, it restores
 * the transpiration curves for the stem and leaf components of the xylem.
 * 
 * @details
 * - For each soil layer in the xylem:
 *   - Restores the historical root curves using `restoreCurves()`.
 *   - Updates the `cavitated` state to the historical `cavitated_t` value.
 *   - Updates the `failure` state to the historical `failure_t` value.
 * - Restores the transpiration curves for the stem and leaf using their respective
 *   `restoreTranspirationCurve()` methods.
 */
void Plant::gethistory() //'restores historical element e(P) curves
{
    for (int z = 1; z <= layers; z++)
    {
        xylem.soils[z]->root.restoreCurves();
        xylem.soils[z]->cavitated = xylem.soils[z]->cavitated_t;
        xylem.soils[z]->failure = xylem.soils[z]->failure_t;
    }
    
    xylem.stem.restoreTranspirationCurve();
    xylem.leaf.restoreTranspirationCurve();
}