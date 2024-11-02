#include "09Plant.h"

Plant::Plant(int stage) {
    stage_id = stage;
    std::cout.precision(FIO_PRECISION);
}

void Plant::setConfig() { // sets up model configuration
    // // Model Configurations
    // Plant Community
    std::cout << "Plant Community ----------------------------------" << std::endl;
    std::cout << "            Multi-Species Mode: ";
    // Are we working with multiple species? Values: y; n
    if (config_data.getColumnValue("i_multipleSP") == "y")
    {
        std::cout << "On" << std::endl;
        config_data.getColumnValue(species_no, "i_speciesN");
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
    config_data.getColumnValue(ground, "i_gWaterEnable");
    if (ground){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    std::cout << "     Soil water redistribution: "; //i_soilRedEnable
    // turns on/off soil redistribution routine. Values: on: y; off: n
    config_data.getColumnValue(soilred, "i_soilRedEnable");
    if(soilred){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    std::cout << "        Soil water evaporation: "; //i_soilEvapEnable
    // turns on/off soil evaporation routine. Values: on: y; off: n
    config_data.getColumnValue(sevap, "i_soilEvapEnable");
    if(sevap){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }

    // // Climate
    std::cout << std::endl;
    std::cout << "Climate ------------------------------------------" << std::endl;
    std::cout << "               Rainfall inputs: "; //i_rainEnable
    config_data.getColumnValue(raining, "i_rainEnable");
    if (!raining) {
        std::cout << "off" << std::endl;
    }
    else {
        //enabled is made the "else" case becasue the "default" state is to process rain, so any input OTHER than //n// enables
        std::cout << "on" << std::endl;
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
        if (config_data.getColumnValue("i_useGSDataOpt") == "y") { // different parameter name
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
    if (config_data.getColumnValue("i_predawnsMode") == "y")
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
    config_data.getColumnValue(refilling, "i_refilling");
    if(refilling){
        std::cout << "On" << std::endl;
    } else {
        std::cout << "Off" << std::endl;
    }
    
    std::cout << "              Xylem hysteresis: "; 
    // turns on/off xylem hysteresis from previous growing season. Values: n(off); y(on) 
    if(config_data.getColumnValue("i_cavitFatigue") == "y"){// "i_cavitFatigue"
        hysteresis = true;
        std::cout << "On" << std::endl;
    } else {
        hysteresis = false;// default
        std::cout << "Off" << std::endl;
    }

    std::cout << "   Cavitation fatigue in roots: ";
    if (config_data.getColumnValue("i_stemOnly") == "n")
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
    if (config_data.getColumnValue("i_iter_gwEnable") == "y")
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
    
    int sp_n = 0;
    param_data.getColumnValue(iter_gwInc, "i_iter_gwInc", sp_n);
    param_data.getColumnValue(iter_gwStart, "i_iter_gwStart", sp_n);
    param_data.getColumnValue(iter_gwEnd, "i_iter_gwEnd", sp_n);

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

    param_data.getColumnValue(iter_ffcInc, "i_iter_ffcInc",sp_n);
    param_data.getColumnValue(iter_ffcStart, "i_iter_ffcStart",sp_n);
    param_data.getColumnValue(iter_ffcEnd, "i_iter_ffcEnd",sp_n);
    
    std::cout << "           Iterate Field BA:GA: ";
    // Iterating BA:GA for each stand. Values: off; on
    if (config_data.getColumnValue("i_iter_bagaEnable") == "y")
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

    param_data.getColumnValue(iter_bagaInc, "i_iter_bagaInc",sp_n);
    param_data.getColumnValue(iter_bagaStart, "i_iter_bagaStart",sp_n);
    param_data.getColumnValue(iter_bagaEnd, "i_iter_bagaEnd",sp_n);
    param_data.getColumnValue(iter_bagaRef, "i_iter_bagaRef",sp_n);
    param_data.getColumnValue(iter_bagaCutoff, "i_iter_bagaCutoff",sp_n);
    
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
    climate_forcing_data_path = config_data.getColumnValue("i_ClimateData");
    std::cout << climate_forcing_data_path << std::endl;
    
    std::cout << "   Path to growing season data: " << std::endl;
    growing_season_limits_data_path = config_data.getColumnValue("i_GSData");
    std::cout << growing_season_limits_data_path << std::endl;

    std::cout << "  Path to time-step header file: " << std::endl;
    data_header_file_path = config_data.getColumnValue("i_dataheader");
    std::cout << data_header_file_path << std::endl;
    
    std::cout << "   Path to annual summary header file: " << std::endl;
    sum_header_file_path = config_data.getColumnValue("i_sumheader");
    std::cout << sum_header_file_path << std::endl;
}

void Plant::initModelVars() {
    // this can be called on the start of every iteration BUT NOT ON NEW YEARS
    isNewYear = true;
    // this is a good place for general initialization too
    gs_yearIndex = 0;
    gs_prevDay = 0;
    gs_inGrowSeason = false;
    gs_doneFirstDay = false;
    year_cur = 0;
    year_start = 0;
    yearVal = 0;
}

void Plant::cleanModelVars() {
    max_plc_x = 0;
}

/* Get Van Genuchten alpha // override if provided */
// This function obtains the Van Genuchten parameters for a given texture
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
    } else {
        std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
        std::cout << "SOURCE: Incorrect soil texture category" << std::endl;
        std::cout << "ACTION: Model stops " << std::endl;
        std::cout << std::endl;
        abort();
    }
    for (int k = 0; k < layers; k++){
        soils[k]->rhizosphere.setVanGenAlpha(a);
        soils[k]->rhizosphere.setVanGenN(n);
        soils[k]->kkmax = soilkmax;
        soils[k]->rhizosphere.setThetasat(thetasat);
    }
}

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
    temp = temp * 0.000001; // ambient co2 in moles per mole
    param.setModelParam(temp, "ca");
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

    for (int k = 0; k < layers; k++) {
        SoilLayer *t = new SoilLayer; // Xylem destructor deletes the allocated SoilLayers
        t->root.setCwb(c_temp);
        t->root.setBwb(b_temp);
        xylem.top_soil.root.setBwb(b_temp);
        xylem.top_soil.root.setCwb(c_temp);
        xylem.soils.push_back(t);

        for (int j = 0; j < sizeof(b_fatigue[0]) / sizeof(b_fatigue[0][0]); j++) {
            xylem.soils[k]->root.setFatigue(b_fatigue[0][j], j);
            xylem.top_soil.root.setFatigue(b_fatigue[0][j], j);
        }
    }
    get_vgparams(param_data.getColumnValue("i_texture", species_no), layers, xylem.soils); // extract VG parameters for a given texture

    // STEMS
    param_data.getColumnValue(c_temp, "i_cs",species_no);
    param_data.getColumnValue(b_temp, "i_bs",species_no);
    std::cout << "Weibull parameter before accounting for previous year PLC: " << std::endl;
    std::cout << "stem b: " << b_temp << " stem c: " << c_temp << std::endl;
    if (hysteresis == true)
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

        temp = param.getModelParam("lsc_input") / (std::exp(-(std::pow((param.getModelParam("lsc_bref") / b_temp), c_temp))));//--//
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
    param.setModelParam(temp, "hav_max");
    param_data.getColumnValue(temp, "i_hdvmax",species_no); //
    param.setModelParam(temp, "hdv_max");
    param_data.getColumnValue(temp, "i_svvmax",species_no); //
    param.setModelParam(temp, "sv_vmax");
    param_data.getColumnValue(temp, "i_hajmax",species_no); //
    param.setModelParam(temp, "haj_max");
    param_data.getColumnValue(temp, "i_hdjmax",species_no); //
    param.setModelParam(temp, "hd_jmax");
    param_data.getColumnValue(temp, "i_svjmax",species_no); //
    param.setModelParam(temp, "svj_max");
    param_data.getColumnValue(temp, "i_lightCurv",species_no); //
    param.setModelParam(temp, "light_curv");

    // initial conditions calculations ---------------------------------------------------------------------------------------------------------
    temp = 1000000;
    param.setModelParam(temp, "g_max");
    temp = temp * (1.0 / param.getModelParam("leaf_per_basal")) * (1 / 3600.0) * 55.56 * 1000.0; //convert to gmax per leaf area in mmol m-2s-1
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
        for (int k = 0; k < layers; k++) {
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
        xylem.stem.setKmax(1.0 / rSatStem);

        // any chance this worked? (later me: It did, once I remembered the unit conversion for lscMax)
        std::cout << "Conductance calculations: kMaxTree = " << param.getModelParam("ksatp") << " lsc = " << param.getModelParam("lsc") << " kMaxStem = " << xylem.stem.getKmax() << " kMaxRoot = " << param.getModelParam("k_sat_root") << std::endl;
    }
    temp = param.getModelParam("height") * 0.01; //pressure drop from gravity in MPa
    param.setModelParam(temp, "p_grav");
    temp = param.getModelParam("ksatp") / 500.0; // e increment in kg hr-1 m-2 basal area for composite curve
    param.setModelParam(temp, "e_inc");
    temp = param.getModelParam("ksatp") / 2000.0; //"instantaneous K" cutoff for global K(P) curves for each element
    param.setModelParam(temp, "k_min");

    xylem.num_layers = layers;
    //for this soil data we want to use the original anchor-offset system
    for (int k = 0; k < layers; k++) { // set layer depths and % root ksat
        xylem.soils[k]->layer_depth = 0.01 * log(1.0 - (k + 1) * 0.995 / layers) / log(param.getModelParam("beta")); // lower depth of each layer converted to m
    }
    double depthmax = xylem.soils[layers - 1]->layer_depth;
    // calculate transport distances
    std::vector<double> temp_vert_distance;
    // first get vertical distance to biomass center of each layer
    for (int k = 0; k < layers * 2.0; k++)
    {
        temp_vert_distance.push_back(0.01 * log(1 - (k + 1) * 0.995 / (layers * 2.0)) / log(param.getModelParam("beta"))); // get half depths
    }
    int i = 0;
    for (int k = 0; k < layers * 2; k += 2) // To layers * 2 Step 2
    {
        xylem.soils[i]->vert_distance = temp_vert_distance[k]; //take every other vertdistance
        i = i + 1;
    }
    // now get radial distances
    for (int k = 0; k < layers; k++) //To layers //get thicknesses
    {
        if (k == 0)
        {
            xylem.soils[k]->depth = xylem.soils[k]->layer_depth; //Cells(8 + k, 11)
        }
        else {
            xylem.soils[k]->depth = xylem.soils[k]->layer_depth - xylem.soils[k - 1]->layer_depth; // depth is layer thickness in meters
        } // endif
    }
    temp = xylem.soils[0]->depth * PI * pow((depthmax * param.getModelParam("aspect")), 2.0); //volume of first, and hence all, layers
    param.setModelParam(temp, "vol");
    double shallow = 0;
    // get radial widths of each layer and transport length
    for (int k = 0; k < layers; k++) {
        xylem.soils[k]->radius = pow(temp / (xylem.soils[k]->depth * PI), 0.5); // width in m
        xylem.soils[k]->length = xylem.soils[k]->radius + xylem.soils[k]->vert_distance; // transport distance
        if (k == 0)
            shallow = xylem.soils[k]->length;
        xylem.soils[k]->length = xylem.soils[k]->length / shallow;    
    }

    unknowns = layers + 1; // number of unknowns to be solved for; also dimensions of matrix
    temp = 1 - param.getModelParam("rock_frac"); // fraction of volume with no rocks
    for (int k = 0; k < layers; k++) { // read in soil properties
        xylem.soils[k]->rhizosphere.setThetasat(xylem.soils[k]->rhizosphere.getThetaSat() * temp); // reduce for actual rock-free fraction of soil
    }
    /* Now add toplayer (layer 0) of rootless soil 2 cm thick w. same properties as layer 1 */
    xylem.top_soil.rhizosphere.setVanGenAlpha(xylem.soils[0]->rhizosphere.getVanGenAlpha());
    xylem.top_soil.rhizosphere.setVanGenN(xylem.soils[0]->rhizosphere.getVanGenN());
    xylem.top_soil.kkmax = xylem.soils[0]->kkmax;
    xylem.top_soil.rhizosphere.setThetasat(xylem.soils[0]->rhizosphere.getThetaSat());
    xylem.top_soil.depth = 0.02; //sets top layer to 2 cm

    /* Now solve for kmax rhizosphere that gives the desired ave % rhizosphere resistance */
    int z = 0; //use layer 1 as stand in for whole root system
    xylem.soils[0]->root.setKmax(param.getModelParam("k_sat_root")); //set to whole root system
    double x = 0.5; //start by finding kmaxrh at 0.5 MPa...a deliberate under-shoot
    double rootr, rstem, rleaf, rplant, rhizor, vp, vgterm, kinc, sum, rrhizofrac;
    rootr = 1.0 / xylem.soils[0]->root.wb(x);
    rstem = 1.0 / xylem.stem.wb(x);
    rleaf = 1.0 / xylem.leaf.wb(x);
    rplant = rootr + rstem + rleaf; //rplant here is just the xylem part
    rhizor = rplant * (param.getModelParam("rhizo_targ") / (1.0 - param.getModelParam("rhizo_targ"))); //solve for what rhizor has to be at the target
    vp = 1.0 / (pow((xylem.soils[z]->rhizosphere.getVanGenAlpha() * x), xylem.soils[z]->rhizosphere.getVanGenN()) + 1); //van genuchten terms // vp = 1 / ((a(z) * x) ^ n(z) + 1) 
    vgterm = pow(vp, ((xylem.soils[z]->rhizosphere.getVanGenN() - 1) / (2.0 * xylem.soils[z]->rhizosphere.getVanGenN()))) * pow((pow((1 - vp), ((xylem.soils[z]->rhizosphere.getVanGenN() - 1) / xylem.soils[z]->rhizosphere.getVanGenN())) - 1), 2.0); //van genuchten terms // vgterm = vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
    xylem.soils[0]->rhizosphere.setKmax((1.0 / rhizor) / vgterm); //solve for kmaxrh[1]
    kinc = xylem.soils[0]->rhizosphere.getKmax() * 0.1;
    do //loop through rhizosphere kmax
    {
        xylem.soils[0]->rhizosphere.setKmax(xylem.soils[0]->rhizosphere.getKmax() + kinc); //increase from deliberate undershoot
        x = 0;
        sum = 0;
        do //loop through pressures
        {
            x = x + 0.1;
            rootr = 1.0 / xylem.soils[0]->root.wb(x);
            rstem = 1.0 / xylem.stem.wb(x);
            rleaf = 1.0 / xylem.leaf.wb(x);
            rhizor = 1.0 / xylem.soils[0]->rhizosphere.vg(x);
            rplant = rootr + rstem + rleaf + rhizor;
            rrhizofrac = rhizor / rplant; //fraction of resistance in rhizosphere
            sum = sum + rrhizofrac; //add up fractions
        } while (!((1.0 / rplant) < param.getModelParam("k_min"))); //Loop Until 1 / rplant < kmin //average over full range
        sum = sum / (x / 0.1); //average fraction
    } while (!(sum < param.getModelParam("rhizo_targ"))); // Until sum < rhizotarg //loop until desired soil limitation is reached
    
    xylem.soils[0]->rhizosphere.setKmax(xylem.soils[0]->rhizosphere.getKmax() / layers);  //divide whole root rhizokmax into equal portions for each layer
    // end of soil limitation adjustment
    // now set soil layer parameters based on aroot of entire root system
    for (int k = 0; k < layers; k++)
    {
        xylem.soils[k]->rhizosphere.setKmax(xylem.soils[0]->rhizosphere.getKmax()); // soil to root MAXIMUM conductance in kg hr-1 MPa-1//re - set for individual layers
                                                                                    // note: this is kMAX...at P=0; not at saturated PD
    }

    // t = 0;
    //loop to find ksatroot for each layer
    double coef = 0.0;
    do
    {
        coef = coef + 0.01;
        sum = 0.0;
        for (int k = 0; k < layers; k++)//k = 0 To layers //soil layers from top to bottom
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

void Plant::componentPCrits() {
    xylem.calc_net_flow(param.getModelParam("p_inc"), param.getModelParam("k_min"));
}

int Plant::modelTimestepIter(int &dd) {
    double transpiration, md, psynact, gcmd, lavpdmd, cinc, lai, laisl, laperba;
    double transpirationsh, mdsh, psynactsh, gcmdsh, lavpdmdsh, cincsh, qsh, laish;
    std::vector<double> eplantl, pleaf, psyn, gcanw, lavpd, cin, leaftemp;
    std::vector<double> psynsh, gcanwsh, lavpdsh, cinsh, leaftempsh;

    bool failure = false;
    int check = 0,
        weird,
        total;
    std::string failspot, // used for debugging?
                night;
    std::vector<double> kminroot;
    int    tod, // time of today
           timestep;
    double kminstem,
           kminleaf,
           kminplant,
           gwflow,
           drainage,
           ca,  // carbon assimilatiomn
           obssolar,
           vpd,
           airtemp,
           maxvpd,
           wind,
           us,
           soiltemp,
           chalk,
           qsl,
           lightcomp,
           transpiration_tree,
           atree,
           kplantold,
           patm;

    if (dd == 0 || isNewYear)
    {
        failure = 0;//'=1 for system failure at midday...terminates run
        failspot = "no failure";
        this->componentPCrits();//'gets pcrits for each component
        failspot = "no failure";

        for (int i = 0; i < layers; i++)// k = 1 To layers ;//'exclude the top layer
        {
            kminroot.push_back(xylem.soils[i]->root.getKmax());
        }

        kminstem = xylem.stem.getKmax();
        kminleaf = xylem.leaf.getKmax();
        kminplant = param.getModelParam("ksatp");

        gwflow = 0; //'inflow to bottom of root zone
        drainage = 0; //'drainage from bottom of root zone

        gs_prevDay = 0;
        gs_inGrowSeason = false;
        gs_doneFirstDay = false; // prevents PLC misses on first year
                                //[/HNT]
    }

    //dd++;
    long testCount = 0;
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
            return 0; // failure
        }
    }
    
    if ( dd == 0 || isNewYear){// Get CO2 for current year
        ca = getCarbonByYear(yearVal); // get current year atmospheric CO2
        std::cout << "Atmospheric CO2 concentration for " << yearVal << ": " << ca << std::endl;
        ca = ca * 0.000001;
        param.setStemBWb(gs_yearIndex, xylem.stem.getBwb());
        param.setRootBWb(gs_yearIndex, xylem.soils[0]->root.getBwb());
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
        } //End if// //'tod if
    } //End if// //'dd>1 if
        //[HNT] multi-year support
    else
    {
        gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's the first data point, test if we're starting in growing season
    }

    gs_inGrowSeason = isInGrowSeasonSimple(jd); // just always call this!
                                                //[/HNT]

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
    } //End if//
    us = wind * 0.1; //'understory windspeed in m s-1
    soiltemp = data.getColumnValue("T-soil", dd); //'surface temp of soil
    if (vpd > maxvpd) { //if//
        vpd = maxvpd;
        data.setColumnValue(maxvpd * param.getModelParam("p_atm"), dd, "D-MD"); //'print out maximum vpd
    } //End if//

    getsoilwetness(dd, timestep, tod, night, lai, laish, laisl, transpiration_tree, laperba, atree, cinc, ca); //'after initializing, start updating water contents of soil layers
//     solarcalc(); //'get radiation for timestep
//     if (qsl > lightcomp /*[HNT] multiyear*/ && gs_inGrowSeason /*[/HNT]*/) { //if//
//         night = "n"; //'it//'s light enough for sun layer to do business
//     }
//     else {
//         night = "y"; //'too dark
//     } //End if// //'end night if

//     gwflow = 0; //'re-set inflow to bottom of root zone
//     drainage = 0; //'re-set drainage from bottom of root zone
//                     //'Call getpredawns //'update soil pressure of each layer
//                     //'if failure = 1 { //if// Exit do
//     chalk = 0; //'einc counter

// twentyMarker:

//     getpredawns(); //'passed initializing...update soil pressure of each layer

//     if (failure == 1)
//         return -1;//break;//Exit do

//     int test = 0; //'=1 if stem or leaf fails
//     int p = -1; //'E(Pleaf) point counter
//     double einc = param.getModelParam("e_inc");
//     double e = -einc; //'total e: einc and e are still in kg hr-1 m-2 basal area, as are conductances
//     double psynmax = -100;
//     double psynmaxsh = -100;
//     int skip = 0; //'this turns off psynthesis
//     double sum = 0;

//     do //'this loop obtains and stores the entire composite curve
//     {
//         p = p + 1;
//         e = e + einc;
//         newtonrhapson(); //'solves for p//'s and e//'s in fingers of chicken foot
//                         //'if check >= 400 { //if// GoTo 60: //'NR can//'t find a solution
//         if (check > 500) { //if//
//             break;//Exit do
//         } //End if//
//         //'test for total failure
//         sum = 0;
//         for (int k = 0; k < layers; k++)//k = 1 To layers
//         {
//             sum = sum + layer[k];
//         } //end for k
//         if (sum == layers) { //if//
//             failspot = "below ground"; //'total failure below ground
//             break;//Exit do
//         } //End if//
//         stem(); //'gets stem and leaf pressures
//         leaf();
//         if (test == 1)
//             break;//Exit do
//         compositecurve(); //'stores the entire composite curve
//                         //'if skip = 0 { //if//
//         leaftemps(); //'gets sun layer leaf temperature from energy balance
//         leaftempsshade(); //'gets shade layer leaf temperature
//         assimilation(); //'gets sun layer photosynthesis
//         assimilationshade(); //'gets shade layer photosynthesis

//     } while (!(sum == layers || test == 1 || night == "y" && (dd > 1 && !isNewYear) || check >= 500)); //'loop to complete failure unless it//'s night

//     if (chalk > 0) { //if//
//         weird = 0; //'done our best
//         failspot = "convergence";
//         for (int z = 0; z < layers; z++)//z = 0 To layers //'restore layers to functioning if they//'ve been turned off by convergence failure
//         {
//             if (kminroot[z] != 0)
//                 layer[z] = 0;
//         } //end for z

//         goto fortyMarker; //'got as much of the composite curve as is going to happen
//     } //End if//
//     if (dd == 0 || isNewYear || night == "n") { //if//

//         if (check >= 500) { //if// //'try once more
//                             //'Stop
//             chalk = chalk + 1;

//             if (ecritsystem == 0)
//             {
//             einc = param.getModelParam("ksatp") / 500.0;
//             std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
//             }

//             goto twentyMarker;
//         } //End if//

//         if (total > 500 || total < 400) { //if//
//             einc = ecritsystem / 450.0; //'re-set Einc
//             if (ecritsystem == 0)
//             {
//             einc = param.getModelParam("ksatp") / 500.0;
//             std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
//             }
//             testCount++; // [DEBUG]
//             if (testCount > 10)
//             {
//             testCount = 0;
//             goto fortyMarker;
//             }
//             goto twentyMarker; //'recalculate the composite curve
//         } //End if// //'total ok

//     } //End if// //'night <>"n" or it//'s not the first round

// fortyMarker:

//     bool isNight = true;
//     if (night == "n")
//         isNight = false;
//     if (night == "n" && psynmax > 0 && psynmaxsh > 0 && weird == 0) { //if//
//                                                                         //DoEvents //'[HNT] this was required to prevent a hard lock -- this portion of the loop is the most intensive, so let Excel take a "breath" by processing system events to prevent lockup
//         canopypressure(); //'returns canopy P and associated output
//                         //'if check >= 2000 { //if// GoTo 60:
//         if (refilling == false)
//             updatecurves(); //'updates element E(P) curves as required for midday exposure for no refilling
//     } //End if// //'night <> "n", psynmax IF

//     if (soilred == true) { //if//
//         soilflow(); //'gets vertical soil flow between layers in m3/m2
//     }
//     else {
//         for (int z = 0; z < layers; z++)//z = 0 To layers
//         {
//             xylem.soils[z]->soilredist = 0;
//         } //end for z
//     } //End if// //'soil red <> y
//     if (ground == true)
//         deepflow(); //'gets groundwater flow into bottom layer in m3/m2
//                     //'} //End if// //'pet <> y or n
//     if (gs_inGrowSeason && sevap == true) { //if//
//         soilevaporation(); //'gets soil evaporation rate
//     }
//     else {
//         param.setModelParam(0, "soil_evap");
//     } //End if//

//     if (failure == 0 || weird == 1) { //if//
//                                         //Debug.Print "DOING A LOOP-8 " & dd

//         if (night == "y" || psynmax == 0 || psynmaxsh == 0) { //if// //'set everything to starting point
//             int k = 0;
//             transpiration = eplantl[k]; //'all gas exchange values are for closed stomata
//             md = pleaf[k];
//             psynact = psyn[k];
//             gcmd = gcanw[k]; //'g for water in mmol
//             lavpdmd = lavpd[k] * patm;
//             cinc = cin[k];
//             halt = k;
//             transpirationsh = eplantl[k]; //'all gas exchange values are from most recent historical values
//             mdsh = pleaf[k];
//             psynactsh = psynsh[k];
//             gcmdsh = gcanwsh[k]; //'g for water in mmol
//             lavpdmdsh = lavpdsh[k] * patm;
//             cincsh = cinsh[k];
//             haltsh = k; //'halt is index of midday datum
//         } //End if// //'night<>y

//         data.setColumnValue(pleaf[0], dd, "P-PD"); //'the predawn
//                                                                     //'SUN LAYER OUTPUT
//         data.setColumnValue(md, dd, "P-MD"); //'the midday
//         data.setColumnValue(transpiration, dd, "E-MD"); //'midday transpiration, mmol s-1 m-2 leaf area
//         data.setColumnValue(gcmd, dd, "GW"); //'midday canopy diffusive conductance to water, mmol s-1m-2
//         data.setColumnValue(lavpdmd, dd, "leaf-air-vpd"); //'leaf-to-air vpd
//         data.setColumnValue(leaftemp[halt], dd, "leaftemp"); //'leaf temp
//         data.setColumnValue(psynact, dd, "Anet-la"); //'net A in umol s-1m-2 leaf area
//         data.setColumnValue(cinc * patm * 1000, dd, "ci"); //'partial pressure of CO2 in Pa
//         data.setColumnValue(qsl, dd, "PPFD"); //'umol s-1m-2 photon flux density
//         data.setColumnValue(mdsh, dd, "S-P-MD"); //'the midday
//         data.setColumnValue(transpirationsh, dd, "S-E-MD"); //'midday transpiration, mmol s-1 m-2 leaf area
//         data.setColumnValue(gcmdsh, dd, "S-GW"); //'midday canopy diffusive conductance to water, mmol s-1m-2
//         data.setColumnValue(lavpdmdsh, dd, "S-leaf-air-vpd"); //'leaf-to-air vpd
//         data.setColumnValue(leaftempsh[haltsh], dd, "S-leaftemp"); //'leaf temp
//         data.setColumnValue(psynactsh, dd, "S-Anet-la"); //'A in umol s-1m-2 leaf area
//         data.setColumnValue(cincsh * patm * 1000, dd, "S-ci"); //'partial pressure of CO2 in Pa
//         data.setColumnValue(qsh, dd, "S-PPFD"); //'umol s-1m-2 photon flux density

//         if (night == "n")
//             transpiration_tree = laisl / lai * transpiration + laish / lai * transpirationsh; //'weighted mean
//         if (night == "y")
//             transpiration_tree = transpiration;
//         if (night == "n")
//             atree = laisl / lai * psynact + laish / lai * psynactsh; //'weighted mean
//         if (night == "y")
//             atree = (psynact + psynactsh) / 2.0; //'simple average at night when there//'s no sun or shade leaves

//         data.setColumnValue(transpiration_tree, dd, "S-E-tree"); //'weighted mean
//         data.setColumnValue(atree, dd, "Anet-tree"); // Anet Tree per Leaf Area (umol s-1m-2)
//         //'Cells(16 + dd, o + 35) = dpamax //'shade leaf dpa
//         //'HYDRAULIC OUTPUT (BASED ON SUN MD)
//         data.setColumnValue(pcritsystem, dd, "Pcrit");
//         data.setColumnValue(ecritsystem * (1 / laperba) * (1.0 / 3600.0) * 55.4 * 1000, dd, "Ecrit");
//         data.setColumnValue(xylem.stem.getPressure(halt), dd, "P-stem");
//         data.setColumnValue(xylem.root_pressure[halt], dd, "P-root");
//         data.setColumnValue(xylem.stem.getKComp(halt), dd, "K-stem");
//         data.setColumnValue(xylem.leaf.getKComp(halt), dd, "K-leaf");

//         if (transpiration > 0) { //if//
//             kplantold = this->e_p[halt] / (pleaf[halt] - pleaf[0]);  //'whole plant k at midday in kg hr-1 m-2 basal area...sun value
//             data.setColumnValue(kplantold, dd, "K-plant");
//         } //End if//
//         if (transpiration == 0)
//             data.setColumnValue(kplantold, dd, "K-plant"); //'use most recent kplant

//         if (kplantold < gs_data.getColumnValue("K-plant", gs_yearIndex) || gs_data.getColumnValue("K-plant", gs_yearIndex) == 0)
//             gs_data.setColumnValue(kplantold, gs_yearIndex, "K-plant");

//         // k = o + dColF_CP_kroot1 - 1;//43;
//         sum = 0;
//         for (int z = 0; z < layers; z++)//z = 1 To layers
//         {
//             std::ostringstream oss;
//             oss << "K-root-" << z;
//             data.setColumnValue(xylem.soils[z]->root.getKComp(halt), dd, oss.str()); //'root k at midday, sun fluxes
//             sum = sum + xylem.soils[z]->root.getKComp(halt);
//         } //end for z
//         data.setColumnValue(sum, dd, "K-root-all"); //'total root k at midday
//         if (failure == 0) { //if//
//             double tempDouble = 0.0;
//             tempDouble = 1 / (1 / kminleaf + 1 / kminstem + 1 / data.getColumnValue("K-root-all")); //'total xylem k
//             data.setColumnValue(tempDouble, dd, "K-xylem"); //1 / (1 / kminleaf + 1 / kminstem + 1 / dSheet.Cells(rowD + dd, colD + k + 1 + layers)); //'total xylem k
//             if (tempDouble < gs_data.getColumnValue("K-xylem", gs_yearIndex) || gs_data.getColumnValue("K-xylem", gs_yearIndex) == 0)
//                 gs_data.setColumnValue(tempDouble, gs_yearIndex, "K-xylem");
//         } //End if//
//         for (int z = 0; z < layers; z++)//z = 1 To layers
//         {
//             std::ostringstream oss;
//             oss << "E-root-" << z;
//             data.setColumnValue(xylem.soils[z]->root.getEcomp(halt) * (1 / laperba) * (1.0 / 3600.0) * 55.56 * 1000, dd, oss.str()); //'uptake in mmol s-1m-2 leaf area...sun rate
//         } //end for z
//         //'Cells(16 + dd, k + 2 + 2 * layers) = failspot //'position of failure at critical point

//         // TODO since all IO is being handled as double currently, cannot add these failure notes
//         // temporarily putting in an obvious number to flag failure
//         for (int z = 1; z <= 1; z++)//z = 1 To 1
//         {
//             std::ostringstream oss;
//             oss << "Layer-" << z << "-failure";
//             if (layer[z] == 1)
//                 data.setColumnValue(-1137, dd, oss.str()); // layer_failure[z]; //'layers failed at critical point
//         } //end for z

//         //Debug.Print "DOING A LOOP-9 " & dd
//     } //End if// //'failure IF (basically...failure can//'t happen!)


//     if (dd == 0 || isNewYear) { //if// //'NOTE: must be sure that pcritsystem is computed for dd=1!!! (i.e., it//'s not computed at night)
//         int x = pcritsystem; //'estimate of "permanent wilting point"
//         for (int z = 0; z < layers; z++)//z = 1 To layers
//         {
//             xylem.soils[z]->swclimit = swc(x); //'theta/thetasat at critical point
//             xylem.soils[z]->swclimit = xylem.soils[z]->swclimit * xylem.soils[z]->rhizosphere.getThetaSat(); //'convert to water content
//             xylem.soils[z]->swclimit = xylem.soils[z]->swclimit * xylem.soils[z]->depth; //'water content left over in m3/m2 ground area
//                                                 //'sumsoil = sumsoil + (fc[z] - xylem.soils[z]->swclimit) //'sum is total m3 water per m2 ground withdrawn
//         } //end for z

//         //Debug.Print "DOING A LOOP-11 " & dd
//     } //End if//

//         //'now...need to reset layer failure status to midday values (not critical point)
//     for (int z = 0; z < layers; z++)//z = 0 To layers
//     {
//         if (layer[z] == 1) { //if// //'check to see if kminroot[z]=0; otherwise, restore it to function
//             if (layer_failure[z] == "root" && kminroot[z] > 0) { //if//
//             layer[z] = 0;
//             layer_failure[z] = "ok";
//             } //End if//
//             if (layer_failure[z] == "rhizosphere" && kminroot[z] > 0) { //if//
//             layer[z] = 0;
//             layer_failure[z] = "ok";
//             } //End if//
//         } //End if//
//     } //end for z

    if (isNewYear)
        isNewYear = false; // always set this

    return -1;
}

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

// void Plant::modelProgramNewYear()
// {   
//     std::cout << std::endl;
//     std::cout << "Starting new year" << std::endl;
//     // save the iterate-able water system states
//     std::string oldGround;
//     std::string oldRaining;
//     bool oldRainEnabled;
//     oldGround = ground;
//     oldRaining = raining;
//     oldRainEnabled = rainEnabled;
    
//     // clear all the model variables
//     cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
//     maxPLCx = gs_ar_PLCx[gs_yearIndex-1];// previous year mean PLC?
//     std::cout << "Previous year Xylem PLC: " << maxPLCx << std::endl;
//     readin(); // get all the parameters again

//     gs_prevDay = 0;
//     gs_inGrowSeason = false;

//     // set up the iteration details based on what we saved
//     if (iter_runSupplyCurve) {
//         ground = "n"; //disable ground water for the first two iterations
//         raining = "n"; //disable rain for the first iteration
//         rainEnabled = false; //disable rain for the first iteration
//         ground = oldGround;
//         raining = oldRaining;
//         rainEnabled = oldRainEnabled;

//         if (iter_bagaEnable) {
//             param.getModelParam("ba_per_ga") = iter_baga * 0.0001;
//             //if we set the BA:GA { also set the LAI
//             lai = (laperba * param.getModelParam("ba_per_ga")) / treeToPhotoLAI;
//             std::cout << "DEBUG: set BA:GA to " << iter_baga << " m2/ha -> " << param.getModelParam("ba_per_ga") << "m2/m2. treeLAI/fotoLAI = " << treeToPhotoLAI << " and LAI set to " << lai << std::endl;
//         } //End if

//     } //End if

//     for (k = 0; k <= layers; k++) // To layers //assign source pressures, set layer participation
//     {
//         layer_failure[k] = "ok";
//         layer[k] = 0; //1 if out of function
//     } //Next k

//         //[HNT] all this is done in C now
//     if (true) { //Not useDLL {
//         failure = 0; //=1 for system failure at midday...terminates run
//         failspot = "no failure";
//         componentpcrits(); //gets pcrits for each component
//         failspot = "no failure";

//         for (k = 1; k <= layers; k++) {//k = 1 To layers //exclude the top layer
//             kminroot[k] = ksatr[k];
//         } //Next k

//         kminstem = ksats;
//         kminleaf = ksatl;
//         kminplant = ksatp;

//         gwflow = 0; //inflow to bottom of root zone
//         drainage = 0; //drainage from bottom of root zone
//     } //End if

//         //if kmaxset = False Or leafpercentset = False { dd = 0 //start with initializing row
//         //if kmaxset = True And leafpercentset = True { dd = 1 //skip initializing row
//         //dd = 0
//         //liveGraphNextUpdate = liveGraphUpdate

//         //cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
//         //readin();
//         //Call loadParamsToC
//         //Call CPP_setIterationCount(iter_Counter)
// }

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

double swc(const double &x, const double &van_gen_alpha, const double &van_gen_n) //'gives soil water content, theta/thetasat from soil water potential in MPa=x
{
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((van_gen_alpha * x), van_gen_n))), (1 - 1 / van_gen_n));
}

//'gets predawns accounting for rain, groundwater flow, transpiration, and redistribution via soil and roots
void Plant::getsoilwetness(const int &dd, 
                           const int &timestep,
                           const int &tod, 
                           const std::string &night,
                           const double &lai,
                           const double &laish,
                           const double &laisl,
                           const double &transpiration_tree,
                           const double &laperba,
                           const double &atree,
                           const double &cinc,
                           const double &ca) 
{
    double tempDouble = 0.0;
    double drainage = 0;
    double runoff = 0;
    double waterold,
           layerflow;

    std::vector<double> thetafracres(layers+1),
                        thetafracfc(layers+1),
                        thetafc(layers+1),
                        water(layers+1),
                        fc(layers+1);
    if (dd == 0 || isNewYear) { //if// //'every layer starts at initial % of field capacity

        if (!(useGSData && gs_yearIndex > 0))
        {
            waterold = 0; //'total root zone water content (mmol m-2 ground area)

            /* Set top layer */
            thetafracres[0] = swc(10, xylem.top_soil.rhizosphere.getVanGenAlpha(), xylem.top_soil.rhizosphere.getVanGenN()); //'residual thetafrac
            thetafracfc[0] = (1 - thetafracres[0]) * param.getModelParam("field_cap_frac") + thetafracres[0]; //'thetafrac at field capacity
            thetafc[0] = xylem.soils[0]->rhizosphere.getThetaSat() * thetafracfc[0]; //'water content at field capacity

            /* Set the rest */
            for (int z = 1; z <= layers; z++)//z = 0 To layers //'
            {
                // int x = 10; //'MPa water potential for getting residual thetafrac
                thetafracres[z] = swc(10, xylem.soils[z - 1]->rhizosphere.getVanGenAlpha(), xylem.soils[z - 1]->rhizosphere.getVanGenN()); //'residual thetafrac
                thetafracfc[z] = (1 - thetafracres[z]) * param.getModelParam("field_cap_frac") + thetafracres[z]; //'thetafrac at field capacity
                thetafc[z] = xylem.soils[z-1]->rhizosphere.getThetaSat() * thetafracfc[z]; //'water content at field capacity
            } //for//z //'

            /* Set top layer */
            water[0] = thetafc[0] * xylem.top_soil.depth; //'field capacity estimated as 1/2 saturated capacity, water content of layer in m3 water per m2 ground area
            fc[0] = water[0]; //'records field capacity in m3 water volume per m2 ground area.
            water[0] = param.getModelParam("ffc") * water[0]; //'start off with initial fraction of field capacity
                                        //'Cells(16 + dd, 42 + z) = water[z]
            waterold = waterold + water[0]; //'in m3/m2 ground area
            for (int z = 1; z <= layers; z++)//z = 1 To layers
            {
                water[z] = thetafc[z] * xylem.soils[z-1]->depth; //'field capacity estimated as 1/2 saturated capacity, water content of layer in m3 water per m2 ground area
                fc[z] = water[z]; //'records field capacity in m3 water volume per m2 ground area.
                water[z] = param.getModelParam("ffc") * water[z]; //'start off with initial fraction of field capacity
                                            //'Cells(16 + dd, 42 + z) = water[z]
                waterold = waterold + water[z]; //'in m3/m2 ground area
            } //for//z
                                                           
            data.setColumnValue(waterold * 1000, dd, "water-content");  //'root zone water content in mm m-2 ground area
                                                                        //'waterold = waterold * 55555556# //'convert m3 water per ground area to mmol water per ground area
                                                                        // [HNT] starting water now counts as an input
            gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "input");
        }

        if (gs_yearIndex == 0) // if it's the first year, store this as the off-season starting water because we don't have a real value
            gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "initial-off");

        // store the initial water content to check how much we consume at the end
        gs_data.setColumnValue(waterold * 1000, gs_yearIndex, "initial");
        // [/HNT]
    } //End if// //'dd=1 if
        //'if pet = "y" Or pet = "n" { //if// //'do the original routine...doesn//'t run for PET scenario
    if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'get flows that happened during previous timestep

        for (int z = 0; z < layers; z++)//z = 0 To layers - 1 //'transpiration, root and soil redistribution
        {
            if (night == "n") { //if// //'it//'s day, must adjust elayer for sun vs. shade weighting
                layerflow = xylem.soils[z]->root.getEcomp(halt) * laisl / lai + xylem.soils[z]->root.getEcomp(haltsh) * laish / lai; //'weighted flow; NOTE: elayer = 0 for layer 0
            }
            else {
                layerflow = xylem.soils[z]->root.getEcomp(halt); //'no adjustment necessary at night; NOTE: elayer = 0 for layer 0
            } //End if// //'night if
            layerflow = layerflow * param.getModelParam("ba_per_ga") * 1.0 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
            layerflow = layerflow + xylem.soils[z]->soilredist * 1.0 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow). NOTE: xylem.soils(0)->soilredist includes soil evaporation for layer 0
            water[z] = water[z] - layerflow; //'subtracts rootflow from layer on per ground area basis
        } //for//z
        //'now do the bottom layer and potential groundwater input
        if (night == "n") { //if// //'it//'s day, must adjust layerflow for sun vs. shade weighting
            layerflow = xylem.soils[layers-1]->root.getEcomp(halt) * laisl / lai + xylem.soils[layers-1]->root.getEcomp(haltsh) * laish / lai; //'weighted flow
        }
        else {
            layerflow = xylem.soils[layers-1]->root.getEcomp(halt); //'no adjustment necessary at night
        } //End if// //'night if
        layerflow = layerflow * param.getModelParam("ba_per_ga") * 1 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
        layerflow = layerflow + xylem.soils[layers - 1]->soilredist * 1 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow)
                                                                            //'water(layers) = water(layers) - layerflow //'subtracts rootflow from layer on per ground area basis
        if (layerflow < 0) { //if// //'water is added
            for (int z = layers; z > 0; z--)//z = layers To 0 Step -1 //'start at bottom, go up
            {
                double deficit = xylem.soils[z - 1]->rhizosphere.getThetaSat() * xylem.soils[z - 1]->depth - water[z]; //'m of water required to wet up layer to SATURATION
                if (-1 * layerflow - deficit >= 0) { //if// //'there//'s enough to wet the layer...remember, negative flow is flow into the layer
                    water[z] = xylem.soils[z - 1]->rhizosphere.getThetaSat() * xylem.soils[z - 1]->depth; //'m water at saturation in layer
                    layerflow = layerflow + deficit; //'reduce what//'s left over for next layers
                }
                else { //'just soak up all the groundwater
                    water[z] = water[z] - layerflow; //'add to bottom layer
                    layerflow = 0; //' all gone
                } //End if// //'wetting if
            } //for//z
            /* Top layer */
            double deficit = xylem.top_soil.rhizosphere.getThetaSat() * xylem.top_soil.depth - water[0]; //'m of water required to wet up layer to SATURATION
            if (-1 * layerflow - deficit >= 0) { //if// //'there//'s enough to wet the layer...remember, negative flow is flow into the layer
                water[0] = xylem.top_soil.rhizosphere.getThetaSat() * xylem.top_soil.depth; //'m water at saturation in layer
                layerflow = layerflow + deficit; //'reduce what//'s left over for next layers
            }
            else { //'just soak up all the groundwater
                water[0] = water[0] - layerflow; //'add to bottom layer
                layerflow = 0; //' all gone
            } //End if// //'wetting if
            runoff = runoff - layerflow; //'add what//'s left to runoff...
        }
        else { //'groundwater is positive...bottom layer is losing water
            water[layers] = water[layers] - layerflow; //'subtract from bottom layer
        } //End if// //'layerflow if

        //'now reset any exhausted layers to extraction limit
        if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            if (water[z] < xylem.soils[z - 1]->swclimit)
                water[z] = xylem.soils[z - 1]->swclimit; //'water at limit
        } //for//z

        bool rainOverride = false;
        // if (iter_Counter == 0 && tod == 23 && (stage_id == STAGE_ID_HIST_OPT || stage_id == STAGE_ID_FUT_OPT))
        //     rainOverride = true;

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
            } //End if// //'wetting up to field capacity "if"
                //}
            } //for//z

            // If there's drainage, and the ground water is on, that should be used to fill up to saturation
            if (rain > 0.0 && ground) // the remaining "drainage" is also still stored in the rain variable
            {
            // this is kind of inefficient, but the rain routine actually drained all the layers to FC even if GW was on and we should have been filling to sat
            // now we start at the bottom and fill the layers to saturation using the drainage

            if (rain > 0) {
                for (int j = layers; j > 0; j--) //j = z - 1 To 0 Step -1 //'go back up to fill profile to saturation
                {
                    double deficit = xylem.soils[j - 1]->rhizosphere.getThetaSat() * xylem.soils[j - 1]->depth - water[j];
                    if (deficit >= 0) { //if// //'got capacity
                        if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                            water[j] = xylem.soils[j - 1]->rhizosphere.getThetaSat() * xylem.soils[j - 1]->depth; //'saturate the layer
                            rain = rain - deficit; //'reduce rain
                        }
                        else { //'rain absorbed by layer
                            water[j] = water[j] + rain;
                            rain = 0; //'use up rain
                        } //End if// //'deficit=>0 "if"
                    }
                    else { //'deficit<0...layer//'s saturated
                        rain = rain - deficit; //'increase rain by super-saturated amount (deficit is negative)
                        water[j] = xylem.soils[j - 1]->rhizosphere.getThetaSat() * xylem.soils[j - 1]->depth; //'reset to saturation
                    } //End if// //'deficit <>0 if
                } //for//j

                /* Top layer */
                double deficit = xylem.top_soil.rhizosphere.getThetaSat() * xylem.top_soil.depth - water[0];
                if (deficit >= 0) { //if// //'got capacity
                    if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                        water[0] = xylem.top_soil.rhizosphere.getThetaSat() * xylem.top_soil.depth; //'saturate the layer
                        rain = rain - deficit; //'reduce rain
                    }
                    else { //'rain absorbed by layer
                        water[0] = water[0] + rain;
                        rain = 0; //'use up rain
                    } //End if// //'deficit=>0 "if"
                }
                else { //'deficit<0...layer//'s saturated
                    rain = rain - deficit; //'increase rain by super-saturated amount (deficit is negative)
                    water[0] = xylem.top_soil.rhizosphere.getThetaSat() * xylem.top_soil.depth; //'reset to saturation
                } //End if// //'deficit <>0 if
            }

            runoff = runoff + rain; //'whatever is left over will run off
            drainage = 0; //'no drainage if any layer is rising above field capacity
            rain = 0; //'reset rain to zero
            }
            //'sumdrain = sumdrain + drainage //'total drainage
        } //End if// //'rain if

        //'now check for exhausted layers
        if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
        for (int z = 1; z <= layers; z++)//z = 1 To layers
        {
            if (water[z] < xylem.soils[z - 1]->swclimit)
            layer[z] = 1; //'water exhausted
        } //for//z
        //'now get water content change over PREVIOUS time step
        if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'now get updated new water content
            double waternew = 0;
            for (int z = 0; z <= layers; z++)//z = 0 To layers //'check for exhausted layers
            {
                waternew = water[z];
            } //for//z

            //'waternew = waternew * 55555556# //'new water content at beginning of timestep
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
            data.setColumnValue(data.getColumnValue("rain", dd) + data.getColumnValue("end-ground-water", dd), dd, "end-total-water-input"); //'total input per timestep in mm
            
            gs_data.setColumnValue(gs_data.getColumnValue("input", gs_yearIndex) + data.getColumnValue("end-total-water-input", dd), gs_yearIndex, "input");

            data.setColumnValue(runoff * 1000, dd, "end-runoff"); //'excess root zone water per timestep in mm
            waterold = waternew; //'reset to beginning of current timestep
        } //End if// //'dd>1 if
    } //End if// //'dd>1 if
        //'} //End if////'pet if
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
    } //End if// //'dd>1 if

    double kpday1, kxday1;
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
    } //End if// //'dd=16 if
    
    if (gs_doneFirstDay) { //if// //'calculate plc relative to midday of day 1
        if (iter_refK < 0.000000001) // || iter_Counter == 0 // no longer appropriate to test for iter_Counter == 0 ... may be doing a seperate stress profile that refers to a saved refK, which will have been loaded at start of modelProgramMain
            tempDouble = 100 * (1 - data.getColumnValue("K-plant", dd - 1) / kpday1); //' done to avoid repeating this calculation
        else // if we haven't loaded a ref K it will be set to zero, and we need to fall back to the old method.. otherwise use refK to calculate PLC
        {
            tempDouble = 100 * (1 - data.getColumnValue("K-plant", dd - 1) / iter_refK); // PLCp calculated from refK if it exists
            if (tempDouble < 0.0) // if we're using the refK, it's possible for this to be negative briefly -- should be considered zero
                tempDouble = 0.0;
        }
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
            // if we've done the first GS day but we're NOT in the growing season, it's the winter following GS.
            // this is considered NEXT YEAR's off-season input! First check if we already have a value for the FINAL water for THIS YEAR, because if it's zero then
            // this is the first timestep of the winter and we need to store it
            if (gs_data.getColumnValue("water-final", gs_yearIndex) <= 0.0 && gs_data.getColumnValue("initial-off", gs_yearIndex + 1) <= 0.0) // ok to use == or <= with double here because this will be memset to 0 if it hasn't been set
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
                gs_data.setColumnValue(gs_data.getColumnValue("water-initial-off", gs_yearIndex + 1) + data.getColumnValue("water-content", dd), 
                                       gs_yearIndex + 1, 
                                       "water-initial-off"); // add the stored input to the input tally, for NEXT YEAR
            }
        }
        // [/HNT]
    } //End if// //'dd>16 if
    else if (!gs_inGrowSeason)// have NOT done first day and are NOT in growing season
    {
        // we must be in the pre-GS winter of what we called the NEXT year above... so gs_yearIndex is now that year, and we add to the off-season input for THIS year
        gs_data.setColumnValue(gs_data.getColumnValue("water-input-off", gs_yearIndex) + data.getColumnValue("end-total-water-input", dd), 
                               gs_yearIndex, 
                               "water-input-off");
    }
}