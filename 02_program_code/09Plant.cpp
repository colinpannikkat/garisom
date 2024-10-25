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

    for (int k = 0; k < param.getModelParam("layers"); k++) {
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
    get_vgparams(param_data.getColumnValue("i_texture", species_no), param.getModelParam("layers"), xylem.soils); // extract VG parameters for a given texture

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
        for (int k = 0; k < param.getModelParam("layers"); k++) {
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

    double layers = param.getModelParam("layers");
    xylem.num_layers = layers;
    //for this soil data we want to use the original anchor-offset system
    for (int k = 0; k < layers; k++) { // set layer depths and % root ksat
        xylem.soils[k]->layer_depth = 0.01 * log(1.0 - ((k + 1) * 0.995) / layers) / log(param.getModelParam("beta")); // lower depth of each layer converted to m
    }
    double depthmax = xylem.soils[layers - 1]->layer_depth;
    // calculate transport distances
    std::vector<double> temp_vert_distance;
    // first get vertical distance to biomass center of each layer
    for (int k = 0; k < layers * 2.0; k++)
    {
        temp_vert_distance.push_back(0.01 * log(1 - ((k + 1) * 0.995) / (layers * 2.0)) / log(param.getModelParam("beta"))); // get half depths
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
        xylem.soils[k]->radius = pow(temp / xylem.soils[k]->depth * PI, 0.5); // width in m
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