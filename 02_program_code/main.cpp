#include "01Utilities.h"
#include "09Plant.h"

// For making output directory
#include <filesystem>

int main(int argc, char *argv[])
{
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|    CARBON GAIN VS HYDRAULIC RISK MODEL V 3.0    |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;

    Plant *plantModel = new Plant(0);
    std::string param_data = CONFIG_FILE_PATH;
    std::string config_data = PARAMETER_FILE_PATH;
    int config_setting = CONFIG_SETTING; // default to just use first row in configuration
    std::string output_dir = OUT_DIR;

    // If arguments are provided, use provided config and param data files, and set config_setting
    if (argc > 1) {
        if (strcmp(argv[1], "--help") == 0) {
            std::cout << "Usage: ./run [parameter_file] [config_file] [config_setting] [output_dir]" << std::endl;
            std::cout << "All arguments are optional:" << std::endl;
            std::cout << "  parameter_file: Path to the parameter data file (default: " << PARAMETER_FILE_PATH << ")" << std::endl;
            std::cout << "  config_file: Path to the configuration data file (default: " << CONFIG_FILE_PATH << ")" << std::endl;
            std::cout << "  config_setting: Configuration setting index (default: " << CONFIG_SETTING + 1 << ")" << std::endl;
            std::cout << "  output_dir: Path to output directory (default: " << OUT_DIR << ")" << std::endl;
            std::cout << "Example: ./run params.csv config.csv 2 ./out" << std::endl;
            exit(0);
        } else {
            param_data = argv[1];
        }

        if (argc > 2) // if config file is provided
            config_data = argv[2];

        if (argc > 3) // config setting provided
            config_setting = std::atoi(argv[3]) - 1;

        if (argc > 4) { // output dir specified
            output_dir = argv[4];

            if (std::filesystem::create_directories(output_dir)) {
                printf("Directory %s made successfully.\n\n", output_dir.c_str());
            } else {
                printf("Error creating directory %s, it either exists already, or some other error has occured.\n\n", output_dir.c_str());
            }
        }
    }

    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|            READING MODEL INPUT FILES            |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;
    bool lrSuccess = locateRanges(plantModel->config_data, config_data, plantModel->param_data, param_data);
    // It also loads parameter and configuration file

    if (!lrSuccess)
    {
        std::cout << "Unrecoverable model failure!" << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        return 1; // failure, unrecoverable
    }
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "|              MODEL CONFIGURATION               |" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << std::endl;

    plantModel->cleanModelVars();
    std::cout << "Model variables cleaned." << std::endl;

    plantModel->setConfig(config_setting);
    std::cout << "Configuration set." << std::endl;

    plantModel->initModelVars();
    std::cout << "Model variables initialized." << std::endl;

    plantModel->readin();
    std::cout << "Model input read." << std::endl;

    if (plantModel->stage_id == STAGE_ID_FUT_STRESS_NOACCLIM) // override some if we're doing the odd "no acclimation stress profile"
    {
        // std::cout << "Stage " << plantModel->stage_id << "; NoAcclim Stress Profile, overriding historical ca " << ca << " -> " << stage_CO2Fut << " and ksatp " << ksatp << " -> " << stage_KmaxFut << std::endl;
    }

    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "|             CLIMATE FORCING FILES              |" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    readDataSheet(plantModel->data, plantModel->climate_forcing_data_path, plantModel->data_header_file_path);
    readGSSheet(plantModel->gs_data, plantModel->growing_season_limits_data_path, plantModel->useGSData);

    readGrowSeasonData(plantModel->param, plantModel->gs_data);
    std::cout << "Growing season data read." << std::endl;

    if ((plantModel->iter_useAreaTable)) {
        readSiteAreaValues();
        std::cout << "Site area values read." << std::endl;
    }

    plantModel->resetLayerStatus();

    int dd = -1,
    successCode = 0;

    std::string species = plantModel->param_data("i_sp", plantModel->species_no);
    std::string region = plantModel->param_data("i_region", plantModel->species_no);
    std::string site = plantModel->param_data("i_site", plantModel->species_no);

    do //loop through time steps
        {

            dd += 1;
            successCode = plantModel->modelTimestepIter(dd);

            if (successCode == -1)
            {
                plantModel->data.output(output_dir + "/timesteps_output_" + species + "_" + region + "_" + site + ".csv");
                plantModel->gs_data.output(output_dir + "/sum_output_" + species + "_" + region + "_" + site + ".csv");
                std::cout << "Unrecoverable model failure!" << std::endl;
                return 1; // failure, unrecoverable
            }
            else if (successCode > 0) // this returns the year if we've incremented it -- not necessary in the full C version (also only supports 1 year right now)
            {
                dd = dd - 1; //we need to repeat this timestep because we bailed early when finding a new year
                plantModel->gs_yearIndex = successCode;

                // if we're running without growing season limits, we need to record the "end of GS" water content now
                // because we did not complete the previous timestep, back up 1 more to grab a value
                if (!plantModel->useGSData && plantModel->gs_yearIndex > 0)
                    plantModel->gs_data.setColumnValue(plantModel->data.getColumnValue("water-content", dd - 1), plantModel->gs_yearIndex - 1, "water-final");  // make sure this goes with the previous year


                plantModel->modelProgramNewYear();
            }
            else // 0 = success, VBA bool convention
            {
                                        // if we're running without growing season limits, we need to record the "end of GS" water content now
                if (!plantModel->useGSData && plantModel->gs_yearIndex > 0)
                    plantModel->gs_data.setColumnValue(plantModel->data.getColumnValue("water-content", dd - 1), plantModel->gs_yearIndex, "water-final"); // if this was the end of the set of years, gs_yearIndex will not have been changed so use as-is
            }
            if (dd % 1000 == 0)
                std::cout << "Timestep " << dd << " completed" << std::endl;
        } while (!(plantModel->data.getColumnValue("julian-day", dd + 1) < 0.01)); // loop until the jd value on next row is zero -- it's an integer, but everything is stored in the array as double

    plantModel->data.output(output_dir + "/timesteps_output_" + species + "_" + region + "_" + site + ".csv");
    plantModel->gs_data.output(output_dir + "/sum_output_" + species + "_" + region + "_" + site + ".csv");

    delete plantModel;

    return 0;
}
