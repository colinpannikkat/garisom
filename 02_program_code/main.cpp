#include "01Utilities.h"
#include "09Plant.h"

int main(int argc, char *argv[])
{
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|    CARBON GAIN VS HYDRAULIC RISK MODEL V 3.0    |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;

    Plant *plantModel = new Plant(0);
    std::string param_data;
    std::string config_data;

    // If arguments are provided, use provided config and param data files
    if (argc > 1) {
        param_data = argv[1];
        if (argc == 3) { // if config file is provided
            config_data = argv[2];
        } else {
            config_data = CONFIG_FILE_PATH;
        }
    } else { // default
        param_data = PARAMETER_FILE_PATH;
        config_data = CONFIG_FILE_PATH;
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
        return 0; // failure, unrecoverable
    }
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << "|              MODEL CONFIGURATION               |" << std::endl;
    std::cout << " ------------------------------------------------" << std::endl;
    std::cout << std::endl;

    plantModel->cleanModelVars();
    std::cout << "Model variables cleaned." << std::endl;

    plantModel->setConfig();
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
    
    readDataSheet(plantModel->data, plantModel->sum_data, plantModel->climate_forcing_data_path, plantModel->data_header_file_path, plantModel->sum_header_file_path);
    readGSSheet(plantModel->gs_data, plantModel->growing_season_limits_data_path, plantModel->useGSData);

    readGrowSeasonData(plantModel->param, plantModel->gs_data);
    std::cout << "Growing season data read." << std::endl;

    if ((plantModel->iter_useAreaTable)) {
        readSiteAreaValues();
        std::cout << "Site area values read." << std::endl;
    }

    plantModel->resetLayerStatus();

    std::cout << "Calculating critical points for components" << std::endl;
    plantModel->componentPCrits();

    for (int k = 1; k <= plantModel->layers; k++) // k = 1 To layers //exclude the top layer
    {
        plantModel->xylem.soils[k]->root.setKmin(plantModel->xylem.soils[k]->root.getKmax());
    }

    plantModel->xylem.stem.setKmin(plantModel->xylem.stem.getKmax());
    plantModel->xylem.leaf.setKmin(plantModel->xylem.leaf.getKmax());
    plantModel->kmin = plantModel->param.getModelParam("ksatp");

    int dd = 0,
    successCode = 0;

    do //loop through time steps
        {

            dd += 1;
            successCode = plantModel->modelTimestepIter(dd);

            if (successCode == -1)
            {
                plantModel->data.output("test_fail_output.csv");
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

    plantModel->data.output("high-stress_test_output.csv");
    plantModel->gs_data.output("high-stress_test_sum_output.csv");

    delete plantModel;

    return 0;
}
