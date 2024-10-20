#include "01Utilities.h"
#include "09Plant.h"

int main()
{
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|    CARBON GAIN VS HYDRAULIC RISK MODEL V 3.0    |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;

    Plant plantModel(0);

    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|            READING MODEL INPUT FILES            |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;
    bool lrSuccess = locateRanges(plantModel.config_data, plantModel.param_data); //Finds all of the input/output sections across the workbook
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
    plantModel.setConfig();

    plantModel.cleanModelVars();
    plantModel.initModelVars();
    plantModel.readin();

    if (plantModel.stage_id == STAGE_ID_FUT_STRESS_NOACCLIM) // override some if we're doing the odd "no acclimation stress profile"
    {
        // std::cout << "Stage " << plantModel.stage_id << "; NoAcclim Stress Profile, overriding historical ca " << ca << " -> " << stage_CO2Fut << " and ksatp " << ksatp << " -> " << stage_KmaxFut << std::endl;
    }

    readDataSheet(plantModel.data, plantModel.sum_data, plantModel.climate_forcing_data_path, plantModel.data_header_file_path, plantModel.sum_header_file_path);
    readGSSheet(plantModel.gs_data, plantModel.growing_season_limits_data_path);
    readGrowSeasonData(plantModel.param, plantModel.gs_data);

    if ((plantModel.iter_useAreaTable)) {
        readSiteAreaValues();
    }

    plantModel.componentPCrits();

    return 0;
}
