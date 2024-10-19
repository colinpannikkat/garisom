#include "01Utilities.h"
#include "09Plant.h"


int main()
{
    Plant plantModel(0);

    locateRanges(plantModel.config_data, plantModel.param_data);
    plantModel.setConfig();
    // cleanModelVars();
    // initModelVars();
    // readin(); //get all global parameters
    // readDataSheet();
    // readGSSheet();
    // readGrowSeasonData(); //hnt todo cleanup - should just put this in readin?
    // if ((iter_useAreaTable)) {
    //     readSiteAreaValues();
    // }

    plantModel.readin();

    readDataSheet(plantModel.data, plantModel.sum_data, plantModel.climate_forcing_data_path, plantModel.data_header_file_path, plantModel.sum_header_file_path);
    readGSSheet(plantModel.gs_data, plantModel.growing_season_limits_data_path);
    readGrowSeasonData(plantModel.param, plantModel.gs_data);

    return 0;
}
