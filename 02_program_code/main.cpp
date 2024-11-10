#include "01Utilities.h"
#include "09Plant.h"

int main()
{
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|    CARBON GAIN VS HYDRAULIC RISK MODEL V 3.0    |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;

    Plant *plantModel = new Plant(0);

    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << "|            READING MODEL INPUT FILES            |" << std::endl;
    std::cout << " -------------------------------------------------" << std::endl;
    std::cout << std::endl;
    bool lrSuccess = locateRanges(plantModel->config_data, plantModel->param_data); //Finds all of the input/output sections across the workbook
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
    plantModel->setConfig();
    std::cout << "Configuration set." << std::endl;

    plantModel->cleanModelVars();
    std::cout << "Model variables cleaned." << std::endl;

    plantModel->initModelVars();
    std::cout << "Model variables initialized." << std::endl;

    plantModel->readin();
    std::cout << "Model input read." << std::endl;

    if (plantModel->stage_id == STAGE_ID_FUT_STRESS_NOACCLIM) // override some if we're doing the odd "no acclimation stress profile"
    {
        // std::cout << "Stage " << plantModel->stage_id << "; NoAcclim Stress Profile, overriding historical ca " << ca << " -> " << stage_CO2Fut << " and ksatp " << ksatp << " -> " << stage_KmaxFut << std::endl;
    }

    readDataSheet(plantModel->data, plantModel->sum_data, plantModel->climate_forcing_data_path, plantModel->data_header_file_path, plantModel->sum_header_file_path);
    readGSSheet(plantModel->gs_data, plantModel->growing_season_limits_data_path);

    readGrowSeasonData(plantModel->param, plantModel->gs_data);
    std::cout << "Growing season data read." << std::endl;

    if ((plantModel->iter_useAreaTable)) {
        readSiteAreaValues();
        std::cout << "Site area values read." << std::endl;
    }

    plantModel->resetLayerStatus();

    std::cout << "Calculating critical points for components" << std::endl;
    plantModel->componentPCrits();

    int dd = 0,
    successCode = 0;

    do //loop through time steps
        {
            successCode = plantModel->modelTimestepIter(dd);

            if (successCode == 1)
            {
                std::cout << "Unrecoverable model failure!" << std::endl;
                return 1; // failure, unrecoverable
            }
            else if (successCode > 0) // this returns the year if we've incremented it -- not necessary in the full C version (also only supports 1 year right now)
            {
                std::cout << "Repeating previous year" << std::endl;
                dd = dd - 1; //we need to repeat this timestep because we bailed early when finding a new year
                plantModel->gs_yearIndex = successCode;

                // if we're running without growing season limits, we need to record the "end of GS" water content now
                // because we did not complete the previous timestep, back up 1 more to grab a value
                // if (!plantModel->useGSData && plantModel->gs_yearIndex > 0)
                //     plantModel->gs_ar_waterFinal_GS[plantModel->gs_yearIndex - 1] = dSheet.Cells(rowD + dd - 1, colD + dColF_End_watercontent); // make sure this goes with the previous year

                // modelProgramNewYear();
            }
            // else // 0 = success, VBA bool convention
            // {
            //     int breakpoint = 1137; // success, in the C version we just continue instead of outputting
            //                             // do all CSV writing at the end

            //                             // if we're running without growing season limits, we need to record the "end of GS" water content now
            //     if (!plantModel->useGSData && plantModel->gs_yearIndex > 0)
            //         plantModel->gs_ar_waterFinal_GS[plantModel->gs_yearIndex] = dSheet.Cells(rowD + dd - 1, colD + dColF_End_watercontent); // if this was the end of the set of years, gs_yearIndex will not have been changed so use as-is
            // }

            if (dd % 1000 == 0)
                std::cout << "Timestep " << dd << " completed" << std::endl;
            dd += 1;
        } while (!(plantModel->data.getColumnValue("julian-day", dd) < 0.01)); // loop until the jd value on next row is zero -- it's an integer, but everything is stored in the array as double

    delete plantModel;

    return 0;
}
