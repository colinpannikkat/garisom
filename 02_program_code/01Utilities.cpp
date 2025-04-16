#include "01Utilities.h"
#include "11Parameters.h"

/**
 * @brief Reads configuration and parameter data from CSV files.
 *
 * This function sets the filenames for the configuration and parameter files,
 * reads the model configuration from the configuration file, and reads the
 * parameters from the parameter file. It outputs the progress to the console.
 *
 * @param config_data Reference to a CSVData object where the configuration data will be stored.
 * @param param_data Reference to a CSVData object where the parameter data will be stored.
 * @return true if the operation is successful.
 */
bool locateRanges(CSVData<std::string> &config_data, std::string config_file_name, CSVData<std::string> &param_data, std::string param_file_name) {
    
    // Reading model controls
    std::cout << std::endl;
    std::cout << "Reading model configuration from " << config_file_name << std::endl;

    config_data = CSVData<std::string>(config_file_name);

    std::cout << "Reading parameters from " << param_file_name << std::endl;

    param_data = CSVData<std::string>(param_file_name);       

    return true;
}

/**
 * @brief Reads the growing season limits from a CSV file and populates the provided CSVData object.
 *
 * This function reads data from a specified CSV file containing growing season limits and stores it in the provided
 * CSVData object. If no file name is provided, an assertion error will occur.
 * 
 * Even if years are being treated as independent (no growing season), a file with years
 * should still be provided to act as the base for summary output.
 *
 * @param gs_data Reference to a CSVData object where the growing season data will be stored.
 * @param gs_file_name Reference to a string containing the name of the CSV file to be read.
 *
 * @note If the file cannot be opened, an error message will be printed to the console.
 * @note If no file name is provided, an assertion error will occur.
 */
void readGSSheet(CSVData<double> &gs_data, std::string &gs_file_name, bool use_gs_data) {

    gs_data.clear();
    assert(gs_file_name != "");

    std::cout << std::endl;
    std::cout << "Reading growing season limits." << std::endl;

    std::ifstream data_file(gs_file_name);
    if (!data_file.is_open())
        std::cout << "FAILED to open GS RANGE file " << gs_file_name << std::endl;

    std::cout << "Opened GS RANGE file " << gs_file_name << std::endl;

    gs_data = CSVData<double>(gs_file_name);
    
    std::cout << "Finished GS RANGE read." << std::endl;
    std::cout << std::endl;

    if (!use_gs_data) {

        std::cout << "Years will be treated as independent." << std::endl;
        std::cout << std::endl;

    }
}

/**
 * @brief Loads climate forcing variables from CSV files into the provided data structures.
 * 
 * This function reads data from the specified CSV files and populates the provided CSVData objects.
 * It first clears the existing data in the `data` object, then initializes both `data` and `sum_data`
 * with the contents of the specified files.
 * 
 * @param data Reference to a CSVData object where the main data will be loaded.
 * @param sum_data Reference to a CSVData object where the summary data will be loaded.
 * @param data_file_name The name of the file containing the main data.
 * @param header_file_name The name of the file containing the header information for the main data.
 * @param sum_header_file_name The name of the file containing the header information for the summary data.
 */
void readDataSheet(CSVData<double> &data, std::string &data_file_name, std::string &header_file_name) {

    data = CSVData<double>(data_file_name, header_file_name);

    std::cout << "Finished climate data read." << std::endl;
}

/**
 * @brief Reads grow season data from a CSVData object and populates the Parameters object.
 *
 * This function iterates over each row in the provided CSVData object and extracts the values
 * for "Year", "Start_day", "End_day", and "ca_ppm". These values are then used to populate
 * the corresponding fields in the Parameters object.
 *
 * @param param Reference to a Parameters object that will be populated with the grow season data.
 * @param gs_data Reference to a CSVData<double> object containing the grow season data.
 */
void readGrowSeasonData(Parameters &param, CSVData<double> &gs_data) {
    int n_rows = gs_data.row_size();
    int temp_int = 0;

    for (int i = 0; i < n_rows; i++) {
        temp_int = gs_data.getColumnValue("Year", i);
        param.setGsArYear(i, temp_int);

        temp_int = gs_data.getColumnValue("Start_day", i);
        param.setGsArStart(i, temp_int);

        temp_int = gs_data.getColumnValue("End_day", i);
        param.setGsArEnd(i, temp_int);

        param.setGsArPpm(i, gs_data.getColumnValue("ca_ppm", i));
    }

    std::cout << "Finished grow season data read." << std::endl;
}

void readSiteAreaValues() {
    std::cout << "Nothing" << std::endl;
    return;
}