#include "01Utilities.h"

/**
 * @brief Converts a string to an integer.
 *
 * This function attempts to convert the given string to an integer using std::stoi.
 * If the conversion fails due to an invalid argument, it throws a runtime_error with
 * an appropriate error message.
 *
 * @param str The string to be converted to an integer.
 * @return The integer value represented by the string.
 * @throws std::runtime_error If the conversion fails due to an invalid argument.
 */
int toInt(const std::string& str) {
    try {
        return std::stoi(str);
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Conversion to int failed: " + str);
    }
}

/**
 * @brief Converts a string to a double.
 *
 * This function attempts to convert the given string to a double using std::stod.
 * If the conversion fails due to an invalid argument, it throws a runtime_error.
 *
 * @param str The string to be converted to a double.
 * @return The double value represented by the string.
 * @throws std::runtime_error if the conversion fails.
 */
double toDouble(const std::string& str) {
    try {
        return std::stod(str);
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Conversion to double failed: " + str);
    }
}

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * This function removes leading spaces and trailing spaces, newlines, and carriage returns
 * from the input string.
 *
 * @param str The input string to be trimmed.
 * @return A new string with leading and trailing whitespace removed.
 */
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(" \n\r"); // trim carriage return or new line
    if (first == std::string::npos) {
        return ""; // String is empty or all whitespace
    }
    return str.substr(first, (last - first + 1));
}

/**
 * @brief Reads the header line from a CSV file and stores the column names with their respective indices.
 * 
 * This function reads the first line of the provided CSV file stream, which is expected to be the header line.
 * It then parses the header line to extract column names, trims any leading or trailing whitespace from each column name,
 * and stores the column names along with their corresponding indices in the `header` map.
 * 
 * @param data_file A reference to an input file stream object representing the CSV file.
 */
template <typename T>
void CSVData<T>::readHeader(std::ifstream &data_file) {
    std::string line;

    // Read header
    if (getline(data_file, line)) {
        std::stringstream ss(line);
        std::string column;
        int index = 0;
        while (getline(ss, column, ',')) {
            column = trim(column);
            if (column.empty() || column == "\0" || column == "\r") continue;
            header[column] = index++;
        }
    }
}

/**
 * @brief Reads data from a CSV file and stores it in the CSVData object as string.
 * 
 * This function reads each line from the provided input file stream, parses the line into cells
 * separated by commas, and stores the non-empty cells in a vector. Each vector of cells (representing
 * a row) is then added to the data member of the CSVData object.
 * 
 * @param data_file The input file stream from which to read the CSV data.
 */
template <>
void CSVData<std::string>::readData(std::ifstream &data_file) {
    std::string line;
    
    // Read data
    while (getline(data_file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (getline(ss, cell, ',')) {
            if (cell.empty() || cell == "\0" || cell == "\r") continue;
            row.push_back(cell);
        }
        if (row.size() == 0) continue;  // for some reason this is some data with 
                                        // invisible character that needs to be ignored
        data.push_back(row);
    }
}

/**
 * @brief Reads data from a CSV file and stores it in the CSVData object as double.
 * 
 * This function reads each line from the provided input file stream, parses the line into cells
 * separated by commas, and stores the non-empty cells in a vector. Each vector of cells (representing
 * a row) is then added to the data member of the CSVData object.
 * 
 * @param data_file The input file stream from which to read the CSV data.
 */
template <>
void CSVData<double>::readData(std::ifstream &data_file) {
    std::string line;
    
    // Read data
    while (getline(data_file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> row;
        while (getline(ss, cell, ',')) {
            if (cell.empty() || cell == "\0" || cell == "\r") continue;
            row.push_back(toDouble(cell));
        }
        if (row.size() == 0) continue;  // for some reason this is some data with 
                                        // invisible character that needs to be ignored
        data.push_back(row);
    }
}

/**
 * @brief Constructs a CSVData object and initializes it by reading data from the specified file.
 * 
 * This constructor attempts to open the specified CSV file and read its contents. If the file
 * cannot be opened, the program will print an error message and terminate.
 * 
 * @param data_file_name The name of the CSV file to be read.
 * 
 * @throws std::runtime_error If the file cannot be opened.
 */
template <typename T>
CSVData<T>::CSVData(const std::string& data_file_name) {
            
    std::ifstream data_file(data_file_name);

    if (!data_file.is_open()) {
        std::cout << "UNRECOVERABLE: Failed to open file " << data_file_name << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    data.clear();
    header.clear();

    readHeader(data_file);
    readData(data_file);
}

/**
 * @brief Constructs a CSVData object by reading data and header files.
 *
 * This constructor initializes a CSVData object by opening and reading the specified
 * data and header files. If the files cannot be opened, the program will terminate
 * with an error message. Optionally, the first line of the data file can be skipped,
 * which is useful if the first line contains headers.
 *
 * @param data_file_name The name of the data file to be read.
 * @param header_file_name The name of the header file to be read.
 * @param skip_first A boolean flag indicating whether to skip the first line of the data file. Default is true.
 */
template <typename T>
CSVData<T>::CSVData(const std::string &data_file_name, const std::string &header_file_name, bool skip_first/*=true*/) {
    std::ifstream data_file(data_file_name);
    
    if (!data_file.is_open()) {
        std::cout << "UNRECOVERABLE: Failed to open file " << data_file_name << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    // Easy hack to skip first line in data file (header)
    if (skip_first) {
        std::string line;
        getline(data_file, line);
    }

    std::ifstream header_file(header_file_name);

    if (!header_file.is_open()) {
        std::cout << "UNRECOVERABLE: Failed to open file " << header_file_name << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    readHeader(header_file);
    readData(data_file);
}

/**
 * @brief Checks if the CSV data is empty.
 * 
 * This function returns true if the CSV data container is empty, 
 * indicating that there are no data entries present.
 * 
 * @return true if the CSV data is empty, false otherwise.
 */
template <typename T>
bool CSVData<T>::empty() {
    return data.empty();
}

/**
 * @brief Prints the CSV data in a paginated format.
 *
 * This function prints the CSV data with a specified number of rows per page.
 * It displays the column headers at the top of each page and clears the terminal
 * screen before printing each page. The user can navigate through the pages by
 * pressing any key to continue or 'q' to quit.
 *
 * @param page_size The number of rows to display per page.
 */
template <typename T>
void CSVData<T>::print(size_t page_size) {
    size_t total_rows = data.size();
    size_t total_pages = (total_rows + page_size - 1) / page_size;

    // Print column headers in order
    std::vector<std::string> ordered_headers(header.size());
    for (const auto& pair : header) {
        ordered_headers[pair.second] = pair.first;
    }
    for (const auto& column : ordered_headers) {
        std::cout << column << "\t";
    }
    std::cout << std::endl;

    for (size_t page = 0; page < total_pages; ++page) {
        size_t start_row = page * page_size;
        size_t end_row = std::min(start_row + page_size, total_rows);

        // Clear the terminal screen
        std::cout << "\033[2J\033[1;1H";

        // Print column headers again for each page
        for (const auto& header : ordered_headers) {
            std::cout << std::setw(5) << header << "\t";
        }
        std::cout << std::endl;

        std::cout << "Page " << (page + 1) << " of " << total_pages << std::endl;

        // Print rows
        for (size_t row = start_row; row < end_row; ++row) {
            for (const auto& cell : data[row]) {
            try {
                double num = std::stod(cell);
                std::cout << std::setw(FIO_PRECISION+5) << std::fixed << std::setprecision(FIO_PRECISION) << num << "\t";
            } catch (const std::invalid_argument&) {
                std::cout << std::setw(FIO_PRECISION+5) << cell << "\t";
            }
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

        if (page < total_pages - 1) {
            std::cout << "Press any key to continue to the next page or 'q' to quit: ";
            char c;
            std::cin >> c;
            if (c == 'q' || c == 'Q') {
                break;
            }
        }
    }
}

/**
 * @brief Prints the CSV data in a paginated format.
 *
 * This function prints the CSV data with a specified number of rows per page.
 * It displays the column headers at the top of each page and clears the terminal
 * screen before printing each page. The user can navigate through the pages by
 * pressing any key to continue or 'q' to quit.
 *
 * @param page_size The number of rows to display per page.
 */
template <>
void CSVData<double>::print(size_t page_size) {
    size_t total_rows = data.size();
    size_t total_pages = (total_rows + page_size - 1) / page_size;

    // Print column headers in order
    std::vector<std::string> ordered_headers(header.size());
    for (const auto& pair : header) {
        ordered_headers[pair.second] = pair.first;
    }
    for (const auto& column : ordered_headers) {
        std::cout << column << "\t";
    }
    std::cout << std::endl;

    for (size_t page = 0; page < total_pages; ++page) {
        size_t start_row = page * page_size;
        size_t end_row = std::min(start_row + page_size, total_rows);

        // Clear the terminal screen
        std::cout << "\033[2J\033[1;1H";

        // Print column headers again for each page
        for (const auto& header : ordered_headers) {
            std::cout << std::setw(5) << header << "\t";
        }
        std::cout << std::endl;

        std::cout << "Page " << (page + 1) << " of " << total_pages << std::endl;

        // Print rows
        for (size_t row = start_row; row < end_row; ++row) {
            for (const auto& cell : data[row])
                std::cout << std::setw(FIO_PRECISION+5) << std::fixed << std::setprecision(FIO_PRECISION) << cell << "\t";
            std::cout << std::endl;
        }

        std::cout << std::endl;

        if (page < total_pages - 1) {
            std::cout << "Press any key to continue to the next page or 'q' to quit: ";
            char c;
            std::cin >> c;
            if (c == 'q' || c == 'Q') {
                break;
            }
        }
    }
}

/**
 * @brief Returns the number of rows in the CSV data.
 * 
 * This function calculates and returns the total number of rows
 * present in the CSV data by determining the size of the data container.
 * 
 * @return int The number of rows in the CSV data.
 */
template <typename T>
int CSVData<T>::row_size() {
    return data.size();
}

/**
 * @brief Returns the number of columns in the CSV data.
 * 
 * This function assumes that the CSV data is stored in a 2D vector
 * and returns the size of the first row, which represents the number
 * of columns in the CSV data.
 * 
 * @return int The number of columns in the CSV data.
 */
template <typename T>
int CSVData<T>::col_size() {
    return data[0].size();
}

/**
 * @brief Retrieves the value of a specified column in a given row.
 * 
 * This function fetches the value from the CSV data for the specified row and column name.
 * It asserts that the row index is within the bounds of the data and that the column name exists.
 * If the column name does not exist, an error message is printed to stderr.
 * 
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The index of the row from which to retrieve the value. Default is 0.
 * @return std::string The value of the specified column in the given row. 
 *         Returns an empty string if the column name does not exist.
 */
template <typename T>
T CSVData<T>::getColumnValue(const std::string &column_name, int row/*=0*/) {
    assert(row < data.size());
    if (header.find(column_name) != header.end()) {
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        return data[row][col_index];
    } else {
        fprintf(stderr, "Column name %s does not exist\n", column_name.c_str());
    }
}

/**
 * @brief Retrieves the value of a specified column in a given row.
 * 
 * This function fetches the value from the CSV data for the specified row and column name.
 * It asserts that the row index is within the bounds of the data and that the column name exists.
 * If the column name does not exist, an error message is printed to stderr.
 * 
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The index of the row from which to retrieve the value. Default is 0.
 * @return std::string The value of the specified column in the given row. 
 *         Returns an empty string if the column name does not exist.
 */
template <>
std::string CSVData<std::string>::getColumnValue(const std::string &column_name, int row/*=0*/) {
    assert(row < data.size());
    if (header.find(column_name) != header.end()) {
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        return data[row][col_index];
    } else {
        fprintf(stderr, "Column name %s does not exist\n", column_name.c_str());
    }
    return "";
}

/**
 * @brief Retrieves the value of a specified column in a given row.
 * 
 * This function fetches the value from the CSV data for the specified row and column name.
 * It asserts that the row index is within the bounds of the data and that the column name exists.
 * If the column name does not exist, an error message is printed to stderr.
 * 
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The index of the row from which to retrieve the value. Default is 0.
 * @return std::string The value of the specified column in the given row. 
 *         Returns an empty string if the column name does not exist.
 */
template <>
double CSVData<double>::getColumnValue(const std::string &column_name, int row/*=0*/) {
    assert(row < data.size());
    if (header.find(column_name) != header.end()) {
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        return data[row][col_index];
    } else {
        fprintf(stderr, "Column name %s does not exist\n", column_name.c_str());
    }
    return 0.0;
}

/**
 * @brief Retrieves the value from a specified column and row in the CSV data.
 *
 * This function searches for the specified column name in the CSV header and retrieves
 * the value from the given row. If the column name is found and the row index is valid,
 * the value is stored in the provided string reference and the function returns true.
 * If the column name does not exist or the row index is invalid, an error message is
 * printed to stderr, the provided string reference is set to an empty string, and the
 * function returns false.
 *
 * @param store A reference to a string where the retrieved value will be stored.
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The index of the row from which to retrieve the value. Default is 0.
 * @return true if the value is successfully retrieved, false otherwise.
 */
template <>
bool CSVData<std::string>::getColumnValue(std::string &store, const std::string &column_name, int row/*=0*/) {
    assert(row < data.size());
    if (header.find(column_name) != header.end()) {
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        store = data[row][col_index];
        return true;
    } else {
        fprintf(stderr, "Column name %s does not exist\n", column_name.c_str());
    }
    store = "";
    return false;
}

/**
 * @brief Retrieves the value of a specified column in a given row and stores it as an integer.
 * 
 * This function fetches the value from the specified column and row in the CSV data,
 * converts it to an integer, and stores it in the provided reference variable.
 * 
 * @param store Reference to an integer where the retrieved value will be stored.
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The row index from which to retrieve the value. Default is 0.
 * @return true if the value was successfully retrieved and converted to an integer, false otherwise.
 */
template <>
bool CSVData<std::string>::getColumnValue(int &store, const std::string &column_name, int row/*=0*/) {
    std::string value;
    getColumnValue(value, column_name, row);
    if (value == "") {
        store = 0;
        return false;
    }
    store = toInt(value);
    return true;
}

/**
 * @brief Retrieves the value of a specified column in a given row and stores it as an bool.
 * 
 * This function fetches the value from the specified column and row in the CSV data,
 * converts it to an bool, and stores it in the provided reference variable.
 * 
 * Boolean conversions
 *      'y' = true
 *      'n' or any other value = false     
 * 
 * @param store Reference to an integer where the retrieved value will be stored.
 * @param column_name The name of the column from which to retrieve the value.
 * @param row The row index from which to retrieve the value. Default is 0.
 * @return true if the value was successfully retrieved and converted to an integer, false otherwise.
 */
template <>
bool CSVData<std::string>::getColumnValue(bool &store, const std::string &column_name, int row/*=0*/) {
    std::string value;
    getColumnValue(value, column_name, row);
    if (value == "") {
        store = false;
        return false;
    }
    if (value == "y") store = true;
    else store = false;
    return true;
}

/**
 * @brief Retrieves the value of a specified column in a given row and converts it to a double.
 * 
 * This function attempts to retrieve the value from the specified column and row in the CSV data.
 * If the value is successfully retrieved and is not an empty string, it converts the value to a double
 * and stores it in the provided reference. If the value is an empty string, it sets the reference to 0.0.
 * 
 * @param store Reference to a double where the retrieved value will be stored.
 * @param row The row number from which to retrieve the value.
 * @param column_name The name of the column from which to retrieve the value.
 * @return true if the value was successfully retrieved and converted to a double, false if the value was an empty string.
 */
template <>
bool CSVData<std::string>::getColumnValue(double &store, const std::string &column_name, int row/*=0*/) {
    std::string value;
    getColumnValue(value, column_name, row);
    if (value == "") {
        store = 0.0;
        return false;
    }
    store = toDouble(value);
    return true;
}

/**
 * @brief Sets the value of a specific column in a given row of the CSV data.
 * 
 * @param value The value to set in the specified column.
 * @param row The row index where the value should be set.
 * @param column_name The name of the column where the value should be set.
 * @return true if the value is successfully set.
 */
template <typename T>
bool CSVData<T>::setColumnValue(T value, int row, const std::string &column_name) {
    int col_index = header[column_name];
    data[row][col_index] = value;
    return true;
}

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
bool locateRanges(CSVData<std::string> &config_data, CSVData<std::string> &param_data) {
    // Set the parameter and nametable filenames
    std::string config_file_name = CONFIG_FILE_PATH;
    std::string param_file_name = PARAMETER_FILE_PATH;
    
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
 * CSVData object. If no file name is provided, the function will skip the loading process and treat years as independent.
 *
 * @param gs_data Reference to a CSVData object where the growing season data will be stored.
 * @param gs_file_name Reference to a string containing the name of the CSV file to be read.
 *
 * @note If the file cannot be opened, an error message will be printed to the console.
 * @note If no file name is provided, a message indicating that the load is being skipped will be printed to the console.
 */
void readGSSheet(CSVData<double> &gs_data, std::string &gs_file_name) {

    gs_data.empty();

    if (gs_file_name != "") { // only if we actually selected a file, and we're using the GS data files in the first place

        std::cout << std::endl;
        std::cout << "Reading growing season limits." << std::endl;

        std::ifstream data_file(gs_file_name);
        if (!data_file.is_open())
            std::cout << "FAILED to open GS RANGE file " << gs_file_name << std::endl;

        std::cout << "Opened GS RANGE file " << gs_file_name << std::endl;

        gs_data = CSVData<double>(gs_file_name);
        
        std::cout << "Finished GS RANGE read." << std::endl;
        std::cout << std::endl;

    } else {

        std::cout << "Skipping GS RANGE load -- Years will be treated as independent." << std::endl;
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
void readDataSheet(CSVData<double> &data, CSVData<double> &sum_data, std::string &data_file_name, std::string &header_file_name, std::string &sum_header_file_name) {

    data = CSVData<double>(data_file_name, header_file_name);
    sum_data = CSVData<double>(data_file_name, sum_header_file_name);

    std::cout << "Finished climate DATA read." << std::endl;
}

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
}

void readSiteAreaValues() {
}