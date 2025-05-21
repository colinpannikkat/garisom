#ifndef UTILITIES_H
#define UTILITIES_H

/* Never let user compile with ffast-math */
#ifdef __FAST_MATH__
#error "-ffast-math is broken, don't use it"
#endif

#include <stdio.h>
#include <cmath>    // math utility functions
#include <string>   // the C++ String Class, easier to deal with than char arrays for this application
#include <ctime>    // timers for performance testing
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <assert.h>
#include <iomanip>
#include <algorithm>

// Input File Names
#define CONFIG_FILE_PATH "../03_test_data/configuration.csv"
#define PARAMETER_FILE_PATH "../03_test_data/parameters.csv"
#define DATA_HEADER_FILE_PATH "./dataheader.csv"
#define OUT_DIR "."

// Config Setting (used for determining which row of configuration file to use)
#define CONFIG_SETTING 0

// Output Precision
#define FIO_PRECISION 12

// staging
#define STAGE_ID_NONE 0 // use this mode to skip all of the stage code and just run based on param sheet settings (a "normal" run, the default config)
#define STAGE_ID_HIST_OPT 1
#define STAGE_ID_HIST_STRESS 2
#define STAGE_ID_FUT_OPT 3
#define STAGE_ID_FUT_STRESS 4
#define STAGE_ID_FUT_STRESS_NOACCLIM 5

// Data Arrays
#define MAX_SUMMARY_COLS 121
#define MAX_SUMMARY_ROWS 2001
#define PARAMFILE_MAXROWS 2001
#define PARAMFILE_MAXCOLS 101
#define DATAFILE_MAXROWS 2000001
#define DATAFILE_MAXCOLS 101
#define CONFIGFILE_MAXROWS 10
#define CONFIGFILE_MAXCOLS 101
#define CURVE_MAX 100001

// Constants
#define PROFT_MAX_RUN_MEAN 1        // running mean for profit maximization
#define DPA_MAX_CUTOFF 1.1          // cutoff for stopping dpamax search
#define MIN_WIND_THRESH 0.4515      // m s-1'minimum wind threshold
#define TRAP_ITER_MAX 70            // itmax for trapzd routine
#define FLOW_ITER_LIMIT 100000      // itmax for calc_flow_rate
#define XYLEM_EPSX 0.0001           // acceptable error for e integral for xylem
#define RHIZO_EPSX 0.001            // acceptable error for e integral for rhizosphere
#define GMAX 1000000                // maximum G, wet soil, vpd = 0, in kg m - 2 hr - 1 basal area
#define PI 3.14159
#define SBC 0.0000000567            // Stefan-Boltzmann constant in W m-2 K-4
#define SHA 29.3                    // Specific heat of air in J mol-1C-1
#define GAS 8.3144598               // Universal gas constant J mol-1K-1
#define OA 0.20999999999999999                     // Mole fraction of O2
#define SOLAR 1362                  // Solar constant W m-2
#define ABS_SOLAR 0.5               // Absorptivity of solar for leaves
#define ABS_PAR 0.8                 // Absorptivity of PAR for leaves
#define ABS_NIR 0.2                 // Absorptivity of near infrared for leaves
#define MAX_YEARS 90

// Iteration Limits

/**
 * @class CSVData
 * @brief A class to handle CSV data, including reading headers and data, and accessing and setting column values.
 * 
 * This class provides functionalities to read CSV files, store the data in a structured format,
 * and access or modify the data using column names and row indices. Can only instantiate with <std::string>
 * or <double> types.
 * 
 * @tparam T The type of data stored in the CSV.
 */
template <typename T>
class CSVData {
    private:
        std::unordered_map<std::string, int> header;
        std::vector<std::vector<T>> data;

        int num_cols;

    public:
        CSVData() {num_cols = 0;};
        CSVData(const std::string &data_file_name);
        CSVData(const std::string &data_file_name, const std::string &header_file_name, bool skip_first=true);

        void readHeader(std::ifstream &dataFile);
        void readData(std::ifstream &dataFile);

        bool empty();
        void clear();
        void print(size_t page_size);
        int row_size();
        int col_size();
        void output(std::string out_file_name);
        
        // Generic getter, used when accessing numerical data
        T& getColumnValue(const std::string& column_name, int row=0);
        T& operator()(const std::string &column_name, int row=0);

        // Overloaded getter for when accessing data that is stored in string format and requires conversion
        // Typically required if class is instantiated with <std::string> type
        bool getColumnValue(std::string &store, const std::string &column_name, int row=0);
        bool getColumnValue(double &store, const std::string &column_name, int row=0);
        bool getColumnValue(bool &store, const std::string& column_name, int row=0);
        bool getColumnValue(int &store, const std::string &column_name, int row=0);

        // Generic setter, will set as template type, therein it is important that the type
        // of data matches the CSVData type
        bool setColumnValue(T value, int row, const std::string &column_name);
};

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
inline int toInt(const std::string& str) {
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
inline double toDouble(const std::string& str) {
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
inline std::string trim(const std::string& str) {
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
            num_cols++;
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
inline void CSVData<std::string>::readData(std::ifstream &data_file) {
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
        if (row.size() < num_cols) {
            while (row.size() != num_cols) {
                row.push_back("");
            }
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
inline void CSVData<double>::readData(std::ifstream &data_file) {
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
        if (row.size() < num_cols) {
            while (row.size() != num_cols) {
                row.push_back(0.0);
            }
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

    num_cols = 0;
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
    
    std::ifstream header_file(header_file_name);
    this->num_cols = 0;

    if (!header_file.is_open()) {
        std::cout << "UNRECOVERABLE: Failed to open file " << header_file_name << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    readHeader(header_file);

    if (!data_file_name.empty()) {

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

        readData(data_file);
    }
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
inline bool CSVData<T>::empty() {
    return data.empty();
}

template<typename T>
inline void CSVData<T>::clear() {
    data.clear();
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
inline void CSVData<T>::print(size_t page_size) {
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
inline void CSVData<double>::print(size_t page_size) {
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
template <typename T>
inline void CSVData<T>::output(std::string out_file_name) {
    std::ofstream out_file(out_file_name);
    if (!out_file.is_open()) {
        std::cerr << "Failed to open output file: " << out_file_name << std::endl;
        return;
    }

    // Write header, in unordered set so must order
    std::vector<std::pair<std::string, int>> ordered_header(header.begin(), header.end());
    std::sort(ordered_header.begin(), ordered_header.end(), [](const auto &a, const auto &b) {
        return a.second < b.second;
    });
    for (const auto &pair : ordered_header) {
        out_file << pair.first << ",";
    }
    out_file.seekp(-1, std::ios_base::cur); // Remove the last comma
    out_file << "\n";

    // Write data
    for (const auto& row : data) {
        for (const auto& cell : row) {
            if constexpr (std::is_same_v<T, double>) {
                if (std::abs(cell) < std::pow(10, -FIO_PRECISION)) {
                    out_file << 0 << ",";
                } else {
                    if (cell == static_cast<int>(cell)) {
                        out_file << static_cast<int>(cell) << ",";
                    } else {
                        out_file << std::fixed << std::setprecision(FIO_PRECISION) << cell << ",";
                    }
                }
            } else {
                out_file << cell << ",";
            }
        }
        out_file.seekp(-1, std::ios_base::cur); // Remove the last comma
        out_file << "\n";
    }

    std::cout << "Succesfully wrote to: " << out_file_name << std::endl;

    out_file.close();
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
inline int CSVData<T>::row_size() {
    return data.size();
}

/**
 * @brief Returns the number of columns in the CSV data.
 * 
 * This function retrieves the size of the first row in the CSV data,
 * which corresponds to the number of columns in the dataset.
 * 
 * @tparam T The type of data stored in the CSV.
 * @return int The number of columns in the CSV data.
 */
template <typename T>
inline int CSVData<T>::col_size() {
    return num_cols;
}

template <typename T>
T& CSVData<T>::operator()(const std::string &column_name, int row) {
    return getColumnValue(column_name, row);
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
T& CSVData<T>::getColumnValue(const std::string &column_name, int row/*=0*/) {
    assert(row < data.size());
    if (header.find(column_name) != header.end()) {
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        return data[row][col_index];
    } else {
        // fprintf(stderr, "Column name %s does not exist, creating...\n", column_name.c_str());
        setColumnValue(T(), row, column_name);
        int col_index = header[column_name];
        assert(col_index < data[row].size());
        return data[row][col_index];
    }
}

// /**
//  * @brief Retrieves the value of a specified column in a given row.
//  * 
//  * This function fetches the value from the CSV data for the specified row and column name.
//  * It asserts that the row index is within the bounds of the data and that the column name exists.
//  * If the column name does not exist, an error message is printed to stderr.
//  * 
//  * @param column_name The name of the column from which to retrieve the value.
//  * @param row The index of the row from which to retrieve the value. Default is 0.
//  * @return std::string The value of the specified column in the given row. 
//  *         Returns an empty string if the column name does not exist.
//  */
// template <>
// inline std::string CSVData<std::string>::getColumnValue(const std::string &column_name, int row/*=0*/) {
//     assert(row < data.size());
//     if (header.find(column_name) != header.end()) {
//         int col_index = header[column_name];
//         assert(col_index < data[row].size());
//         return data[row][col_index];
//     } else {
//         fprintf(stderr, "Column name %s does not exist\n", column_name.c_str());
//     }
//     return "";
// }

// /**
//  * @brief Retrieves the value of a specified column in a given row.
//  * 
//  * This function fetches the value from the CSV data for the specified row and column name.
//  * It asserts that the row index is within the bounds of the data and that the column name exists.
//  * If the column name does not exist, an error message is printed to stderr.
//  * 
//  * @param column_name The name of the column from which to retrieve the value.
//  * @param row The index of the row from which to retrieve the value. Default is 0.
//  * @return std::string The value of the specified column in the given row. 
//  *         Returns an empty string if the column name does not exist.
//  */
// template <>
// inline double CSVData<double>::getColumnValue(const std::string &column_name, int row/*=0*/) {
//     assert(row < data.size());
//     if (header.find(column_name) != header.end()) {
//         int col_index = header[column_name];
//         assert(col_index < data[row].size());
//         return data[row][col_index];
//     } else {
//         fprintf(stderr, "Column name %s does not exist, creating...\n", column_name.c_str());
//         setColumnValue(0.0, row, column_name);
//     }
//     return 0.0;
// }

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
inline bool CSVData<std::string>::getColumnValue(std::string &store, const std::string &column_name, int row/*=0*/) {
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
inline bool CSVData<std::string>::getColumnValue(int &store, const std::string &column_name, int row/*=0*/) {
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
inline bool CSVData<std::string>::getColumnValue(bool &store, const std::string &column_name, int row/*=0*/) {
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
inline bool CSVData<std::string>::getColumnValue(double &store, const std::string &column_name, int row/*=0*/) {
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
inline bool CSVData<T>::setColumnValue(T value, int row, const std::string &column_name) {
    if (header.find(column_name) == header.end()) {
        // Add new column to header
        header[column_name] = num_cols++;
        // Add new column to each row in data
        for (auto& row_data : data) {
            row_data.push_back(T());
        }
    }
    int col_index = header[column_name];

    // Resize!
    if (row >= data.size()) {
        data.resize(row + 1, std::vector<T>(num_cols, T()));
    }
    data[row][col_index] = value;
    return true;
}

class Parameters; // moved header file include to .cpp, in future remove below functions as there is circular dependency

bool locateRanges(CSVData<std::string> &config_data, std::string config_data_file, CSVData<std::string> &param_data, std::string param_data_file);
void readGSSheet(CSVData<double> &gs_data, std::string &gs_file_name, bool use_gs_data);
void readDataSheet(CSVData<double> &data, std::string &data_file_name, std::string &header_file_name);
void readGrowSeasonData(Parameters &param, CSVData<double> &gs_data);
void readSiteAreaValues();

#endif