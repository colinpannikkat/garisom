#include "11Parameters.h"

Parameters::Parameters() {

    /* Reserve some initial capacity for vectors to avoid frequent reallocations */
    gs_ar_year.reserve(100);
    gs_ar_start.reserve(100);
    gs_ar_end.reserve(100);
    gs_ar_ppm.reserve(100);
    // layer_depth.reserve(100);
    // vert_distance.reserve(100);
    // depth.reserve(100);
    // radius.reserve(100);
    // length.reserve(100);
}

double& Parameters::getModelParam(const std::string &param_name) {
    if (model_parameters.find(param_name) != model_parameters.end()) {
        return model_parameters[param_name];
    } else {
        fprintf(stderr, "Column name %s does not exist\n", param_name.c_str());
    }
    throw std::invalid_argument("Column name does not exist");
}

void Parameters::setModelParam(double value, const std::string &param_name) {
    model_parameters[param_name] = value;
}

int& Parameters::getGsArYear(int index) {
    if (index < 0 || index >= gs_ar_year.size()) {
        throw std::out_of_range("Index out of range");
    }
    return gs_ar_year[index];
}

void Parameters::setGsArYear(int index, int value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= gs_ar_year.size()) {
        gs_ar_year.resize(index + 1);
    }
    gs_ar_year[index] = value;
}

int& Parameters::getGsArStart(int index) {
    if (index < 0 || index >= gs_ar_start.size()) {
        throw std::out_of_range("Index out of range");
    }
    return gs_ar_start[index];
}

void Parameters::setGsArStart(int index, int value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= gs_ar_start.size()) {
        gs_ar_start.resize(index + 1);
    }
    gs_ar_start[index] = value;
}

int& Parameters::getGsArEnd(int index) {
    if (index < 0 || index >= gs_ar_end.size()) {
        throw std::out_of_range("Index out of range");
    }
    return gs_ar_end[index];
}

void Parameters::setGsArEnd(int index, int value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= gs_ar_end.size()) {
        gs_ar_end.resize(index + 1);
    }
    gs_ar_end[index] = value;
}

double& Parameters::getGsArPpm(int index) {
    if (index < 0 || index >= gs_ar_ppm.size()) {
        throw std::out_of_range("Index out of range");
    }
    return gs_ar_ppm[index];
}

void Parameters::setGsArPpm(int index, double value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= gs_ar_ppm.size()) {
        gs_ar_ppm.resize(index + 1);
    }
    gs_ar_ppm[index] = value;
}

// double& Parameters::getLayerDepth(int index) {
//     if (index < 0 || index >= layer_depth.size()) {
//         throw std::out_of_range("Index out of range");
//     }
//     return layer_depth[index];
// }

// void Parameters::setLayerDepth(int index, double value) {
//     if (index < 0) {
//         throw std::out_of_range("Index out of range");
//     }
//     if (index >= layer_depth.size()) {
//         layer_depth.resize(index + 1);
//     }
//     layer_depth[index] = value;
// }

// double& Parameters::getVertDistance(int index) {
//     if (index < 0 || index >= vert_distance.size()) {
//         throw std::out_of_range("Index out of range");
//     }
//     return vert_distance[index];
// }

// void Parameters::setVertDistance(int index, double value) {
//     if (index < 0) {
//         throw std::out_of_range("Index out of range");
//     }
//     if (index >= vert_distance.size()) {
//         vert_distance.resize(index + 1);
//     }
//     vert_distance[index] = value;
// }

// double& Parameters::getDepth(int index) {
//     if (index < 0 || index >= depth.size()) {
//         throw std::out_of_range("Index out of range");
//     }
//     return depth[index];
// }

// void Parameters::setDepth(int index, double value) {
//     if (index < 0) {
//         throw std::out_of_range("Index out of range");
//     }
//     if (index >= depth.size()) {
//         depth.resize(index + 1);
//     }
//     depth[index] = value;
// }

// double& Parameters::getRadius(int index) {
//     if (index < 0 || index >= radius.size()) {
//         throw std::out_of_range("Index out of range");
//     }
//     return radius[index];
// }

// void Parameters::setRadius(int index, double value) {
//     if (index < 0) {
//         throw std::out_of_range("Index out of range");
//     }
//     if (index >= radius.size()) {
//         radius.resize(index + 1);
//     }
//     radius[index] = value;
// }

// double& Parameters::getLength(int index) {
//     if (index < 0 || index >= length.size()) {
//         throw std::out_of_range("Index out of range");
//     }
//     return length[index];
// }

// void Parameters::setLength(int index, double value) {
//     if (index < 0) {
//         throw std::out_of_range("Index out of range");
//     }
//     if (index >= length.size()) {
//         length.resize(index + 1);
//     }
//     length[index] = value;
// }