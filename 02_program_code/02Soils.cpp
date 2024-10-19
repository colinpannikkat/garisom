// Functions to calculate soil parameters at start
#include "02Soils.h"

void RhizosphereComponent::setVanGenAlpha(double alpha) {
    this->van_gen_alpha = alpha;
}
double RhizosphereComponent::getVanGenAlpha() { return van_gen_alpha; }

void RhizosphereComponent::setVanGenN(double n) {
    this->van_gen_n = n;
}
double RhizosphereComponent::getVanGenN() { return van_gen_n; }

void RhizosphereComponent::setThetasat(double thetasat) {
    this->thetasat = thetasat;
}
double RhizosphereComponent::getThetaSat() { return thetasat; }

double RhizosphereComponent::vg(double pressure) {
    double vp = 1 / (pow((van_gen_alpha * pressure), van_gen_n) + 1);
    return this->getKmax() * pow(vp, ((van_gen_n - 1) / (2 * van_gen_n))) * pow((pow((1 - vp), ((van_gen_n - 1) / van_gen_n)) - 1), 2);
}

/*Get Van Genuchten alpha // override if provided*/
// This function obtains the Van Genuchten parameters for a given texture
void soils::get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[2001][101], long &rowLR, long &colLR) {
    std::string a, n, soilkmax, thetasat;
    if (texture == "sand"){
        a = "1479.5945"; 
        n = "2.68";
        soilkmax = "30305.88"; 
        thetasat = "0.43";
    } else if (texture == "loamy sand"){
        a = "1265.3084";
        n = "2.28";
        soilkmax = "14897.84";
        thetasat = "0.41";
    } else if (texture == "sandy loam") {
        a = "765.3075";
        n = "1.89";
        soilkmax = "4510.168";
        thetasat = "0.41";
    } else if (texture == "loam") {
        a = "367.3476";
        n = "1.56";
        soilkmax = "1061.216";
        thetasat = "0.43";
    } else if (texture == "silt") {
        a = "163.2656";
        n = "1.37";
        soilkmax = "255.1";
        thetasat = "0.46";
    } else if (texture == "silt loam") {
        a = "204.082";
        n = "1.41";
        soilkmax = "459.18";
        thetasat = "0.45";
    } else if (texture == "sandy clay loam") {
        a = "602.0419";
        n = "1.48";
        soilkmax = "1336.724";
        thetasat = "0.39";
    } else if (texture == "clay loam") {
        a = "193.8779";
        n = "1.31";
        soilkmax = "265.304";
        thetasat = "0.41";
    } else if (texture == "silty clay loam") {
        a = "102.041";
        n = "1.23";
        soilkmax = "71.428";
        thetasat = "0.43";
    } else if (texture == "sandy clay") {
        a = "275.5107";
        n = "1.23";
        soilkmax = "122.448";
        thetasat = "0.38";
    } else if (texture == "silty clay") {
        a = "51.0205";
        n = "1.09";
        soilkmax = "20.408";
        thetasat = "0.36";
    } else if (texture == "clay") {
        a = "81.6328";
        n = "1.09";
        soilkmax = "204.08";
        thetasat = "0.38";
    } else {
        std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
        std::cout << "SOURCE: Incorrect soil texture category" << std::endl;
        std::cout << "ACTION: Model stops " << std::endl;
        std::cout << std::endl;
        abort();
    }
    for(long k=1;k<=layers;k++){
        soillayersTable[rowLR + k][colLR + 3] = a; 
        soillayersTable[rowLR + k][colLR + 4] = n;
        soillayersTable[rowLR + k][colLR + 5] = soilkmax;
        soillayersTable[rowLR + k][colLR + 6] = thetasat;
    }
}

