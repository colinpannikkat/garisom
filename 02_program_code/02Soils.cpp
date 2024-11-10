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

/**
 * @brief Calculates the van Genuchten (vg) function for a given pressure.
 *
 * This function computes the van Genuchten function, which is commonly used 
 * in soil science to describe the relationship between soil water content 
 * and soil water potential (pressure). The function uses the van Genuchten 
 * parameters (alpha and n) and the maximum hydraulic conductivity (k_max) 
 * to perform the calculation.
 *
 * @param pressure The soil water potential (pressure) for which the vg function is calculated.
 * @return The calculated value of the van Genuchten function for the given pressure.
 */
double RhizosphereComponent::vg(double pressure) {
    double vp = 1 / (pow((van_gen_alpha * pressure), van_gen_n) + 1);
    return k_max * pow(vp, ((van_gen_n - 1) / (2 * van_gen_n))) * pow((pow((1 - vp), ((van_gen_n - 1) / van_gen_n)) - 1), 2);
}

double RhizosphereComponent::rvg(double pressure) { //'gives soil Y in MPa from soil theta/thetasat=x
    double aa = pow((pow(pressure, (1 / (1 - 1 / van_gen_n))) + 1), (1 / van_gen_n));
    double bb = (pow(pressure, (1 / (van_gen_n - 1))) * van_gen_alpha);
    return aa / bb;
}

double RhizosphereComponent::swc(double pressure) { //'gives soil water content, theta/thetasat from soil water potential in MPa=x
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((van_gen_alpha * pressure), van_gen_n))), (1 - 1 / van_gen_n));
}

void RhizosphereComponent::trapzd(const double &p1, const double &p2, double &s, const int &t, int &it) { //integrates root element z weibull

    double sum = 0, x = 0, del = 0;

    if (t == 0)
    {
        s = 0.5 * (p2 - p1) * (vg(p1) + vg(p2));
        it = 1;
    }
    else
    {
        del = (p2 - p1) / it;
        x = p1 + 0.5 * del;
        sum = 0;
        for (int j = 0; j < it; j++)
        {
            sum = sum + vg(x);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / it);
        it = 2 * it;
    }
}

void RhizosphereComponent::qtrap(double &p1, double &p2, double &s) { //'evaluates accuracy of root element z integration
    int it = 0;    // keeping track of iterations
    double olds = -1; //'starting point unlikely to satisfy if statement below
    for (int t = 0; t < TRAP_ITER_MAX; t++)
    {
        trapzd(p1, p2, s, t, it);
        if (std::abs(s - olds) < (0.001 * std::abs(olds))) // changing EPSX to 0.001 is the only reason for redeclaration
            return;
        olds = s;
    }
}

void RhizosphereComponent::calc_flow_rate(const double &p_inc, const double &k_min) {
    memset(e_p, 0, sizeof(e_p));
    memset(k, 0, sizeof(k));

    double p1 = 0, p2 = 0, s = 0, e = 0;
    int i = 1;
    this->e_p[0] = 0;
    this->k[0] = k_max;
    do {
        p2 = p1 + p_inc;
        qtrap(p1, p2, s);
        e += s;
        this->e_p[i] = e;
        this->k[i] = vg(p2); //weibull k
        p1 = p2; //reset p1 for next increment
        i += 1;
        if (i == 100000)
            break;
    } while (!(this->k[i - 1] < k_min));
    this->p_crit = p2;
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

