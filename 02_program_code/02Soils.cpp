// Functions to calculate soil parameters at start
#include "02Soils.h"

/* Getters and setters for the Rhizosphere component attributes. */

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
 * @brief Calculates the conductivity for a given pressure using the van 
 * Genuchten (vg) function.
 *
 * This function computes the van Genuchten function, which is used to determine
 * rhizosphere hydraulic conductance. The function uses the van Genuchten 
 * parameters (alpha and n) and the maximum hydraulic conductivity (k_max) 
 * to perform the calculation.
 *
 * @param pressure The soil water potential (pressure) for which the vg function is calculated.
 * @return The rhizosphere hydraulic conductance at a certain pressure.
 */
double RhizosphereComponent::vg(double pressure) {
    double vp = 1 / (pow((van_gen_alpha * pressure), van_gen_n) + 1);
    return k_max * pow(vp, ((van_gen_n - 1) / (2 * van_gen_n))) * pow((pow((1 - vp), ((van_gen_n - 1) / van_gen_n)) - 1), 2);
}

/**
 * @brief Computes the soil water retention function (Y) in MPa based on the soil 
 *        water content ratio (theta/thetasat).
 * 
 * This function calculates the soil water potential using the van Genuchten 
 * model parameters. It takes the soil pressure as input and computes the 
 * corresponding soil water potential.
 * 
 * @param pressure The soil pressure (in MPa) used in the van Genuchten model.
 *                 This value represents the matric potential of the soil.
 * 
 * @return The computed soil water potential (Y) in MPa.
 * 
 * @note The function assumes that the van Genuchten parameters (van_gen_n and 
 *       van_gen_alpha) are properly initialized and valid for the soil type 
 *       being modeled.
 */
double RhizosphereComponent::rvg(double pressure) {
    double aa = pow((pow(pressure, (1 / (1 - 1 / van_gen_n))) + 1), (1 / van_gen_n));
    double bb = (pow(pressure, (1 / (van_gen_n - 1))) * van_gen_alpha);
    return aa / bb;
}

/**
 * @brief Calculates the soil water content (SWC) based on soil water potential.
 * 
 * This function computes the soil water content (theta/thetasat) using the 
 * van Genuchten equation, which relates soil water potential (pressure) to 
 * soil water content. The equation is parameterized by the van Genuchten 
 * parameters `van_gen_alpha` and `van_gen_n`.
 * 
 * @param pressure The soil water potential in MPa.
 * @return The soil water content as a fraction of saturation (theta/thetasat).
 */
double RhizosphereComponent::swc(double pressure) {
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((van_gen_alpha * pressure), van_gen_n))), (1 - 1 / van_gen_n));
}

/**
 * @brief Performs numerical integration using the trapezoidal rule.
 * 
 * This function calculates the integral of rhizosphere conductivity over the 
 * interval [p1, p2] to determine the steady state flow rate.
 * 
 * This function uses a trapezoidal approximation of the integral, refining the 
 * approximation with each call by doubling the number of subintervals.
 * 
 * @param p1 The lower bound of the integration interval.
 * @param p2 The upper bound of the integration interval.
 * @param s Reference to the variable where the computed integral value will be stored.
 * @param t An integer flag indicating whether this is the first call (t == 0) or a subsequent call.
 * @param it Reference to the number of subintervals used in the integration. This value is updated
 *           with each call to refine the approximation.
 */
void RhizosphereComponent::trapzd(const double &p1, const double &p2, double &s, const int &t, int &it) {
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

/**
 * @brief Evaluates the accuracy of root element z integration using the trapezoidal rule.
 * 
 * This function iteratively refines the integration result using the trapezoidal rule
 * until the desired accuracy is achieved or the maximum number of iterations is reached.
 * 
 * @param p1 Reference to the first parameter for the integration range.
 * @param p2 Reference to the second parameter for the integration range.
 * @param s Reference to the variable where the integration result is stored.
 * 
 * The function uses a helper function `trapzd` to perform the integral calculations.
 * The iteration stops when the difference between the current and previous integration
 * results is less than 0.1% of the previous result, or when the maximum number of iterations
 * (`TRAP_ITER_MAX`) is reached.
 */
void RhizosphereComponent::qtrap(double &p1, double &p2, double &s) {
    int it = 0;
    double olds = -1;
    for (int t = 0; t < TRAP_ITER_MAX; t++)
    {
        trapzd(p1, p2, s, t, it);
        if (std::abs(s - olds) < (RHIZO_EPSX * std::abs(olds)))
            return;
        olds = s;
    }
}

/**
 * @brief Calculates the flow rate in the rhizosphere component based on 
 * pressure increments and a minimum conductivity threshold.
 * 
 * This function computes the flow rate by iteratively integrating over pressure 
 * increments and updating the cumulative flow and conductivity values. The process 
 * stops when the conductivity drops below the specified minimum threshold or 
 * when the iteration limit is reached.
 * 
 * @param p_inc The pressure increment used for integration.
 * @param k_min The minimum conductivity threshold to terminate the calculation.
 * 
 * @details
 * - The function initializes the cumulative flow (`e_p`) and conductivity (`k`) arrays to zero.
 * - It uses a loop to incrementally calculate the flow and update the conductivity using the `vg` function.
 * - The integration is performed using the `qtrap` function over the pressure range.
 * - The loop terminates when the conductivity falls below `k_min` or after 100,000 iterations.
 * - The critical pressure (`p_crit`) is set to the last pressure value (`p2`) when the loop exits.
 * 
 * @note
 * - The `vg` Van Genuchten function is assumed to calculate the conductivity.
 * - The `qtrap` function is assumed to perform numerical integration over the pressure range.
 * - The maximum number of iterations is hardcoded to 100,000 to prevent infinite loops.
 */
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
        this->k[i] = vg(p2); // vg k
        p1 = p2; //reset p1 for next increment
        i += 1;
        if (i == FLOW_ITER_LIMIT)
            break;
    } while (!(this->k[i - 1] < k_min));
    this->p_crit = p2;
}