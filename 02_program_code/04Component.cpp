#include "04Component.h"

/* Getters */
double Component::getBwb() { return b_wb; }
double Component::getCwb() { return c_wb; }
double Component::getResPercent() { return res_percent; }
double Component::getKmax() { return k_max; }
double Component::getPcrit() { return p_crit; }
double Component::getKmin() { return k_min; }
double Component::getPressure() { return pressure; }

double& Component::getFatigue(int index) {
    if (index < 0 || index >= wb_fatigue.size()) {
        throw std::out_of_range("Index out of range");
    }
    return wb_fatigue[index];
}
double& Component::getPressureComp(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return pressure_comp[index];
}
double& Component::getPressureVirgin(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return pressure_v[index];
}

double& Component::getEp(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return e_p[index];
}

double& Component::getEpVirgin(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return e_pv[index];
}

double& Component::getEComp(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return e_comp[index];
}

double& Component::getK(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return k[index];
}

double& Component::getKVirgin(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return k_v[index];
}

double& Component::getKComp(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return k_comp[index];
}

/* Setters */
void Component::setBwb(double value) { this->b_wb = value; }
void Component::setCwb(double value) { this->c_wb = value; }
void Component::setResPercent(double value) { this->res_percent = value; }
void Component::setKmax(double value) { this->k_max = value; }
void Component::setPcrit(double value) { this->p_crit = value; }
void Component::setKmin(double value) { this->k_min = value; }
void Component::setPressure(double value) { this->pressure = value; }
void Component::setFatigue(int index, double value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= wb_fatigue.size()) {
        wb_fatigue.resize(index + 1);
    }
    wb_fatigue[index] = value;
}
void Component::setPressureComp(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    pressure_comp[index] = value;
}

void Component::setPressureVirgin(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    pressure_v[index] = value;
}

void Component::setEp(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    e_p[index] = value;
}

void Component::setEpVirgin(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    e_pv[index] = value;
}

void Component::setEComp(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    e_comp[index] = value;
}

void Component::setK(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    k[index] = value;
}

void Component::setKVirgin(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    k_v[index] = value;
}

void Component::setKComp(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    k_comp[index] = value;
}

/* Storage functions */

/* Stores transpiration curve */
void Component::storeTranspirationCurve() {
    int i = 0;
    do
    {
        e_pt[i] = e_p[i];

        i = i + 1;
        if (i == (CURVE_MAX - 1))
            break;
    } while (e_p[i] != 0 || e_pt[i] != 0);
}

/* Stores conductivity curve */
void Component::storeConductivityCurve() {
    int i = 0;
    do
    {
        k_t[i] = k[i];

        i++;
        if (i == (CURVE_MAX - 1))
            break;
    } while (k[i] != 0 || k_t[i] != 0);
}

/* Stores both conductivity and transpiration curves */
void Component::storeCurves() {
    int i = 0;
    do
    {
        e_pt[i] = e_p[i];
        k_t[i] = k[i];

        i++;
        if (i == (CURVE_MAX - 1))
            break;
    } while (k[i] != 0 || e_p[i] != 0 || e_pt[i] != 0 || k_t[i] != 0);
}

/* Stores transpiration curve and uses virgin curve */
void Component::storeTranspirationCurveAndUseVirgin() {
    int i = 0;
    do
    {
        e_pt[i] = e_p[i];
        e_p[i] = e_pv[i];

        i = i + 1;
        if (i == (CURVE_MAX - 1))
            break;
    } while (e_p[i] != 0 || e_pv[i] != 0 || e_pt[i] != 0); // I don't believe these checks to be right
}

/* Stores conductivity curve and uses virgin curve */
void Component::storeConductivityCurveAndUseVirgin() {
    int i = 0;
    do
    {
        k_t[i] = k[i];
        k[i] = k_v[i];

        i++;
        if (i == (CURVE_MAX - 1))
            break;
    } while (k[i] != 0 || k_v[i] != 0 || k_t[i] != 0); // I don't believe these checks to be right
}

/* Stores both conductivity and transpiration curves and uses virgin */
void Component::storeCurvesAndUseVirgin() {
    int i = 0;
    do
    {
        e_pt[i] = e_p[i];
        k_t[i] = k[i];

        e_p[i] = e_pv[i];
        k[i] = k_v[i];

        i++;
        if (i == (CURVE_MAX - 1))
            break;
    } while (e_p[i] != 0 ||  e_pv[i] != 0 || e_pt[i] != 0); // I don't believe these checks to be right
}

/* Restores transpiration curve from historical */
void Component::restoreTranspirationCurve() {
    int i = 0;
    do
    {
        e_p[i] = e_pt[i];

        i = i + 1;
        if (i == (CURVE_MAX - 1))
            break;
    } while (e_pt[i] != 0 || e_p[i] != 0);
}

/* Restores conductivity curve from historical */
void Component::restoreConductivityCurve() {
    int i = 0;
    do
    {
        k[i] = k_t[i];

        i = i + 1;
        if (i == (CURVE_MAX - 1))
            break;
    } while (k_t[i] != 0 || k[i] != 0);
}

/* Restores conductivity and transpiration curve from historical */
void Component::restoreCurves() {
    int i = 0;
    do
    {
        k[i] = k_t[i];
        e_p[i] = e_pt[i];

        i = i + 1;
        if (i == (CURVE_MAX - 1))
            break;
    } while (e_pt[i] != 0 || e_p[i] != 0); // I believe these checks are incomplete
}

/* Clears historical curve storage */
void Component::clearHistoricalCurves() {
    memset(this->e_pt, 0, sizeof(this->e_pt));
    memset(this->k_t, 0, sizeof(this->k_t));
}

/* Clean parameters */
void Component::cleanParameters() {
    b_wb = 0.0;
    c_wb = 0.0;
    res_percent = 0.0;
    k_max = 0.0;
    p_crit = 0.0;
    k_min = 0.0;
    pressure = 0.0;
    wb_fatigue.clear();
    std::fill(std::begin(e_p), std::end(e_p), 0.0);
    std::fill(std::begin(e_pv), std::end(e_pv), 0.0);
    std::fill(std::begin(e_comp), std::end(e_comp), 0.0);
    std::fill(std::begin(e_pt), std::end(e_pt), 0.0);
    std::fill(std::begin(k), std::end(k), 0.0);
    std::fill(std::begin(k_v), std::end(k_v), 0.0);
    std::fill(std::begin(k_comp), std::end(k_comp), 0.0);
    std::fill(std::begin(k_t), std::end(k_t), 0.0);
    std::fill(std::begin(pressure_comp), std::end(pressure_comp), 0.0);
    std::fill(std::begin(pressure_v), std::end(pressure_v), 0.0);
}

/* Hydraulic functions */

/**
 * @brief Calculates the conductance of a component for a given pressure.
 * 
 * This function computes the conductance value using a Weibull function that 
 * incorporates the maximum coefficient (k_max), and fitted Weibull parameters 
 * b_wb and c_wb. The formula is as follows:
 * 
 * wb = k_max * exp(-(pow((pressure / b_wb), c_wb)))
 * 
 * @param pressure The input pressure value (-MPa)
 * @return The computed hydraulic conductance value as a double.
 */
double Component::wb(const double &pressure) {
    return this->k_max * exp(-(pow((pressure / this->b_wb), this->c_wb)));
}

/**
 * @brief Performs numerical integration using the trapezoidal rule.
 * 
 * This function calculates the integral of component conductivity over the 
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
void Component::trapzd(const double &p1, const double &p2, double &s, const int &t, int &it) {
    double sum = 0, x = 0, del = 0;

    if (t == 0)
    {
        s = 0.5 * (p2 - p1) * (wb(p1) + wb(p2));
        it = 1;
    }
    else
    {
        del = (p2 - p1) / it;
        x = p1 + 0.5 * del;
        sum = 0;
        for (int j = 0; j < it; j++)
        {
            sum = sum + wb(x);
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
void Component::qtrap(double &p1, double &p2, double &s) { //'evaluates accuracy of root element z integration
    int it = 0;    // keeping track of iterations
    double olds = -1; //'starting point unlikely to satisfy if statement below
    for (int t = 0; t < TRAP_ITER_MAX; t++)
    {
        trapzd(p1, p2, s, t, it);
        if (std::abs(s - olds) <= (XYLEM_EPSX * std::abs(olds)))
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
 * @param virgin Indicates whether to fill virgin conductivity and transpiration.
 * 
 * @details
 * - The function initializes the cumulative flow (`e_p`) and conductivity (`k`) arrays to zero.
 * - It uses a loop to incrementally calculate the flow and update the conductivity using the `vg` function.
 * - The integration is performed using the `qtrap` function over the pressure range.
 * - The loop terminates when the conductivity falls below `k_min` or after 100,000 iterations.
 * - The critical pressure (`p_crit`) is set to the last pressure value (`p2`) when the loop exits.
 * - The virgin curves are not updated by updatecurves().
 * 
 * @note
 * - The `vg` Van Genuchten function is assumed to calculate the conductivity.
 * - The `qtrap` function is assumed to perform numerical integration over the pressure range.
 * - The maximum number of iterations is hardcoded to 100,000 to prevent infinite loops.
 */
void Component::calc_flow_rate(const double &p_inc, const double &k_min, bool virgin) {
    double *e_ptr, *k_ptr;
    if (virgin) {
        e_ptr = this->e_pv;
        k_ptr = this->k_v;
    } else {
        e_ptr = this->e_p;
        k_ptr = this->k;
    }

    memset(e_ptr, 0, sizeof(*e_ptr) * CURVE_MAX);

    double p1 = 0, p2 = 0, s = 0, e = 0;
    int i = 1;
    e_ptr[0] = 0;
    k_ptr[0] = k_max;
    do {
        p2 = p1 + p_inc;
        qtrap(p1, p2, s);
        e += s;
        e_ptr[i] = e;
        k_ptr[i] = wb(p2); //weibull k
        p1 = p2; //reset p1 for next increment

        i += 1;
        if (i == FLOW_ITER_LIMIT)
            break;
    } while (!(k_ptr[i - 1] < k_min));
    this->p_crit = p2;
}

/**
 * @brief Calculates through-flow and interpolates lower and upper limits.
 * 
 * Computes flow between pressures (p1, p2) and interpolates `klower` and 
 * `kupper` values, representing K(P) at lower and upper integration limits.
 * 
 * @param p1 Starting pressure value.
 * @param p2 Ending pressure value.
 * @param p_inc Pressure increment for indexing.
 * @param flow Reference to store calculated flow (difference in e values).
 * @param klower Reference to store interpolated K(P) at lower limit.
 * @param kupper Reference to store interpolated K(P) at upper limit.
 * 
 * @details
 * - Uses linear interpolation for intermediate `e` and `k` values.
 * - `e_p` stores `e` values for pressure indices.
 * - `k` stores K(P) values for pressure indices.
 * - Flow is the difference between interpolated `e` values at `p2` and `p1`.
 */
void Component::calc_through_flow(const double &p1, 
                                  const double &p2, 
                                  const double &p_inc,
                                  double &flow,
                                  double &klower,
                                  double &kupper) {
    double plow = int(p1 / p_inc); //'pressure index below target
    int i = int(p1 / p_inc);
    double elow = e_p[i]; //'e below target
    double klow = k[i];
    double ehigh = e_p[i + 1]; //'e above target
    double khigh = k[i + 1];
    plow = plow * p_inc; //'convert index to pressure
    double estart = (p1 - plow) / p_inc * (ehigh - elow) + elow; //'linear interpolation of starting e
    klower = (p1 - plow) / p_inc * (khigh - klow) + klow; //'linear interpolation of K(P)at lower limit of integration
    plow = int(p2 / p_inc); //'pressure index below target
    i = int(p2 / p_inc);
    elow = e_p[i]; //'e below target
    klow = k[i];
    ehigh = e_p[i + 1]; //'e above target
    khigh = k[i + 1];
    plow = plow * p_inc; //'convert index to pressure
    double efinish = (p2 - plow) / p_inc * (ehigh - elow) + elow; //'linear interpolation of finishing e
    kupper = (p2 - plow) / p_inc * (khigh - klow) + klow; //'linear interpolation of K(P) at upper limit of flow integration
    flow = efinish - estart; //'e upstream flow
}

/**
 * @brief Calculates the downstream component pressure based on the provided 
 * transpiration rate.
 * 
 * This function computes the downstream pressure of a component given the 
 * transpiration rate `e`, the bottom pressure, gravitational pressure drop, 
 * and pressure increment. It uses linear interpolation and iterates through 
 * the pressure array `e_p` to find the appropriate pressure corresponding to 
 * the given transpiration rate.
 * 
 * @param e The transpiration rate (rate of fluid flow).
 * @param bottom_pressure The pressure at the bottom of the component.
 * @param p_grav The gravitational pressure drop.
 * @param p_inc The pressure increment used for indexing.
 * @return int Returns 0 if the calculation is successful and the pressure is 
 *             below the critical pressure (`p_crit`). Returns 1 if the 
 *             pressure exceeds the critical pressure or if the end of the 
 *             `e_p` array is reached.
 * 
 * @note The function assumes that the `e_p` array contains valid pressure 
 *       values and that the array is properly indexed. If `e_p[j]` is zero, 
 *       the function terminates early with a return value of 1.
 */
int Component::calc_pressure(const double &e, const double &bottom_pressure, const double &p_grav, const double &p_inc) {

    //'start with bottom pressure
    double p1 = bottom_pressure + p_grav; //add gravity drop before integration
    double plow = int(p1 / p_inc); //pressure index below target
    int i = int(p1 / p_inc);
    double elow = e_p[i]; //e below target
    double ehigh = e_p[i + 1]; //e above target
    plow = plow * p_inc; //convert index to pressure
    double estart = (p1 - plow) / p_inc * (ehigh - elow) + elow; //linear interpolation of starting e
    double efinish = estart + e;
    int j = i;
    do //find efinish
    {
        j = j + 1;
        if ((j == FLOW_ITER_LIMIT) || (e_p[j] == 0)) {
            return 1;
        }
    } while (!(e_p[j] > efinish));
    //Loop Until es(j) > efinish
    ehigh = e_p[j];
    elow = e_p[j - 1];
    double p2 = ((efinish - elow) / (ehigh - elow)) * p_inc + p_inc * (j - 1);
    this->pressure = p2; // p2 is downstream component pressure at the provided transpiration rate, e
    if (this->pressure >= this->p_crit) {
        return 1;
    }
    return 0;
}