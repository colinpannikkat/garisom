#include "04Component.h"

/* Getters */
double Component::getBwb() { return b_wb; }
double Component::getCwb() { return c_wb; }
double Component::getResPercent() { return res_percent; }
double Component::getResWb() { return res_wb; }
double Component::getCondWb() { return cond_wb; }
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
void Component::setResWb(double value) { this->res_wb = value; }
void Component::setCondWb(double value) { this->cond_wb = value; }
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

double Component::wb(const double &pressure) {
    return this->k_max * exp(-(pow((pressure / this->b_wb), this->c_wb)));
}

void Component::trapzd(const double &p1, const double &p2, double &s, const int &t, int &it) { //integrates root element z weibull

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

void Component::qtrap(double &p1, double &p2, double &s) { //'evaluates accuracy of root element z integration
    int it = 0;    // keeping track of iterations
    double olds = -1; //'starting point unlikely to satisfy if statement below
    for (int t = 0; t < TRAP_ITER_MAX; t++)
    {
        trapzd(p1, p2, s, t, it);
        if (std::abs(s - olds) <= (EPSX * std::abs(olds)))
            return;
        olds = s;
    }
}

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
    memset(k_ptr, 0, sizeof(*k_ptr) * CURVE_MAX);

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
        if (i == 100000)
            break;
    } while (!(k_ptr[i - 1] < k_min));
    this->p_crit = p2;
}

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

/* Gets pressure of component based on previous component (bottom_pressure), does not apply to root and rhizosphere since these are solved via Newton Rhaphson approximation. */
/* Used iteratively over every e from zero to determine specific component pressure */
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
        if (e_p[j] == 0) {
            return 1;
        }
    } while (!(e_p[j] > efinish));
    //Loop Until es(j) > efinish
    ehigh = e_p[j];
    elow = e_p[j - 1];
    double p2 = (((efinish - elow) / (ehigh - elow)) * p_inc) + (p_inc * (j - 1));
    this->pressure = p2; // p2 is downstream component pressure at the provided transpiration rate, e
    if (this->pressure >= this->p_crit) {
        return 1;
    }
    return 0;
}