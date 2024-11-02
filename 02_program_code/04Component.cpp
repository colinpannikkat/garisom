#include "04Component.h"

/* Getters */
double Component::getBwb() { return b_wb; }
double Component::getCwb() { return c_wb; }
double Component::getResPercent() { return res_percent; }
double Component::getResWb() { return res_wb; }
double Component::getCondWb() { return cond_wb; }
double Component::getKmax() { return k_max; }
double& Component::getFatigue(int index) {
    if (index < 0 || index >= wb_fatigue.size()) {
        throw std::out_of_range("Index out of range");
    }
    return wb_fatigue[index];
}
double& Component::getPressure(int index) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    return pressure[index];
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

double& Component::getEcomp(int index) {
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
void Component::setFatigue(int index, double value) {
    if (index < 0) {
        throw std::out_of_range("Index out of range");
    }
    if (index >= wb_fatigue.size()) {
        wb_fatigue.resize(index + 1);
    }
    wb_fatigue[index] = value;
}
void Component::setPressure(int index, double value) {
    if (index < 0 || index >= CURVE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    pressure[index] = value;
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

void Component::setEcomp(int index, double value) {
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