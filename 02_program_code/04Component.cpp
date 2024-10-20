#include "04Component.h"

double Component::wb(const double &pressure) {
    return k_max * exp(-(pow((pressure / b_wb), c_wb)));
}

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

void Component::trapzdwb(const double &p1, const double &p2, double &s, const int &t, int &it) { //integrates root element z weibull

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

void Component::qtrapwb(double &p1, double &p2, double &s) { //'evaluates accuracy of root element z integration
    int it = 0;    // keeping track of iterations
    double olds = -1; //'starting point unlikely to satisfy if statement below
    for (int t = 0; t < TRAP_ITER_MAX; t++)
    {
        trapzdwb(p1, p2, s, t, it);
        if (std::abs(s - olds) <= (EPSX * std::abs(olds)))
            return;
        olds = s;
    }
}

void Component::calc_flow_rate(const double &p_inc, const double &k_min) {
    memset(e_p, 0, sizeof(e_p));

    double p1 = 0, p2 = 0, s = 0, e = 0, k = 0;
    int i = 1;
    e_p[0] = 0;
    do {
        p2 = p1 + p_inc;
        qtrapwb(p1, p2, s);
        e += s;
        e_p[i] = e;
        k = wb(p2); //weibull k
        p1 = p2; //reset p1 for next increment
        i += 1;
        if (i == 100000)
            break;
    } while (!(k < k_min));
    p_crit = p2;
}