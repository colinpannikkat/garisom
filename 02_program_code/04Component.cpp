#include "04Component.h"

double Component::wb(double pressure) {
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