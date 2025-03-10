#include "08Xylem.h"

XylemComponent::~XylemComponent() {
    int n = soils.size();
    for (int k = 0; k < n; k++) {
        delete soils[k];
    }
}

void XylemComponent::cleanParameters() {
    leaf.cleanParameters();
    stem.cleanParameters();
    
    for (int k = 0; k < soils.size(); k++) {
        delete soils[k];
    }
    soils.clear();

    std::fill(std::begin(e_p), std::end(e_p), 0.0);
    std::fill(std::begin(k), std::end(k), 0.0);
    std::fill(std::begin(root_pressure), std::end(root_pressure), 0.0);

    rough = 0.0;
    zdispl = 0.0;
    zh = 0.0;
}

double XylemComponent::fatigue(double &b_wb, const double &sapwood_t, const double &conduit_d, const double &max_plc_x) {
    double recoveryRate = sapwood_t * (conduit_d/10); // conduitD has to be in the right units! ()
    return (b_wb * (100 - max_plc_x)) / (recoveryRate + (100 - max_plc_x));
}

void XylemComponent::calc_net_flow(const double &p_inc, const double &k_min) {
    for (int z = 1; z <= this->num_layers; z++) {
        soils[z]->rhizosphere.calc_flow_rate(p_inc, k_min);
        soils[z]->root.calc_flow_rate(p_inc, k_min);
        soils[z]->root.calc_flow_rate(p_inc, k_min, true);
        soils[z]->root.clearHistoricalCurves();
    }
    stem.calc_flow_rate(p_inc, k_min);
    leaf.calc_flow_rate(p_inc, k_min);
    stem.calc_flow_rate(p_inc, k_min, true);
    leaf.calc_flow_rate(p_inc, k_min, true);

    stem.clearHistoricalCurves();
    leaf.clearHistoricalCurves();
}