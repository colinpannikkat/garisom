#include "08Xylem.h"

XylemComponent::~XylemComponent() {
    int n = soils.size();
    for (int k = 0; k < n; k++) {
        delete soils[k];
    }
}

double XylemComponent::fatigue(double &b_wb, const double &sapwood_t, const double &conduit_d, const double &max_plc_x) {
    double recoveryRate = sapwood_t * (conduit_d/10); // conduitD has to be in the right units! ()
    return (b_wb * (100 - max_plc_x)) / (recoveryRate + (100 - max_plc_x));
}