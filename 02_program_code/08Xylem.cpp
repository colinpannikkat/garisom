#include "08Xylem.h"

/**
 * @brief Destructor for the XylemComponent class.
 * 
 * This destructor is responsible for cleaning up dynamically allocated memory
 * associated with the `soils` vector. It iterates through all elements in the
 * `soils` vector, deletes each dynamically allocated object, and ensures that
 * no memory leaks occur when the XylemComponent object is destroyed.
 * 
 * @note Assumes that all pointers in the `soils` vector were allocated using
 *       `new` and need to be deallocated using `delete`.
 */
XylemComponent::~XylemComponent() {
    int n = soils.size();
    for (int k = 0; k < n; k++) {
        delete soils[k];
    }
}

/**
 * @brief Cleans and resets the parameters of the XylemComponent object.
 * 
 * This function performs the following actions:
 * - Calls the `cleanParameters` method for the `leaf` and `stem` components.
 * - Deletes all dynamically allocated soil objects in the `soils` vector and clears the vector.
 * - Resets the arrays `e_p`, `k`, and `root_pressure` to zero.
 * - Resets the scalar parameters `rough`, `zdispl`, and `zh` to zero.
 * 
 * This ensures that the XylemComponent is in a clean and default state.
 */
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

/**
 * @brief Calculates the fatigue of a xylem component based on various parameters.
 * 
 * This function computes the fatigue of a xylem component using the provided
 * parameters related to sapwood thickness, conduit diameter, and maximum percent
 * loss of conductivity (PLC). The recovery rate is calculated internally based
 * on the sapwood thickness and conduit diameter.
 * 
 * @param b_wb Reference to a double representing the baseline water balance.
 *             This value is modified during the calculation.
 * @param sapwood_t Constant reference to a double representing the change in sapwood 
 *                  per change in diameter at breastheight.
 * @param conduit_d Constant reference to a double representing the conduit diameter.
 *                  Note: The conduit diameter must be in the correct units (um).
 * @param max_plc_x Constant reference to a double representing the maximum percent
 *                  loss of conductivity (PLC).
 * @return A double representing the calculated fatigue of the xylem component.
 */
double XylemComponent::fatigue(double &b_wb, const double &sapwood_t, const double &conduit_d, const double &max_plc_x) {
    double recoveryRate = sapwood_t * (conduit_d/10); // conduitD has to be in the right units! ()
    return (b_wb * (100 - max_plc_x)) / (recoveryRate + (100 - max_plc_x));
}

/**
 * @brief Calculates the flow of water and critical pressure for each xylem 
 * component.
 * 
 * This function computes the flow rates for various parts of the xylem system,
 * including the rhizosphere, root, stem, and leaf. It also clears historical
 * flow curves after performing the calculations.
 * 
 * @param p_inc The pressure increment used in the flow rate calculations.
 * @param k_min The minimum hydraulic conductivity used in the flow rate calculations.
 * 
 * The function performs the following steps:
 * - Iterates through all soil layers and calculates flow rates for the rhizosphere
 *   and root components. It also clears historical flow curves for the root.
 * - Calculates flow rates for the stem and leaf components, including optional
 *   calculations with additional parameters.
 * - Clears historical flow curves for the stem and leaf components.
 */
void XylemComponent::calc_flow(const double &p_inc, const double &k_min) {
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