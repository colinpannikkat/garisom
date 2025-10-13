/**
 * Author: Colin Pannikkat
 * Date: 10/13/25
 * Modified from R script provided by Justine Rojas.
 * 
 * This script provides the function for calculating the biophysical potential 
 * of the cambium, this is the "potential cambial cell division rate if it were 
 * only limited by turgor-driven cell expansion and temperature-dependent cell 
 * metabolism"
 * 
 * Equations from "Probing the interplay of biophysical constraints and 
 * photosynthesis to model tree growth", Cabon et al. (2024).
 */

#include <iostream>

/**
 * @brief This function calculates the biophysical potential of the cambium.
 * 
 * @param psi           Tree water potential at breast height (MPa)
 * @param pi            Cambium osmotic potential (MPa)
 * @param gamma_p       Pressure yield threshold (MPa)
 * @param gamma_thresh  Pressure yield threshold (MPa) for conditional value
 * @param phi_T         Cell wall extensibility
 */
double calc_biophysical_potential(
    double psi, 
    double pi,
    double gamma_p,
    double gamma_thresh,
    double phi_T
){
    double delta = psi - pi;

    if (delta > gamma_thresh)
        return ((1 / log(2)) * phi_T * (delta - gamma_p));
    
    return 0.;
}

/**
 * @brief This function calculates the cell wall extensibility of a plant given
 * the air temperature, T and other environmental constants.
 * 
 * @param T         Air temperature (K)
 * @param phi_max   Maximum cell wall extensibility
 * @param lambda    Scaling parameter, sensitivty to T1
 * @param T1        Threshold temperature (K) for sigmoid (logistic) function to
 *                  equal 0.5
 * @param k         Boltzmann constant (J/K)
 * @param delta_HA  Enthalpy of activation (J/mol)
 * @param R         Universal gas constant (J/mol * K)
 * @param delta_SD  Entropy difference of active and inactive states of enzyme
 *                  system (J/mol * K)
 * @param delta_HD  Enthalpy difference of active and inactive states of enzyme 
 *                  system (J/mol)
 */
double calc_phi_T(
    double T,
    double phi_max,
    double lambda,
    double T1,
    double k,
    double delta_HA,
    double R,
    double delta_SD,
    double delta_HD
){

    // Logistic activation term
    double logistic_term = 1 / (1 + exp(lambda * (T1 - T)));

    // Exponential activation term
    double numerator = k * T * exp(delta_HA / (R * T));
    double denominator = 1 + exp((delta_SD / R) * (1 - (delta_HD / delta_SD)));

    // Final phi(T)
    double phi_T = phi_max * logistic_term * (numerator / denominator);

    return phi_T;
}

int main() {

    double res = 0;

    // Output: 0
    res = calc_biophysical_potential(-1.5, -0.5, 0.2, 0.1, 0.8);
    std::cout << "biophysical potential = " << res << std::endl;
    // Output: 0.519
    res = calc_biophysical_potential(-0.5, -1.2, 0.3, 0.1, 0.9);
    std::cout << "biophysical potential = " << res << std::endl;

    // Example values, output: 1.04418e-12
    double phi_max = 1.0,
           lambda = 0.2,
           T1 = 298.,
           T = 310.,
           k = 1.38e-23,
           delta_HA = 50000.,
           R = 8.314,
           delta_SD = 100.,
           delta_HD = 40000.;

    res = calc_phi_T(T, phi_max, lambda, T1, k, delta_HA, R, delta_SD, delta_HD);
    std::cout << "phi(T) = " << res << std::endl;

    return 0;
}