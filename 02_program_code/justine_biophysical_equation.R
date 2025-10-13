# Author: Justine Rojas
# Date: 9/26/25
# Biophysical equation



# Define the biophysical potential function
biophysical_potential <- function(psi, pi, gamma_p, gamma_thresh, phi_T) {
  delta <- psi - pi

  if (delta > gamma_thresh) {
    return((1 / log(2)) * phi_T * (delta - gamma_p))
  } else {
    return(0)
  }
}

# Example
biophysical_potential(psi = -1.5, pi = -0.5, gamma_p = 0.2, gamma_thresh = 0.1, phi_T = 0.8)
# Output: 0
biophysical_potential(psi = -0.5, pi = -1.2, gamma_p = 0.3, gamma_thresh = 0.1, phi_T = 0.9)
# Output: 0.519



# Cell wall extensibility for biophysical equation
# Define the phi(T) function
phi_T_function <- function(T, phi_max, lambda, T1, k, delta_HA, R, delta_SD, delta_HD) {

  # Logistic activation term
  logistic_term <- 1 / (1 + exp(lambda * (T1 - T)))

  # Exponential activation term
  numerator <- k * T * exp(delta_HA / (R * T))
  denominator <- 1 + exp((delta_SD / R) * (1 - (delta_HD / delta_SD)))

  # Final phi(T)
  phi_T <- phi_max * logistic_term * (numerator / denominator)

  return(phi_T)
}


# Constants (example values)
phi_max <- 1.0         # Max scaling factor
lambda <- 0.2          # Sensitivity to T1
T1 <- 298              # Reference temperature (K)
T <- 310               # Current temperature (K)
k <- 1.38e-23          # Boltzmann constant (J/K)
delta_HA <- 50000      # Enthalpy A (J/mol)
R <- 8.314             # Universal gas constant (J/mol·K)
delta_SD <- 100        # Entropy D (J/mol·K)
delta_HD <- 40000      # Enthalpy D (J/mol)

# Compute phi(T)
phi_T_function(T, phi_max, lambda, T1, k, delta_HA, R, delta_SD, delta_HD)
