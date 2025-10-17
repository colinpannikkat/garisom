#include "07Leaf.h"

/**
 * @brief Calculates leaf temperature and leaf-to-air vapor pressure deficit (VPD).
 *
 * This function computes the leaf temperature and leaf-to-air vapor pressure deficit
 * based on various environmental and plant parameters. It uses energy balance equations
 * and other physiological relationships to determine these values.
 *
 * @param p Index of the plant or leaf for which calculations are performed.
 * @param airtemp Ambient air temperature (in degrees Celsius).
 * @param eplant Array of plant transpiration rates (in kg hr-1 m-2 basal area).
 * @param vpd Vapor pressure deficit of the air (in kPa).
 * @param wind Wind speed (in m s^-1).
 * @param laperba Leaf area per unit ground area (dimensionless).
 * @param leafwidth Width of the leaf (in meters).
 * @param patm Atmospheric pressure (in kPa).
 *
 * @details
 * - The function calculates the total absorbed radiation, heat of vaporization,
 *   radiative conductance, and heat conductance.
 * - It computes the energy balance to determine the leaf temperature.
 * - The leaf-to-air vapor pressure deficit (lavpd) is calculated based on the
 *   saturated mole fraction of vapor in air and the vapor pressure deficit.
 * - Negative values of lavpd are set to zero to ensure physical consistency.
 *
 * @note
 * - The energy balance assumes a two-sided leaf.
 * - Constants such as emissivity and Stefan-Boltzmann constant (SBC) are used
 *   in the calculations.
 * - The function modifies the `eplantl`, `leaftemp`, and `lavpd` arrays for the
 *   given index `p`.
 * - This must be called before tempShade, as tempShade will the convert eplantlp 
 *   into mmol m-2 s-1.
 */
void LeafComponent::temp(const int &tod, const int p,
                         const double &airtemp,
                         const double eplant[],
                         const double &vpd,
                         const double &wind,
                         const double &laperba,
                         const double &leafwidth,
                         const double &patm)  
{
    double rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
    lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
    grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
    gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
    eplantl[p] = eplant[p] * (1.0 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
    double numerator = rabs - emiss * SBC * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
    double denominator = SHA * (grad + gha);
    leaftemp[p] = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftemp[p])) - lavpd[p] + vpd; //'leaf-to-air vpd
    if (lavpd[p] < 0)
        lavpd[p] = 0; //'don//'t allow negative lavpd
}

/**
 * @brief Calculates the temperature and vapor pressure deficit (VPD) for shaded leaves.
 *
 * This function computes the leaf temperature and leaf-to-air vapor pressure deficit (VPD)
 * for shaded leaves based on absorbed radiation, air temperature, atmospheric pressure,
 * and vapor pressure deficit. It also ensures that the calculated VPD is non-negative.
 *
 * @param p Index of the plant or leaf for which the calculations are performed.
 * @param airtemp Air temperature in degrees Celsius.
 * @param patm Atmospheric pressure in kPa.
 * @param vpd Vapor pressure deficit in kPa.
 *
 * The function performs the following calculations:
 * - Computes the total absorbed radiation for shaded leaves.
 * - Calculates the numerator and denominator for the energy balance equation.
 * - Determines the leaf temperature for shaded leaves.
 * - Computes the saturated mole fraction of vapor in air and the leaf-to-air VPD.
 * - Ensures that the leaf-to-air VPD is non-negative.
 * - Converts the transpiration rate to mmol m⁻² s⁻¹.
 */
void LeafComponent::tempShade(const int &tod, const int &p,
                              const double &airtemp,
                              const double &patm,
                              const double &vpd) 
{
    double rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
                                                            //'lambda = -42.9143 * airtemp + 45064.3 //'heat of vaporization for water at air temp in J mol-1
                                                            //'grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * airtemp ^ 2 //'radiative conductance (long wave) at air temp in mol m-2 s-1
                                                            //'gha = 1.4 * 0.135 * (wind / leafwidth) ^ 0.5 //'heat conductance in mol m-2s-1
                                                            //'eplantl[p] = eplant[p] * (1 / laperba) * (1 / 3600) * 55.4 //'convert to E per leaf area in mol m-2s-1
    double numerator = rabs - emiss * SBC * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
    double denominator = SHA * (grad + gha);
    leaftempsh[p] = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftempsh[p])) - lavpdsh[p] + vpd; //'leaf-to-air vpd
    if (lavpdsh[p] < 0)
        lavpdsh[p] = 0; //'don//'t allow negative lavpd
    eplantl[p] = eplantl[p] * 1000; //'convert to mmol m-2 s-1
}

 /**
 * @brief Calculates virgin leaf temperature and leaf-to-air vapor pressure deficit (VPD).
 *
 * This function computes the midday leaf temperature and leaf-to-air vapor pressure 
 * deficit based on various environmental and plant parameters. It uses energy balance 
 * equations and other physiological relationships to determine these values.
 *
 * @param p Index of the plant or leaf for which calculations are performed.
 * @param airtemp Ambient air temperature (in degrees Celsius).
 * @param eplant Array of plant transpiration rates (in kg hr-1 m-2 basal area).
 * @param vpd Vapor pressure deficit of the air (in kPa).
 * @param wind Wind speed (in m s^-1).
 * @param laperba Leaf area per unit ground area (dimensionless).
 * @param leafwidth Width of the leaf (in meters).
 * @param patm Atmospheric pressure (in kPa).
 *
 * @details
 * - The function calculates the total absorbed radiation, heat of vaporization,
 *   radiative conductance, and heat conductance.
 * - It computes the energy balance to determine the leaf temperature.
 * - The leaf-to-air vapor pressure deficit (lavpd) is calculated based on the
 *   saturated mole fraction of vapor in air and the vapor pressure deficit.
 * - Negative values of lavpd are set to zero to ensure physical consistency.
 *
 * @note
 * - The energy balance assumes a two-sided leaf.
 * - Constants such as emissivity and Stefan-Boltzmann constant (SBC) are used
 *   in the calculations.
 * - The function modifies the `emd`, `leaftmd`, and `lavpdmd` arrays for the
 *   given index `p`.
 * - This must be called before tempShade, as tempShade will the convert eplantlp 
 *   into mmol m-2 s-1.
 */
void LeafComponent::tempMd(const int &tod,
                           const int &p,
                           const double &e,
                           const double &airtemp,
                           const double &vpd,
                           const double &wind,
                           const double &laperba,
                           const double &leafwidth,
                           const double &patm)  
{
    double rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
    lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
    grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
    gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
    emd = e * (1.0 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
    double numerator = rabs - emiss * SBC * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
    double denominator = SHA * (grad + gha);
    leaftmd = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftmd)) - lavpdmd + vpd; //'leaf-to-air vpd
    if (lavpdmd < 0)
        lavpdmd = 0; //'don//'t allow negative lavpd
}

/**
 * @brief Calculates the virgin temperature and vapor pressure deficit (VPD) for 
 * shaded leaves.
 *
 * This function computes the midday leaf temperature and leaf-to-air vapor pressure 
 * deficit (VPD) for shaded leaves based on absorbed radiation, air temperature, 
 * atmospheric pressure, and vapor pressure deficit. It also ensures that the 
 * calculated VPD is non-negative.
 *
 * @param p Index of the plant or leaf for which the calculations are performed.
 * @param airtemp Air temperature in degrees Celsius.
 * @param patm Atmospheric pressure in kPa.
 * @param vpd Vapor pressure deficit in kPa.
 *
 * The function performs the following calculations:
 * - Computes the total absorbed radiation for shaded leaves.
 * - Calculates the numerator and denominator for the energy balance equation.
 * - Determines the leaf temperature for shaded leaves.
 * - Computes the saturated mole fraction of vapor in air and the leaf-to-air VPD.
 * - Ensures that the leaf-to-air VPD is non-negative.
 * - Converts the transpiration rate to mmol m⁻² s⁻¹.
 */
void LeafComponent::tempShadeMd(const int &tod, const int &p,
                           const double &e,
                           const double &airtemp,
                           const double &vpd,
                           const double &patm)  
{
    double rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
                                                            //'lambda = -42.9143 * airtemp + 45064.3 //'heat of vaporization for water at air temp in J mol-1
                                                            //'grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * airtemp ^ 2 //'radiative conductance (long wave) at air temp in mol m-2 s-1
                                                            //'gha = 1.4 * 0.135 * (wind / leafwidth) ^ 0.5 //'heat conductance in mol m-2s-1
                                                            //'eplantl(p) = eplant(p) * (1 / laperba) * (1 / 3600) * 55.4 //'convert to E per leaf area in mol m-2s-1
    double numerator = rabs - emiss * SBC * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
    double denominator = SHA * (grad + gha);
    leaftshmd = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftshmd)) - lavpdshmd + vpd; //'leaf-to-air vpd
    if (lavpdshmd < 0)
        lavpdshmd = 0; //'don't allow negative lavpd
    emd = emd * 1000; //'convert to mmol m-2 s-1
}

void LeafComponent::cleanParameters() {
    Component::cleanParameters();
    std::fill(std::begin(lavpd), std::end(lavpd), 0.0);
    std::fill(std::begin(lavpdsh), std::end(lavpdsh), 0.0);
    std::fill(std::begin(leaftemp), std::end(leaftemp), 0.0);
    std::fill(std::begin(leaftempsh), std::end(leaftempsh), 0.0);
    std::fill(std::begin(eplantl), std::end(eplantl), 0.0);

    emd = 0.0;
    lavpdmd = 0.0;
    lavpdshmd = 0.0;
    leaftmd = 0.0;
    leaftshmd = 0.0;

    ssun = 0.0;
    sshade = 0.0;
    sref = 0.0;
    sbottom = 0.0;
    emiss = 0.0;
    la = 0.0;
    lg = 0.0;
    lai = 0.0;
    laisl = 0.0;
    laish = 0.0;
    lambda = 0.0;
    grad = 0.0;
    gha = 0.0;
}