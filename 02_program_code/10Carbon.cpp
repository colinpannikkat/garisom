#include "10Carbon.h"

 //'gets assimilation for sun leaves
void CarbonAssimilationModel::assimilation(const int &p,
                                           const double &gmax,
                                           const double &qmax,
                                           const double &comp25,
                                           const double &thetac,
                                           const double &vmax25,
                                           const double &jmax25,
                                           const double &kc25,
                                           const double &ko25,
                                           const double &svvmax,
                                           const double &svjmax,
                                           const double &havmax,
                                           const double &hdvmax,
                                           const double &hdjmax,
                                           const double &hajmax,
                                           const double &lightcurv,
                                           const std::string &night,
                                           const double eplantl[],
                                           const double lavpd[],
                                           const double leaftemp[]) {
    //'get g from D and e
    if (lavpd[p] == 0) { //if//
        gcanw[p] = gmax; //'maximum g if no lavpd
    }
    else {
        gcanw[p] = eplantl[p] / lavpd[p]; //'gcanopy in mmol m-2s-1 (leaf area)
       
        // std::cout << lavpd[p] << std::endl;
    } //End if//
    // if (p == 1)
    //     exit(1);
        //'gcanw[p] = eplantl[p] / lavpd[p] //'gcanopy in mmol m-2s-1 (leaf area)
    gcanc[p] = (gcanw[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                        //'adjust photosynthetic inputs for Tleaf
    double comp = comp25 * exp((37830 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    double numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (GAS * 298.2))) * exp((havmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
    double denominator = 1 + exp((svvmax * (273.2 + leaftemp[p]) - hdvmax) / (GAS * (273.2 + leaftemp[p])));
    double vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (GAS * 298.2))) * exp((hajmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
    denominator = 1 + exp((svjmax * (273.2 + leaftemp[p]) - hdjmax) / (GAS * (273.2 + leaftemp[p])));
    double jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    double kc = kc25 * exp((79430 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    double ko = ko25 * exp((36380 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    double rday25 = vmax25 * 0.01; //'from Medlyn 2002
    double rday = rday25 * pow(2, ((leaftemp[p] - 25) / 10.0));
    rday = rday * pow((1 + exp(1.3 * (leaftemp[p] - 55))), -1); //'high temp inhibition, collatz // never used outside of this function
    if (night == "n" && gcanc[p] > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            double ci = comp - 0.00000001; //'start at photorespiratory compensation point
            double var = -INFINITY;
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
                ci = ci + 0.00000001;
                double jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psyn[p] = var; //'always start with zero or above
            cin[p] = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=1 if
        if (p > 1) { //if//
            double ci = cin[p - 1]; //'start from previous ci
            double var = -INFINITY;
            double marker = 0.0;
            do //'loop through A-ci curve to find ci
            {
                ci = ci + 0.00000001;
                marker = gcanc[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                double jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psyn[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psyn[p] > psynmax)
                psynmax = psyn[p];
            cin[p] = ci; //'store ci
        } //End if// //'p>1 if
    }
    else { //'it//'s night or stomata are closed
        if (night == "y")
        {
            psyn[p] = 0 - rday;
            psynmax = 0;
            cin[p] = ca; //'respiration accounted for at night
        }
        if (night == "n") {
            psyn[p] = 0;
            psynmax = 0;
            cin[p] = ca; //'it//'s day and p=0, g=0
        }
    } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
}

 //'gets virgin assimilation for sun leaves
void CarbonAssimilationModel::assimilationMd(const int &p,
                                           const double &gmax,
                                           const double &qmax,
                                           const double &comp25,
                                           const double &thetac,
                                           const double &vmax25,
                                           const double &jmax25,
                                           const double &kc25,
                                           const double &ko25,
                                           const double &svvmax,
                                           const double &svjmax,
                                           const double &havmax,
                                           const double &hdvmax,
                                           const double &hdjmax,
                                           const double &hajmax,
                                           const double &lightcurv,
                                           const std::string &night,
                                           const double emd,
                                           const double lavpdmd,
                                           const double leaftmd) {
    //'get g from D and e
    if (lavpdmd == 0) { //if//
        gcanwmd = gmax; //'maximum g if no lavpd
    }
    else {
        gcanwmd = emd / lavpdmd; //'gcanopy in mmol m-2s-1 (leaf area)
    } //End if//
        //'gcanw[p] = eplantl[p] / lavpd[p] //'gcanopy in mmol m-2s-1 (leaf area)
    gcancmd = (gcanwmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                        //'adjust photosynthetic inputs for Tleaf
    double comp = comp25 * exp((37830 * ((leaftmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftmd + 273.15))); //'Bernacchi via Medlyn
    double numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (GAS * 298.2))) * exp((havmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
    double denominator = 1 + exp((svvmax * (273.2 + leaftmd) - hdvmax) / (GAS * (273.2 + leaftmd)));
    double vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (GAS * 298.2))) * exp((hajmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
    denominator = 1 + exp((svjmax * (273.2 + leaftmd) - hdjmax) / (GAS * (273.2 + leaftmd)));
    double jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    double kc = kc25 * exp((79430 * ((leaftmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftmd + 273.15))); //'Bernacchi via Medlyn
    double ko = ko25 * exp((36380 * ((leaftmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftmd + 273.15))); //'Bernacchi via Medlyn
    double rday25 = vmax25 * 0.01; //'from Medlyn 2002
    double rdaymd = rday25 * pow(2, ((leaftmd - 25) / 10.0));
    rdaymd = rdaymd * pow((1 + exp(1.3 * (leaftmd - 55))), -1); //'high temp inhibition, collatz // never used outside of this function
    if (night == "n" && gcancmd > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            double ci = comp - 0.00000001; //'start at photorespiratory compensation point
            double var = -INFINITY;
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
                ci = ci + 0.00000001;
                double jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psynmd[p] = var; //'always start with zero or above
            cinmd = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=1 if
        if (p > 1) { //if//
            double ci = cinmd - 0.00000001; //'cin(p - 1) //'start from previous ci
            double var = -INFINITY;
            double marker = 0.0;
            do //'loop through A-ci curve to find ci
            {
                ci = ci + 0.00000001;
                marker = gcancmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                double jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psynmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynmd[p] > psynmaxmd)
                psynmaxmd = psynmd[p];
            cinmd = ci; //'store ci
        } //End if// //'p>1 if
    }
    else { //'it//'s night or stomata are closed
        if (night == "y")
        {
            psynmd[p] = 0 - rdaymd;
            psynmaxmd = 0;
        }
        if (night == "n") {
            psynmd[p] = 0;
            psynmaxmd = 0;
        }
    } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
}

// Literally same function as above just using shade parameters
//'gets assimilation for shade leaves
void CarbonAssimilationModel::assimilationShade(const int &p,
                                                const double &gmax,
                                                const double &qmax,
                                                const double &comp25,
                                                const double &thetac,
                                                const double &vmax25,
                                                const double &jmax25,
                                                const double &kc25,
                                                const double &ko25,
                                                const double &svvmax,
                                                const double &svjmax,
                                                const double &havmax,
                                                const double &hdvmax,
                                                const double &hdjmax,
                                                const double &hajmax,
                                                const double &lightcurv,
                                                const std::string &night,
                                                const double eplantl[],
                                                const double lavpdsh[],
                                                const double leaftempsh[]) {
    //'get g from D and e
    if (lavpdsh[p] == 0) { //if//
        gcanwsh[p] = gmax; //'maximum g if no lavpd
    }
    else {
        gcanwsh[p] = eplantl[p] / lavpdsh[p]; //'gcanopy in mmol m-2s-1 (leaf area)
    } //End if//
        //'gcanw[p] = eplantl[p] / lavpd[p] //'gcanopy in mmol m-2s-1 (leaf area)
    gcancsh[p] = (gcanwsh[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                        //'adjust photosynthetic inputs for Tleaf
    double comp = comp25 * exp((37830 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    double numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (GAS * 298.2))) * exp((havmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
    double denominator = 1 + exp((svvmax * (273.2 + leaftempsh[p]) - hdvmax) / (GAS * (273.2 + leaftempsh[p])));
    double vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (GAS * 298.2))) * exp((hajmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
    denominator = 1 + exp((svjmax * (273.2 + leaftempsh[p]) - hdjmax) / (GAS * (273.2 + leaftempsh[p])));
    double jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    double kc = kc25 * exp((79430 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    double ko = ko25 * exp((36380 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * GAS * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    double rday25 = vmax25 * 0.01; //'from Medlyn 2002
    double rday = rday25 * pow(2, ((leaftempsh[p] - 25) / 10.0));
    rday = rday * pow((1 + exp(1.3 * (leaftempsh[p] - 55))), -1); //'high temp inhibition, collatz // never used outside of this function
    if (night == "n" && gcancsh[p] > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            double ci = comp - 0.00000001; //'start at photorespiratory compensation point
            double var = -INFINITY;
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
                ci = ci + 0.00000001;
                double jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psynsh[p] = var; //'always start with zero or above
            cinsh[p] = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=1 if
        if (p > 1) { //if//
            double ci = cinsh[p - 1]; //'start from previous ci
            double var = -INFINITY;
            double marker = 0.0;
            do //'loop through A-ci curve to find ci
            {
                ci = ci + 0.00000001;
                marker = gcancsh[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                double jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psynsh[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynsh[p] > psynmaxsh)
                psynmaxsh = psynsh[p];
            cinsh[p] = ci; //'store ci
        } //End if// //'p>1 if
    }
    else { //'it//'s night or stomata are closed
        if (night == "y")
        {
            psynsh[p] = 0 - rday;
            psynmaxsh = 0;
            cinsh[p] = ca; //'respiration accounted for at night
        }
        if (night == "n") {
            psynsh[p] = 0;
            psynmaxsh = 0;
            cinsh[p] = ca; //'it//'s day and p=0, g=0
        }
    } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
}

//'gets virgin assimilation for shade leaves for midday solution
void CarbonAssimilationModel::assimilationShadeMd(const int &p,
                                                const double &gmax,
                                                const double &qmax,
                                                const double &comp25,
                                                const double &thetac,
                                                const double &vmax25,
                                                const double &jmax25,
                                                const double &kc25,
                                                const double &ko25,
                                                const double &svvmax,
                                                const double &svjmax,
                                                const double &havmax,
                                                const double &hdvmax,
                                                const double &hdjmax,
                                                const double &hajmax,
                                                const double &lightcurv,
                                                const std::string &night,
                                                const double emd,
                                                const double lavpdshmd,
                                                const double leaftshmd) {
    //'get g from D and e
    if (lavpdshmd == 0)
    {
        gcanwshmd = gmax; //'set to maximum if no vpd
    }
    else
    {
        gcanwshmd = emd / lavpdshmd; //'gcanopy in mmol m-2s-1 (leaf area)
    }
    //'gcanwsh(p) = eplantl(p) / lavpdsh(p) //'gcanopy in mmol m-2s-1 (leaf area)
    gcancshmd = (gcanwshmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                            //'adjust photosynthetic inputs for Tleaf
    double comp = comp25 * exp((37830 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    double numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (GAS * 298.2))) * exp((havmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
    double denominator = 1 + exp((svvmax * (273.2 + leaftshmd) - hdvmax) / (GAS * (273.2 + leaftshmd)));
    double vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (GAS * 298.2))) * exp((hajmax / (GAS * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
    denominator = 1 + exp((svjmax * (273.2 + leaftshmd) - hdjmax) / (GAS * (273.2 + leaftshmd)));
    double jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    double kc = kc25 * exp((79430 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    double ko = ko25 * exp((36380 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * GAS * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    double rday25 = vmax25 * 0.01; //'from Medlyn 2002
    double rdaymdsh = rday25 * pow(2, ((leaftshmd - 25) / 10.0));
    rdaymdsh = rdaymdsh * pow((1 + exp(1.3 * (leaftshmd - 55))), -1); //'high temp inhibition, collatz // never used outside of this function
    if (night == "n" && gcancshmd > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            double ci = comp - 0.00000001; //'start at photorespiratory compensation point
            double var = -INFINITY;
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
                ci = ci + 0.00000001;
                double jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymdsh; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psynshmd[p] = var; //'always start with zero or above
            cinshmd = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=1 if
        if (p > 1) { //if//
            double ci = cinshmd - 0.00000001; //'start from previous ci backed off a bit
            double var = -INFINITY;
            double marker = 0.0;
            do //'loop through A-ci curve to find ci
            {
                ci = ci + 0.00000001;
                marker = gcancshmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                double jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                double je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + OA / ko);
                double jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymdsh; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psynshmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynshmd[p] > psynmaxshmd)
                psynmaxshmd = psynshmd[p];
            cinshmd = ci; //'store ci
        } //End if// //'p>1 if
    }
    else { //'it//'s night or stomata are closed
        if (night == "y")
        {
            psynshmd[p] = 0 - rdaymdsh;
            psynmaxshmd = 0;
        }
        if (night == "n") {
            psynshmd[p] = 0;
            psynmaxshmd = 0;
        }
    } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
}

/* Sets all parameters associated with assimilation to zero */
void CarbonAssimilationModel::clearParameters() {
    ca = 0;
    cinc = 0;
    cincsh = 0;
    atree = 0;
    qsl = 0;
    qsh = 0;
    psynmax = 0;
    psynmaxsh = 0;
    psynact = 0;
    psynactsh = 0;

    gcanwmd = 0;
    gcanwshmd = 0;
    gcancmd = 0;
    gcancshmd = 0;
    cinmd = 0;
    cinshmd = 0;
    psynmaxmd = 0;
    psynmaxshmd = 0;
    gcmd = 0;
    gcmdsh = 0;

    for (int i = 0; i < CURVE_MAX; ++i) {
        cin[i] = 0;
        cinsh[i] = 0;
        psyn[i] = 0;
        psynsh[i] = 0;
        psynmd[i] = 0;
        psynshmd[i] = 0;
        psync[i] = 0;
        gcanw[i] = 0;
        gcanc[i] = 0;
        gcanwsh[i] = 0;
        gcancsh[i] = 0;
    }
}