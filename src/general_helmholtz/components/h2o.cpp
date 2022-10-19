/*-------------------------------------------------------------------------------+
| The Institute for the Design of Advanced Energy Systems Integrated Platform    |
| Framework (IDAES IP) was produced under the DOE Institute for the              |
| Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021      |
| by the software owners: The Regents of the University of California, through   |
| Lawrence Berkeley National Laboratory,  National Technology & Engineering      |
| Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University |
| Research Corporation, et al.  All rights reserved.                             |
|                                                                                |
| Please see the files COPYRIGHT.md and LICENSE.md for full copyright and        |
| license information.                                                           |
+-------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------
 Specific water functions from:

 International Association for the Properties of Water and Steam (2016).
     IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
     the Properties of Ordinary Water Substance for General Scientific Use,"
     URL: http://iapws.org/relguide/IAPWS95-2016.pdf
 Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
     Thermodynamic Properties of Ordinary Water Substance for General and
     Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.
 Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
     Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
     and Power, 122, 150-182.

 Author: John Eslick
 File: h2o.cpp
--------------------------------------------------------------------------------*/
#include<math.h>
#include<asl.h>
#undef real
#include"../config.h"
#include<iostream>
#include<string> 


double melting_tau_h2o(double pr){
  /*
    Estimate the melting temperature at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid tange of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor temperature.
  */
  double Tn, Pn;
  // Ice I Sublimation, Max error 0.15 to 251 K, 1.4 K to 235 K
  //   fit from 273.16 to 251
  if(pr < param::Pt[h2o]){
    Tn = 273.16;
    Pn = 0.611657;
    return param::T_star[h2o]/(Tn * 0.9995047*pow(pr/Pn, 0.04264942));
  }
  // Ice I Melting, Max error 0.05 K
  if(pr <= 206207.0){
    Tn = 273.16;
    Pn = 0.611657;
    return param::T_star[h2o]/(Tn * (
      -0.0000000000002087276*(pr/Pn)*(pr/Pn) -
      0.0000001637791*pr/Pn +
      1.000026
    ));
  }
  // Ice III Melting, Max error 0.065 K
  if(pr <= 350110.0){
    Tn = 251.165;
    Pn = 209900;
    return param::T_star[h2o]/(Tn * (
      -0.02487839*(pr/Pn)*(pr/Pn) -
      0.09535237*pr/Pn +
      0.9298592
    ));
  }
  // Ice V Melting, Max error 0.055 K
  if(pr <= 632400.0){
    Tn = 256.164;
    Pn = 350100.00;
    return param::T_star[h2o]/(Tn * (
      -0.02304521*(pr/Pn)*(pr/Pn) +
      0.1472252*pr/Pn +
      0.8760317
    ));
  }
  // Ice VI Melting, Max error 0.9 K
  //   shouldn't get anywhere near the max pressure on this, so we'll let this
  //   pick up the rest (goes to about 2,000 MPa)
  Tn = 273.31;
  Pn = 632400.00;
  return param::T_star[h2o]/(Tn * (
    -0.02197019*(pr/Pn)*(pr/Pn) +
    0.2161217*pr/Pn +
    0.8086136
  ));
}

double melting_liquid_delta_h2o(double pr){
  /*
    Estimate the melting liquid density at a given pressure.  This doesn't need
    to be highly accurate, it is just used to partly define the valid range of
    temperatures at a given pressure (kPa). If there is no good metling curve,
    data just supply a resonable upper limit on vapor density.
  */
  if(pr >= 400000){
    return (-9.025000E-11*pr*pr + 2.802900E-04*pr + 1.047067E+03)/param::rho_star[h2o];
  }
  if(pr >= 22500){
    return (-3.908565E-10*pr*pr + 5.195933E-04*pr + 9.992365E+02)/param::rho_star[h2o];
  }
  if(pr >= 7000){
    return (4.954471E-04*pr + 9.998203E+02)/param::rho_star[h2o];
  }
  if(pr >= param::Pt[h2o]){
    return (4.974967E-04*pr + 9.997973E+02)/param::rho_star[h2o];
  }
}


double delta_sat_v_approx_h2o(double tau){
/*
  Approximate saturated vapor density
  This equation is from the original IAPWS-95 paper
*/
  double XX = 1 - 1.0/tau;
  return exp(
    -2.03150240*pow(XX, 2.0/6.0)
    - 2.68302940*pow(XX, 4.0/6.0)
    - 5.38626492*pow(XX, 8.0/6.0)
    - 17.2991605*pow(XX, 18.0/6.0)
    - 44.7586581*pow(XX, 37.0/6.0)
    - 63.9201063*pow(XX, 71.0/6.0)
  );
}

double delta_sat_l_approx_h2o(double tau){
/*
  Approximate saturated vapor liquid
  This equation is from the original IAPWS-95 paper
*/
  double XX = 1 - 1.0/tau;
  return
    1.001
    + 1.99274064*pow(XX, 1.0/3.0)
    + 1.09965342*pow(XX, 2.0/3.0)
    - 0.510839303*pow(XX, 5.0/3.0)
    - 1.75493479*pow(XX, 16.0/3.0)
    - 45.5170352*pow(XX, 43.0/3.0)
    - 6.74694450e5*pow(XX, 110.0/3.0);
}

void phi_h2o_read_ampl(void){
    ASL *asl;
    if(cdata[PHI_AMPL_MODEL(h2o)].asl == nullptr){
        asl = ASL_alloc(ASL_read_pfgh);
        pfgh_read(jac0dim("h2o_expressions.nl", 30), 0);
        cdata[PHI_AMPL_MODEL(h2o)].asl = (void*)asl;
    }
}