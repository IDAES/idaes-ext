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

/*------------------------------------------------------------------------------
Author: John Eslick
File sat.h

 This file contains the functions to solve the saturation curve.  The method
 is described in the following paper:

 Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
     State from Helmholtz Energy Equations of State." Journal of Thermal
     Science and Technology, 3(3), 442-451.
------------------------------------------------------------------------------*/

#ifndef _INCLUDE_SAT_H_
#define _INCLUDE_SAT_H_

f12_struct sat_tau(uint comp, double pr);  // tau_sat as a function of pressure [kPa]
f12_struct sat_t(uint comp, double pr); 
f12_struct sat_p(uint comp, double tau); // p_sat as a function of tau
f12_struct sat_delta_v(uint comp, double tau); // vapor delta_sat as a function of tau
f12_struct sat_delta_l(uint comp, double tau); // liquid delta_sat as a function of tau


f12_struct sat_p_t(uint comp, double t); // pressure [kPa] as a function of temperature [K]
f12_struct sat_h_liq_t(uint comp, double t); // sat. liquid enthalpy [kJ/kg] as a function of temperature [K]
f12_struct sat_h_vap_t(uint comp, double t); // sat. vapor enthalpy [kJ/kg] as a function of temperature [K]
f12_struct sat_s_liq_t(uint comp, double t); // sat. liquid entropy [kJ/kg/K] as a function of temperature [K]
f12_struct sat_s_vap_t(uint comp, double t); // sat. vapor entropy [kJ/kg/K] as a function of temperature [K]
f12_struct sat_u_liq_t(uint comp, double t); // sat. liquid internal energy [kJ/kg] as a function of  temperature [K]
f12_struct sat_u_vap_t(uint comp, double t); // sat. vapor internal energy [kJ/kg] as a function of  temperature [K]
f12_struct sat_v_liq_t(uint comp, double t); // sat. liquid specific volume [m3/kg] as a function of temperature [K]
f12_struct sat_v_vap_t(uint comp, double t); // sat. vapor specific volume [m3/kg] as a function of temperature [K]

f12_struct sat_h_liq_p(uint comp, double pr); // sat. liquid enthalpy [kJ/kg] as a function of pressure [kPa]
f12_struct sat_h_vap_p(uint comp, double pr); // sat. vapor enthalpy [kJ/kg] as a function of pressure [kPa]
f12_struct sat_s_liq_p(uint comp, double pr); // sat. liquid entropy [kJ/kg/K] as a function of pressure [kPa]
f12_struct sat_s_vap_p(uint comp, double pr); // sat. vapor entropy [kJ/kg/K] as a function of pressure [kPa]
f12_struct sat_u_liq_p(uint comp, double pr); // sat. liquid internal energy [kJ/kg] as a function of pressure [kPa]
f12_struct sat_u_vap_p(uint comp, double pr); // sat. vapor internal energy [kJ/kg] as a function of pressure [kPa]
f12_struct sat_v_liq_p(uint comp, double pr); // sat. liquid specific volume [m3/kg] as a function of pressure [kPa]
f12_struct sat_v_vap_p(uint comp, double pr); // sat. vapor specific volume [m3/kg] as a function of pressure [kPa]


#define SAT_PROP_FUNCTION_OF_T(NEW_FUNC, TAU_FUNC) \
f12_struct NEW_FUNC(uint comp, double t){ \
  double Tc = cdata[comp].T_star; \
  f12_struct res, p_res = TAU_FUNC(comp, Tc/t); \
  res.f = p_res.f; \
  res.f_1 = -p_res.f_1 * Tc / t / t; \
  res.f_11 = p_res.f_11 * Tc / t / t * Tc / t / t + 2*p_res.f_1 * Tc / t / t / t; \
  return res; \
}

#define SAT_FUNCTION_OF_TAU(NEW_FUNC, PROP_FUNC, DELTA_FUNC) \
f12_struct NEW_FUNC(uint comp, double tau){ \
  f12_struct res, delta = DELTA_FUNC(comp, tau); \
  f22_struct h = PROP_FUNC(comp, delta.f, tau); \
  res.f = h.f; \
  res.f_1 = h.f_1 * delta.f_1 + h.f_2; \
  res.f_11 = h.f_11 * delta.f_1 * delta.f_1 + 2*h.f_12 * delta.f_1 + h.f_1 * delta.f_11 + h.f_22; \
  return res; \
}

#define SAT_PROP_FUNCTION_OF_PRESSURE(NEW_FUNC, PROP_FUNC) \
f12_struct NEW_FUNC(uint comp, double pr){ \
  f12_struct res; \
  f12_struct tau = sat_tau(comp, pr); \
  f12_struct h = PROP_FUNC(comp, tau.f); \
  res.f = h.f; \
  res.f_1 = h.f_1 * tau.f_1; \
  res.f_11 = h.f_11 * tau.f_1 * tau.f_1 + h.f_1 * tau.f_11; \
  return res; \
}

#endif
