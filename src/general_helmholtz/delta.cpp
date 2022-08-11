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
 Provide functions to solve for liquid and vapor density from T and P

 Author: John Eslick
 File: delta.cpp
--------------------------------------------------------------------------------*/

#include<unordered_map>
#include<boost/functional/hash.hpp>
#include "components/function_pointers.h"
#include "props.h"
#include "solver.h"
#include "delta.h"
#include "math.h"
#include "sat.h"
#include <iostream>

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_delta_liquid2;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_delta_vapor2;

static std::vector<double> nan_vec2 = {
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan(""),
  nan("")
};

double pwrap(double delta, void *dat){
  pressure_wrap_state *d = (pressure_wrap_state*)dat;
  return (pressure(d->comp, delta, d->tau) - d->p)/param::Pc[d->comp];
}

void pwrap_gh(double delta, std::vector<double> *out, void *dat){
  pressure_wrap_state *d = (pressure_wrap_state*)dat;
  pressure2(d->comp, delta, d->tau, out);
  out->at(0) = ((*out)[0] - d->p)/param::Pc[d->comp];
  out->at(1) = (*out)[1]/param::Pc[d->comp];
  out->at(2) = (*out)[2]/param::Pc[d->comp];
}

double delta_vapor(comp_enum comp, double pr, double tau){
  double delta;
  double delta_sat;
  double p_sat;
  pressure_wrap_state ps;
  ps.comp = comp;
  ps.p = pr;
  ps.tau = tau;

  // case 0 super close to the critical point
  if(fabs(tau - 1) < 1e-9 and fabs(pr/param::Pc[comp] - 1) < 1e-9){
    return 1.0;
  }
  // case 1 P > Pc (don't really need to worry about phase change)
  //   This could really be ice, liquid or vapor, but for liquid/vapor there is
  //   no phase change, and for ice, I'll try to pretend it's still liquid and
  //   give a reasonable number anyway for math reasons
  if(pr > param::Pc[comp]){
    std::vector<double> out;
    if(tau > 1.0){
      // This is liquid just use the liquid density
      return delta_liquid(comp, pr, tau);
    }
    else{
      bracket(pwrap, 0.0001, 1.001, &delta, 40, 1e-4, 1e-4, &ps);
      halley(pwrap_gh, delta, &delta, &out, 50, 1e-9, &ps);
    }
    return delta;
  }

  // case 2 P < Psat, this is vapor or ice, if ice, I'll pretend its
  //   vapor and try to return a reasonable number anyway for math reasons
  delta_sat = sat_delta_v(comp, tau)->at(0);
  p_sat = sat_p(comp, tau)->at(0);
  std::vector<double> out;
  if(pr <= p_sat){
    bracket(pwrap, 1e-8, delta_sat, &delta, 3, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-9, &ps);
    return delta;
  }

  // case 3, you're in the liquid region, I'll still try to pretend to have vapor
  //   and see if I can give a good answer by looking between the saturated
  //   liquid density and the vapor density.  There may be multiple roots here,
  //   so I'll start from the sat density and hope to pick up the closest
  halley(pwrap_gh, delta_sat, &delta, &out, 50, 1e-9, &ps);
  if(delta < 0 || std::isnan(delta)){
      return 0.1;
  }
  return delta;
}

double delta_liquid(comp_enum comp, double pr, double tau){
  double delta;
  double delta_sat;
  double p_sat;
  pressure_wrap_state ps;
  ps.comp = comp;
  ps.p = pr;
  ps.tau = tau;

  // case 0 super close to the critical point
  if(fabs(tau - 1) < 1e-9 and fabs(pr/param::Pc[comp] - 1) < 1e-9){
    return 1.0;
  }
  // case 1 P > Pc (don't really need to worry about phase change)
  //   This could really be ice, liquid or vapor, but for liquid/vapor there is
  //   no phase change, and for ice, I'll try to pretend it's still liquid and
  //   give a reasonable number anyway for math reasons
  if(pr >= param::Pc[comp] && tau < 1.0){
    std::vector<double> out;
    bracket(pwrap, 0.0001, 1.001, &delta, 40, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-9, &ps);
    return delta;
  }
  // case 2 Psat < P < Pc, this is liquid or ice, if ice, I'll pretend its
  //   liquid and try to return a reasonable number anyway for math reasons
  delta_sat = sat_delta_l(comp, tau)->at(0);
  p_sat = sat_p(comp, tau)->at(0);
  std::vector<double> out;

  if(pr >= p_sat){
    bracket(pwrap, delta_sat, melting_liquid_density_func[comp](pr)/param::rhoc[comp], &delta, 6, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-9, &ps);
    return delta;
  }

  // case 3, you're in the vapor region, I'll still try to pretend to have liquid
  //   and see if I can give a good answer by looking between the saturated
  //   liquid density and the vapor density.  There may be multiple roots here,
  //   so I'll start from the sat density and hope to pick up the closest
  halley(pwrap_gh, delta_sat, &delta, &out, 50, 1e-9, &ps);

  if(delta < 0 || std::isnan(delta)){
      return 2.0;
  }
  return delta;
}

void delta_liquid2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  double delta_l = delta_liquid(comp, pr, tau);
  std::vector<double> *pr_vec = memo2_pressure(comp, delta_l, tau); // get derivatives
  out->resize(6);
  out->at(0) = delta_l;
  out->at(f2_1) = 1.0/pr_vec->at(f2_1);
  out->at(f2_2) = -pr_vec->at(f2_2)/pr_vec->at(f2_1);
  out->at(f2_11) = -pr_vec->at(f2_11)*out->at(f2_1)*out->at(f2_1)*out->at(f2_1);
  out->at(f2_12) = -(pr_vec->at(f2_12) + pr_vec->at(f2_11)*out->at(f2_2))*out->at(f2_1)*out->at(f2_1);
  out->at(f2_22) = -(out->at(f2_1)*(pr_vec->at(f2_22) + out->at(f2_2)*pr_vec->at(f2_12)) + pr_vec->at(f2_2)*out->at(f2_12));
}

void delta_vapor2(comp_enum comp, double pr, double tau, std::vector<double> *out){
  double delta_v = delta_vapor(comp, pr, tau);
  std::vector<double> *pr_vec = memo2_pressure(comp, delta_v, tau); // get derivatives
  out->resize(6);
  out->at(0) = delta_v;
  out->at(f2_1) = 1.0/pr_vec->at(f2_1);
  out->at(f2_2) = -pr_vec->at(f2_2)/pr_vec->at(f2_1);
  out->at(f2_11) = -pr_vec->at(f2_11)*out->at(f2_1)*out->at(f2_1)*out->at(f2_1);
  out->at(f2_12) = -(pr_vec->at(f2_12) + pr_vec->at(f2_11)*out->at(f2_2))*out->at(f2_1)*out->at(f2_1);
  out->at(f2_22) = -(out->at(f2_1)*(pr_vec->at(f2_22) + out->at(f2_2)*pr_vec->at(f2_12)) + pr_vec->at(f2_2)*out->at(f2_12));
}

std::vector<double> *memo2_delta_liquid(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_delta_liquid2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range const&){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_delta_liquid2.size() > MAX_MEMO_PROP) memo_table_delta_liquid2.clear();
  yvec_ptr = &memo_table_delta_liquid2[std::make_tuple(comp, pr, tau)];
  delta_liquid2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}

std::vector<double> *memo2_delta_vapor(comp_enum comp, double pr, double tau){
  try{
    return &memo_table_delta_vapor2.at(std::make_tuple(comp, pr, tau));
  }
  catch(std::out_of_range const&){
  }
  std::vector<double> *yvec_ptr;
  if(memo_table_delta_vapor2.size() > MAX_MEMO_PROP) memo_table_delta_vapor2.clear();
  yvec_ptr = &memo_table_delta_vapor2[std::make_tuple(comp, pr, tau)];
  delta_vapor2(comp, pr, tau, yvec_ptr);
  return yvec_ptr;
}
