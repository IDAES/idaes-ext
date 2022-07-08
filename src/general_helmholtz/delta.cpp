#include "props.h"
#include "solver.h"
#include "param.h"
#include "delta.h"
#include "math.h"
#include "sat.h"

double pwrap(double delta, void *dat){
  pressure_wrap_state *d = (pressure_wrap_state*)dat;
  return (pressure(d->comp, delta, d->tau) - d->p)/Pc[d->comp];
}

void pwrap_gh(double delta, std::vector<double> *out, void *dat){
  pressure_wrap_state *d = (pressure_wrap_state*)dat;
  pressure2(d->comp, delta, d->tau, out);
  out->at(0) = ((*out)[0] - d->p)/Pc[d->comp];
  out->at(1) = (*out)[1]/Pc[d->comp];
  out->at(2) = (*out)[2]/Pc[d->comp];
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
  if(fabs(tau - 1) < 1e-8 and fabs(pr/Pc[comp] - 1) < 1e-8){
    return 1.0;
  }
  // case 1 P > Pc (don't really need to worry about phase change)
  //   This could really be ice, liquid or vapor, but for liquid/vapor there is
  //   no phase change, and for ice, I'll try to pretend it's still liquid and
  //   give a reasonable number anyway for math reasons
  if(pr > Pc[comp]){
    std::vector<double> out;
    bracket(pwrap, 0, rho_max[comp], &delta, 20, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-10, &ps);
    return delta;
  }

  // case 2 P < Psat, this is vapor or ice, if ice, I'll pretend its
  //   vapor and try to return a reasonable number anyway for math reasons
  delta_sat = sat_delta_v(comp_enum::h2o, tau)->at(0);
  p_sat = sat_p(comp_enum::h2o, tau)->at(0);
  std::vector<double> out;
  if(pr <= p_sat){
    bracket(pwrap, 0, delta_sat, &delta, 3, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-10, &ps);
    return delta;
  }

  // case 3, you're in the liquid region, I'll still try to pretend to have vapor
  //   and see if I can give a good answer by looking between the saturated
  //   liquid density and the vapor density.  There may be multiple roots here,
  //   so I'll start from the sat density and hope to pick up the closest
  halley(pwrap_gh, delta_sat, &delta, &out, 50, 1e-10, &ps);
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
  if(fabs(tau - 1) < 1e-8 and fabs(pr/Pc[comp] - 1) < 1e-8){
    return 1.0;
  }
  // case 1 P > Pc (don't really need to worry about phase change)
  //   This could really be ice, liquid or vapor, but for liquid/vapor there is
  //   no phase change, and for ice, I'll try to pretend it's still liquid and
  //   give a reasonable number anyway for math reasons
  if(pr > Pc[comp]){
    std::vector<double> out;
    bracket(pwrap, 0, rho_max[comp], &delta, 20, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-10, &ps);
    return delta;
  }

  // case 2 Psat < P < Pc, this is liquid or ice, if ice, I'll pretend its
  //   liquid and try to return a reasonable number anyway for math reasons
  delta_sat = sat_delta_l(comp_enum::h2o, tau)->at(0);
  p_sat = sat_p(comp_enum::h2o, tau)->at(0);
  std::vector<double> out;
  if(pr >= p_sat){
    bracket(pwrap, delta_sat, rho_max[comp], &delta, 3, 1e-4, 1e-4, &ps);
    halley(pwrap_gh, delta, &delta, &out, 50, 1e-10, &ps);
    return delta;
  }

  // case 3, you're in the vapor region, I'll still try to pretend to have liquid
  //   and see if I can give a good answer by looking between the saturated
  //   liquid density and the vapor density.  There may be multiple roots here,
  //   so I'll start from the sat density and hope to pick up the closest
  halley(pwrap_gh, delta_sat, &delta, &out, 50, 1e-10, &ps);
  return delta;
}
