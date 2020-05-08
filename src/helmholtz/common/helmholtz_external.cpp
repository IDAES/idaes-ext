/*------------------------------------------------------------------------------
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 Main Helmholtz calculations. For now, the non-analytic terms are not included,
 because they cause a singularity at the critical point.  Well behaved
 calculations are prefered to extreme accuracy in near the critical point.

 Author: John Eslick
 File helmholtz_external.cpp
------------------------------------------------------------------------------*/
#include<stdio.h>
#include<cmath>
#include<iostream>

#include"helmholtz_memo.h"
#include"helmholtz_external.h"
#include"helmholtz_config.h"
#include"helmholtz_deriv_parts.h"

#define TAU_LOW 0.15
#define TAU_HIGH 4.0

/*------------------------------------------------------------------------------
Basic Property Calculations

Your guide to variables in this section:
  T_c: critical temperature (K) defined in iapws95_param.h
  rho_c: critical density (kg/m3) defined in iapws95_param.h
  R: ideal gas constant (kJ/kg/K) defined in iapws95_param.h
  delta: reduced density rho/rho_c
  tau: 1/reduced temperature T_c/T
*-----------------------------------------------------------------------------*/

s_real p(s_real delta, s_real tau){ //pressure (kPa)
  return R*(delta*rho_c)*(T_c/tau)*(1 + delta*phir_delta(delta, tau));
}

s_real u(s_real delta, s_real tau){ // internal energy (kJ/kg)
  return R*T_c*(phi0_tau(tau) + phir_tau(delta, tau));
}

s_real s(s_real delta, s_real tau){ // entropy (kJ/kg/K)
  return R*(tau*(phi0_tau(tau) + phir_tau(delta, tau)) -
         phi0(delta, tau) - phir(delta, tau));
}

s_real h(s_real delta, s_real tau){ // enthalpy (kJ/kg)
  return R*(T_c/tau)*(1 + tau*(phi0_tau(tau) + phir_tau(delta, tau)) +
         delta*phir_delta(delta, tau));
}

s_real f(s_real delta, s_real tau){ // Helmholtz free energy (kJ/kg)
  return R*(T_c/tau)*(phi0(delta, tau) + phir(delta, tau));
}

s_real g(s_real delta, s_real tau){ // Gibbs free energy (kJ/kg)
  return h(delta, tau) - T_c/tau*s(delta, tau);
}

s_real cv(s_real delta, s_real tau){ // Constant volume heat capacity (kJ/kg/K)
  return -R*tau*tau*(phi0_tau2(tau) + phir_tau2(delta, tau));
}

s_real cp(s_real delta, s_real tau){ // Constant pressure heat capacity (kJ/kg/K)
  return cv(delta, tau) + R*XD/XE;
}

s_real w(s_real delta, s_real tau){ // Speed of sound (m/s)
  return s_sqrt(XB*1000*(XE + R*XD/cv(delta, tau)));
}

/*------------------------------------------------------------------------------
Basic Properties Calculations, 1st and 2nd derivatives
------------------------------------------------------------------------------*/

s_real p_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XA_d + XA_d*XC + XA*XC_d;
    grad[1] = XA_t + XA_t*XC + XA*XC_t;
    if(hes!=NULL){
      hes[0] = XA_dd + XA_dd*XC + 2*XA_d*XC_d + XA*XC_dd;
      hes[1] = XA_dt + XA_dt*XC + XA_d*XC_t + XA_t*XC_d + XA*XC_dt;
      hes[2] = XA_tt + XA_tt*XC + 2*XA_t*XC_t + XA*XC_tt;}}
  return p(delta, tau);
}

s_real u_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*XH_d;
    grad[1] = XB_t*XH + XB*XH_t;
    if(hes!=NULL){
      hes[0] = XB*XH_dd;
      hes[1] = XB_t*XH_d + XB*XH_dt;
      hes[2] = XB_tt*XH + 2*XB_t*XH_t + XB*XH_tt;}}
  return u(delta, tau);
}

s_real s_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = R*(XH_d - phi0_delta(delta) - phir_delta(delta, tau));
    grad[1] = R*(XH_t - phi0_tau(tau) - phir_tau(delta, tau));
    if(hes!=NULL){
      hes[0] = R*(XH_dd - phi0_delta2(delta) - phir_delta2(delta, tau));
      hes[1] = R*(XH_dt - phir_delta_tau(delta, tau));
      hes[2] = R*(XH_tt - phi0_tau2(tau) - phir_tau2(delta, tau));}}
  return s(delta, tau);
}

s_real h_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*(XH_d + XC_d);
    grad[1] = XB_t*(1 + XH + XC) + XB*(XH_t + XC_t);
    if(hes!=NULL){
      hes[0] = XB*(XH_dd + XC_dd);
      hes[1] = XB_t*(XH_d + XC_d) + XB*(XH_dt + XC_dt);
      hes[2] = XB_tt*(1 + XH + XC) + 2*XB_t*(XH_t + XC_t) + XB*(XH_tt + XC_tt);}}
  return h(delta, tau);
}


s_real hvpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::HVPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, h;
  delta = delta_vap(pr, tau, gradd, hesd);
  h = h_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::HVPT_FUNC, pr, tau, h, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return h;
}

s_real hlpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::HLPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, h;
  delta = delta_liq(pr, tau, gradd, hesd);
  h = h_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::HLPT_FUNC, pr, tau, h, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return h;
}

s_real svpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::SVPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, s;
  delta = delta_vap(pr, tau, gradd, hesd);
  s = s_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::SVPT_FUNC, pr, tau, s, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return s;
}

s_real slpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::SLPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, s;
  delta = delta_liq(pr, tau, gradd, hesd);
  s = s_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::SLPT_FUNC, pr, tau, s, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return s;
}

s_real uvpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::UVPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, u;
  delta = delta_vap(pr, tau, gradd, hesd);
  u = u_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::UVPT_FUNC, pr, tau, s, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return u;
}

s_real ulpt_with_derivs(s_real pr, s_real tau, s_real *grad, s_real *hes){
  s_real val = memoize::get_bin(memoize::ULPT_FUNC, pr, tau, grad, hes);
  if(!std::isnan(val)) return val;
  bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
  // Since I'm going to cache results, grad and hes will get calculated
  // whether requested or not.  If they were NULL allocate space.
  if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

  s_real gradh[2], hesh[3]; // derivatives of h(delta, tau)
  s_real gradd[2], hesd[3]; // derivatives of delta(p, tau)
  s_real delta, u;
  delta = delta_liq(pr, tau, gradd, hesd);
  u = u_with_derivs(delta, tau, gradh, hesh);
  grad[0] = gradh[0]*gradd[0];
  grad[1] = gradh[1] + gradh[0]*gradd[1];
  hes[0] = hesh[0]*gradd[0]*gradd[0] + gradh[0]*hesd[0];
  hes[1] = hesh[1]*gradd[0] + hesh[0]*gradd[0]*gradd[1] + gradh[0]*hesd[1];
  hes[2] = hesh[2] + 2*hesh[1]*gradd[1] + hesh[0]*gradd[1]*gradd[1] + gradh[0]*hesd[2];

  memoize::add_bin(memoize::ULPT_FUNC, pr, tau, s, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes; //   function
  return u;
}

s_real vf_with_derivs(s_real ht, s_real pr, s_real *grad, s_real *hes){
    s_real val = memoize::get_bin(memoize::VF_FUNC, ht, pr, grad, hes);
    if(!std::isnan(val)) return val;
    bool free_grad = 0, free_hes = 0; // grad and/or hes not provided so allocate
    // Since I'm going to cache results, grad and hes will get calculated
    // whether requested or not.  If they were NULL allocate space.
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

    s_real tau, gradt[1], hest[1], gradhv[2], heshv[3], gradhl[2], heshl[3];
    tau = sat_tau_with_derivs(pr, gradt, hest);
    s_real hv = hvpt_with_derivs(pr, tau, gradhv, heshv);
    s_real hl = hlpt_with_derivs(pr, tau, gradhl, heshl);

    if(pr >= P_c){zero_derivs2(grad, hes); return 0.0;} //Classify supercritical as liquid
    else if(pr <= P_t){zero_derivs2(grad, hes); return 1.0;} //below tripple point going to assume vapor, because it's that or a solid and I'm not dealing with solids.
    else if(hl > ht){ zero_derivs2(grad, hes); return 0.0;}
    else if(hv < ht){ zero_derivs2(grad, hes); return 1.0;}

    s_real dhvdp = gradhv[0] + gradhv[1]*gradt[0];
    s_real dhldp = gradhl[0] + gradhl[1]*gradt[0];
    grad[0] = 1.0/(hv - hl);
    grad[1] = -dhldp/(hv - hl) - (ht - hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
    s_real d2hvdp2 = heshv[0] + 2*heshv[1]*gradt[0] + heshv[2]*gradt[0]*gradt[0] + gradhv[1]*hest[0];
    s_real d2hldp2 = heshl[0] + 2*heshl[1]*gradt[0] + heshl[2]*gradt[0]*gradt[0] + gradhl[1]*hest[0];
    hes[0] = 0;
    hes[1] = -1.0/(hv-hl)/(hv-hl)*(dhvdp - dhldp);
    hes[2] = -d2hldp2/(hv-hl) + 2*dhldp/(hv-hl)/(hv-hl)*(dhvdp - dhldp) +
              2*(ht-hl)/(hv-hl)/(hv-hl)/(hv-hl)*(dhvdp - dhldp)*(dhvdp - dhldp) -
              (ht-hl)/(hv-hl)/(hv-hl)*(d2hvdp2 - d2hldp2);
    s_real vf = (ht - hl)/(hv - hl);
    memoize::add_bin(memoize::VF_FUNC, ht, pr, vf, grad, hes);
    if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
    if(free_hes) delete[] hes; //   function
    return vf;
}

s_real vfs_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes){
    s_real val = memoize::get_bin(memoize::VFS_FUNC, st, pr, grad, hes);
    if(!std::isnan(val)) return val;
    bool free_grad = 0, free_hes = 0; // if 1 free grad and hes at end
    // Since I'm going to cache results, grad and hes will get calculated
    // whether requested or not.  If they were NULL allocate space.
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

    s_real tau, gradt[1], hest[1], gradsv[2], hessv[3], gradsl[2], hessl[3];
    tau = sat_tau_with_derivs(pr, gradt, hest);
    s_real sv = svpt_with_derivs(pr, tau, gradsv, hessv);
    s_real sl = slpt_with_derivs(pr, tau, gradsl, hessl);

    if(pr >= P_c){zero_derivs2(grad, hes); return 0.0;} //Classify supercritical as liquid
    else if(pr <= P_t){zero_derivs2(grad, hes); return 1.0;} //below tripple point going to assume vapor, because it's that or a solid and I'm not dealing with solids.
    else if(sl > st){ zero_derivs2(grad, hes); return 0.0;} //saturated liquid entropy > total entropy so liquid
    else if(sv < st){ zero_derivs2(grad, hes); return 1.0;} //saturated vapor enthalpy < total entropy so vapor

    s_real dsvdp = gradsv[0] + gradsv[1]*gradt[0];
    s_real dsldp = gradsl[0] + gradsl[1]*gradt[0];
    grad[0] = 1.0/(sv - sl);
    grad[1] = -dsldp/(sv - sl) - (st - sl)/(sv-sl)/(sv-sl)*(dsvdp - dsldp);
    s_real d2svdp2 = hessv[0] + 2*hessv[1]*gradt[0] + hessv[2]*gradt[0]*gradt[0] + gradsv[1]*hest[0];
    s_real d2sldp2 = hessl[0] + 2*hessl[1]*gradt[0] + hessl[2]*gradt[0]*gradt[0] + gradsl[1]*hest[0];
    hes[0] = 0;
    hes[1] = -1.0/(sv-sl)/(sv-sl)*(dsvdp - dsldp);
    hes[2] = -d2sldp2/(sv-sl) + 2*dsldp/(sv-sl)/(sv-sl)*(dsvdp - dsldp) +
              2*(st-sl)/(sv-sl)/(sv-sl)/(sv-sl)*(dsvdp - dsldp)*(dsvdp - dsldp) -
              (st-sl)/(sv-sl)/(sv-sl)*(d2svdp2 - d2sldp2);
    s_real vf = (st - sl)/(sv - sl);
    memoize::add_bin(memoize::VFS_FUNC, st, pr, vf, grad, hes);
    if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
    if(free_hes) delete[] hes; //   function
    return vf;
}

s_real vfu_with_derivs(s_real st, s_real pr, s_real *grad, s_real *hes){
    s_real val = memoize::get_bin(memoize::VFU_FUNC, st, pr, grad, hes);
    if(!std::isnan(val)) return val;
    bool free_grad = 0, free_hes = 0; // if 1 free grad and hes at end
    // Since I'm going to cache results, grad and hes will get calculated
    // whether requested or not.  If they were NULL allocate space.
    if(grad==NULL){grad = new s_real[2]; free_grad = 1;}
    if(hes==NULL){hes = new s_real[3]; free_hes = 1;}

    s_real tau, gradt[1], hest[1], gradsv[2], hessv[3], gradsl[2], hessl[3];
    tau = sat_tau_with_derivs(pr, gradt, hest);
    s_real uv = uvpt_with_derivs(pr, tau, gradsv, hessv);
    s_real ul = ulpt_with_derivs(pr, tau, gradsl, hessl);

    if(pr >= P_c){zero_derivs2(grad, hes); return 0.0;} //Classify supercritical as liquid
    else if(pr <= P_t){zero_derivs2(grad, hes); return 1.0;} //below tripple point going to assume vapor, because it's that or a solid and I'm not dealing with solids.
    else if(ul > ut){ zero_derivs2(grad, hes); return 0.0;} //saturated liquid entropy > total entropy so liquid
    else if(uv < ut){ zero_derivs2(grad, hes); return 1.0;} //saturated vapor enthalpy < total entropy so vapor

    s_real dsvdp = gradsv[0] + gradsv[1]*gradt[0];
    s_real dsldp = gradsl[0] + gradsl[1]*gradt[0];
    grad[0] = 1.0/(uv - ul);
    grad[1] = -dsldp/(uv - ul) - (ut - ul)/(uv-ul)/(uv-ul)*(dsvdp - dsldp);
    s_real d2svdp2 = hessv[0] + 2*hessv[1]*gradt[0] + hessv[2]*gradt[0]*gradt[0] + gradsv[1]*hest[0];
    s_real d2sldp2 = hessl[0] + 2*hessl[1]*gradt[0] + hessl[2]*gradt[0]*gradt[0] + gradsl[1]*hest[0];
    hes[0] = 0;
    hes[1] = -1.0/(uv-ul)/(uv-ul)*(dsvdp - dsldp);
    hes[2] = -d2sldp2/(uv-ul) + 2*dsldp/(uv-ul)/(uv-ul)*(dsvdp - dsldp) +
              2*(ut-ul)/(uv-sl)/(uv-ul)/(uv-ul)*(dsvdp - dsldp)*(dsvdp - dsldp) -
              (st-sl)/(sv-sl)/(sv-sl)*(d2svdp2 - d2sldp2);
    s_real vf = (ut - ul)/(uv - ul);
    memoize::add_bin(memoize::VFU_FUNC, ut, pr, vf, grad, hes);
    if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
    if(free_hes) delete[] hes; //   function
    return vf;
}

s_real g_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  s_real h_grad[2], h_hes[3], s_grad[2], s_hes[3], s_val;
  s_real *h_grad_ptr=NULL, *h_hes_ptr=NULL, *s_grad_ptr=NULL, *s_hes_ptr=NULL;

  if(grad!=NULL){
    h_grad_ptr=h_grad;
    s_grad_ptr=s_grad;
    if(hes!=NULL){
      h_hes_ptr=h_hes;
      s_hes_ptr=s_hes;}}

  if(grad!=NULL){
    h_with_derivs(delta, tau, h_grad_ptr, h_hes_ptr);
    s_val = s_with_derivs(delta, tau, s_grad_ptr, s_hes_ptr);
    grad[0] = h_grad[0] - T_c/tau*s_grad[0];
    grad[1] = h_grad[1] + T_c/tau/tau*s_val - T_c/tau*s_grad[1];
    if(hes!=NULL){
      hes[0] = h_hes[0] - T_c/tau*s_hes[0];
      hes[1] = h_hes[1] + T_c/tau/tau*s_grad[0] - T_c/tau*s_hes[1];
      hes[2] = h_hes[2] - 2.0*T_c/tau/tau/tau*s_val + 2.0*T_c/tau/tau*s_grad[1]
        - T_c/tau*s_hes[2];}}
  return g(delta, tau);
}

s_real f_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = XB*(phi0_delta(delta) + phir_delta(delta,tau));
    grad[1] = XB_t*(phi0(delta, tau) + phir(delta,tau)) +
      XB*(phi0_tau(tau) + phir_tau(delta,tau));
    if(hes!=NULL){
      hes[0] = XB*(phi0_delta2(delta) + phir_delta2(delta,tau));
      hes[1] = XB_t*(phi0_delta(delta) + phir_delta(delta,tau)) +
        XB*phir_delta_tau(delta,tau);
      hes[2] = XB_tt*(phi0(delta, tau) + phir(delta,tau)) +
        2.0*XB_t*(phi0_tau(tau) + phir_tau(delta,tau)) +
        XB*(phi0_tau2(tau) + phir_tau2(delta,tau));}}
  return f(delta, tau);
}

s_real cv_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  if(grad!=NULL){
    grad[0] = -tau*tau*R*(phi0_delta_tau2 + phir_delta_tau2(delta, tau));
    grad[1] = -2.0*tau*R*(phi0_tau2(tau) + phir_tau2(delta, tau)) -
      tau*tau*R*(phi0_tau3(tau) + phir_tau3(delta, tau));
    if(hes!=NULL){
      hes[0] = -tau*tau*R*(phi0_delta2_tau2 + phir_delta2_tau2(delta, tau));
      hes[1] = -2.0*tau*R*(phi0_delta_tau2 + phir_delta_tau2(delta, tau)) -
        tau*tau*R*(phi0_delta_tau3 + phir_delta_tau3(delta, tau));
      hes[2] = -2.0*R*(phi0_tau2(tau) + phir_tau2(delta, tau)) -
        4.0*tau*R*(phi0_tau3(tau) + phir_tau3(delta, tau)) -
        tau*tau*R*(phi0_tau4(tau) + phir_tau4(delta, tau));}}
  return cv(delta, tau);
}

s_real cp_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  s_real cv_grad[2], cv_hes[3];
  s_real *cv_grad_ptr=NULL, *cv_hes_ptr=NULL;
  if(grad!=NULL){
    cv_grad_ptr = cv_grad;
    if(hes!=NULL){
      cv_hes_ptr = cv_hes;}}
  cv_with_derivs(delta, tau, cv_grad_ptr, cv_hes_ptr);
  if(grad!=NULL){
    grad[0] = cv_grad[0] + R*XD_d/XE - R*XD*XE_d/XE/XE;
    grad[1] = cv_grad[1] + R*XD_t/XE - R*XD*XE_t/XE/XE;
    if(hes!=NULL){
      hes[0] = cv_hes[0] + R*XD_dd/XE - 2.0*R*XD_d*XE_d/XE/XE -
        R*XD*XE_dd/XE/XE + 2.0*R*XD*XE_d*XE_d/XE/XE/XE;
      hes[1] = cv_hes[1] + R*XD_dt/XE - R*XD_d*XE_t/XE/XE - R*XD_t*XE_d/XE/XE -
        R*XD*XE_dt/XE/XE + 2.0*R*XD*XE_d*XE_t/XE/XE/XE;
      hes[2] = cv_hes[2] + R*XD_tt/XE - 2.0*R*XD_t*XE_t/XE/XE -
        R*XD*XE_tt/XE/XE + 2.0*R*XD*XE_t*XE_t/XE/XE/XE;}}
  return cp(delta, tau);
}

s_real w_with_derivs(s_real delta, s_real tau, s_real *grad, s_real *hes){
  s_real ww = w(delta, tau), cv;
  s_real cv_grad[2], cv_hes[3];
  s_real *cv_grad_ptr=NULL, *cv_hes_ptr=NULL;
  if(grad!=NULL){
    cv_grad_ptr = cv_grad;
    if(hes!=NULL){
      cv_hes_ptr = cv_hes;}}
  cv = cv_with_derivs(delta, tau, cv_grad_ptr, cv_hes_ptr);
  if(grad!=NULL){
    s_real A_t, A_tt, A_d, A_dd, A_dt;
    A_t = -R*T_c/tau/tau*(XE + XD*R/cv) +
           R*T_c/tau*(XE_t + R*XD_t/cv - cv_grad[1]*R*XD/(cv*cv));
    A_d =  R*T_c/tau*(XE_d + R*XD_d/cv - cv_grad[0]*R*XD/(cv*cv));
    A_tt = 2*R*T_c/(tau*tau*tau)*(XE + XD*R/cv) -
           2*R*T_c/(tau*tau)*(XE_t + R*XD_t/cv - cv_grad[1]*R*XD/(cv*cv)) +
           R*T_c/tau*(XE_tt +
                      R*XD_tt/cv -
                      2*cv_grad[1]*R*XD_t/(cv*cv) -
                      R*XD*cv_hes[2]/(cv*cv) +
                      2*cv_grad[1]*cv_grad[1]*R*XD/(cv*cv*cv));
    A_dd = R*T_c/tau*(XE_dd +
                      R*XD_dd/cv -
                      2*cv_grad[0]*R*XD_d/(cv*cv) -
                      R*XD*cv_hes[0]/(cv*cv) +
                      2*cv_grad[0]*cv_grad[0]*R*XD/(cv*cv*cv));
    A_dt = -R*T_c/tau/tau*(XE_d + R*XD_d/cv - cv_grad[0]*R*XD/(cv*cv)) +
            R*T_c/tau*(XE_dt +
                       R*XD_dt/cv -
                       cv_grad[1]*R*XD_d/(cv*cv) -
                       cv_grad[0]*R*XD_t/(cv*cv) -
                       R*XD*cv_hes[1]/(cv*cv) +
                       2*cv_grad[0]*cv_grad[1]*R*XD/(cv*cv*cv));
    grad[0] = 1000*(0.5/ww) * A_d; //the 1000 if from km to m
    grad[1] = 1000*(0.5/ww) * A_t; //the 1000 if from km to m
    if(hes!=NULL){
      hes[0] = -1e6*(0.25/(ww*ww*ww))*A_d*A_d + 1000*(0.5/ww)*A_dd;
      hes[1] = -1e6*(0.25/(ww*ww*ww))*A_t*A_d + 1000*(0.5/ww)*A_dt;
      hes[2] = -1e6*(0.25/(ww*ww*ww))*A_t*A_t + 1000*(0.5/ww)*A_tt;}}
  return ww;
}

s_real sat_delta_liq_with_derivs(s_real tau, s_real *grad, s_real *hes){
  // For efficiency I'm requiring memoization here.  Delta_v and delta_l are
  // calculated together so should just remember both.  This will break if you
  // turn off memoization.
  s_real delta_l, delta_v;
  s_real val = memoize::get_un(memoize::DL_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  sat(tau, &delta_l, &delta_v);
  // Should be there now
  val = memoize::get_un(memoize::DL_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  return (s_real)NAN;
}

s_real sat_delta_vap_with_derivs(s_real tau, s_real *grad, s_real *hes){
  // For efficiency I'm requiring memoization here.  Delta_v and delta_l are
  // calculated together so should just remember both.  This will break if you
  // turn off memoization.
  s_real delta_l, delta_v;
  s_real val = memoize::get_un(memoize::DV_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  sat(tau, &delta_l, &delta_v);
  // Should be there now
  val = memoize::get_un(memoize::DV_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  return (s_real)NAN;
}

s_real sat_p_with_derivs(s_real tau, s_real *grad, s_real *hes, bool limit){
  if(tau < 1.0 && limit){ // above critical T
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return P_c;
  }
  if(tau > T_c/T_t && limit){ // below tripple
    if(grad != NULL) grad[0] = 0;
    if(hes != NULL) hes[0] = 0;
    return P_t;
  }
  //Now check if there is a stored value
  s_real val = memoize::get_un(memoize::P_SAT_FUNC, tau, grad, hes);
  if(!std::isnan(val)) return val;
  // grad and/or hes not provided so allocate for memoization
  bool free_grad = 0, free_hes = 0;
  if(grad==NULL){grad = new s_real[1]; free_grad = 1;}
  if(hes==NULL){hes = new s_real[1]; free_hes = 1;}
  s_real grad_delta[1], hes_delta[1], grad_p[2], hes_p[3];
  s_real delta = sat_delta_vap_with_derivs(tau, grad_delta, hes_delta);
  s_real Psat = p_with_derivs(delta, tau, grad_p, hes_p);
  grad[0] = grad_p[1] + grad_p[0]*grad_delta[0];
  hes[0] = hes_p[2] + 2*hes_p[1]*grad_delta[0] + grad_p[0]*hes_delta[0]
           + hes_p[0]*grad_delta[0]*grad_delta[0];
  memoize::add_un(memoize::P_SAT_FUNC, tau, Psat, grad, hes);
  if(free_grad) delete[] grad; // free grad and hes if not allocated by calling
  if(free_hes) delete[] hes;   // function
  return Psat;
}
