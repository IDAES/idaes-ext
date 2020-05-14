/*------------------------------------------------------------------------------
 Institute for the Design of Advanced Energy Systems Process Systems
 Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
 software owners: The Regents of the University of California, through
 Lawrence Berkeley National Laboratory,  National Technology & Engineering
 Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
 University Research Corporation, et al. All rights reserved.

 Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
 license information, respectively. Both files are also available online
 at the URL "https://github.com/IDAES/idaes".
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 This provides the ASL function interface for Helmholtz equations of state.

 Author: John Eslick
 File: helmholtz_asl_funcs.cpp
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include"helmholtz_external.h"
#include"helmholtz_phi.h"
#include"helmholtz_asl_funcs.h"
#include"helmholtz_config.h"

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info */
    int typ = FUNCADD_REAL_VALUED;
    addfunc("p_EOS_TAG", (rfunc)p_EOS_TAG, typ, 2, NULL);
    addfunc("u_EOS_TAG", (rfunc)u_EOS_TAG, typ, 2, NULL);
    addfunc("s_EOS_TAG", (rfunc)s_EOS_TAG, typ, 2, NULL);
    addfunc("h_EOS_TAG", (rfunc)h_EOS_TAG, typ, 2, NULL);
    addfunc("g_EOS_TAG", (rfunc)g_EOS_TAG, typ, 2, NULL);
    addfunc("f_EOS_TAG", (rfunc)f_EOS_TAG, typ, 2, NULL);
    addfunc("cv_EOS_TAG", (rfunc)cv_EOS_TAG, typ, 2, NULL);
    addfunc("cp_EOS_TAG", (rfunc)cp_EOS_TAG, typ, 2, NULL);
    addfunc("w_EOS_TAG", (rfunc)w_EOS_TAG, typ, 2, NULL);
    addfunc("hvpt_EOS_TAG", (rfunc)hvpt_EOS_TAG, typ, 2, NULL);
    addfunc("hlpt_EOS_TAG", (rfunc)hlpt_EOS_TAG, typ, 2, NULL);
    addfunc("svpt_EOS_TAG", (rfunc)svpt_EOS_TAG, typ, 2, NULL);
    addfunc("slpt_EOS_TAG", (rfunc)slpt_EOS_TAG, typ, 2, NULL);
    addfunc("uvpt_EOS_TAG", (rfunc)uvpt_EOS_TAG, typ, 2, NULL);
    addfunc("ulpt_EOS_TAG", (rfunc)ulpt_EOS_TAG, typ, 2, NULL);
    addfunc("tau_EOS_TAG", (rfunc)tau_EOS_TAG, typ, 2, NULL);
    addfunc("memo_test_tau_EOS_TAG", (rfunc)memo_test_tau_EOS_TAG, typ, 2, NULL);
    addfunc("tau_sp_EOS_TAG", (rfunc)tau_sp_EOS_TAG, typ, 2, NULL);
    addfunc("tau_up_EOS_TAG", (rfunc)tau_up_EOS_TAG, typ, 2, NULL);
    addfunc("p_stau_EOS_TAG", (rfunc)p_stau_EOS_TAG, typ, 2, NULL);
    addfunc("vf_EOS_TAG", (rfunc)vf_EOS_TAG, typ, 2, NULL);
    addfunc("vfs_EOS_TAG", (rfunc)vfs_EOS_TAG, typ, 2, NULL);
    addfunc("vfu_EOS_TAG", (rfunc)vfu_EOS_TAG, typ, 2, NULL);
    addfunc("delta_liq_EOS_TAG", (rfunc)delta_liq_EOS_TAG, typ, 2, NULL);
    addfunc("delta_vap_EOS_TAG", (rfunc)delta_vap_EOS_TAG, typ, 2, NULL);
    addfunc("delta_sat_l_EOS_TAG", (rfunc)delta_sat_l_EOS_TAG, typ, 1, NULL);
    addfunc("delta_sat_v_EOS_TAG", (rfunc)delta_sat_v_EOS_TAG, typ, 1, NULL);
    addfunc("p_sat_EOS_TAG", (rfunc)p_sat_EOS_TAG, typ, 1, NULL);
    addfunc("tau_sat_EOS_TAG", (rfunc)tau_sat_EOS_TAG, typ, 1, NULL);
    addfunc("phi0_EOS_TAG", (rfunc)phi0_EOS_TAG, typ, 2, NULL);
    addfunc("phi0_delta_EOS_TAG", (rfunc)phi0_delta_EOS_TAG, typ, 1, NULL);
    addfunc("phi0_delta2_EOS_TAG", (rfunc)phi0_delta2_EOS_TAG, typ, 1, NULL);
    addfunc("phi0_tau_EOS_TAG", (rfunc)phi0_tau_EOS_TAG, typ, 1, NULL);
    addfunc("phi0_tau2_EOS_TAG", (rfunc)phi0_tau2_EOS_TAG, typ, 1, NULL);
    addfunc("phir_EOS_TAG", (rfunc)phir_EOS_TAG, typ, 2, NULL);
    addfunc("phir_delta_EOS_TAG", (rfunc)phir_delta_EOS_TAG, typ, 2, NULL);
    addfunc("phir_delta2_EOS_TAG", (rfunc)phir_delta2_EOS_TAG, typ, 2, NULL);
    addfunc("phir_tau_EOS_TAG", (rfunc)phir_tau_EOS_TAG, typ, 2, NULL);
    addfunc("phir_tau2_EOS_TAG", (rfunc)phir_tau2_EOS_TAG, typ, 2, NULL);
    addfunc("phir_delta_tau_EOS_TAG", (rfunc)phir_delta_tau_EOS_TAG, typ, 2, NULL);
}

inline void cast_deriv2(s_real *g1, double *g2, s_real *h1, double *h2){
  if(g2 != NULL){
    g2[0] = (double)g1[0];
    g2[1] = (double)g1[1];
  }
  if(h2 != NULL){
    h2[0] = (double)h1[0];
    h2[1] = (double)h1[1];
    h2[2] = (double)h1[2];
  }
}

inline void cast_deriv1(s_real *g1, double *g2, s_real *h1, double *h2){
  if(g2 != NULL) g2[0] = (double)g1[0];
  if(h2 != NULL) h2[0] = (double)h1[0];
}

double p_EOS_TAG(arglist *al){
  return p_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double u_EOS_TAG(arglist *al){
  return u_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double s_EOS_TAG(arglist *al){
  return s_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double h_EOS_TAG(arglist *al){
  return h_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double g_EOS_TAG(arglist *al){
  return g_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double f_EOS_TAG(arglist *al){
  return f_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double cv_EOS_TAG(arglist *al){
  return cv_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double cp_EOS_TAG(arglist *al){
  return cp_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double w_EOS_TAG(arglist *al){
  return w_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double hvpt_EOS_TAG(arglist *al){
  return hvpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double hlpt_EOS_TAG(arglist *al){
  return hlpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double svpt_EOS_TAG(arglist *al){
  return svpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double slpt_EOS_TAG(arglist *al){
  return slpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double uvpt_EOS_TAG(arglist *al){
  return uvpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double ulpt_EOS_TAG(arglist *al){
  return ulpt_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double tau_EOS_TAG(arglist *al){
    return tau_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double memo_test_tau_EOS_TAG(arglist *al){
  return mem_tau_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double tau_sp_EOS_TAG(arglist *al){
  return tau_from_sp_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double tau_up_EOS_TAG(arglist *al){
  return tau_from_up_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double p_stau_EOS_TAG(arglist *al){
  return p_from_stau_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double vf_EOS_TAG(arglist *al){
  return vf_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double vfs_EOS_TAG(arglist *al){
  return vfs_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double vfu_EOS_TAG(arglist *al){
  return vfu_with_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double delta_sat_l_EOS_TAG(arglist *al){
  return sat_delta_liq_with_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double delta_sat_v_EOS_TAG(arglist *al){
  return sat_delta_vap_with_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double p_sat_EOS_TAG(arglist *al){
  return sat_p_with_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double tau_sat_EOS_TAG(arglist *al){
  return sat_tau_with_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double delta_liq_EOS_TAG(arglist *al){
  return delta_liq(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double delta_vap_EOS_TAG(arglist *al){
  return delta_vap(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

//with deriv functions for phi, only provided for testing through ASL interface
double phi0_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta(delta);
    grad[1] = phi0_tau(tau);
    if(hes != NULL){
      hes[0] = phi0_delta2(delta);
      hes[1] = 0;
      hes[2] = phi0_tau2(tau);
    }
  }
  return phi0(delta, tau);
}

double phi0_delta_derivs(double delta, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta2(delta);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_delta3(delta);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_delta(delta);
}

double phi0_delta2_derivs(double delta, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_delta3(delta);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_delta4(delta);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_delta2(delta);
}

double phi0_tau_derivs(double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_tau2(tau);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_tau3(tau);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_tau(tau);
}

double phi0_tau2_derivs(double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phi0_tau3(tau);
    grad[1] = 0;
    if(hes != NULL){
      hes[0] = phi0_tau4(tau);
      hes[1] = 0;
      hes[2] = 0;
    }
  }
  return phi0_tau2(tau);
}

double phir_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta(delta, tau);
    grad[1] = phir_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2(delta, tau);
      hes[1] = phir_delta_tau(delta, tau);
      hes[2] = phir_tau2(delta, tau);
    }
  }
  return phir(delta, tau);
}

double phir_delta_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta2(delta, tau);
    grad[1] = phir_delta_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta3(delta, tau);
      hes[1] = phir_delta2_tau(delta, tau);
      hes[2] = phir_delta_tau2(delta, tau);
    }
  }
  return phir_delta(delta, tau);
}

double phir_delta2_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta3(delta, tau);
    grad[1] = phir_delta2_tau(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta4(delta, tau);
      hes[1] = phir_delta3_tau(delta, tau);
      hes[2] = phir_delta2_tau2(delta, tau);
    }
  }
  return phir_delta2(delta, tau);
}

double phir_tau_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta_tau(delta, tau);
    grad[1] = phir_tau2(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2_tau(delta, tau);
      hes[1] = phir_delta_tau2(delta, tau);
      hes[2] = phir_tau3(delta, tau);
    }
  }
  return phir_tau(delta, tau);
}

double phir_tau2_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta_tau2(delta, tau);
    grad[1] = phir_tau3(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta2_tau2(delta, tau);
      hes[1] = phir_delta_tau3(delta, tau);
      hes[2] = phir_tau4(delta, tau);
    }
  }
  return phir_tau2(delta, tau);
}

double phir_delta_tau_derivs(double delta, double tau, double *grad, double *hes){
  if(grad != NULL){
    grad[0] = phir_delta2_tau(delta, tau);
    grad[1] = phir_delta_tau2(delta, tau);
    if(hes != NULL){
      hes[0] = phir_delta3_tau(delta, tau);
      hes[1] = phir_delta2_tau2(delta, tau);
      hes[2] = phir_delta_tau3(delta, tau);
    }
  }
  return phir_delta_tau(delta, tau);
}

double phi0_EOS_TAG(arglist *al){
  return phi0_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phi0_delta_EOS_TAG(arglist *al){
  return phi0_delta_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double phi0_delta2_EOS_TAG(arglist *al){
  return phi0_delta2_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double phi0_tau_EOS_TAG(arglist *al){
  return phi0_tau_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double phi0_tau2_EOS_TAG(arglist *al){
  return phi0_tau2_derivs(al->ra[al->at[0]], al->derivs, al->hes);
}

double phir_EOS_TAG(arglist *al){
  return phir_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phir_delta_EOS_TAG(arglist *al){
  return phir_delta_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phir_delta2_EOS_TAG(arglist *al){
  return phir_delta2_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phir_tau_EOS_TAG(arglist *al){
  return phir_tau_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phir_tau2_EOS_TAG(arglist *al){
  return phir_tau2_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}

double phir_delta_tau_EOS_TAG(arglist *al){
  return phir_delta_tau_derivs(al->ra[al->at[0]], al->ra[al->at[1]], al->derivs, al->hes);
}
