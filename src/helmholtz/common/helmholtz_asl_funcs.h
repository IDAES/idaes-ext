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
 This provides the ASL function interface for Helmholtz equation of state
 property functions.

 Author: John Eslick
 File: helmholtz_asl_funcs.h
------------------------------------------------------------------------------*/

#include <funcadd.h>

#ifndef _INCLUDE_HELMHOLTZ_ASL_FUNCS_H_
#define _INCLUDE_HELMHOLTZ_ASL_FUNCS_H_

double p_EOS_TAG(arglist *al);
double u_EOS_TAG(arglist *al);
double s_EOS_TAG(arglist *al);
double h_EOS_TAG(arglist *al);
double g_EOS_TAG(arglist *al);
double f_EOS_TAG(arglist *al);
double cv_EOS_TAG(arglist *al);
double cp_EOS_TAG(arglist *al);
double w_EOS_TAG(arglist *al);

double hvpt_EOS_TAG(arglist *al);
double hlpt_EOS_TAG(arglist *al);
double svpt_EOS_TAG(arglist *al);
double slpt_EOS_TAG(arglist *al);
double uvpt_EOS_TAG(arglist *al);
double ulpt_EOS_TAG(arglist *al);
double vf_EOS_TAG(arglist *al);
double vfs_EOS_TAG(arglist *al);
double vfu_EOS_TAG(arglist *al);
double tau_EOS_TAG(arglist *al);
double tau_sp_EOS_TAG(arglist *al);
double tau_up_EOS_TAG(arglist *al);

double delta_sat_l_EOS_TAG(arglist *al);
double delta_sat_v_EOS_TAG(arglist *al);
double p_sat_EOS_TAG(arglist *al);
double tau_sat_EOS_TAG(arglist *al);

double delta_liq_EOS_TAG(arglist *al);
double delta_vap_EOS_TAG(arglist *al);

double phi0_EOS_TAG(arglist *al);
double phi0_delta_EOS_TAG(arglist *al);
double phi0_delta2_EOS_TAG(arglist *al);
double phi0_tau_EOS_TAG(arglist *al);
double phi0_tau2_EOS_TAG(arglist *al);
double phir_EOS_TAG(arglist *al);
double phir_delta_EOS_TAG(arglist *al);
double phir_delta2_EOS_TAG(arglist *al);
double phir_tau_EOS_TAG(arglist *al);
double phir_tau2_EOS_TAG(arglist *al);
double phir_delta_tau_EOS_TAG(arglist *al);

#endif
