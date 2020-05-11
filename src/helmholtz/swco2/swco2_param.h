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
 This file provides the Span-Wagner parameters.

 References:
   Span, R., and W. Wanger (1996). "A New Equation of State for Carbon Dioxide
       Covering the Fluid Region from the Triple-Point Temperature to 1100 K as
       Pressures up to 800 MPa." Journal of Physical and Chemical Reference Data,
       25, 1509.
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
       State from Helmholtz Energy Equations of State." Journal of Thermal
       Science and Technology, 3(3), 442-451.

 Author: John Eslick
 File: swco2_param.h
------------------------------------------------------------------------------*/

#ifndef _INCLUDE_SWCO2_PARAM_H_
#define _INCLUDE_SWCO2_PARAM_H_

#include "helmholtz_config.h"
#include "swco2_guess.h"

#define LIQUID_DELTA_GUESS delta_p_tau_liq_guess_swco2(p, tau)
#define VAPOR_DELTA_GUESS delta_p_tau_vap_guess_swco2(p, tau)
#define DELTA_LIQ_SAT_GUESS delta_sat_l_approx(tau)
#define DELTA_VAP_SAT_GUESS delta_sat_v_approx(tau)

#define TAU_LOW 0.15
#define TAU_HIGH 4.0
#define P_LOW 0.0001
#define P_HIGH 1e10

const s_real R = 0.1889241;  // Specific gas constant (kJ/kg/K)

// Critiacal point for water and R
const s_real T_c = 304.128;   // Critical T (K)
const s_real rho_c = 467.6;   // Critical density (kg/m^3)
const s_real P_c = 7377.3;    // Critical Pressure (kPa)
const s_real T_t = 216.592;   // Triple point temperature (K)
const s_real P_t = 517.95;    // Triple point pressure (kPa)
const s_real P_max = 100*P_c;  // Max pressure where answer is sure to be right
const s_real P_min = 1;       // Min pressure where answer is sure to be right
const s_real T_max = 800;     // Max temp where answer is sure to be right
const s_real T_min = T_t;     // Min temp where answer is sure to be right

// To generalize the equation of state there are parameters to set the number
// of terms in each summation.  So far we are looking at IAPWS95 and Span-Wagner
// They have the same types of terms but different numbers of term.  For both,
// we are currently dropping the nonanalytic terms, which we've labed S4 there
// See docs for the form of the terms in S1, S2, S3, and S4.

const unsigned char S1_set[2] = {1, 7};
const unsigned char S2_set[2] = {8, 34};
const unsigned char S3_set[2] = {35, 39};
const unsigned char S4_set[2] = {40, 42};  // we don't currently use these

//
// Constants
//

const s_real param[] = {
  0,	//padding	0
  8.37304456,	//n0_1	1	0 index at 0
  -3.70454304,	//n0_2	2
  2.5,	//n0_3	3
  1.99427042,	//n0_4	4
  0.62105248,	//n0_5	5
  0.41195293,	//n0_6	6
  1.04028922,	//n0_7	7
  0.08327678,	//n0_8	8
  3.15163,	//gamma0_4	9	0 index at 5
  6.1119,	//gamma0_5	10
  6.77708,	//gamma0_6	11
  11.32384,	//gamma0_7	12
  27.08792,	//gamma0_8	13
  3.88568232031610E-01,	//n_1	14 0 index at 13
  2.93854759427400E+00,	//n_2	15
  -5.58671885349340E+00,	//n_3	16
  -7.67531995924770E-01,	//n_4	17
  3.17290055804160E-01,	//n_5	18
  5.48033158977670E-01,	//n_6	19
  1.22794112203350E-01,	//n_7	20
  2.16589615432200E+00,	//n_8	21
  1.58417351097240E+00,	//n_9	22
  -2.31327054055030E-01,	//n_10	23
  5.81169164314360E-02,	//n_11	24
  -5.53691372053820E-01,	//n_12	25
  4.89466159094220E-01,	//n_13	26
  -2.42757398435010E-02,	//n_14	27
  6.24947905016780E-02,	//n_15	28
  -1.21758602252460E-01,	//n_16	29
  -3.70556852700860E-01,	//n_17	30
  -1.67758797004260E-02,	//n_18	31
  -1.19607366379870E-01,	//n_19	32
  -4.56193625087780E-02,	//n_20	33
  3.56127892703460E-02,	//n_21	34
  -7.44277271320520E-03,	//n_22	35
  -1.73957049024320E-03,	//n_23	36
  -2.18101212895270E-02,	//n_24	37
  2.43321665592360E-02,	//n_25	38
  -3.74401334234630E-02,	//n_26	39
  1.43387157568780E-01,	//n_27	40
  -1.34919690832860E-01,	//n_28	41
  -2.31512250534800E-02,	//n_29	42
  1.23631254929010E-02,	//n_30	43
  2.10583219729400E-03,	//n_31	44
  -3.39585190263680E-04,	//n_32	45
  5.59936517715920E-03,	//n_33	46
  -3.03351180556460E-04,	//n_34	47
  -2.13654886883200E+02,	//n_35	48
  2.66415691492720E+04,	//n_36	49
  -2.40272122045570E+04,	//n_37	50
  -2.83416034239990E+02,	//n_38	51
  2.12472844001790E+02,	//n_39	52
  -6.66422765407510E-01,	//n_40	53
  7.26086323498970E-01,	//n_41	54
  5.50686686128420E-02,	//n_42	55
  1,	//d1	56	0 index at 55
  1,	//d2	57
  1,	//d3	58
  1,	//d4	59
  2,	//d5	60
  2,	//d6	61
  3,	//d7	62
  1,	//d8	63
  2,	//d9	64
  4,	//d10	65
  5,	//d11	66
  5,	//d12	67
  5,	//d13	68
  6,	//d14	69
  6,	//d15	70
  6,	//d16	71
  1,	//d17	72
  1,	//d18	73
  4,	//d19	74
  4,	//d20	75
  4,	//d21	76
  7,	//d22	77
  8,	//d23	78
  2,	//d24	79
  3,	//d25	80
  3,	//d26	81
  5,	//d27	82
  5,	//d28	83
  6,	//d29	84
  7,	//d30	85
  8,	//d31	86
  10,	//d32	87
  4,	//d33	88
  8,	//d34	89
  2,	//d35	90
  2,	//d36	91
  2,	//d37	92
  3,	//d38	93
  3,	//d39	94
  0.00,	//t1	95	0 index at 94
  0.75,	//t2	96
  1.00,	//t3	97
  2.00,	//t4	98
  0.75,	//t5	99
  2.00,	//t6	100
  0.75,	//t7	101
  1.50,	//t8	102
  1.50,	//t9	103
  2.50,	//t10	104
  0.00,	//t11	105
  1.50,	//t12	106
  2.00,	//t13	107
  0.00,	//t14	108
  1.00,	//t15	109
  2.00,	//t16	110
  3.00,	//t17	111
  6.00,	//t18	112
  3.00,	//t19	113
  6.00,	//t20	114
  8.00,	//t21	115
  6.00,	//t22	116
  0.00,	//t23	117
  7.00,	//t24	118
  12.00,	//t25	119
  16.00,	//t26	120
  22.00,	//t27	121
  24.00,	//t28	122
  16.00,	//t29	123
  24.00,	//t30	124
  8.00,	//t31	125
  2.00,	//t32	126
  28.00,	//t33	127
  14.00,	//t34	128
  1.00,	//t35	129
  0.00,	//t36	130
  1.00,	//t37	131
  3.00,	//t38	132
  3.00,	//t39	133
  1,	//c8	134	0 index at 126
  1,	//c9	135
  1,	//c10	136
  1,	//c11	137
  1,	//c12	138
  1,	//c13	139
  1,	//c14	140
  1,	//c15	141
  1,	//c16	142
  2,	//c17	143
  2,	//c18	144
  2,	//c19	145
  2,	//c20	146
  2,	//c21	147
  2,	//c22	148
  2,	//c23	149
  3,	//c24	150
  3,	//c25	151
  3,	//c26	152
  4,	//c27	153
  4,	//c28	154
  4,	//c29	155
  4,	//c30	156
  4,	//c31	157
  4,	//c32	158
  5,	//c33	159
  6,	//c34	160
  25,	//alpha35	161	0 index at 126
  25,	//alpha36	162
  25,	//alpha37	163
  15,	//alpha38	164
  20,	//alpha39	165
  325,	//beta35	166	0 index at 131
  300,	//beta36	167
  300,	//beta37	168
  275,	//beta38	169
  275,	//beta39	170
  1.16,	//theta35	171	0 index at 136
  1.19,	//theta36	172
  1.19,	//theta37	173
  1.25,	//theta38	174
  1.22,	//theta39	175
  1.00,	//eps35	176	0 index at 141
  1.00,	//eps36	177
  1.00,	//eps37	178
  1.00,	//eps38	179
  1.00	//eps39	180
};

// Offset for the zero index of each parameter in the param array
const unsigned int n0_offset = 0;
const unsigned int gamma0_offset = 5;
const unsigned int c_offset = 126;
const unsigned int d_offset = 55;
const unsigned int t_offset = 94;
const unsigned int n_offset = 13;
const unsigned int alpha_offset = 126;
const unsigned int theta_offset = 136;
const unsigned int eps_offset = 141;
const unsigned int beta_offset = 131;

#endif
