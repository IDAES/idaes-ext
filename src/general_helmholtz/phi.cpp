//#include <funcadd.h>
#include <adolc/adolc.h>
#include<unordered_map>
#include<boost/functional/hash.hpp>
#include"h2o.h"
#include"config.h"
#include"param.h"
#include"phi.h"
#include"function_pointers.h"

unsigned int taped_ideal[NCOMPS] = {0};
unsigned int taped_real[NCOMPS] = {0};

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_phi_real;

std::unordered_map<
  std::tuple<comp_enum, double, double>,
  std::vector<double>,
  boost::hash<std::tuple<comp_enum, double, double>>
> memo_table_phi_ideal;

std::vector<double> *phi_real(comp_enum comp, double delta, double tau){
  /*
    Calculate the dimensionless real part of Helmholtz free energy (phi^r), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */
  try{
    return &memo_table_phi_real.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  double x[2] = {delta, tau};
  double *y[1];
  double *S[2];
  static double S0[2] = {1, 0};
  static double S1[2] = {0, 1};
  std::vector<double> *yvec_ptr;
  static double y0[15];
  S[0] = S0;
  S[1] = S1;
  y[0] = y0;

  if(taped_real[comp] == 0){
    phi_real_tape_func[comp]();
  }
  tensor_eval(taped_real[comp], 1, 2, 4, 2, x, y, S);
  if(memo_table_phi_real.size() > MAX_MEMO_PHI) memo_table_phi_real.clear();
  yvec_ptr = &memo_table_phi_real[std::make_tuple(comp, delta, tau)];
  yvec_ptr->assign(y[0], y[0] + 15);
  return yvec_ptr;
}

std::vector<double> *phi_ideal(comp_enum comp, double delta, double tau){
  /*
    Calculate the dimensionless ideal part of Helmholtz free energy (phi^0), and
    derivatives up to fourth order. The second order derivatives are need to
    thermodynamic relations for the properties.  Then the 3rd and 4th order
    derivatives are needed to calculate first and second order derivatives of
    the properties.

    The results of this function are cached, so they don't need to be repeated
    when calculating multiple properties.  After the first calculation the AD
    tape is also reused for subsequent calculations.
  */
  try{
    return &memo_table_phi_ideal.at(std::make_tuple(comp, delta, tau));
  }
  catch(std::out_of_range){
  }
  double x[2] = {delta, tau};
  double *y[1];
  double *S[2];
  static double S0[2] = {1, 0};
  static double S1[2] = {0, 1};
  std::vector<double> *yvec_ptr;
  static double y0[15];
  S[0] = S0;
  S[1] = S1;
  y[0] = y0;

  if(taped_ideal[comp] == 0){
    phi_ideal_tape_func[comp]();
  }
  tensor_eval(taped_ideal[comp], 1, 2, 4, 2, x, y, S);
  if(memo_table_phi_ideal.size() > MAX_MEMO_PHI) memo_table_phi_ideal.clear();
  yvec_ptr = &memo_table_phi_ideal[std::make_tuple(comp, delta, tau)];
  yvec_ptr->assign(y[0], y[0] + 15);
  return yvec_ptr;
}

void phi_real_for_sat(comp_enum comp, double delta, double tau, std::vector<double> *yvec_ptr){
  /*
    Calculate the dimensionless real part of Helmholtz free energy (phi^r), and
    derivatives second order with respect to delta. After the first calculation
    the AD tape is reused for subsequent calculations. This is not cached, since
    it is used specifically in solving for the saturation curve calculations.
  */
  double x[2] = {delta, tau};
  double *y[1];
  double *S[2];
  static double S0[2] = {1, 0}; // partials with respect to delta
  static double S1[2] = {0, 0}; // don't need the tau derivaties for this
  static double y0[3];
  S[0] = S0;
  S[1] = S1;
  y[0] = y0;
  if(taped_real[comp] == 0){
    phi_real_tape_func[comp]();
  }
  tensor_eval(taped_real[comp], 1, 2, 2, 1, x, y, S);
  yvec_ptr->assign(y[0], y[0] + 3);
}
