#ifndef _INCLUDE_SAT_H_
#define _INCLUDE_SAT_H_

std::vector<double> *sat_p(comp_enum comp, double tau);
std::vector<double> *sat_delta_v(comp_enum comp, double tau);
std::vector<double> *sat_delta_l(comp_enum comp, double tau);

#endif
