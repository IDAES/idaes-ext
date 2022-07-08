
#ifndef _INCLUDE_DELTA_H_
#define _INCLUDE_DELTA_H_

struct pressure_wrap_state {
  comp_enum comp;
  double p;
  double tau;
};

double delta_liquid(comp_enum comp, double pr, double tau);
double delta_vapor(comp_enum comp, double pr, double tau);

#endif
