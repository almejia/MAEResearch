//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DOPANT_H
#define DOPANT_H

#include "Kernel.h"

class Dopant;

template <>
InputParameters validParams<Dopant>();

class Dopant : public Kernel
{
public:
  Dopant(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _vac_var;
  const VariableValue & _phi_var;

  const PostprocessorValue & _pps_value;

  //  const VariableValue & _l_d;

  unsigned int _v_var;
  unsigned int _p_var;

  //  unsigned int _ld;

  Real _R; /** Universal constant [J/mol*K] */
  Real _T; /** Temperature [K] */

  Real _Z;   /** Dopant chemcial formula index */
  Real _ny;  /** Cation site density [mol/m^3] */
  Real _nyy; /** Cation-cation density [mol/m^3] */
  Real _nyv; /** Anion-cation density [mol/m^3] */
  Real _F;   /** 96485.0 */
  Real _fy;  /** Dopant self-interaction [J/mol] */
  Real _fyv; /** Dopant-vacancy interaction [J/mol] */
  Real _cd;  /** Dopant gradient energy coefficient [J/mol-m] */

  Real _lambd; /** Lagrange multiplier */

  Real _Disc_flag;
  Real _By1;
  Real _By2;
  Real _Bvy1;
  Real _Bvy2;
  Real _Bv1;
  Real _Bv2;
  Real _Byp;
  Real _Bvp;

};

#endif /* DOPANT_H */
