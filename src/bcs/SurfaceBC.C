//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceBC.h"

registerMooseObject("MooseApp", SurfaceBC);

template <>
InputParameters
validParams<SurfaceBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Specify input parameters that we want users to be able to set:
  params.addParam<Real>("alpha", 1.0, "Value multiplied by on the boundary");
  params.addParam<Real>("R", 8.314, "Universal constant (R = 8.314)");
  params.addParam<Real>("T", 1673, "Temperature (T).");
  params.addParam<Real>("Pco", 1, "CO Partial pressure");
  params.addParam<Real>("eta2", 1, "Pre-exponential parameter inc");
  params.addParam<Real>("dHdd2", 1, "Enthalpy** for inc");
  params.addParam<Real>("dS2", 1, "Entropy for inc");
  params.addParam<Real>("dH2", 1, "Enthalpy for inc");
  params.addParam<Real>("u", 1, "Mobility coefficient");
  params.addParam<Real>("n", 1, "Surface site density");
  
  params.addRequiredCoupledVar("theta", "[THETA] value at boundary");
  params.addRequiredCoupledVar("ele", "[ELEC] value at boundary");
  params.addRequiredCoupledVar("vac", "[VAC] value at boundary");
  params.addRequiredCoupledVar("var", "[VAR] value at boundary");
  
  return params;
}

SurfaceBC::SurfaceBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    // store the user-specified parameters from the input file...
    _alpha(getParam<Real>("alpha")),        // for the multiplier
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    _Pco(getParam<Real>("Pco")),
    _eta2(getParam<Real>("eta2")),
    _dHdd2(getParam<Real>("dHdd2")),
    _dS2(getParam<Real>("dS2")),
    _dH2(getParam<Real>("dH2")),
    _u(getParam<Real>("u")),
    _n(getParam<Real>("n")),

    _theta_var(coupledValue("theta")), // for the coupled variable
    _ele_var(coupledValue("ele")),
    _vac_var(coupledValue("vac")),
    _some_var(coupledValue("var"))
{
}

Real
SurfaceBC::computeQpResidual()
{
  // For this Neumann BC grad(u)=alpha * v on the boundary.
  // We use the term produced from integrating the diffusion operator by parts.
  // return -_test[_i][_qp] * _alpha * _some_var_val[_qp];
  return -_eta2 * exp(-_dHdd2/(_R*_T)) * (_theta_var[_qp]*_vac_var[_qp]*_ele_var[_qp]*_ele_var[_qp] - _Pco*(1-_theta_var[_qp]-_vac_var[_qp])*(1-_theta_var[_qp]-_vac_var[_qp])*(1-_ele_var[_qp])*(1-_ele_var[_qp]) / (exp(_dS2/_R)*exp(-_dH2/(_R*_T))) ) / (_alpha*_u*_n*_some_var[_qp]) * _test[_i][_qp];
}
