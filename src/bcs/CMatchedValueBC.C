//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CMatchedValueBC.h"

registerMooseObject("MooseApp", CMatchedValueBC);

template <>
InputParameters
validParams<CMatchedValueBC>()
{
  InputParameters params = validParams<NodalBC>();
  params.addRequiredCoupledVar("v", "The variable whose value we are to match.");
  params.addClassDescription("Implements a NodalBC which equates the values of 2 different coupled Variables "
                             "on a specified boundary.");

  params.addParam<Real>("R", 8.314, "Gas constant (R = 8.314)");
  params.addParam<Real>("T", 773, "Temperature (T).");

  params.addParam<Real>("Vac", 0.04, "8% YSZ oxygen vacancy concentration in electrolyte [0.04]");
  params.addParam<Real>("nv_ysz", 1, "Vacancy site density (nv) in YSZ electrolyte: 9.7656e4 = 8/(avnum*latpar^3).");

  params.addParam<Real>("F", 96485.0, "Energy constant (F). Default is 96485.");
  params.addParam<Real>("Mu_YSZ", 8.4, "Standard chemical potential for oxygen vacancies in YSZ. Default is 8.4 [eV].");

  return params;
}

CMatchedValueBC::CMatchedValueBC(const InputParameters & parameters)
  : NodalBC(parameters), _v(coupledValue("v")), _v_num(coupled("v")), _R(getParam<Real>("R")), _T(getParam<Real>("T")), _vac(getParam<Real>("Vac")), _nv_ysz(getParam<Real>("nv_ysz")), _F(getParam<Real>("F")), _mu_ysz(getParam<Real>("Mu_YSZ"))

{
}

Real
CMatchedValueBC::computeQpResidual()
{
  return _u[_qp] - _mu_ysz * _F - 2 * _F * _v[_qp] - _R * _T * _nv_ysz * log(_vac / (1 - _vac));
}

Real
CMatchedValueBC::computeQpJacobian()
{
  return 1;
}

Real
CMatchedValueBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_num)
    return -2 * _F;
  else
    return 0.;
}
