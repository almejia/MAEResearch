//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CNeumannBC.h"

registerMooseObject("MooseApp", CNeumannBC);

template <>
InputParameters
validParams<CNeumannBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addParam<Real>("var_1", 1.0, "The value of the gradient on the boundary.");
  params.addParam<Real>("var_2", 1.0, "The value of the gradient on the boundary.");

  return params;
}

CNeumannBC::CNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _value_1(getParam<Real>("var_1")), _value_2(getParam<Real>("var_2"))
{
}

Real
CNeumannBC::computeQpResidual()
{
  return -_test[_i][_qp] * _value_1 * _value_2;
}
