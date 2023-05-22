//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SplitPCRes2.h"

registerMooseObject("MooseApp", SplitPCRes2);

template <>
InputParameters
validParams<SplitPCRes2>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");

  params.addParam<Real>("M", 1, "Mobility constant (M = 1)");

  return params;
}

SplitPCRes2::SplitPCRes2(const InputParameters & parameters)
  : Kernel(parameters),

    _M(getParam<Real>("M"))

{
}

Real
SplitPCRes2::computeQpResidual()
{
  //*
    return _M * _grad_u[_qp] * _grad_test[_i][_qp]; //*/ /*_l_d[_qp]/**/
}

//*
Real
SplitPCRes2::computeQpJacobian()
{
  //*
  return _M * _grad_phi[_j][_qp] * _grad_test[_i][_qp]; /**/
}
/**/

Real
SplitPCRes2::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
