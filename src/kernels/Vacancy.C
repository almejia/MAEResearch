//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Vacancy.h"
#include "func_bssanova_2.hpp"


registerMooseObject("MooseApp", Vacancy);

template <>
InputParameters
validParams<Vacancy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");

  params.addRequiredCoupledVar("dop", "Coupled variable [dopant concentration].");
  params.addRequiredCoupledVar("phi", "Coupled variable Phi (electrostatic potential).");

  params.addParam<PostprocessorName>("pps_name", 0.1, "Postprocessor name");


  params.addParam<Real>("R", 8.314, "Universal constant (R = 8.314)");
  params.addParam<Real>("T", 1673, "Temperature (T).");

  params.addParam<Real>("Z", 2, "Dopant chemical formula index");
  params.addParam<Real>(
      "nv", 8.2975e4, "Cation site density (ny). Default is 41488 = 4/(avnum*latpar^3).");
  params.addParam<Real>("nvv",
                        2.4893e5,
                        "Anion-anion next nearest-neighbor bond density (nvv). Default is 165952 = "
                        "4/(avnum*latpar^3).");
  params.addParam<Real>(
      "nyv",
      3.3190e5,
      "Anion-cation next-nearest-neighbor bond density (nyv). Default is 331900 = 8*ny.");
  params.addParam<Real>("F", 96485.0, "Energy constant (F). Default is 96485.");
  params.addParam<Real>("fv", 9e4, "Dopant self interaction energy (fy). Default is 9e4.");
  params.addParam<Real>("fyv", -7e3, "Dopant-vacancy interaction energy (fyv). Default is -7e3.");
  params.addParam<Real>("cv", 1e-9, "Dopant gradient energy coefficient (cd). Default is 8e-9.");

  params.addParam<Real>("lambv", 1.3934712e9, "Lagrange multiplier"); // 1.3903e9

  params.addParam<Real>("Disc", 0, "Discrepancy flag");
  params.addParam<Real>("By1", 0, "By1");
  params.addParam<Real>("By2", 0, "By2");
  params.addParam<Real>("Bvy1", 0, "Bvy1");
  params.addParam<Real>("Bvy2", 0, "Bvy2");
  params.addParam<Real>("Bv1", 0, "Bv1");
  params.addParam<Real>("Bv2", 0, "Bv2");
  params.addParam<Real>("Byp", 0, "Byp");
  params.addParam<Real>("Bvp", 0, "Bvp");

  return params;
}

Vacancy::Vacancy(const InputParameters & parameters)
  : Kernel(parameters),

    _dop_var(coupledValue("dop")),
    _phi_var(coupledValue("phi")),

    _pps_value(getPostprocessorValue("pps_name")),


    _d_var(coupled("dop")),
    _p_var(coupled("phi")),

    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),

    _Z(getParam<Real>("Z")),
    _nv(getParam<Real>("nv")),
    _nvv(getParam<Real>("nvv")),
    _nyv(getParam<Real>("nyv")),
    _F(getParam<Real>("F")),
    _fv(getParam<Real>("fv")),
    _fyv(getParam<Real>("fyv")),
    _cv(getParam<Real>("cv")),

    _lambv(getParam<Real>("lambv")),

    _Disc_flag(getParam<Real>("Disc")),
    _By1(getParam<Real>("By1")),
    _By2(getParam<Real>("By2")),
    _Bvy1(getParam<Real>("Bvy1")),
    _Bvy2(getParam<Real>("Bvy2")),
    _Bv1(getParam<Real>("Bv1")),
    _Bv2(getParam<Real>("Bv2")),
    _Byp(getParam<Real>("Byp")),
    _Bvp(getParam<Real>("Bvp"))

{
}

Real
Vacancy::computeQpResidual()
{

  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Bvp;

  Betas.resize(5);
  Betas(0) = _Bv1;  /*1;  1.9002e10;*/  // Bv1
  Betas(1) = _Bv2;  /*1;  0.6241e10;*/  // Bv2
  Betas(2) = _Bvy1;  /*1;  -2.0418e9;*/  // Bvy
  Betas(3) = _Bvp;  /*1;  6.3635e-9;*/  // Bvp
  Betas(4) = 0;         // Bvvp

  bool diskconfsuccess = false;
  diskconfsuccess = DiskFunc1.SplineConfig(DiskConfile);
  if (!diskconfsuccess) {
    std::cout << "There is an issue with the Discrepancy Configuration File. Press Enter to exit." << std::endl;
    return 0;
  }

  Betas1.resize(DiskFunc1.betsize());

  for (int i = 0; i < Betas1.size(); i++) {
    Betas1[i] = Betas[i];
  } 

  Bvp = _Bvp;  /*0;  6.3635e-9;*/

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _dop_var[_qp];

  DiskFunc1.SplineDeriv2(0, Par1, Betas1, DiskRes);

  //std::cout << "\nParam (v u): " << Par1[0] << " " << Par1[1] << "\t Betas: " << Betas1(0) << " " << Betas(1) << "\t Disc: " << DiskRes << "\t Res:" << ((_nv * 2 * _F * _phi_var[_qp]) + _R * _T * _nv * log(_u[_qp] / (1 - _u[_qp])) + _lambv + _pps_value + DiskRes) * _test[_i][_qp] << std::endl;
  //*

  if (_Disc_flag == 1) {
    return ((_nv * 2 * _F * _phi_var[_qp]) + _R * _T * _nv * log(_u[_qp] / (1 - _u[_qp])) + _lambv + _pps_value + DiskRes) * _test[_i][_qp] /*+ _test[_i][_qp] * Byp * _grad_u[_qp]*/ ; //*/
  }

  else {
  //*
  return _cv * _grad_u[_qp] * _grad_test[_i][_qp] +
         ((_nv * 2 * _F * _phi_var[_qp]) + (_nyv * _fyv * _dop_var[_qp]) +
          2 * _nvv * _fv * _u[_qp] + _R * _T * _nv * log(_u[_qp] / (1 - _u[_qp])) + _lambv) *
	  _test[_i][_qp];/**/
  }
}

//*
Real
Vacancy::computeQpJacobian()
{

  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc2.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Bvp;

  Betas.resize(5);
  Betas(0) = _Bv1;  /*1;  1.9002e10; */ // Bv1
  Betas(1) = _Bv2;  /*1;  0.6241e10; */ // Bv2
  Betas(2) = _Bvy1; /*1;  -2.0418e9; */ // Bvy
  Betas(3) = _Bvp;  /*1;  6.3635e-9; */ // Bvp
  Betas(4) = 0;         // Bvvp

  bool diskconfsuccess = false;
  diskconfsuccess = DiskFunc1.SplineConfig(DiskConfile);
  if (!diskconfsuccess) {
    std::cout << "There is an issue with the Discrepancy Configuration File. Press Enter to exit." << std::endl;
    return 0;
  }

  Betas1.resize(DiskFunc1.betsize());

  for (int i = 0; i < Betas1.size(); i++) {
    Betas1[i] = Betas[i];
  } 

  Bvp = _Bvp;  /*0;  6.3635e-9; */

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _dop_var[_qp];

  DiskFunc1.Spline2ndDeriv(0, Par1, Betas1, DiskRes);

  //*
  if (_Disc_flag == 1) {
    return /*Byp * _grad_phi[_j][_qp] * _test[_i][_qp] +*/ ( _R * _T * _nv * (1 / (_u[_qp] * (1 - _u[_qp]) ) ) + DiskRes) * _test[_i][_qp] * _phi[_j][_qp]; /**/
  }

  else {
  //*
  return +_cv * _grad_phi[_j][_qp] * _grad_test[_i][_qp] +
         2 * _nvv * _fv * _phi[_j][_qp] * _test[_i][_qp] +
         _R * _T * _nv * (1 / (_u[_qp] * (1 - _u[_qp]))) * _test[_i][_qp] * _phi[_j][_qp]; /**/
  }
}
/**/

Real
Vacancy::computeQpOffDiagJacobian(unsigned int jvar)
{

  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc3.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Bvp;

  Betas.resize(1);
  Betas(0) = _Bvy1;  /*1; -2.0418e9; */ // Bvy

  bool diskconfsuccess = false;
  diskconfsuccess = DiskFunc1.SplineConfig(DiskConfile);
  if (!diskconfsuccess) {
    std::cout << "There is an issue with the Discrepancy Configuration File. Press Enter to exit." << std::endl;
    return 0;
  }

  Betas1.resize(DiskFunc1.betsize());

  for (int i = 0; i < Betas1.size(); i++) {
    Betas1[i] = Betas[i];
  } 

  Bvp = _Bvp;  /*0;  6.3635e-9; */

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _dop_var[_qp];


  DiskFunc1.SplineDeriv3(0, Par1, Betas1, DiskRes);


  if (jvar == _d_var)
  {
    if (_Disc_flag == 1) {
      return DiskRes *_test[_i][_qp] * _phi[_j][_qp];
    }
    else {
      return _nyv * _fyv * _test[_i][_qp] * _phi[_j][_qp];
    }
  }

  if (jvar == _p_var)
  {
    return _nv * 2 * _F * _test[_i][_qp] * _phi[_j][_qp];
  }

  return 0.0;
}
