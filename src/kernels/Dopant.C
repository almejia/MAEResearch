//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Dopant.h"
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include "func_bssanova_2.hpp"


registerMooseObject("MooseApp", Dopant);

template <>
InputParameters
validParams<Dopant>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");

  params.addRequiredCoupledVar("vac", "Coupled variable [vacancy concentration].");
  params.addRequiredCoupledVar("phi", "Coupled variable Phi (electrostatic potential).");

  params.addParam<PostprocessorName>("pps_name", 0.1, "Postprocessor name");

  //  params.addRequiredCoupledVar("ld", "Coupled Lagrange Multiplier Lambda_d.");

  params.addParam<Real>("R", 8.314, "Universal constant (R = 8.314)");
  params.addParam<Real>("T", 1673, "Temperature (T).");

  params.addParam<Real>("Z", 2, "Dopant chemical formula index");
  params.addParam<Real>(
      "ny", 4.1488e4, "Cation site density (ny). Default is 41488 = 4/(avnum*latpar^3).");
  params.addParam<Real>(
      "nyy",
      1.6595e5,
      "Cation-cation next-nearest-neighbor bond density (nyy). Default is 165952 = 4*ny.");
  params.addParam<Real>(
      "nyv",
      3.3190e5,
      "Anion-cation next-nearest-neighbor bond density (nyv). Default is 331900 = 8*ny.");
  params.addParam<Real>("F", 96485.0, "Energy constant (F). Default is 96485.");
  params.addParam<Real>("fy", 90000.0, "Dopant self interaction energy (fy). Default is 9e4.");
  params.addParam<Real>(
      "fyv", -7000.0, "Dopant-vacancy interaction energy (fyv). Default is -7e3.");
  params.addParam<Real>("cd", 8e-9, "Dopant gradient energy coefficient (cd). Default is 8e-9.");

  params.addParam<Real>("lambd", -1.599515214e9, "Lagrange multiplier"); //-1.603e9

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

Dopant::Dopant(const InputParameters & parameters)
  : Kernel(parameters),

    _vac_var(coupledValue("vac")),
    _phi_var(coupledValue("phi")),
    
    _pps_value(getPostprocessorValue("pps_name")),

    //    _l_d(coupledValue("ld")),

    _v_var(coupled("vac")),
    _p_var(coupled("phi")),

    //    _ld(coupled("ld")),

    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),

    _Z(getParam<Real>("Z")),
    _ny(getParam<Real>("ny")),
    _nyy(getParam<Real>("nyy")),
    _nyv(getParam<Real>("nyv")),
    _F(getParam<Real>("F")),
    _fy(getParam<Real>("fy")),
    _fyv(getParam<Real>("fyv")),
    _cd(getParam<Real>("cd")),

    _lambd(getParam<Real>("lambd")),

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
Dopant::computeQpResidual()
{
  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Byp;

  Betas.resize(5);
  Betas(0) = _By1;  /*1; 0.8422e10;*/ // By1
  Betas(1) = _By2;  /*1; 2.9252e10;*/ // By2
  Betas(2) = _Bvy1;  /*1; -2.0418e9;*/ // Bvy
  Betas(3) = _Byp;  /*1; 1.1804e-9;*/ // Bvp
  Betas(4) = 0;         // Byyp

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

  Byp = _Byp;  /*0; 1.1804e-9*/

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _vac_var[_qp];

  DiskFunc1.SplineDeriv2(0, Par1, Betas1, DiskRes);

  //  std::cout << "\nParam (u v): " << Par1[0] << " " << Par1[1] << "\t Betas: " << Betas1(0) << " " << Betas(1) << "\t Disc: " << DiskRes << "\t Res:" << ((-_ny * _Z * _F * _phi_var[_qp]) + _R * _T * _ny * log(_u[_qp] / (1 - _u[_qp])) + _lambd + _pps_value + DiskRes) * _test[_i][_qp] << std::endl;
  //*
  
  if (_Disc_flag == 1) {
    return ((-_ny * _Z * _F * _phi_var[_qp]) + _R * _T * _ny * log(_u[_qp] / (1 - _u[_qp])) + _lambd + _pps_value + DiskRes) * _test[_i][_qp] /*+ _test[_i][_qp] * Byp * _grad_u[_qp]*/ ; //*/
  }

  else {
  //*
    return _cd * _grad_u[_qp] * _grad_test[_i][_qp] +
         ((-_ny * _Z * _F * _phi_var[_qp]) + (_nyv * _fyv * _vac_var[_qp]) +
          2 * _nyy * _fy * _u[_qp] + _R * _T * _ny * log(_u[_qp] / (1 - _u[_qp])) + _lambd /*+ _pps_value*/) *
             _test[_i][_qp]; //*/ /*_l_d[_qp]/**/
  }
}

//*
Real
Dopant::computeQpJacobian()
{

  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc2.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Byp;

  Betas.resize(5);
  Betas(0) = _By1;  /*1;  0.8422e10; */// By1
  Betas(1) = _By2;  /*1;  2.9252e10; */// By2
  Betas(2) = _Bvy1;  /*1;  -2.0418e9; */// Bvy
  Betas(3) = _Bvp;  /*1;  1.1804e-9; */// Bvp
  Betas(4) = 0;         // Byyp

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

  Byp = _Byp;  /*0;  1.1804e-9;*/

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _vac_var[_qp];

  DiskFunc1.Spline2ndDeriv(0, Par1, Betas1, DiskRes);

  //*
  if (_Disc_flag == 1) {
    return /*Byp * _grad_phi[_j][_qp] * _test[_i][_qp] +*/ ( _R * _T * _ny * (1 / (_u[_qp] * (1 - _u[_qp]) ) ) + DiskRes) * _test[_i][_qp] * _phi[_j][_qp]; /**/
  }

  else {
  //*
  return +_cd * _grad_phi[_j][_qp] * _grad_test[_i][_qp] +
         2 * _nyy * _fy * _test[_i][_qp] * _phi[_j][_qp] +
         _R * _T * _ny * (1 / (_u[_qp] * (1 - _u[_qp]))) * _test[_i][_qp] * _phi[_j][_qp]; /**/
  }
}
/**/

Real
Dopant::computeQpOffDiagJacobian(unsigned int jvar)
{

  func_bssanova DiskFunc1;
  std::string DiskConfile = "config/config_Disc3.txt";
  boost::numeric::ublas::vector<double> Betas1, Betas;
  boost::numeric::ublas::vector<double> Par1; 
  double DiskRes;
  double DiskRes2;
  double Byp;

  Betas.resize(1);
  Betas(0) = _Bvy1;  /*1;  -2.0418e9; */ // Bvy

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

  Byp = _Byp;  /*0;  1.1804e-9;*/

  Par1.resize(2);
  Par1[0] = _u[_qp];
  Par1[1] = _vac_var[_qp];

  DiskFunc1.SplineDeriv3(0, Par1, Betas1, DiskRes);

  if (jvar == _v_var)
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
    return -_ny * _Z * _F * _test[_i][_qp] * _phi[_j][_qp];
  }

  //  if(jvar == _ld)
  //  {
  //    return _test[_i][_qp] * _phi[_j][_qp];
  //  }

  return 0.0;
}
