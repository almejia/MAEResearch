## Input file: ElePhaNT MOOSE app  

[Mesh]
  type = GeneratedMesh
  dim = 2         # Mesh dimensions: 1D, 2D, 3D
  nx = 20         # 1D: elements X-direction 
  ny = 30         # 2D: elements Y-direction
  nz = 1          # 3D: elements Z-direction
  xmax = 6e-3     # Cell base length [m] = 6.5 [mm]
  ymax = 0.2e-3   # 2D cell height [m] = 2 [mm] (0.2 mm)
  zmax = 2.5e-3   # 3D cell width 2.5 [mm]
  bias_x = 1.0    # bias 0.5 - 2.0 (left - right element-size refinement)
  bias_y = 0.65   # y-bias 0.5 - 2.0 (top - bottom refinement)
[]

[GlobalParameters]

[]

[MeshModifiers]
  [./subdomain_id]
    type = AssignElementSubdomainID
    subdomain_ids = '1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

3 3 3 3 2 2 2 2 2 2 0 0 0 0 0 0 4 4 4 4
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 4 4 4 4
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0

2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0
2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0'
  [../]

  [./side_set1]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = '0 1 2'
    paired_block = '2'
    new_boundary = 'new'
  [../]
  [./side_set2]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = '1 3'
    paired_block = '3'
    new_boundary = 'new2'
  [../]
  [./side_set3]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = '0 1 4'
    paired_block = '4'
    new_boundary = 'new3'
  [../]

  [./delete_subdomain1]
    type = BlockDeleter
    depends_on = side_set1
    block_id = 2
  [../]
  [./delete_subdomain2]
    type = BlockDeleter
    depends_on = side_set2
    block_id = 3
  [../]
  [./delete_subdomain3]
    type = BlockDeleter
    depends_on = side_set3
    block_id = 4
  [../]

  [./break_boundary]
    type = BreakBoundaryOnSubdomain
    depends_on = subdomain_id
  [../]
  [./CeO2_YSZ]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = 0
    paired_block = 1
    new_boundary = 'ceo2_ysz'
  [../]
  [./Au_CeO2]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = 0
    paired_block = 4
    new_boundary = 'ceo2_au'
  [../]
  [./Au_YSZ]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = 1
    paired_block = 4
    new_boundary = 'ysz_au'
  [../]
  [./Pt_YSZ]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = 1
    paired_block = 3
    new_boundary = 'ysz_pt'
  [../]

[]
  
[Variables]
#active = 'q v phi muq muv'
  [./q]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
    scaling = 1
    block = 0 #'0 5'
  [../]
  [./v]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
    scaling = 1
    block = 0 #'0 5'
  [../]
  [./phi]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
    scaling = 1
    block = 0 #'0 5'
  [../]
  
  [./muq]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
    scaling = 1
    block = 0 #'0 5'
  [../]
  [./muv]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
    scaling = 1
    block = 0 #'0 5'
  [../]

  [./phi_ysz]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
    scaling = 1
    block = 1
  [../]

[]

[Kernels]
#active = 'phi_kernel electron_kernel vacancy_kernel mu_electron mu_vacancy'
  [./phi_kernel]
    type = EPotential
    variable = phi
    ele = q
    vac = v
    block = 0 #'0 5'
  [../]
  [./electron_kernel]
    type = Electron
    variable = q
    vac = v
    phi = phi
    mu = muq
    lambd = 0 # -1.599515251458632e9 # fixed Lagrange multiplier value
    block = 0 #'0 5'
  [../]
  [./vacancy_kernel]
    type = NewVacancy
    variable = v
    ele = q
    phi = phi
    mu = muv
    lambv = 0 # 1.393471200379236e9 # fixed Lagrange multiplier value
    block = 0 #'0 5'
  [../]

  [./mu_electron]
    type = SplitPCRes
    variable = muq
    con = q
    block = 0 #'0 5'
  [../]
  [./mu_vacancy]
    type = SplitPCRes
    variable = muv
    con = v
    block = 0 #'0 5'
  [../]

  [./phi_ysz]
    type = Diffusion
    variable = phi_ysz
    block = 1
  [../]
		       
#  [./u_dot]
#    type = TimeDerivative
#    variable = u
#  [../]
#  [./v_dot]
#    type = TimeDerivative
#    variable = v
#  [../]
[]

#[InterfaceKernels]
#  [./interface]
#    type = PenaltyInterfaceDiffusion
#    variable = muv
#    neighbor_var = phi_ysz
#    boundary = ceo2_ysz
#    penalty = 1
#  [../]
#[]
		       
[BCs]
#active = 'Neumann_inlet_u Neumann_outlet_u Neumann_inlet_v Neumann_outlet_v Neumann_inlet_phi Neumann_outlet_phi'

  [./Neumann_inlet_u]
    type = NeumannBC
    variable = q
    boundary = top
    value = 2.5
  [../]
  [./Neumann_outlet_u]
    type = NeumannBC
    variable = q
    boundary = 'new ceo2_ysz' # top
    value = 0
  [../]

  [./Neumann_inlet_v]
    type = NeumannBC
    variable = v
    boundary = top
    value = 4.2746 # e9 using the larger and negative value for this BC leads to solver failure?
  [../]
  [./Neumann_outlet_v]
    type = NeumannBC
    variable = v
    boundary = 'new ceo2_ysz' # top
    value = 0
  [../]

  [./Neumann_inlet_phi]
    type = NeumannBC
    variable = phi
    boundary = top
    value = 0
  [../]
  [./Neumann_outlet_phi]
    type = NeumannBC # DirichletBC
    variable = phi
    boundary = 'new ceo2_ysz' # top
    value = 0
  [../]

  [./muq_ce_au]
    type = CDirichletBC
    variable = muq
    boundary = ceo2_au
    var1 = 0.01   # VOLTAGE applied
    var2 = 94685  # Energy coefficient
  [../]

  [./ysz_au]
    type = NeumannBC
    variable = phi_ysz
    boundary = ysz_au
    value = 0
  [../]
  [./phi_ysz_pt]
    type = DirichletBC
    variable = phi_ysz
    boundary = ysz_pt
    value = 0
  [../]

  [./interface]
    type = CMatchedValueBC
    variable = muv
    boundary = ceo2_ysz
    v = phi_ysz
  [../]
[]

[Postprocessors]
active = ''
  [./average_concentration]
    type = ElementAverageValue
    variable = q
  [../]
[]

[VectorPostprocessors]
  [./line_x]
    type = LineValueSampler
    start_point = '4e-3 0.2e-3 0'
    end_point = '4e-3 0.19995e-3 0'
    num_points = 256
    variable = 'q v phi muq muv'
    sort_by = 'id'
  [../]
  [./line_y]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '0 5e-9 0'
    num_points = 125
    variable = 'q v phi muq muv'
    sort_by = 'id'
  [../]
[]


[Executioner]
#  type = Transient
#  solve_type = NEWTON
#  l_max_its = 30
#  l_tol = 1e-6
#  nl_max_its = 50
#  nl_abs_tol = 1e-9
#  end_time = 1800
#  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
#  petsc_options_value = 'asm  31  preonly  ilu  1'

  
  type = Steady

  solve_type = NEWTON
  nl_max_its = 30
  l_max_its = 100
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_tol = 1.0e-8
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-8
		       
#  petsc_options = '-ksp_converged_reason -snes_converged_reason'

#  solve_type = PJFNK
#  nl_max_its = 30
#  l_max_its = 1000
#  petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre boomeramg'

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Adaptivity]
  marker = errorfrac
  steps = 2
  [./Indicators]
    [./error]
    type = GradientJumpIndicator
    variable = phi
    [../]
  [../]
  [./Markers]
    [./errorfrac]
    type = ErrorFractionMarker
    refine = 0.25
    coarsen = 0.1
    indicator = error
    [../]
  [../]
[]

[Outputs]
#  file_base = out_test
#  console = true
#  exodus = true
#  csv = true
#  [./outtest]
#    type = Exodus
#  [../]

  execute_on = 'nonlinear'
#  [./csv]
#    type = csv
#    execute_on = 'initial nonlinear final'
#  [../]
  [./exodus]
    type = Exodus
    execute_on = 'initial nonlinear final'
  [../]

[]
