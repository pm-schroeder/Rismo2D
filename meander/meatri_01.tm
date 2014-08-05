$RISMO2D 40300

# ==================================================================================================
# RISMO2D TIMESTEP FILE CREATED BY TIEGRIS ON Tue Aug  5 17:35:54 2014
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# TIME STEPPING (first,last,repetition)

#    first      : number of first time step (must be >= 1)
#    last       : number of last time step
#    repetition : number of repetition time step

$TM_STEPS        1        5        4

# --------------------------------------------------------------------------------------------------
# TIME INTERVAL (deltaTime,relaxFlow,relaxTurb)

#    deltaTime  : length of time step in seconds

#    relaxFlow  : minimum length of relaxed time step in flow equations
#    relaxTurb  : minimum length of relaxed time step in turbulence equations

$TM_INTERVAL   3600.0000    0.100000    0.100000

# --------------------------------------------------------------------------------------------------
# TIME WEIGHTING (thetaFlow, thetaTurb, thetaSedi)

#    thetaFlow  : time weighting in flow equations
#    thetaTurb  : time weighting in turbulence equations
#    thetaSedi  : time weighting in sediment equations

$TM_WEIGHT   0.50   0.50   0.50

# --------------------------------------------------------------------------------------------------
# OUTPUT FOR TIME STEPS

# write result files at the following time steps

#    1) enter a line with comma-separated numbers of time steps

#    2) an expression of the following form
#          <first number> to <last number> step <n>
#       e.g.: $TM_OUTPUT   5  to  20  step  5
#             is equivalent to
#             $TM_OUTPUT   5,10,15,20

$TM_OUTPUT  1,2,3,5

$TM_OUTPUT       50  to     1000  step       50

# --------------------------------------------------------------------------------------------------
# OPERATION SET FOR TIME STEP

$TM_STEP_NO  1

# --------------------------------------------------------------------------------------------------
#  operation cycles                                                                 $TM_CYCLE

#       no : variables
#     ------------------ 2D shallow water equations ------------------------------------------------
#        1 :  UV,S       velocities UV quadratic / water elevation S linear interpolation
#        2 :  UV,S       velocities UV quadratic / water elevation S linear interpolation / different time integration
#        4 :  UV,S       velocities UV linear / water elevation S constant

#     ------------------ dispersion models ---------------------------------------------------------
#       20 :  Dxx        secondary flow due to streamline curvature

#     ------------------ initialisations -----------------------------------------------------------
#       30 :  S          initialisation of water elevation, interpolation with section file
#       35 :             compute a elementwise constant water elevation

#       40 :             divergence free flow, quadratic velocities UV / linear water elevation S

#     ------------------ sediment transport models -------------------------------------------------
#       50 :  Qb         bed load sediment transport rate without bottom evolution
#       51 :  C          suspended load sediment transport, advection-diffusion equation
#       52 :  Qb         bed load sediment transport
#       53 :  C,Qb       bed load  + suspended load sediment transport
#       55 :  dz         bottem evolution computed from changes between Qb and initial Qb

#     ------------------ turbulence models ---------------------------------------------------------

#                        one-equation model
#       71 :  K,D        k differential equation, epsion: algebraic equation

#                        two equation model
#       80 :  K,D        initialisation of parameters K,D,vt from bottom shear stress
#       81 :  K,D        k-epsilon model - quadratic interpolation
#       84 :  K,D        k-epsilon model - linear interpolation
#       87 :  K,D        k-epsilon model - linear interpolation on quartered elements

#     ------------------ further operation cycles --------------------------------------------------
#       90 :             dry-rewet algorithm
#       95 :             reordering of elements
#       98 :             output of scaled model
#       99 :             output of results

# --------------------------------------------------------------------------------------------------
#  time gradients:  0 = instationary, 1 = stationary computation                    $TM_STATIONARY

# --------------------------------------------------------------------------------------------------
#  turbulence modeling                                                              $TM_TURBULENCE

#      ##1 :  constant eddy viscosity
#      ##2 :  algebraic bottom friction model for eddy viscosity
#             vt = cn * H * Utau
#      ##3 :  algebraic shear stress model for eddy viscosity
#             experimental version
#      ##4 :  Prandtls mixing length model
#             vt = lm * lm * sqrt(flow gradients)
#      ##5 :  Large Eddy Simulation by Smagorinsky
#             vt = ls * ls * sqrt(flow gradients)
#      ##6 :  Prandtl-Kolmogorov relation
#             eddy viscosity determined from previously computed K,D values
#             vt = cm * cd * K * K / D
#      ##7 :  LES by Smagorinsky + algebraic bottom friction
#             eddy viscosity computed as the sum of models ##2 and ##5
#   ----------------------------------------------------------------------------------------------
#      #1# :  anisotropic eddy viscosity (ELDERs model)
#   ----------------------------------------------------------------------------------------------
#      1## :  apply lower limit Vtmin to the eddy viscosity
#      2## :  iterative adaption of eddy viscosity in cycle
#      3## :  1## + 2##

# --------------------------------------------------------------------------------------------------
#  dispersion modeling:  0 = deactivated, 1 = activated                             $TM_DISPERSION

# --------------------------------------------------------------------------------------------------
#  maximum number of iterations                                                     $TM_MAXITER

# --------------------------------------------------------------------------------------------------
#  equation solver (number of solver in Rismo startup file (*.ris)                  $TM_SOLVER

# OPERATION CYCLES

$TM_CYCLE       40
$TM_STATIONARY   1
$TM_TURBULENCE   1
$TM_DISPERSION   0
$TM_MAXITER      1
$TM_SOLVER       1

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS

#  You may set multiple boundary conditions for lines
#  and nodes, if these are compatible.

#  type of conditions: D = Dirichlet, N=Neumann, C=Cauchy

#  kind of conditions:  $TM_LINE <id> <line number> <parameter>     // lines
#                       $TM_NODE <id> <node number> <parameter>     // nodes, x,y ignored

#   id           type  kind   description                                          parameter
#  ------------  ----  ----  ---------------------------------------------------  -----------
#                             inlet boundary conditions (Dirichlet)
#   INLET          D   L+N    specific discharge q=(U*h,V*h) [m3/s/m]              U,V,x,y
#   QINLET         D    L     discharge Q [m3/s]                                   Q,S
#   QTINLET        D    L     discharge hydrograph Q(t) [m3/s]                     Q,S,t

#                             natural outlet boundary conditions
#   OUTLET         N   L+N    water elevation So [mNN]                             S,x,y
#   OPEN                      open boundary (no conditions at all)
#   SQOUTLET       N    L     S(Q) relation [mNN]                                  S,Q
#   STOUTLET       N   L+N    S(t) relation [mNN       ]                           S,t

#                             further boundary conditions (flow model)
#   SLIP           D   L+N    flow direction tg()=U/V [-]                          U,V,x,y
#   LOGLAW         N   L+N    logarithmic law of the wall; distance and roughness  dw,ks
#   SET_UV         D   L+N    set flow velocities (U,V) [m/s]                      U,V,x,y
#   SET_S          D   L+N    water elevation S [mNN] (Dirichlet)                  S,x,y
#   SOURCE              N     sink / sourcs Qs [m3/s]                              Q

#                             boundary conditions for k-epsilon-Modell
#   SET_KD         D   L+N    set k und epsilon K [m2/s2], D [m2/s2]               K,D,x,y

#                             boundary conditions for sediment transport models
#   SET_C          D   L+N    sediment conzentration Cs [kg/m3]                    Cs,x,y
#   RATE_C         C   L+N    sediment rate Qb [kg/s]                              Qb,x,y

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT LINES

$TM_LINE  INLET           1     0.00450     0.00000     0.00000     0.00000
$TM_LINE  INLET           1     0.00450     0.00000     0.00000     1.44440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.48440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.84440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     1.88440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     2.98000
$TM_LINE  OUTLET          2     0.22330     0.00000     0.00000
$TM_LINE  LOGLAW          7     0.01000     0.00100

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT NODES


# --------------------------------------------------------------------------------------------------
# OPERATION SET FOR TIME STEP

$TM_STEP_NO  2


# OPERATION CYCLES

$TM_CYCLE       20   1
$TM_STATIONARY   1   1
$TM_TURBULENCE   1   1
$TM_DISPERSION   0   1
$TM_MAXITER      1  20
$TM_SOLVER       1   5

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS


# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT LINES

$TM_LINE  INLET           1     0.00450     0.00000     0.00000     0.00000
$TM_LINE  INLET           1     0.00450     0.00000     0.00000     1.44440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.48440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.84440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     1.88440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     2.98000
$TM_LINE  OUTLET          2     0.22330     0.00000     0.00000
$TM_LINE  LOGLAW          7     0.01000     0.00100

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT NODES


# --------------------------------------------------------------------------------------------------
# OPERATION SET FOR TIME STEP

$TM_STEP_NO  3


# OPERATION CYCLES

$TM_CYCLE       20   1  80
$TM_STATIONARY   1   1   1
$TM_TURBULENCE 112 112 112
$TM_DISPERSION   0   1   0
$TM_MAXITER     20  20  20
$TM_SOLVER       1   5   1

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS


# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT LINES

$TM_LINE  INLET           1     0.00450     0.00000     0.00000     0.00000
$TM_LINE  INLET           1     0.00450     0.00000     0.00000     1.44440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.48440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.84440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     1.88440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     2.98000
$TM_LINE  OUTLET          2     0.22330     0.00000     0.00000
$TM_LINE  LOGLAW          7     0.01000     0.00100

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT NODES


# --------------------------------------------------------------------------------------------------
# OPERATION SET FOR TIME STEP

$TM_STEP_NO  4


# OPERATION CYCLES

$TM_CYCLE       87  20   1
$TM_STATIONARY   1   1   1
$TM_TURBULENCE 112 112 112
$TM_DISPERSION   0   0   1
$TM_MAXITER     50  30  30
$TM_SOLVER       5   1   5

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS


# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT LINES

$TM_LINE  INLET           1     0.00450     0.00000     0.00000     0.00000
$TM_LINE  INLET           1     0.00450     0.00000     0.00000     1.44440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.48440
$TM_LINE  INLET           1     0.02190     0.00000     0.00000     1.84440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     1.88440
$TM_LINE  INLET           1     0.00410     0.00000     0.00000     2.98000
$TM_LINE  OUTLET          2     0.22330     0.00000     0.00000
$TM_LINE  LOGLAW          7     0.01000     0.00100

# --------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS AT NODES

