% Author: Danie Ludick (dludick@sun.ac.za)
% Project: PEC plate example (2018.07.31)
%
% Note: Each project directory / example directory needs to have a sunem_initialise.m
% script, that is used to setup the correct environment for that example.
%
% Refer to the /doc folder for more information

% --------------------------------------------------------------------------------------------------
% Initialise the environment
% --------------------------------------------------------------------------------------------------
% Project output directory: './dipoles/'
% Debug: True/False
Const = sunem_initialise('pec_plate',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;
Const.useCppEngine          = true;
Const.useMPI                = true;
Const.cppBuildDir           = '../../src/CPPMoMinterface/content/CPPMoMsolver/';
Const.num_procs             = 4;
Const.num_threads           = 2;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'pec_plate.mat';
Const.FEKOstrfilename          = 'pec_plate.str';
Const.FEKOrhsfilename          = 'pec_plate.rhs';
Const.FEKOoutfilename          = 'pec_plate.out';
Const.FEKOefefilename          = 'pec_plate.efe';
Const.FEKOffefilename          = 'pec_plate.ffe';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      = 'sunem_mom_pec_plate.str';

% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO] = extractFEKOMoMmatrixEq(Const);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectorsFEKO);

[CSOL] = runCppMoMSolver(Const, Solver_setup, yVectorsFEKO, xVectorsFEKO);



% 2018.06.10: If we are going to run the SUNEM MoM solver, then we need to extract our own internal
% MoM matrix equation. Note: we can only do this after the solver data (i.e. geometry, etc. is setup)
[Const, zMatricesSUNEM, yVectorsSUNEM] = extractSUNEMMoMmatrixEq(Const, Solver_setup);

% Compare now the above with that calculated using FEKO
compareMatrices(Const, zMatricesFEKO.values, zMatricesSUNEM.values);

% For the RHS vectors, we calculate the rel. error norm % (2-norm)
yVectorErr = calculateErrorNormPercentage(yVectorsFEKO.values, yVectorsSUNEM.values);
message_fc(Const,sprintf('Rel. error norm. for V(RHS) compared to FEKO sol. %f percent',yVectorErr));

% --------------------------------------------------------------------------------------------------
% Run the EM solver 
% --------------------------------------------------------------------------------------------------
% (Note: We either pass our own (internal) matrices, or those read from FEKO)
[Solution] = runEMsolvers(Const, Solver_setup, zMatricesSUNEM, yVectorsSUNEM, xVectorsFEKO);

