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
Const = sunem_initialise('sphere',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;

Const.use_mpi_processes     = 3;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'sphere.mat';
Const.FEKOstrfilename          = 'sphere.str';
Const.FEKOrhsfilename          = 'sphere.rhs';
Const.FEKOoutfilename          = 'sphere.out';
Const.FEKOefefilename          = 'sphere.efe';
Const.FEKOffefilename          = 'sphere.ffe';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      = 'sunem_mom_sphere.str';

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

% TEST
%[CISOL] = readFileFromCpp(Const);
%error = calculateErrorNormPercentage(xVectorsFEKO.Isol, CISOL);
%disp(error);
writeFileForCpp(Const, Solver_setup, yVectorsFEKO);

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

