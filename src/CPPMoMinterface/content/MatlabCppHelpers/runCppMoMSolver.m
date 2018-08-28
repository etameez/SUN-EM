function [mom] = runCppMoMsolver(Const, Solver_setup, yVectors, refIsol)
    %runCppMoMsolver
    %   Usage:
    %       [mom] = runCppMoMsolver(Const, Solver_setup, yVectors, refIsol)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       yVectors
    %           The Yrhs-vector data
    %       refIsol
    %           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
    %
    %   Output Arguments:
    %       mom
    %           Structs containing MoM solution and timing data
    %
    %   Description:
    %       Runs the MoM solution based on the Z and Y data that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files using a C++
    %       engine
    %
    %   =======================
    %   Written by Tameez Ebrahim on August 28, 2018.
    %   Stellenbosch University
    narginchk(4,4);
    
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running C++ MoM solver'));
    
    % Initialise the return values
    mom  = [];
    mom.name = 'mom';
    mom.numSols = 1;
    
    % Lets write the .mom file for C++ to read
    writeFileForCpp(Const, Solver_setup, yVectors);
    
    % Lets run the C++ solver
    
    if isunix
       tic;
       if Const.useMPI
           cpp_exec_path = [Const.cppBuildDir 'build_mpi/mom_mpi '];
           dot_mom_path = [pwd '/' Const.OutputDirName '.mom'];
           mpi_cmd = ['mpiexec -np ' num2str(Const.num_procs) ' '];
           omp_cmd = ['-genv OMP_NUM_THREADS=' num2str(Const.num_threads) ' '];
           full_cmd = [mpi_cmd omp_cmd cpp_exec_path dot_mom_path];
           [status, cmdout] = unix(full_cmd);
           
       else
           cpp_exec_path = [Const.cppBuildDir 'build/mom'];
           dot_mom_path = [pwd '/' Const.OutputDirName '.mom'];
           full_cmd = [cpp_exec_path ' ' dot_mom_path]; 
           [status, cmdout] = unix(full_cmd);
       end
       if status ~= 0
           error(cmdout);
       else
           mom.totsolTime = toc;
           message_fc(Const,sprintf('Finished MoM solver in %f sec.',mom.totsolTime));
       end
    else
        message_fc(Const,sprintf('C++ engine currently only supports Linux'));
        error('C++ engine currently only supports Linux');
    end

    
    
    % Lets read the currents from the .sol file
    mom.Isol = readFileFromCpp(Const); 
    mom.relError = calculateErrorNormPercentage(refIsol.Isol, mom.Isol);
    message_fc(Const,sprintf('Rel. error norm. for Sol. to reference sol. %f percent', mom.relError));
end

























