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
    Nmom = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    mom.numSols = 1;
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    numRHSperFreq = 1;         % The number of solutions per frequency point.
                                                   % For now, should be 1 (TO-DO: Update)
                                                   
    % Some info about the solution configurations
    message_fc(Const,sprintf('  numSols : %d', mom.numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));
    
    % Calculate the solution vector (all frequency points, all RHSes)
    mom.Isol = complex(zeros(Nmom,mom.numSols));
    
    % The timing calculations also need to take into account that there is a
    % frequency loop
    mom.setupTime = zeros(1,numFreq);
    mom.factorisationTime = zeros(1,numFreq);
    % Zero also the total times (for all frequency iterations)
    mom.totsolTime = 0.0;
    
    for freq=1:numFreq
    
        % Lets write the .mom file for C++ to read
        Const.cppFreq = Const.freqData(freq);
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
                mom.totsolTime = mom.totsolTime + toc;
                freqTime = toc;
                message_fc(Const,sprintf('Finished MoM solver for frequency %f in %f sec.',Const.cppFreq, freqTime));
            end
        else
            message_fc(Const,sprintf('C++ engine currently only supports Linux'));
            error('C++ engine currently only supports Linux');
        end 
    
        % Lets read the currents from the .sol file
        mom.Isol(:,freq) = readFileFromCpp(Const); 

    end
    message_fc(Const,sprintf('Finished MoM solver in %f sec.',mom.totsolTime));
    
    mom.relError = zeros(1,mom.numSols);
    for freq=1:numFreq
        for solNum=1:numRHSperFreq
            index = solNum + (freq-1)*numRHSperFreq;
            mom.relError(index) = calculateErrorNormPercentage(refIsol.Isol(:,index), mom.Isol(:,index));
            message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d of freq. %d of %d compared to reference sol. %f percent',solNum, ...
                numRHSperFreq, freq, numFreq, mom.relError(index)));
        end
    end
    
    % Write sol to .str file
    
end

























