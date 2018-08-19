function [IsolCPP] = readFileFromCpp(Const)
    %calculate_directivity
    %   Date: 2018.08.19
    %   Usage:
    %       IsolCPP = readFileFromCpp(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing simulation parameters
    %       
    %       Output Arguments:
    %           IsolCPP
    %       
    %       Output:
    %           A 1d complex vector containing the solution from Cpp
    %
    %   Description:
    %       Reads the solution from a .sol file written by C++ mom solver.
    %
    %   =======================
    %   Written by Tameez Ebrahim on 2018.07.25
       
    file_name = [Const.OutputDirName '.sol'];
    fid = fopen(file_name, 'r');
    
    if fid == -1
        disp(['Error: cannot read file ' file_name]);
    else
        IsolCPP = cell2mat( textscan(fid, '%f') );
    end