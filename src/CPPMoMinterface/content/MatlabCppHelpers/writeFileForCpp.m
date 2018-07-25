function [] = writeFileForCpp(Const)
    %calculate_directivity
    %   Date: 2018.07.25
    %   Usage:
    %       writeFileForCpp(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing simulation parameters
    %
    %   Output Arguments:
    %       None
    %
    %   Description:
    %       Prints relevent information to a .mom file for reading in a C++
    %       program
    %
    %   =======================
    %   Written by Tameez Ebrahim on 2018.07.25

    % --------------------------------
    % Open file
    % --------------------------------
    file_name = [Const.OutputDirName '.mom'];
    fid = fopen(file_name, 'a');
    
    if fid == -1
        disp(['Error: cannot open or create file ' file_name]);
    else
        % --------------------------------
        % Write file header
        % --------------------------------
        C = fix(clock);
        
        fprintf(fid, 'MOM MATLAB OUTPUT FILE\n');
        fprintf(fid, 'Created On: %d-%d-%d, %d:%d:%d\n\n', C(3),C(2),C(1),C(4),C(5),C(6));
        
        % --------------------------------
        % Write Const
        % --------------------------------
        fprintf(fid, 'CONST START\n');
        
        field_names = fieldnames(Const);
              
        for i = 1:length(field_names)
            fprintf(fid, '%-30s\t\t%s\n' ,string(field_names(i)) ,string(Const.(string(field_names(i)))));
        end
        fprintf(fid, 'CONST END\n');
    
    end
    % --------------------------------
    % Write file header
    % --------------------------------
    fclose(fid);

end