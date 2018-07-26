function [] = writeFileForCpp(Const, FEKO_data, yVectors)
    %calculate_directivity
    %   Date: 2018.07.25
    %   Usage:
    %       writeFileForCpp(Const, FEKO_data, yVectors)
    %
    %   Input Arguments:
    %       Const: A global struct containing simulation parameters
    %       FEKO_data 
    %           The struct containing the frequency list data (e.g. number
    %           of points, frequency sample values, etc.), triangle information and
    %           basis function setup information.
    %       yVectors
    %           The Yrhs-vector data
    %       
    %       Output Arguments:
    %           None
    %       
    %       Output:
    %           A .mom file containing all of the relevan data from the
    %           input arguments.
    %
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
        
        fprintf(fid, '%-30s\t\t%d\n\n', 'NUM_FIELDS', length(field_names));
              
        for i = 1:length(field_names)
            fprintf(fid, '%-30s\t\t%s\n' ,string(field_names(i)) ,string(Const.(string(field_names(i)))));
        end
        fprintf(fid, 'CONST END\n\n');
        
        
        
        % --------------------------------
        % Write FEKO_data
        % --------------------------------
        fprintf(fid, 'FEKO_DATA START\n');
        
        
        % --------------------------------
        % Write Co-ordinates of the Nodes
        % --------------------------------
        
        fprintf(fid, 'NODES START\n');
        fprintf(fid, '%-20s\t\t%d\n', 'NUM_NODES', length(FEKO_data.nodes_xyz));
        fprintf(fid, '%s\t\t%s\t\t%s\n', 'X_COORD', 'Y_COORD', 'Z_COORD');
        
        for i = 1:length(FEKO_data.nodes_xyz)
            fprintf(fid, '%f\t%f\t%f\n', FEKO_data.nodes_xyz(i, 1),FEKO_data.nodes_xyz(i, 2),FEKO_data.nodes_xyz(i, 3));
        end
        
        fprintf(fid, 'NODES END\n\n');
        
        % --------------------------------
        % Write Triangle Information
        % --------------------------------
        fprintf(fid, 'TRIANGLES START\n');
        fprintf(fid, '%-20s\t\t%d\n\n', 'NUM_TRIANGLES', FEKO_data.num_metallic_triangles);
        fprintf(fid, '%s\t%s\t%s\t%s\t\t%s\n', 'V', 'V', 'V', 'CENTRE', 'AREA');
        
        for i = 1:FEKO_data.num_metallic_triangles
            fprintf(fid, '%d\t%d\t%d\t%f\t%f\n', FEKO_data.triangle_vertices(i, 1), ...
                                                 FEKO_data.triangle_vertices(i, 2), ... 
                                                 FEKO_data.triangle_vertices(i, 3), ...
                                                 FEKO_data.triangle_centre_point(i),...
                                                 FEKO_data.triangle_area_m2(i));
        end
        fprintf(fid, 'TRIANGLES END\n\n');
        
        % --------------------------------
        % Write Edge Information
        % --------------------------------
        fprintf(fid, 'EDGES START\n');
        fprintf(fid, '%-20s\t\t%d\n\n', 'NUM_EDGES', FEKO_data.num_metallic_edges);
        fprintf(fid, '%s\t%s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s', ...
                        'V', ...
                        'V', ...
                        'CENTRE_X', ...
                        'CENTRE_Y', ...
                        'CENTRE_Z', ...
                        'LENGTH', ...
                        'TMINUS', ...
                        'TPLUS', ...
                        'TMFREEV', ...
                        'TPFREEV', ...
                        'RHOCMX', ...
                        'RHOCMY', ...
                        'RHOCMZ', ...
                        'RHOCPX', ...
                        'RHOCPY', ...
                        'RHOCPZ');
        fprintf(fid, '\n'); % For some reason a newline doesn't get created in prev fprintf
        
        for i = 1:FEKO_data.num_metallic_edges
            fprintf(fid, '%d\t%d\t%8f\t%8f\t%8f\t%8f\t%8d\t%8d\t%8d\t%8d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\n', ...
                          FEKO_data.rwg_basis_functions_shared_edge_nodes(i, 1), ...
                          FEKO_data.rwg_basis_functions_shared_edge_nodes(i, 2), ...
                          FEKO_data.rwg_basis_functions_shared_edge_centre(i, 1), ...
                          FEKO_data.rwg_basis_functions_shared_edge_centre(i, 2), ...
                          FEKO_data.rwg_basis_functions_shared_edge_centre(i, 3), ...
                          FEKO_data.rwg_basis_functions_length_m(i), ...
                          FEKO_data.rwg_basis_functions_triangleMinus(i), ...
                          FEKO_data.rwg_basis_functions_trianglePlus(i), ...
                          FEKO_data.rwg_basis_functions_triangleMinusFreeVertex(i), ...
                          FEKO_data.rwg_basis_functions_trianglePlusFreeVertex(i), ...
                          FEKO_data.rho_c_mns(i, 1), ...
                          FEKO_data.rho_c_mns(i, 2), ...
                          FEKO_data.rho_c_mns(i, 3), ...
                          FEKO_data.rho_c_pls(i, 1), ...
                          FEKO_data.rho_c_pls(i, 2), ...
                          FEKO_data.rho_c_pls(i, 3));
        end
        fprintf(fid, 'EDGES END\n');  
        fprintf(fid, 'FEKO_DATA END\n\n');
        
        
        
        % --------------------------------
        % Write yVectors (Vrhs)
        % --------------------------------
        fprintf(fid, 'VRHS START\n');
        fprintf(fid, '%-20s\t\t%d\n\n', 'NUM_VALUES', length(yVectors.values));
        fprintf(fid, '%s\n', 'VALUE');
        
        for i = 1:length(yVectors.values)
            fprintf(fid, '%f\n', yVectors.values(i));
        end
        
        fprintf(fid, 'VRHS END');
        
    
    end

    fclose(fid);

end