function [yVectors] = calculateVrhsProper(Const, Solver_setup)
    
    % For an incident plane wave
    % Only for x-directed plane wave
    freq = Solver_setup.frequencies.samples(1);
    k = 2 * pi / (Const.C0 / freq);
    %jk = [0 0 1i*k];
    jk = 1i * k;
    
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;

    yVectors.values = zeros(num_dofs,1); % RHS matrix
    E_pls = zeros(num_dofs,3); % Sampled incident electric field
    E_mns = zeros(num_dofs,3);
    % x-polarized field, unity magnitude
    %E(:,1) = EMag; % Special case
    %E(:,2) = 1; % Special case
    
    for m = 1:num_dofs
        t_pls = Solver_setup.rwg_basis_functions_trianglePlus(m);
        t_mns = Solver_setup.rwg_basis_functions_triangleMinus(m);
        %E_pls(m) = exp(dot(jk, Solver_setup.triangle_centre_point(t_pls,:)));
        %E_mns(m) = exp(dot(jk, Solver_setup.triangle_centre_point(t_mns,:)));
        
        % MATLAB has a bug in their dot product
        % If you do a dot product between a complex and real number, 
        % it incorrectly uses the complex conjugate on the complex number
        % when it should be taking the conjugate of the real number for no
        % net change
        E_pls(m) = exp(jk * Solver_setup.triangle_centre_point(t_pls, 3));
        E_mns(m) = exp(jk * Solver_setup.triangle_centre_point(t_mns, 3));               
    end
    
    for m = 1:num_dofs              
        yVectors.values(m) = E_pls(m,1) * rho_c_pls(m,1)/2 + E_mns(m,1) * rho_c_mns(m,1)/2;
        yVectors.values(m) = ell(m)*yVectors.values(m);
    end