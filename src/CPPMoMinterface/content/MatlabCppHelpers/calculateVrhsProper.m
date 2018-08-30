function [yVectors] = calculateVrhsProper(Const, Solver_setup, EMag, theta_0, phi_0, prop_direction)
    
    narginchk(2,2);
    % For an incident plane wave
    % Only for x-directed plane wave
    freq = Solver_setup.frequencies.samples(1);
    k = 2 * pi / (Const.C0 / freq);
    %jk = [0 0 1i*k];
    %jk = 1i * k;
    
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;

    yVectors.values = zeros(num_dofs,1); % RHS matrix
    E_pls = zeros(num_dofs,3); % Sampled incident electric field
    E_mns = zeros(num_dofs,3);
    
    theta_rad = Const.theta_0 * Const.DEG2RAD;
    phi_rad = Const.phi_0 * Const.DEG2RAD;
    
    if Const.prop_direction == 0
        Ex = Const.EMag * cos(phi_rad) * cos(theta_rad);
        Ey = Const.EMag * sin(phi_rad) * cos(theta_rad);
        Ez = Const.EMag * -sin(theta_rad);
        
        if theta_rad == pi
            Ez = 0;
        end
        
        E = [Ex Ey Ez]; 
    end
    
    kx = sin(theta_rad) * cos(phi_rad);
    ky = sin(theta_rad) * sin(phi_rad);
    kz = cos(theta_rad);
    jk = [1i * k * kx 1i * k * ky 1i * k * kz];
    
    for m = 1:num_dofs
        t_pls = Solver_setup.rwg_basis_functions_trianglePlus(m);
        t_mns = Solver_setup.rwg_basis_functions_triangleMinus(m);
        
        E_pls(m,:) = E.* exp(dot(Solver_setup.triangle_centre_point(t_pls,:), jk));
        E_mns(m,:) = E.* exp(dot(Solver_setup.triangle_centre_point(t_mns,:), jk));               
    end
    
    for m = 1:num_dofs              
        yVectors.values(m) = 0.5 * dot(rho_c_pls(m,:), E_pls(m,:)) + ...
                             0.5 * dot(rho_c_mns(m,:), E_mns(m,:));
        yVectors.values(m) = ell(m)*yVectors.values(m);
    end