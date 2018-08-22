#include "mom_solver.h"

/**
 *   \file mom_solver.cpp
 *   \brief Solve the method of moments(MoM) of an electromagnetic problem
 *
 *  Detailed description
 *  This class will first solve for Zmn using a face pair approach. Then Vm = ZmnIn will be solved
 *  for In using LU decomposition on Zmn. 
 *
 *  Author:  Tameez Ebrahim
 *  Created: 27 July 2018
 *
 */

#include <cfenv>

MoMSolver::MoMSolver(std::vector<Node> nodes,
                     std::vector<Triangle> triangles,
                     std::vector<Edge> edges,
                     std::vector<double> vrhs,
                     std::map<std::string, std::string> const_map)
{
    // This is the constructor for the MoMSolver class
    // It assigns the arguments to internal variables and also declares some constants
  
    this->nodes = nodes;
    this->triangles = triangles;
    this->edges = edges;
    this->vrhs = vrhs;
    this->const_map = const_map;

    // Lets define some constants
    // All the constants are given, or calculated by those given in the const_map
    // Remember that this data is from the .mom file so it would be wise to rather
    // change the initial variables there than here in the code
    // Also note that const_map is of type std::map<std::string, std::string>
    // so all values need to be converted 
    // Lets start by getting the speed of light
    this->c = std::stod(this->const_map["C0"]);

    // Then lets get the frequency
    // There are three frequency values in the const_map
    // freqStart, freqEnd and freqData
    // The one needed is freqData
    this->frequency = std::stod(this->const_map["freqData"]);

    // Lets get omega(w)
    // This is used in equations 17 and 33 of RWG80
    // The formula for w is w = 2 * pi * frequency
    // Note that pi is called as M_PI from the <cmath> library
    // The math library naming is quite strange so it is necessary
    // to remember that <cmath> is for C++ while <math> is for C
    // If M_PI happens to not work add #define _USE_MATH_DEFINES before <cmath>
    // in mom_solver.h
    this->omega = 2 * M_PI * this->frequency;

    // Lets get lambda
    // This lambda is not to be confused with the lambdas used in the quadrature step
    // The formula for lambda is lambda = speed of light / frequency
    this->lambda = this->c / this->frequency;

    // Then lets get k(propagation constant)
    // The formula for k is k = 2 * pi / lambda
    // If there is an issue with M_PI, see the comment on omega for details
    this->k = 2 * M_PI / this->lambda;

    // Lets finally define j(sqrt(-1)), the imaginary constant
    // j is synonymous with i in MATLAB
    // Currently, I can't find if <complex> defines it
    // So, lets set j as 0 + 1i
    std::complex<double> complex_constant(0, 1);
    this->j = complex_constant;
}

void MoMSolver::calculateZmnByFace()
{
    // Lets time how long it takes to get Zmn
    // For reference on a problem of around mxn = 304x304 takes MATLAB 10 seconds
    // As of 04 AUgust 2018, this takes 1.3 seconds
    // Lets use the high precision clock from <chrono> to be accurate
    // std::chrono::high_resolution_clock::time_point z_mn_time_start = std::chrono::high_resolution_clock::now();
    // TIME PROFILE
    // After 1000 iterations
    // ---- Total Times ----
    // Time to fill Zmn = 1.2562s
    // Time to calculate A and Phi = 0.5537114s
    // Time to calculate 4 I's = 0.180459s
    //
    // ---- Breakdown ----
    // Zmn : 0.719081
    // A_P : 0.356655
    // _I_ : 0.180459
    //
    // The profile shows that the most time is spent in Zmn

    // Lets calculate Zmn by face pair combinations
    // This will be done according to RWG80

    // Before starting, lets assign the quadrature weights and points
    // First lets get the number of quadrature points from const_map
    // Remember that all values in const_map are strings
    int num_quadrature_points = std::stoi(this->const_map["QUAD_PTS"]); 

    // Then lets assign the weights and values
    this->quadrature_weights_values = getQuadratureWeightsAndValues(num_quadrature_points); 

    // Lets start by looping over the faces twice to get faces p and q
    // p -> observation triangle 
    // q -> source triangle

    // Declare Zmn TODO: change
    // Remember that Zmn will be complex
    // It is of size mxn where m == n == number of edges
    // The number of edges is defined in const_map
    //std::complex<double> z_mn[this->edges.size()][this->edges.size()] = {std::complex<double>(0,0)}; 
    // std::vector<std::complex<double>> row_vector(this->edges.size(), 0);
    // this->z_mn = std::vector<std::vector<std::complex<double>>>(this->edges.size(), row_vector);
    this->z_mn = new std::complex<double>[this->edges.size() * this->edges.size()];
    for(int i = 0; i < this->edges.size(); i++)
    {
        for(int j = 0; j < this->edges.size(); j++)
        {
        this->z_mn[j * this->edges.size() + i] = std::complex<double>(0, 0);
        }
    }

    for(int p = 0; p < this->triangles.size(); p++)
    {
        for(int q = 0; q < this->triangles.size(); q++)
        {
        // Lets calculate Apq and Phipq
	    // Apq is the magnetic vector potential and Phipq is the scalar potential
	    // For brevity, they will be referred to as A and Phi
	    // They are found using equations 32 and 33 in RWG80
	    // and are used in equation 17 in RWG80
        // To get them, the four I values need to be calculated
	    // These I values are noted in equations 34a-d in RWG80
        // Lets put them in a node vector
        // The form is [A_1 A_2 A_3 Phi]
	    // There are three A values returned. Each corresponds to a triangle
	    // edge. In equation 32 of RWG80 there exists the ri coordinate.
	    // 1 O\
	    //   | O 3  <-- This is a triangle
	    // 2 O/
	    // If we imagine the the above is a triangle with the O's as it's vertices
	    // the ri coordinate will refer to each of the vertices in turn.
	    // Now, how does one choose which A value of the three to use?
	    // This is where the free vertices come into play
	    // The Triangle class has three vertices [v1 v2 v3] >> [1 2 3] in the illustration
	    // A_1 has ri = v1, A_2 has ri = v2 and A_3 has ri = v3
	    // The edge is checked for the free vertex and then the appropriate A is used
	    // If the horizontal edge, | in the illustration needs to be assigned an A value,
	    // then it is trivial to see that A_3 will be used since the free vertex is vertex 3
	    // Remember that in the mesh, triangles will rarely look like the illustration so it
	    // is necessary to check for the free vertex in code rather than pre-empting what it will be
        std::vector<Node> a_and_phi = this->calculateAAndPhi(p, q);

	        // Now let loop over the source triangle edges.
	        // It is important to note that since the only edges associated with a
	        // triangle are the interior(not boundary) edges
	        // This means that the number of edges varies and is not a constant three
            for(int e = 0; e < this->triangles[q].getEdges().size(); e++)
            {
	            // Let us get the right A value for the edge
	            // The explanation on how to choose the correct A value is ^^ (two comments up)
                // Remember that A is agnostic of the triangles positivity so it is fine
	            // to check both the minus and plus free vertices
	            // Also remember that the length of the edge was not multiplied in the calculateAAndphi function.
	            // Only now is the pairing between the edge and the A value so it is imperative not to forget
	            // to multiply by the length
	            // The if statement conditions look quite messy, but remember that both the Edge and Triangle
	            // classes only have indices to the other
                Node a_pq;
                if(this->triangles[q].getVertex1() == this->edges[this->triangles[q].getEdges()[e]].getPlusFreeVertex() ||
                   this->triangles[q].getVertex1() == this->edges[this->triangles[q].getEdges()[e]].getMinusFreeVertex())
                {
                    a_pq = a_and_phi[0].getScalarMultiply(this->edges[this->triangles[q].getEdges()[e]].getLength());   
                }

                if(this->triangles[q].getVertex2() == this->edges[this->triangles[q].getEdges()[e]].getPlusFreeVertex() ||
                   this->triangles[q].getVertex2() == this->edges[this->triangles[q].getEdges()[e]].getMinusFreeVertex())
                {
                    a_pq = a_and_phi[1].getScalarMultiply(this->edges[this->triangles[q].getEdges()[e]].getLength());   
                }

                if(this->triangles[q].getVertex3() == this->edges[this->triangles[q].getEdges()[e]].getPlusFreeVertex() ||
                   this->triangles[q].getVertex3() == this->edges[this->triangles[q].getEdges()[e]].getMinusFreeVertex())
                {
                    a_pq = a_and_phi[2].getScalarMultiply(this->edges[this->triangles[q].getEdges()[e]].getLength());   
                }

		        // Lets create the variable Phi
		        // Since C++ only supports the returning of one entity, Phi was inserted into a node
		        // as the x coordinate
		        // This is just a little workaround to get the value easily
		        // As with A above, now that the edge is known it is important to remember to multiply
		        // the recieved phi value by the edge length
                std::complex<double> phi = a_and_phi[3].getXComplexCoord() * this->edges[this->triangles[q].getEdges()[e]].getLength();

		        // Now lets check if the source triangle is a minus or plus triangle for the edge
		        // This is checked using the plus/minus triangle index variable contained in the Edge
		        // class
                if(this->edges[this->triangles[q].getEdges()[e]].getMinusTriangleIndex() == q)
                {
		            // The triangle is a minus triangle for the edge
		            // Lets multiply phi by -1 as shown in equation 33 in RWG80
		            // Remember that phi is a complex number so it has to be multiplied
		            // with another complex number
		            // This does not always have to be done, but it is safer to do so
		            // The compiler will sometimes throw a type exeption error if forgotten
                    phi = phi * std::complex<double>(-1.0, 0);
                }
                else
                {
		            // The triangle is a plus triangle for the edge
		            // Lets multiply a_pq with -1 as shown in equation 32 in RWG80
		            // The getScalarMultiply member of the Node class is used for
		            // easy scalar multiplication of a vector
                    a_pq = a_pq.getScalarMultiply(-1.0);
                }

                for(int r = 0; r < this->triangles[p].getEdges().size(); r++)
                {
		            // Lets loop over the edges of the observation triangle
		            // These edges will be the m index of the Zmn matrix
                    // Lets first create two variables, one for rho_c
		            // and another for the phi sign
		            Node rho_c;
                    double phi_sign;

		            // Let us check whether the triangle is a plus or minus triangle for the edge
		            // As noted above, remember that the classes Edge and Triangle only contain
		            // indices to the other
                    if(this->edges[this->triangles[p].getEdges()[r]].getMinusTriangleIndex() == p)
                    {
		                // The triangle is a minus triangle for the edge
		                // Lets get the minus rho c node and assign it to rho_c
		                // Lets also set the phi_sign to -1
                        rho_c = this->edges[this->triangles[p].getEdges()[r]].getRhoCMinus();
                        phi_sign = -1;
                    }
                    else
                    {
		                // The triangle is a plus triangle for the edge
		                // Lets get the plus rho c node and assign it to rho_c
		                // Lets also set the phi_sign to 1
                        rho_c = this->edges[this->triangles[p].getEdges()[r]].getRhoCPlus();
                        phi_sign = 1;
                    }

                    // Now that all the necessary values are assembled, lets fill in some
		            // of the relevant entries of Zmn
		            // The index n is from the source triangle and
		            // the index m is from the observation triangle
		            // It is necessary to remember that this is not the final Zmn value being
		            // assigned, but only a piece of it contributed by the source and observation
		            // triangles
		            // The formula for the final Zmn value is found in equation 17 of RWG80
		            // As just noted, since this is just a contribution to the final value
		            // the formula needs to be tweaked as such
		            // Zmn = Zmn + edge_length_of_observation_triangle *
		            //             ((j * omega * A * rho_c / 2) - phi * phi_sign)
                    

                    this->z_mn[this->triangles[p].getEdges()[r] + this->triangles[q].getEdges()[e] * this->edges.size()] +=
                        this->edges[this->triangles[p].getEdges()[r]].getLength() *
                            (this->j * 
                            this->omega *
                            a_pq.getDot(rho_c) / 
                            std::complex<double>(2, 0) -
                            phi * 
                            phi_sign);  
                }                   
            } 
        }
    }
}

void MoMSolver::calculateVrhsInternally()
{
    // Lets calculate the Vrhs data internally
    // This wil just be for a nomally incident x-directed plane wave
    // TODO: make more general

    Node E(1, 0, 0);
    this->vrhs_internal.resize(this->edges.size());

    for(int i = 0; i < this->edges.size(); i++)
    {
        this->vrhs_internal[i] = E.getDotNoComplex(this->edges[i].getRhoCPlus()) / 2 +
                                 E.getDotNoComplex(this->edges[i].getRhoCMinus()) / 2;
        this->vrhs_internal[i] = this->vrhs_internal[i] * this->edges[i].getLength();  
    }
}

void MoMSolver::calculateJMatrixLAPACK()
{
    // Lets calculate the Ilhs(current)

    // Lets start by defining the matix size
    // Zmn will always be square so n_rows = n_cols
    int matrix_size = this->edges.size();

    // Lets define the lda for Zmn
    // This is equal to matrix_size due to Zmn being square
    // This can just be replaced by matrix_size, but it is left
    // just for completeness sake.
    int zmn_lda = std::max(1, matrix_size);
    
    // Lets define the pivot array
    // LAPACK will use this array to store the pivot
    // information from the LU decomposition
    int piv[matrix_size];
    
    // Info is the variable needed to gauge the success of the LAPACK
    // routines. If equal to 0, then everything went fine.
    // Lets set it to an arbitrary number
    int info = 256;

    // Lets get the LU decomposition of Zmn first
    // zgetrf_ is a LAPACK routine to calculate the LU decomposition
    // of a matrix. See the LAPACK website for full function details
    // Remember that LAPACK overwrites this->z_mn with it's LU
    // decomposition so the Zmn matrix will not be available after.
    // If Zmn is needed, then just make a copy before calling the function
    zgetrf_(&matrix_size, &matrix_size, this->z_mn, &zmn_lda, piv, &info);

    // Lets now calculate Ilhs(current)
    // Lets first define the matrix orientation
    char tran = 'N';

    // Now lets define the number of columns in the B matrix
    // In this case, the B matric is Vrhs with one column
    // It is important to remember that LAPACK routines are written
    // in Fortran which only supports pass by reference
    int one = 1;

    // Now lets call zgetrs_ which solves for Ilhs(current) using
    // the LU decomposition from zgetrf_
    // The full function details can be found on the LAPACK website
    zgetrs_(&tran, &matrix_size, &one, this->z_mn, &matrix_size, piv, &this->vrhs_internal[0], &matrix_size, &info);
}

std::vector<Node> MoMSolver::calculateAAndPhi(int p, int q)
{
    // Lets calculate the magnetic vector potential and the scalar potential
    // The formulae are equations 32 and 33 in RWG80
    // First, lets get the four I's calculated in calculateIpq 
    // Remember, the vector is ordered in [Ipq Ipq_xi Ipq_eta Ipq_zeta]
    std::vector<std::complex<double>> i_vector = this->calculateIpq(p, q);

    // There are three values for the magnetic vector potential (Apq) and one for the scalar
    // potential(Phipq)
    // Remember, from calculateIpq, r1 -> vertex_1, r2 -> vertex_2 and r3 -> vertex_3
    // Lets first create the return vector
    std::vector<Node> a_and_phi_vector;

    // Lets loop over the three triangle vertices
    for(int i = 0; i < 3; i++)
    {
        // Lets do the multiplication of the r's by the I's
        // It's a bit messy with needing to create a new node for each multiplication 
        // r1 * I_xi        -> a_pq_node_1 
        // r2 * I_eta       -> a_pq_node_2
        // r3 * I_zeta      -> a_pq_node_3
        // ri * I           -> a_pq_node_4
        Node a_pq_node_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(i_vector[1]);
        Node a_pq_node_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(i_vector[2]);
        Node a_pq_node_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(i_vector[3]);
        Node a_pq_node_4 = this->nodes[this->triangles[q].getVertices()[i]].getScalarMultiply(i_vector[0]);

        // Now lets add the four terms together
        // Again, it is quite messy
        Node a_pq_node = a_pq_node_1.getAddNode(a_pq_node_2);
        a_pq_node = a_pq_node.getAddNode(a_pq_node_3);
        a_pq_node = a_pq_node.getSubtractComplexNode(a_pq_node_4);

        
        // Now lets multiply the added nodes by mu / 4pi
        // Remember that the length is not multiplied yet
        double a_pq_multiplier = std::stod(this->const_map["MU_0"]) / (4 * M_PI);                
        a_pq_node = a_pq_node.getScalarMultiply(a_pq_multiplier);
        a_and_phi_vector.push_back(a_pq_node);
    }

    // Lets create the multiplier for phi_pq
    // The formula is in equation 33 in RWG80
    // Remember that the edge length is not multiplied yet since the edge is unknown
    // The formula is 1 / (j * 2 * pi * omega * epsilon_0)
    // The compiler sometimes throws type exeption errors so it is necessary to
    // cnvert pi and 2 into complex values
    // If there is a problem with using M_PI, the solution is in the constructor
    std::complex<double> phi_pq_multiplier = std::complex<double>(1, 0) / 
                                            (this->j * std::complex<double>(2.0, 0) * 
                                            std::complex<double>(M_PI, 0) *
                                            this->omega *
                                            std::stod(this->const_map["EPS_0"])); 
    
    // Now lets multiply the multiplier with Ipq as shown in RWG80
    std::complex<double> phi_pq = i_vector[0] * phi_pq_multiplier;

    // C++ does not support multiple return values
    // A solution to this is to pass a pointer to the function and change its value
    // This seems a bit out of place considering the flow of the rest of the code
    // so a different approach will be used
    // The return type of the function is std::vector<Node>
    // It is therefore easy to decide to just convert the phi_pq value to a node
    // This is done below
    // It is important to remember that phi_pq is NOT a node, but just a complex
    // value masquerading as one to get out of the function
    // It is therefore necessary to account for this when calling the function
    // to avoid any nasty suprises
    Node phi_pq_node(phi_pq, 0, 0);
    a_and_phi_vector.push_back(phi_pq_node);    

    return  a_and_phi_vector;
}

std::vector<std::complex<double>> MoMSolver::calculateIpq(int p, int q)
{
    // Lets calculate the 4 Ipq integrals(equations 34a-d in RWG80)
    // The first three(a-c) will be calculated using numerical quadrature
    // The final one(d) is calculated from the first three
    // p and q, the triangle indices are the only arguments needed
    // All the other parameters are specified in const_map
    // If p == q, then singularity treatment may be applied
    // This is defined in const_map whether to use the singularity treatment or not
    std::vector<std::complex<double>> i_vector;


    // Lets start by checking for a singularity(p == q)
    if(p == 19823) // TODO: change to p == q && SING == True
    {
        // TODO: Add singularity treatment
        int x = 0;
    }
    else
    {
        // Now lets calculate the four I's using numerical quadrature
        // Remember that the weights and values for quadrature are already stored in 
        // quadrature_weights_values
        // It is of the form [weight lambda_1 lambda_2 lambda_3]
        // with the number of rows equalling the number of quadrature points

        // First lets create the four I's
        // Remember that they are complex
        std::complex<double> Ipq;        // RWG80 34a
        std::complex<double> Ipq_xi;     // RWG80 34b
        std::complex<double> Ipq_eta;    // RWG80 34c
        std::complex<double> Ipq_zeta;   // RWG80 34d

        // Now lets loop over the quadrature points
        for(int i = 0; i < this->quadrature_weights_values.size(); i++)
        {
            // Lets first get r' as noted in equation 30 in RWG80
            // So, it is [xi * vertex_1, eta * vertex_2, zeta * vertex_3]  
            // xi_r_1       -> xi * r1
            // eta_r_2      -> eta * r2
            // zeta_r_3     -> zeta * r3
            Node xi_r_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(this->quadrature_weights_values[i][1]); 
            Node eta_r_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(this->quadrature_weights_values[i][2]); 
            Node zeta_r_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(this->quadrature_weights_values[i][3]);

            Node r_prime = xi_r_1.getAddNode(eta_r_2);
            r_prime = r_prime.getAddNode(zeta_r_3); 

            // Now, lets get Rp as noted in equation 27 in RWG80
            // C++ doesn't have an easy way to take the norm of a vector
            // Lets do it manually from the node class
            // Here is Rp = norm(triangle_centre - r_prime)
            double R_p = this->triangles[p].getCentre().getDifferenceBetween(r_prime).getNorm();

            // Now lets get Greens function
            // That is, Green = (e^(-j*k*R_p)) / R_p
            // Remeber that the imaginary unit 'j' is defined as this->j
	        // This is noted in as the common term in equations 34a-d in RWG80
            std::complex<double> greens_function = std::exp(std::complex<double>(-1.0, 0) * this->j * this->k * R_p) / R_p;

            // Finally, lets calculate three  of the four I's
	        // These equations are the same as equations 34a-d in RWG80
	        // The numerical quadrature rule for a triangular domain is
	        // I = area * sum<num_quadrature_components>(wi (lambda_1, lambda_2, lambda_3))
	        // It is important not to confuse theses lambdas with the single lambda defined
	        // in the constructor
	        // It is easily noticed that no area value is being multiplied in the instances below.
	        // The local coordinate system is being used as noted in equation 28 in RWG80 so the
	        // area referred to in the quadrature formula is equal to 1
	        // It is also important to note that while not defined in the formula, w needs to be normalised
	        // This is done by w = 0.5 * w
            Ipq = Ipq + (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * greens_function);  
               
            Ipq_xi = Ipq_xi + (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                               this->quadrature_weights_values[i][1] * 
                               greens_function);

            Ipq_eta = Ipq_eta + (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                                 this->quadrature_weights_values[i][2] * 
                                 greens_function);  
        }

        // Lets calculate the final I and push to the vector
        Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
        i_vector.push_back(Ipq);
        i_vector.push_back(Ipq_xi);
        i_vector.push_back(Ipq_eta);
        i_vector.push_back(Ipq_zeta);
    }
    return i_vector;
}




























