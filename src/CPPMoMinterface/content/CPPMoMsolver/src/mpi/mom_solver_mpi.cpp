#include "mom_solver_mpi.h"

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

MoMSolverMPI::MoMSolverMPI(std::vector<Node> nodes,
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

void MoMSolverMPI::calculateVrhsInternally()
{
    // Lets calculate the Vrhs data internally
    // This wil just be for a nomally incident x-directed plane wave
    // TODO make more general
    // TODO change to complex for ease of use in calculations

    Node E(1, 0, 0);
    this->vrhs_internal = std::vector<double>(this->edges.size(), 0);

    for(int i = 0; i < this->edges.size(); i++)
    {
        this->vrhs_internal[i] = E.getDotNoComplex(this->edges[i].getRhoCPlus()) / 2 +
                                 E.getDotNoComplex(this->edges[i].getRhoCMinus()) / 2;
        this->vrhs_internal[i] = this->vrhs_internal[i] * this->edges[i].getLength();  
    }
}

void MoMSolverMPI::calculateJMatrix()
{
    // Lets calcualte the I vector
    // LU-decomposition /w partial pivot
    // Using Eigen3

    // First lets put the values into relevant Eigen datatypes
    // TODO After OpenMP switch all to Matrices to Eigen Datatypes
    // TODO change function name
    Eigen::MatrixXcd m(this->edges.size(), this->edges.size()); 

    for(int i = 0; i < this->edges.size(); i++)
    {
        for(int j = 0; j < this->edges.size(); j++)
        {
            m(i, j) = this->z_mn[i][j];
        }
    }

    Eigen::VectorXcd v(this->edges.size());

    for(int i = 0; i < this->vrhs_internal.size(); i++)
    {
        std::complex<double> temp(this->vrhs_internal[i]);
        v(i) = temp;
    }

    // Now lets solve for I
    Eigen::VectorXcd i_lhs = m.partialPivLu().solve(v);
}

void MoMSolverMPI::calculateZmnByFaceMPI()
{


    // Lets get the mpi rank and size
    // rank -> which process is running
    // size -> total number of processes 
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Lets load the quadrature weights and values for all processes
    int num_quadrature_points = std::stoi(this->const_map["QUAD_PTS"]); 
    this->quadrature_weights_values = getQuadratureWeightsAndValues(num_quadrature_points); 
 

    // Lets assign the p values per process
    // p is from 0 -> number_of_triangles
    // Lets split up p as eqully as possible to the processes
    // First lets declare a vector to store the p values
    std::vector<int> sub_p_values;
    int index;

    // Now lets split up p
    // The variable index is used as a starting point for th values to be split
    // e.g. p is:
    // 0 1 2 3 4 5 6 7 8 9
    // and there are two processes
    // so process 0 will get from 0 -> 4
    // and process 1 will get from 5 -> 9
    if(rank == 0)
    {
        for(int i = 0; i < this->numValuesMPI(size, rank, this->triangles.size()); i++)
        {
            sub_p_values.push_back(i);
        }
    }
    else
    {   
        index = 0;
        for(int i = 0; i < rank; i++)
        {
            index += this->numValuesMPI(size, i, this->triangles.size());
        }
        for(int i = index; i < this->numValuesMPI(size, rank, this->triangles.size()) + index; i++)
        {
            sub_p_values.push_back(i);
        }
    }

    // Lets calculate all the portions of Zmn specified by the processes p values
    std::vector<double> sub_zmn = this->workMPI(sub_p_values);
    
    // Lets gather all the vector sizes
    // Lets first create a vector to store the sizes  
    std::vector<int> proc_vector_size;
    proc_vector_size.resize(size);

    // Now lets get the size of sub_zmn for each process
    int z_mn_size = sub_zmn.size();

    // Now lets send them to the main process
    MPI_Gather(&z_mn_size, 1, MPI_INT, &proc_vector_size[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now that the size of the data is available to the root process
    // Lets gather all the data
    // First lets define two vectors
    // The first is to store the data
    // The seconds is to store the displacements of the data
    std::vector<double> all_zmn_data;
    std::vector<int> disps;

    // Lets resize the disps
    disps.resize(size);

    // Now lets fill in disps
    for(int i = 0; i < size; i++)
    {
        if(i == 0)
        {
            // The first displacement is always 0
            disps[i] = 0;
        }
        else
        {
            disps[i] = disps[i - 1] + proc_vector_size[i - 1];
        }
    }

    // Lets resize all_zmn_data but just for the root process
    if(rank == 0)
    {
        all_zmn_data.resize(disps[size - 1] + proc_vector_size[size - 1]);
    }

    // Lets now gather all the Zmn sub data into the main process
    MPI_Gatherv(&sub_zmn[0], sub_zmn.size(), MPI_DOUBLE, &all_zmn_data[0],
                 &proc_vector_size[0], &disps[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Now lets fill all of Zmn
    if(rank == 0)
    {
        // Lets resize the internal Zmn
        std::vector<std::complex<double>> row_vector(this->edges.size(), 0);
        this->z_mn = std::vector<std::vector<std::complex<double>>>(this->edges.size(), row_vector);

        int index;
        for(int i = 0; i < all_zmn_data.size() / 4; i++)
        {
            index = i * 4;
            std::complex<double> temp(all_zmn_data[index + 2], all_zmn_data[index + 3]);
            this->z_mn[(int)all_zmn_data[index]][(int)all_zmn_data[index + 1]] += temp; 
        }
        // for(int i = 0; i < this->edges.size(); i++)
        // {
        //     for(int j = 0; j < this->edges.size(); j++)
        //     {
        //         std::cout << this->z_mn[i][j];
        //     }
        //     std::cout << std::endl;
        // }
    }
}

int MoMSolverMPI::numValuesMPI(int num_procs, int rank, int data_length)
{
    // Efficiently distribute work between processes
    // https://stackoverflow.com/questions/5657158/how-to-distribute-a-vector-of-n-elements-across-p-processors

    return (data_length + rank) / num_procs;
}

std::vector<double> MoMSolverMPI::workMPI(std::vector<int> p_values) // RENAME
{
    // See MoMSolverMPI::calculateZmnByFace() for full commentary

    std::vector<double> partial_zmn;
    for(int i = 0; i < p_values.size(); i++)
    {
        int p = p_values[i];

        for(int q = 0; q < this->triangles.size(); q++)
        {
            std::vector<Node> a_and_phi = this->calculateAAndPhi(p, q);

            for(int e = 0; e < this->triangles[q].getEdges().size(); e++)
            {
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

                std::complex<double> phi = a_and_phi[3].getXComplexCoord() * this->edges[this->triangles[q].getEdges()[e]].getLength();

                if(this->edges[this->triangles[q].getEdges()[e]].getMinusTriangleIndex() == q)
                {
                    phi = phi * std::complex<double>(-1.0, 0);
                }
                else
                {
                    a_pq = a_pq.getScalarMultiply(-1.0);
                }

                for(int r = 0; r < this->triangles[p].getEdges().size(); r++)
                {
                    Node rho_c;
                    double phi_sign;

                    if(this->edges[this->triangles[p].getEdges()[r]].getMinusTriangleIndex() == p)
                    {
                        rho_c = this->edges[this->triangles[p].getEdges()[r]].getRhoCMinus();
                        phi_sign = -1;
                    }
                    else
                    {
                        rho_c = this->edges[this->triangles[p].getEdges()[r]].getRhoCPlus();
                        phi_sign = 1;
                    }

                    partial_zmn.push_back(this->triangles[p].getEdges()[r]);
                    partial_zmn.push_back(this->triangles[q].getEdges()[e]);

                    std::complex<double> temp_zmn_value = 
                    this->edges[this->triangles[p].getEdges()[r]].getLength() *
                    (this->j * 
                        this->omega *
                        a_pq.getDot(rho_c) / 
                        std::complex<double>(2, 0) -
                        phi * 
                        phi_sign);
                    partial_zmn.push_back(temp_zmn_value.real());  
                    partial_zmn.push_back(temp_zmn_value.imag());  
                }                   
            } 
        }
    }

    return partial_zmn;
}

std::vector<Node> MoMSolverMPI::calculateAAndPhi(int p, int q)
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

std::vector<std::complex<double>> MoMSolverMPI::calculateIpq(int p, int q)
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
    if(p == 19823) // TODO change to p == q && SING == True
    {
        // TODO Add singularity treatment
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




























