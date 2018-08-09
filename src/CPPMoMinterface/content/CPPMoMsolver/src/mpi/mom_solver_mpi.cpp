#include "mom_solver_mpi.h"

/**
 *   \file mom_solver.cpp
 *   \brief Solve the method of moments(MoM) of an electromagnetic problem
 *
 *  Detailed description
 *  This class will first solve for Zmn using a face pair approach. Then Vm = ZmnIn will be solved
 *  for In using LU decomposition on Zmn. This is the parallel implementation. For more detailed 
 *  comments see the serial implementation. The comments in this file will mostly focus on the 
 *  parallel commands. 
 *
 *  Author:  Tameez Ebrahim
 *  Created: 09 August 2018 
 *
 */

MoMSolverMPI::MoMSolverMPI(std::vector<Node> nodes,
                     std::vector<Triangle> triangles,
                     std::vector<Edge> edges,
                     std::vector<double> vrhs,
                     std::map<std::string, std::string> const_map)
{
    this->nodes = nodes;
    this->triangles = triangles;
    this->edges = edges;
    this->vrhs = vrhs;
    this->const_map = const_map;

    // Lets define some constants
    this->c = std::stod(this->const_map["C0"]);
    this->frequency = std::stod(this->const_map["freqData"]);
    this->omega = 2 * M_PI * this->frequency;
    this->lambda = this->c / this->frequency;
    this->k = 2 * M_PI / this->lambda;
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
    // A quick explanation of MPI
    // Say there are 4 processes numbered 0, 1, 2, 3
    // They all run in parallel, i.e. at the same time
    // Process 0 is also called the root or master process
    // Each of these processes are independent of the other
    // Imagine it as each process is running on a different computer
    // This means that nothing is explicitly shared between them
    // They however can communicate using MPI
    // Hence the name Message Passing Interface
    // There are two methods to pass messages, blocking and non-blocking
    // This function uses blocking and so it will be explained here
    // Blocking means that the process will wait for the message to be sent or received
    // before continuing. Lets see an example
    // Say MPI_Recv(some args here) is used. The process that is receiving the message will not
    // do anything until it gets the message. This allows for easier program flow
    // For improved performance, non-blocking can be used but it is not viable for this
    // application.
    // The integral send/receive functions that will be used here are MPI_Gather() and MPI_Gatherv
    // They will be explained below as used

    // Lets get the mpi rank and size
    // rank -> which process is running
    // size -> total number of processes 
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Lets load the quadrature weights and values for all processes
    // Remember that each process has its own copy of this and most following
    // declarations. If something is restricted, it will be enclosed by an
    // if-statement. So, for example, 
    // if(rank == some_process)
    // {
    //      Do something only in process some_process
    // }
    // It is important to remember this as it sometimes can be confusing
    // when trying to access something not available in the process
    int num_quadrature_points = std::stoi(this->const_map["QUAD_PTS"]); 
    this->quadrature_weights_values = getQuadratureWeightsAndValues(num_quadrature_points); 
 
    // Lets assign the p values per process
    // p is from 0 -> number_of_triangles
    // Lets split up p as eqully as possible to the processes
    // First lets declare a vector to store a portion of the p values
    // relevant to the process
    std::vector<int> sub_p_values;

    // Lets also declare an index
    // This will help divide the total p values into it's smaller
    // sub_p_values
    int index;

    // Now lets split up p
    // The variable index is used as a starting point for th values to be split
    // e.g. p is:
    // 0 1 2 3 4 5 6 7 8 9
    // and there are two processes
    // so process 0 will get from 0 -> 4
    // and process 1 will get from 5 -> 9
    // In the above example index == 0 for process 0 and index == 5 for process 1
    // One might be wondering why the following is executed in all processes and not just
    // root and then sent to the relevant proccesses.
    // Besides the memory requirement, MPI's main overhead comes from the sending and receiving
    // of messages. None of the processes can do anything until they receive the information so they
    // will be idle due to blocking. It is therefore faster for each process to calculate the 
    // relevant p_values they will have to use.
    // This is made possible by the fact the p is linear and sequential.
    // If p was, for example, a list of random values (0, 82, 3, 4, ...) then it would need
    // to be sent
    if(rank == 0)
    {
        // The index does not need used here since process 0 will always start at index == 0
        // The number of values allocated to the process is determined by numValuesMPI
        for(int i = 0; i < this->numValuesMPI(size, rank, this->triangles.size()); i++)
        {
            sub_p_values.push_back(i);
        }
    }
    else
    {   
        // It gets a bit more complex here since non root processes need to calculate
        // how many p_values were sent to the preceeding processes and set the index accordingly
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
    // It is important to remember to resize the vector
    // The space needs to be allocated for MPI to write the value into
    // Only root is going to receive the data so it only needs to be allocated there
    // TODO the resize nedds to only happen in root. Check to make sure 
    std::vector<int> proc_vector_size;
    if(rank == 0)
    {
        proc_vector_size.resize(size);
    }

    // Now lets get the size of sub_zmn for each process
    int z_mn_size = sub_zmn.size();

    // Now lets send them to the main process
    // This is important because to receive the actual data, MPI needs to be told
    // how many values it has to receive.
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

                    // This is the only difference between the serial implentation
                    // The indices of the partial Zmn value needs to be returned aswell
                    // so the main process knows where to put it
                    // Therefore, lets first push the indices to the vector
                    partial_zmn.push_back(this->triangles[p].getEdges()[r]);
                    partial_zmn.push_back(this->triangles[q].getEdges()[e]);

                    // Now lets calculate the partial Zmn value
                    std::complex<double> temp_zmn_value = 
                    this->edges[this->triangles[p].getEdges()[r]].getLength() *
                    (this->j * 
                        this->omega *
                        a_pq.getDot(rho_c) / 
                        std::complex<double>(2, 0) -
                        phi * 
                        phi_sign);

                    // It is important to remember that MPI does not natively support complex values
                    // Therefore, it is necessary to seperate them into their real and imaginary
                    // parts and push them individually
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
    // The full commentary can be found in the serial implementation
    // The file is mom_solver.cpp in the src/ directory

    std::vector<std::complex<double>> i_vector = this->calculateIpq(p, q);

    // Lets make this compatible with OpenMP
    // To do this, lets make the for loop independent
    // Lets first resize the vector to accomodate(sp?) the values
    std::vector<Node> a_and_phi_vector;
    a_and_phi_vector.resize(4);


    for(int i = 0; i < 3; i++)
    {
        Node a_pq_node_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(i_vector[1]);
        Node a_pq_node_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(i_vector[2]);
        Node a_pq_node_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(i_vector[3]);
        Node a_pq_node_4 = this->nodes[this->triangles[q].getVertices()[i]].getScalarMultiply(i_vector[0]);

        Node a_pq_node = a_pq_node_1.getAddNode(a_pq_node_2);
        a_pq_node = a_pq_node.getAddNode(a_pq_node_3);
        a_pq_node = a_pq_node.getSubtractComplexNode(a_pq_node_4);

        double a_pq_multiplier = std::stod(this->const_map["MU_0"]) / (4 * M_PI);                
        a_pq_node = a_pq_node.getScalarMultiply(a_pq_multiplier);
        a_and_phi_vector[i] = (a_pq_node);
    }

    std::complex<double> phi_pq_multiplier = std::complex<double>(1, 0) / 
                                            (this->j * std::complex<double>(2.0, 0) * 
                                            std::complex<double>(M_PI, 0) *
                                            this->omega *
                                            std::stod(this->const_map["EPS_0"])); 
    
    std::complex<double> phi_pq = i_vector[0] * phi_pq_multiplier;

    Node phi_pq_node(phi_pq, 0, 0);
    a_and_phi_vector[3] = phi_pq_node;    

    return  a_and_phi_vector;
}

std::vector<std::complex<double>> MoMSolverMPI::calculateIpq(int p, int q)
{
    // The full commentary can be found in the serial implementation
    // The file is mom_solver.cpp in the src/ directory

    std::vector<std::complex<double>> i_vector;


    if(p == 19823) // TODO change to p == q && SING == True
    {
        // TODO Add singularity treatment
        int x = 0;
    }
    else
    {
        std::complex<double> Ipq;        // RWG80 34a
        std::complex<double> Ipq_xi;     // RWG80 34b
        std::complex<double> Ipq_eta;    // RWG80 34c
        std::complex<double> Ipq_zeta;   // RWG80 34d

        // Lets make this independant for OpenMP
        std::vector<std::complex<double>> Ipq_vector(this->quadrature_weights_values.size());
        std::vector<std::complex<double>> Ipq_xi_vector(this->quadrature_weights_values.size());
        std::vector<std::complex<double>> Ipq_eta_vector(this->quadrature_weights_values.size());

        for(int i = 0; i < this->quadrature_weights_values.size(); i++)
        {
            Node xi_r_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(this->quadrature_weights_values[i][1]); 
            Node eta_r_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(this->quadrature_weights_values[i][2]); 
            Node zeta_r_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(this->quadrature_weights_values[i][3]);

            Node r_prime = xi_r_1.getAddNode(eta_r_2);
            r_prime = r_prime.getAddNode(zeta_r_3); 

            double R_p = this->triangles[p].getCentre().getDifferenceBetween(r_prime).getNorm();

            std::complex<double> greens_function = std::exp(std::complex<double>(-1.0, 0) * this->j * this->k * R_p) / R_p;

            Ipq_vector[i] = (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * greens_function);  
               
            Ipq_xi_vector[i] = (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                                    this->quadrature_weights_values[i][1] * 
                                    greens_function);

            Ipq_eta_vector[i] = (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                                    this->quadrature_weights_values[i][2] * 
                                    greens_function);  
        }

        for(int i = 0; i < this->quadrature_weights_values.size(); i++)
        {
            Ipq     += Ipq_vector[i];
            Ipq_xi  += Ipq_xi_vector[i];
            Ipq_eta += Ipq_eta_vector[i];
        }

        Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
        i_vector.push_back(Ipq);
        i_vector.push_back(Ipq_xi);
        i_vector.push_back(Ipq_eta);
        i_vector.push_back(Ipq_zeta);
    }
    return i_vector;
}




























