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
    this->frequency = std::stod(this->const_map["cppFreq"]);
    this->omega = 2 * M_PI * this->frequency;
    this->lambda = this->c / this->frequency;
    this->k = 2 * M_PI / this->lambda;
    std::complex<double> complex_constant(0, 1);
    this->j = complex_constant;

    // int size;
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void MoMSolverMPI::calculateVrhsInternally()
{
    // Lets calculate the Vrhs data internally
    // This wil just be for a nomally incident x-directed plane wave

    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(std::stoi(this->const_map["edge_feed"]) == 1)
    {
        if(rank == 0)
        {
            this->vrhs_internal.resize(this->edges.size());
            int edge_index = std::stoi(this->const_map["feed_edge"]) - 1;

            this->vrhs_internal[edge_index] = -1 * this->edges[edge_index].getLength() * 
                                               std::stod(this->const_map["EMag"]);
        }
    }
    else
    {
        

        std::vector<int> sub_edge_values;
        int index;

        if(rank == 0)
        {
            // The index does not need used here since process 0 will always start at index == 0
            // The number of values allocated to the process is determined by numValuesMPI
            for(int i = 0; i < this->numValuesMPI(size, rank, this->edges.size()); i++)
            {
                sub_edge_values.push_back(i);
            }
        }
        else
        {   
            // It gets a bit more complex here since non root processes need to calculate
            // how many edge_values were sent to the preceeding processes and set the index accordingly
            index = 0;
            for(int i = 0; i < rank; i++)
            {
                index += this->numValuesMPI(size, i, this->edges.size());
            }
            for(int i = index; i < this->numValuesMPI(size, rank, this->edges.size()) + index; i++)
            {
                sub_edge_values.push_back(i);
            }
        }

        std::vector<double> sub_vrhs;

        std::complex<double> tmp_vrhs_value;

        Node E_plus;
        Node E_minus;
        int m;

        double theta = std::stod(this->const_map["theta_0"]) * std::stod(this->const_map["DEG2RAD"]);
        double phi = std::stod(this->const_map["phi_0"]) * std::stod(this->const_map["DEG2RAD"]);
        double propagation_direction = std::stoi(this->const_map["prop_direction"]);
        double e_mag = std::stod(this->const_map["EMag"]);

        double e_x;
        double e_y;
        double e_z;    

        if(propagation_direction == 0)
        {
            e_x = e_mag * std::cos(phi) * std::cos(theta);
            e_y = e_mag * std::sin(phi) * std::cos(theta);
            e_z = e_mag * -std::sin(theta);

            if(theta == M_PI)
            {
                e_z = 0.0;
            }
        }
        else if(propagation_direction == 1)
        {

        }

        Node e(e_x, e_y, e_z);

        Node k(this->k * std::sin(theta) * std::cos(phi),
                this->k * std::sin(theta) * std::sin(phi),
                this->k * cos(theta)); 

        Node e_plus;
        Node e_minus;

        for(int i = 0; i < sub_edge_values.size(); i++)
        {
            m = sub_edge_values[i];
            int triangle_plus = this->edges[m].getPlusTriangleIndex();
            int triangle_minus = this->edges[m].getMinusTriangleIndex();

            e_plus = e.getScalarMultiply(std::exp(this->j * k.getDotNoComplex(this->triangles[triangle_plus].getCentre())));
            e_minus = e.getScalarMultiply(std::exp(this->j * k.getDotNoComplex(this->triangles[triangle_minus].getCentre())));
            

            tmp_vrhs_value = 0.5 * e_plus.getDot(this->edges[m].getRhoCPlus()) + 
                             0.5 * e_minus.getDot(this->edges[m].getRhoCMinus());
            tmp_vrhs_value *= this->edges[m].getLength();

            sub_vrhs.push_back(tmp_vrhs_value.real()); 
            sub_vrhs.push_back(tmp_vrhs_value.imag()); 
        }

        std::vector<int> receive_count;
        std::vector<int> displs;
        std::vector<double> tmp_vrhs;
        
        if(rank == 0)
        {
            for(int i = 0; i < size; i++)
            {
                receive_count.push_back(this->numValuesMPI(size, i, this->edges.size()) * 2);
            }

            displs.resize(size);
            displs[0] = 0;
            for(int i = 1; i < size; i++)
            {
                displs[i] = displs[i - 1] + receive_count[i - 1];
            }
            tmp_vrhs.resize(this->edges.size() * 2);
        }

        MPI_Gatherv(&sub_vrhs[0], sub_vrhs.size(), MPI_DOUBLE, &tmp_vrhs[0], 
                        &receive_count[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(rank == 0)
        {
            for(int i = 0; i < tmp_vrhs.size(); i+=2)
            {
                this->vrhs_internal.push_back(std::complex<double>(tmp_vrhs[i], tmp_vrhs[i + 1]));
            }
        }
    }
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
    this->num_points_j = 3;
    this->num_points_i = 4;
    this->quadrature_weights_values_j = getGaussLegendreQuadratureWeightsAndValues(num_points_j); 
    this->quadrature_weights_values_i = getGaussLegendreQuadratureWeightsAndValues(num_points_i);  
    
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
    //std::vector<double> sub_zmn = this->workMPI(sub_p_values);
    std::vector<double> sub_zmn = this->workMPIMP(sub_p_values);
    
    // Lets gather all the vector sizes
    // Lets first create a vector to store the sizes 
    // It is important to remember to resize the vector
    // The space needs to be allocated for MPI to write the value into
    std::vector<int> proc_vector_size;
    proc_vector_size.resize(size);

    // Now lets get the size of sub_zmn for each process
    int z_mn_size = sub_zmn.size();

    // Now lets send them to the main process
    // This is important because to receive the actual data, MPI needs to be told
    // how many values it has to receive.
    MPI_Gather(&z_mn_size, 1, MPI_INT, &proc_vector_size[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now that the size of the data is available to the root process
    // Lets gather all the data
    std::vector<double> all_zmn_data;

    // Now lets send the data from the processes to the root process
    // A very easy way to do this is MPI_Gatherv but it runs into
    // memory issues for large problems(around zmn = 5000x5000)
    // The root process does not need to send it's calculated data
    // to itself, so lets send from all the other processes
    if(rank != 0)
    {
        MPI_Send(&sub_zmn[0], sub_zmn.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // Now at the root process
    // Lets receive the data
    if(rank == 0)
    {
        // First lets allocate some memory for Zmn
        // Remember that this must only be done in the root process
        this->zmn = new std::complex<double>[this->edges.size() * this->edges.size()];

        // Lets initialize all the entries in this->zmn to 0
        // Since this->zmn is dynamically allocated a different
        // indexing method will be used.
        // For indices i and j, the equivalent of zmn[i][j] is,
        // zmn[i + zmn_length * j]
        // A vector of vectos is not used due to SCALAPACK needing
        // contiguos(sp?) data both in the row and column directions
        // A vector of vectors is only contiguous in its rows
        for(int i = 0; i < this->edges.size(); i++)
        {
            for(int j = 0; j < this->edges.size(); j++)
            {
                this->zmn[j * this->edges.size() + i] = std::complex<double>(0, 0); 
            }
        }
        
        // Lets receive the data from the sub processes
        // First check if there are more than one process
        // If there is only one process, then nothing needs to
        // received since nothing was sent
        if(size > 1)
        {
            // Lets loop over the non-root processes
            for(int i = 1; i < size; i++)
            {
                // Lets clear the data from the previous process and allocate 
                // the necessary space for the received data
                all_zmn_data.clear();
                all_zmn_data.resize(proc_vector_size[i]);

                // Lets receive the data
                MPI_Recv(&all_zmn_data[0], proc_vector_size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Now lets start filling in Zmn
                // Remember that the data received is of the form
                // row_index, col_index, real_part, imag_part ...
                int index;
                for(int j = 0; j < all_zmn_data.size() / 4; j++)
                {
                    index = j * 4;
                    std::complex<double> temp(all_zmn_data[index + 2], all_zmn_data[index + 3]);
                    this->zmn[(int)all_zmn_data[index] * this->edges.size() + (int)all_zmn_data[index+1]] += temp;
                }
            }
        }

        // Lets add the Zmn contributions from those calculated by the root process
        // Remember that the roor process cannot send to itself, so this is a necessary step
        int index;
        for(int j = 0; j < sub_zmn.size() / 4; j++)
        {
            index = j * 4;
            std::complex<double> temp(sub_zmn[index + 2], sub_zmn[index + 3]);
            this->zmn[(int)sub_zmn[index] * this->edges.size() + (int)sub_zmn[index+1]] += temp;
        }
    }
}

int MoMSolverMPI::numValuesMPI(int num_procs, int rank, int data_length)
{
    // Efficiently distribute work between processes
    // https://stackoverflow.com/questions/5657158/how-to-distribute-a-vector-of-n-elements-across-p-processors
    return (data_length + rank) / num_procs;
}

std::vector<double> MoMSolverMPI::workMPIMP(std::vector<int> p) // RENAME
{
    // See MoMSolverMPI::calculateZmnByFace() for full commentary

    // Lets create a vector to store the processes section of Zmn
    // All the threads will have access to this
    std::vector<double> zmn;

    // Lets declare an OpenMP parallel region
    // This will allow the code to be split between threads
    // It is important to remember to create a parallel region
    // before calling other OpenMP parallel commands
    // It is also important to remember that creating a parallel
    // region is quite expensive so it is better to create one 
    // large one as opposed to many smaller ones
    // This is the reasin that OpenMP is not used to assist with
    // A_pq and Phi_pq calculations 
    #pragma omp parallel
    {
        // Lets create a vector to store a small portion of 
        // the processes total selection of Zmn
        // This vector is unique to each thread
        std::vector<double> partial_zmn;

        // Lets make the two top for loops parallel
        #pragma omp for nowait collapse(2)
        for(int i = 0; i < p.size(); i++)
        {
            for(int q = 0; q < this->triangles.size(); q++)
            {
                std::vector<Node> a_and_phi = this->calculateAAndPhi(p[i], q);

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

                    for(int r = 0; r < this->triangles[p[i]].getEdges().size(); r++)
                    {
                        Node rho_c;
                        double phi_sign;

                        if(this->edges[this->triangles[p[i]].getEdges()[r]].getMinusTriangleIndex() == p[i])
                        {
                            rho_c = this->edges[this->triangles[p[i]].getEdges()[r]].getRhoCMinus();
                            phi_sign = -1;
                        }
                        else
                        {
                            rho_c = this->edges[this->triangles[p[i]].getEdges()[r]].getRhoCPlus();
                            phi_sign = 1;
                        }

                        // This is the only difference between the serial implentation
                        // The indices of the partial Zmn value needs to be returned aswell
                        // so the main process knows where to put it
                        // Therefore, lets first push the indices to the vector
                        partial_zmn.push_back(this->triangles[p[i]].getEdges()[r]);
                        partial_zmn.push_back(this->triangles[q].getEdges()[e]);

                        // Now lets calculate the partial Zmn value
                        std::complex<double> temp_zmn_value = 
                        this->edges[this->triangles[p[i]].getEdges()[r]].getLength() *
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

        // Lets now gather the smaller portions of Zmn from each thread
        // and place them in the common Zmn vector
        // Rememeber that this is not the full Zmn, but just the portion
        // allocated to the specific process
        #pragma omp critical
        zmn.insert(zmn.end(), partial_zmn.begin(), partial_zmn.end()); 

    }
    return zmn;
}

void MoMSolverMPI::calculateJMatrixSCALAPACK()
{   
    // Lets calculate Ilhs(current) using SCALAPACK
    // First an LU decompostion of Zmn and then solving
    // for Ilhs using it

    // Lets get the rank of the current process and  
    // the number of processes
    // This can also be done using vanilla MPI
    // functions, but to be consistent in the function,
    // lets use BLACS
    // BLACS is a layer on top of vanilla MPI and provides
    // the needed framework for SCALAPACK to function
    // It is important to differentiate between
    // BLAS and BLACS
    // BLAS     -> Low level linear algebra operations
    //          -> Used by SCALAPACK to do the major matix
    //              operations
    //          -> Won't be mentioned again
    //
    // BLACS    -> Layer over MPI
    //          -> Framework for the parallel piece of 
    //              SCALAPACK
    //          -> Used to do all the sending and receiving
    //              and any MPI stuff
    int rank;
    int size;
    Cblacs_pinfo(&rank, &size);

    // Lets get the BLACS context
    // This is similar to tags in vanilla MPI
    // There can be many simultaneous BLACS contexts
    int context;
    Cblacs_get(0, 0, &context);

    // Lets define some constants
    // SCALAPACK is written in Fortran, which
    // only supports pass by reference
    int zero = 0;
    int one = 1;
    int matrix_size = this->edges.size(); 

    // A quick overview of SCALAPACK usage
    // SCALAPACK uses LAPACK routines to do linear algebra operations
    // The LAPACK operations are used on small subsets of data to eventually
    // find the solution to the overall data
    // These smaller subsets can be done simultaneously
    // MPI facilitates these simultaneous operations.
    // This is where BLACS comes in
    // It proveds higher level MPI functionallity for ease of use
    //
    // Lets discuss the process
    // First lets create a process grid
    // This determines how the data will by distributed among the processes
    // An example for four processes in a 2x2 square grid:
    //    0 1
    //    - -
    // 0 |0|1|  -> The numbers in the blocks are the processor ranks
    // 1 |2|3|  -> The 0's and 1's outside are the grid co-ordinated
    //    - -   
    // The shape of the grid affects the placement of the data so
    // different grid shapes influence the speed of the calculations
    // There is no hard and fast rule on the optimal grid size and shape,
    // but there is merit to having the grid as square as possible
    // According to the SCALAPACK website, a grid with more columns
    // than rows is better for LU decompositon
    // Thus, a good grid to aim for is rows x (rows + 1) 
    //
    // Lets now create the process grid
    std::vector<int> proc_grid = this->getProcessGrid(size);
    int proc_rows = proc_grid[0];
    int proc_cols = proc_grid[1];
    Cblacs_gridinit(&context, "Row-major", proc_rows, proc_cols); 

    // Lets get the amount of rows and columns in the grid
    // Lets also get the position in the grid of the current process
    int procrows;
    int proccols;
    int my_row;
    int my_col;
    Cblacs_gridinfo(context, &procrows, &proccols, &my_row, &my_col);

    // Continuing the above overview,
    // Lets discuss how the data is distributed on the process grid
    // SCALAPACK uses a block cyclic distribution pattern
    // The data is divided into blocks and cyclically distribited to each process
    // See http://www.netlib.org/scalapack/slug/node75.html for details 
    // Lets define the block size
    int block_size = 64;
    if(block_size > this->edges.size())
    {
        block_size = this->edges.size() / 4;
    }

    // Lets get the number or local rows and columns
    // Zmn
    int local_rows = numroc_(&matrix_size, &block_size, &my_row, &zero, &proc_rows);
    int local_cols = numroc_(&matrix_size, &block_size, &my_col, &zero, &proc_cols);

    // Vrhs
    int v_local_rows = numroc_(&matrix_size, &block_size, &my_row, &zero, &procrows);
    int v_local_cols = numroc_(&one, &block_size, &my_col, &zero, &proccols);

    // Create the local part of Zmn per process
    std::complex<double> *A_local;
    A_local = new std::complex<double>[local_rows * local_cols];

    // Lets send the global Zmn to the local Zmn's in a block cyclic manner
    int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
    for (int r = 0; r < matrix_size; r += block_size, sendr=(sendr+1)%proc_rows) {
        sendc = 0;
        // Number of rows to be sent
        // Is this the last row block?
        int nr = block_size;
        if (matrix_size-r < block_size)
            nr = matrix_size-r;
 
        for (int c = 0; c < matrix_size; c += block_size, sendc=(sendc+1)%proccols) {
            // Number of cols to be sent
            // Is this the last col block?
            int nc = block_size;
            if (matrix_size-c < block_size)
                nc = matrix_size-c;
 
            if (rank == 0) {
                // Send a nr-by-nc submatrix to process (sendr, sendc)
                Czgesd2d(context, nr, nc, this->zmn+matrix_size*c+r, matrix_size, sendr, sendc);
            }
 
            if (my_row == sendr && my_col == sendc) {
                // Receive the same data
                // The leading dimension of the local matrix is nrows!
                Czgerv2d(context, nr, nc, A_local+local_rows*recvc+recvr, local_rows, 0, 0);
                recvc = (recvc+nc)%local_cols;
            }
 
        }
 
        if (my_row == sendr)
            recvr = (recvr+nr)%local_rows;
    }

    std::complex<double> B_local[v_local_rows];

    // Now lets distribute Vrhs
    int send_row = 0;
    int send_col = 0;
    int rec_row = 0;
    int rec_col = 0;
    int v_local_index = 0;

    for(int i = 0; i < matrix_size; i += block_size)
    {
        // Check if these are the last rows
        
        int num_rows = block_size;
        if(matrix_size - i < block_size)
        {
            num_rows = matrix_size - i;
        }

        // Switch back to sending to (0,0)
        if(send_row >= proc_rows)
        {
            send_row = 0;
        }

        if(rank == 0)
        {
            // SEND
            Czgesd2d(context, num_rows, 1, &this->vrhs_internal[i], matrix_size, send_row, send_col);

        }

        if(my_row == send_row && my_col == send_col)
        {
            // RECEIVE
            Czgerv2d(context, num_rows, 1, &B_local[v_local_index], local_rows, 0, 0);
            v_local_index += block_size;
        }
        send_row++;
    }

    // Lets solve the equation 
    // Lets first get the descriptions
    // First for Zmn
    int info = 0;
    int lda = std::max(1, local_rows);
    int desc[9]; 

    descinit_(desc, &matrix_size, &matrix_size, &block_size, &block_size, &zero, &zero, &context, &lda, &info);

    int v_lda = std::max(1, v_local_rows);
    int v_desc[9];

    descinit_(v_desc, &matrix_size, &one, &block_size, &block_size, &zero, &zero, &context, &v_lda, &info);

    // Now the LU decomposition
    int local_pivot[local_rows * block_size];
    pzgetrf_(&matrix_size, &matrix_size, A_local, &one, &one, desc, local_pivot, &info); 

    // Lets solve for I
    pzgetrs_("T", &matrix_size, &one, A_local, &one, &one, desc, local_pivot, B_local, &one, &one, v_desc, &info);

    // Lets gather the solved vector
    // get from each process in block cyclic manner
    send_row = 0;
    v_local_index = 0;

    if(rank == 0)
    {
        this->ilhs.resize(matrix_size);
    }

    for(int i = 0; i < matrix_size; i += block_size)
    {
        // Check if these are the last rows
        
        int num_rows = block_size;
        if(matrix_size - i < block_size)
        {
            num_rows = matrix_size - i;
        }

        // Switch back to sending to (0,0)
        if(send_row >= proc_rows)
        {
            send_row = 0;
        }

        if(my_row == send_row && my_col == send_col)
        {
            // SEND
            Czgesd2d(context, num_rows, 1, &B_local[v_local_index], matrix_size, 0, 0);
            v_local_index += block_size;
        }

        if(rank == 0)
        {
            // RECEIVE
            Czgerv2d(context, num_rows, 1, &this->ilhs[i], local_rows, send_row, send_col);
        }
        send_row++;
    }

    // Lets print the solved vector
    // if(rank == 0)
    // {
    //     for(int i = 0; i < matrix_size; i++)
    //     {
    //         std::cout << this->ilhs[i] << std::endl;
    //     }
    // }
}

std::vector<Node> MoMSolverMPI::calculateAAndPhi(int p, int q)
{
    // The full commentary can be found in the serial implementation
    // The file is mom_solver.cpp in the src/ directory

    std::vector<std::complex<double>> i_vector = this->calculateIpq(p, q);

    std::vector<Node> a_and_phi_vector;

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
        a_and_phi_vector.push_back(a_pq_node);
    }

    std::complex<double> phi_pq_multiplier = std::complex<double>(1, 0) / 
                                            (this->j * std::complex<double>(2.0, 0) * 
                                            std::complex<double>(M_PI, 0) *
                                            this->omega *
                                            std::stod(this->const_map["EPS_0"])); 
    
    std::complex<double> phi_pq = i_vector[0] * phi_pq_multiplier;

    Node phi_pq_node(phi_pq, 0, 0);
    a_and_phi_vector.push_back(phi_pq_node);    

    return  a_and_phi_vector;
}

std::vector<std::complex<double>> MoMSolverMPI::calculateIpq(int p, int q)
{
    // The full commentary can be found in the serial implementation
    // The file is mom_solver.cpp in the src/ directory

    std::vector<std::complex<double>> i_vector;

    bool sing = true;
    if(p == q && sing) // TODO: change to p == q && SING == True
    {
        i_vector = this->getIpqSING(p);
    }
    else
    {
        std::complex<double> Ipq;        // RWG80 34a
        std::complex<double> Ipq_xi;     // RWG80 34b
        std::complex<double> Ipq_eta;    // RWG80 34c
        std::complex<double> Ipq_zeta;   // RWG80 34d

        for(int i = 0; i < this->quadrature_weights_values.size(); i++)
        {
            Node xi_r_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(this->quadrature_weights_values[i][1]); 
            Node eta_r_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(this->quadrature_weights_values[i][2]); 
            Node zeta_r_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(this->quadrature_weights_values[i][3]);

            Node r_prime = xi_r_1.getAddNode(eta_r_2);
            r_prime = r_prime.getAddNode(zeta_r_3); 

            double R_p = this->triangles[p].getCentre().getDifferenceBetween(r_prime).getNorm();

            std::complex<double> greens_function = std::exp(std::complex<double>(-1.0, 0) * this->j * this->k * R_p) / R_p;

            Ipq = Ipq + (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * greens_function);  
               
            Ipq_xi = Ipq_xi + (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                                    this->quadrature_weights_values[i][1] * 
                                    greens_function);

            Ipq_eta = Ipq_eta+ (std::complex<double>(0.5,0) * this->quadrature_weights_values[i][0] * 
                                    this->quadrature_weights_values[i][2] * 
                                    greens_function);  
        }

        Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
        i_vector.push_back(Ipq);
        i_vector.push_back(Ipq_xi);
        i_vector.push_back(Ipq_eta);
        i_vector.push_back(Ipq_zeta);
    }
    return i_vector;
}

std::vector<int> MoMSolverMPI::getProcessGrid(int num_procs)
{
    // THIS IS FOR ONLY UP TIL 12 PROCS
    // TODO: optimize by hand
    int proc_rows;
    int proc_cols;

    if(num_procs%2 == 0)
    {
        proc_cols = std::sqrt(num_procs);
        proc_rows = num_procs / proc_cols;

        if((proc_rows * proc_cols) != num_procs)
        {
            proc_cols--;
            proc_rows += proc_cols;
        }
    }
    else
    {
        proc_rows = num_procs;
        proc_cols = 1;
    }
      
    std::vector<int> return_vector;
    return_vector.push_back(proc_rows);
    return_vector.push_back(proc_cols);
    return return_vector;
}


std::vector<std::complex<double>> MoMSolverMPI::getIlhs()
{
    return this->ilhs;
}

std::vector<std::complex<double>> MoMSolverMPI::getIpqSING(int p)
{
    // Lets declare and initialise the four I's to zero
    std::complex<double> Ipq = std::complex<double>(0, 0);        // RWG80 34a
    std::complex<double> Ipq_xi = std::complex<double>(0, 0);     // RWG80 34b
    std::complex<double> Ipq_eta = std::complex<double>(0, 0);    // RWG80 34c
    std::complex<double> Ipq_zeta = std::complex<double>(0, 0);   // RWG80 34d
    std::complex<double> Ipq_zt = std::complex<double>(0, 0);   // RWG80 34d

    // Lets declare the vectors that will store the I values for the three sub triangles
    // Don't forget to initialize them to zero
    // std::vector<std::complex<double>> sub_ipq(3, std::complex<double>(0, 0));
    // std::vector<std::complex<double>> sub_ipq_xi(3, std::complex<double>(0, 0));
    // std::vector<std::complex<double>> sub_ipq_eta(3, std::complex<double>(0, 0));

    // Lets loop over the three as of yet declared sub triangles
    for(int sub_triangle_index = 0; sub_triangle_index < 3; sub_triangle_index++)
    {
        // Lets declare the sub triangles now
        // The point of division is the centre of the triangle
        // The centre point will therefore be common
        Node r1_prime = this->triangles[p].getCentre();
        Node r2_prime;
        Node r3_prime;

        // Lets declare the 3x3 matrix T noted in Equation 17
        std::vector<std::array<double, 3>> T;

        // Lets then allocate the other points to the sub triangle
        // This corresponds to Equation 4 in KW05
        // Lets also fill in T
        // Remember that the observation point is the centre of the 
        // observed triangle.
        // This means that xi values needed in T are the simplex
        // coordinates of the triangles centre.
        // The simplex centre thus equals [1/3 1/3 1/3]
        switch(sub_triangle_index)
        {
            case 0:
            {
                r2_prime = this->nodes[this->triangles[p].getVertex2()]; 
                r3_prime = this->nodes[this->triangles[p].getVertex3()];
                  
                std::array<double, 3> row_1 = {1.0 / 3.0, 0, 0}; 
                std::array<double, 3> row_2 = {1.0 / 3.0, 1, 0}; 
                std::array<double, 3> row_3 = {1.0 / 3.0, 0, 1};
                    
                T.push_back(row_1); 
                T.push_back(row_2); 
                T.push_back(row_3); 
                break;
            }
            case 1:
            {
                r2_prime = this->nodes[this->triangles[p].getVertex3()];
                r3_prime = this->nodes[this->triangles[p].getVertex1()];

                std::array<double, 3> row_1 = {1.0 / 3.0, 0, 1}; 
                std::array<double, 3> row_2 = {1.0 / 3.0, 0, 0}; 
                std::array<double, 3> row_3 = {1.0 / 3.0, 1, 0};
                    
                T.push_back(row_1); 
                T.push_back(row_2); 
                T.push_back(row_3); 
                break;
            }

            case 2:
            {
                r2_prime = this->nodes[this->triangles[p].getVertex1()];
                r3_prime = this->nodes[this->triangles[p].getVertex2()];

                std::array<double, 3> row_1 = {1.0 / 3.0, 1, 0}; 
                std::array<double, 3> row_2 = {1.0 / 3.0, 0, 1}; 
                std::array<double, 3> row_3 = {1.0 / 3.0, 0, 0};
                    
                T.push_back(row_1); 
                T.push_back(row_2); 
                T.push_back(row_3); 
                break;
            }
        }

        // Lets get some of the geometric quantities needed
        // Lets get l1, l2 and l3 prime
        // This corresponds to Equation 4 in KW05
        Node l1_prime = r3_prime.getSubtractComplexNode(r2_prime);
        Node l2_prime = r1_prime.getSubtractComplexNode(r3_prime);
        Node l3_prime = r2_prime.getSubtractComplexNode(r1_prime);
        Node l1_prime_hat = l1_prime.getHat(); 

        // Lets get n_hat'
        Node l1_cross_l2 = l1_prime.getCrossProduct(l2_prime);
        Node n_hat_prime = l1_cross_l2.getScalarDivide(l1_cross_l2.getNorm());

        // Lets get A'
        double a_prime = 0.5 * n_hat_prime.getDotNoComplex(l1_cross_l2); 

        // Finally, lets get h_1'
        Node h1_prime = l1_prime.getCrossProduct(n_hat_prime);
        h1_prime = h1_prime.getScalarMultiply(2 * a_prime / std::pow(l1_prime.getNorm(), 2)); 
        Node h1_prime_hat = h1_prime.getHat();

        // // TEST
        // std::cout << "Sub => " << sub_triangle_index << std::endl;
        // std::cout << "V1 => ";
        // this->nodes[this->triangles[p].getVertex1()].printNode();
        // std::cout << "V2 => ";
        // this->nodes[this->triangles[p].getVertex2()].printNode();
        // std::cout << "V3 => ";
        // this->nodes[this->triangles[p].getVertex3()].printNode();
        // std::cout << "Centre => ";
        // this->triangles[p].getCentre().printNode();
        // std::cout << "R1p => ";
        // r1_prime.printNode();
        // std::cout << "R2p => ";
        // r2_prime.printNode();
        // std::cout << "R3p => ";
        // r3_prime.printNode();
        // std::cout << "l1p => ";
        // l1_prime.printNode();
        // std::cout << "l2p => ";
        // l2_prime.printNode();
        // std::cout << "l3p => ";
        // l3_prime.printNode();
        // std::cout << "l1pH => ";
        // l1_prime_hat.printNode();
        // std::cout << "LXL => ";
        // l1_cross_l2.printNode();
        // std::cout << "n_p => ";
        // n_hat_prime.printNode();
        // std::cout << "A => " << a_prime << std::endl;
        // std::cout << "h1p => ";
        // h1_prime.printNode();
        // std::cout << std::endl;
        // std::cout << std::endl;
        
        // Now lets loop over num_points_j
        // This corresponds to Equation 10 in KW05
        for(int j = 0; j < this->num_points_j; j++)
        {
            // Lets get u_l and u_v
            // This corresponds to Equation 9 in KW05
            // To do this we need x_l     > Equation 6
            //                    x_u     > Equation 7
            //                    y'      > Equation 13
            // from KW05
            // Lets first get x_l and x_v
            double x_l = h1_prime_hat.getCrossProduct(l2_prime).getDotNoComplex(n_hat_prime) * 
                         (1 - this->quadrature_weights_values_j[j][1]); 
            double x_u = h1_prime_hat.getCrossProduct(l2_prime).getDotNoComplex(n_hat_prime.getScalarMultiply((double)-1)) * 
                         (1 - this->quadrature_weights_values_j[j][1]);

            // Now lets get y'
            double y_prime = h1_prime.getNorm() * (1 - this->quadrature_weights_values_j[j][1]);

            // Finally, lets get u_l and u_u
            // Remember that the point of observation is the centre of the triangle
            // This means that variable z is equal to zero as the observation point
            // is on the triangle and not elevated
            double u_l = std::asinh(x_l / std::sqrt(std::pow(y_prime, 2)));             
            double u_u = std::asinh(x_u / std::sqrt(std::pow(y_prime, 2)));             

            // if(sub_triangle_index == 0 && j == 0)
            // {
            //     std::cout << "xl => " << x_l << std::endl;
            //     std::cout << "xu => " << x_u << std::endl;
            //     std::cout << "yp => " << y_prime << std::endl;
            //     std::cout << "ul => " << u_l << std::endl;
            //     std::cout << "uu => " << u_u << std::endl;
            // }
            // Now lets loop over num_points_i
            // This corresponds to Equation 10 in KW05
            for(int i = 0; i < this->num_points_i; i++)
            {
                // Continuing with Equation 10, lets get R^(i,j)
                // This corresponds to Equation 18 in KW05

                // Lets first get xi_1^k, xi_2^k and xi_3^k
                // as shown in Equation 17 in KW05

                // We already have T, so all thats left is to 
                // get xi_1'^(i,j), xi_2'^(i,j) and xi_3'^(i,j)

                // Lets start with xi_1'^(i,j) as shown in Equation 13 in KW05
                // It is simply equal to the Gauss-Legendre quadrature point
                double xi1_ij = this->quadrature_weights_values_j[j][1];

                // Now lets get xi_3'^(i,j) as shown in Equation 16 in KW05
                // First, lets get u^ij as shown in Equation 14 in KW05
                double u_ij = u_l * (1 - this->quadrature_weights_values_i[i][1]) +
                              u_u * this->quadrature_weights_values_i[i][1];

                // Then lets get x'^(i,j) as shown in Equation 15 in KW05
                // Remember that since the point of observation of the observed
                // triangle is the triangle's centre, z is equal to zero
                double x_prime_ij = std::sqrt(std::pow(y_prime, 2)) * std::sinh(u_ij);

                // Since all the pieces are available, lets get xi_3'^(i,j)
                // Lets break Equation 13 up to have an easier time
                Node sub_numerator = h1_prime_hat.getScalarMultiply(y_prime);
                sub_numerator = sub_numerator.getSubtractComplexNode(l1_prime_hat.getScalarMultiply(x_prime_ij));
                sub_numerator = l3_prime.getCrossProduct(sub_numerator);
                double xi3_ij = n_hat_prime.getDotNoComplex(sub_numerator) / (2 * a_prime);

                // Lets now get the final piece, xi_2^ij shown in Equation 16 in KW05
                double xi2_ij = 1 - xi3_ij - xi1_ij;

                // Now lets continue with Equation 17
                // Lets first create a vector of the three xi_ij's
                std::vector<double> xi_ij_vector;
                xi_ij_vector.push_back(xi1_ij);
                xi_ij_vector.push_back(xi2_ij);
                xi_ij_vector.push_back(xi3_ij);

                // Now lets get the three xi_k's
                std::vector<double> xi_k_vector;
                xi_k_vector = this->getMatrixXVector(T, xi_ij_vector);

                // Now lets finally get R^(i,j) from Equation 18 in KW05
                // Remember that r1, r2 and r3 are the original triangle vertices
                // r is the observation point which is the triangle centre
                // Lets break it up a bit for clarity
                Node R_term_1 = this->triangles[p].getCentre();
                Node R_term_2 = this->nodes[this->triangles[p].getVertex1()].getScalarMultiply(xi_k_vector[0]);
                Node R_term_3 = this->nodes[this->triangles[p].getVertex2()].getScalarMultiply(xi_k_vector[1]);
                Node R_term_4 = this->nodes[this->triangles[p].getVertex3()].getScalarMultiply(xi_k_vector[2]);

                Node R_inner = R_term_1.getSubtractComplexNode(R_term_2);
                R_inner = R_inner.getSubtractComplexNode(R_term_3);
                R_inner = R_inner.getSubtractComplexNode(R_term_4);

                double R_ij = R_inner.getNorm();

                // if(sub_triangle_index == 0 && j == 0 && i == 0)
                // {
                //     std::cout << "xi1ij => " << xi1_ij << std::endl;
                //     std::cout << "xi2ij => " << xi2_ij << std::endl;
                //     std::cout << "xi3ij => " << xi3_ij << std::endl;
                //     std::cout << "R_ij => " << R_ij << std::endl;
                //     std::cout << "xik1 => " << xi_k_vector[0] << std::endl;
                //     std::cout << "xik2 => " << xi_k_vector[1] << std::endl;
                //     std::cout << "xik3 => " << xi_k_vector[2] << std::endl;
                //     std::cout << "uij => " << u_ij << std::endl;
                //     std::cout << "xprimeji => " << x_prime_ij << std::endl;

                //     std::cout << std::endl;
                //     for(int pp = 0; pp  < 3;pp++)
                //     {
                //         std::cout << T[pp][0] << " "
                //                   << T[pp][1] << " "
                //                   << T[pp][2] << std::endl;
                //     }
                // }


                // Now, lets finish off Equation 10
                // Lets calculate the three of the four I's
                Ipq += 
                       this->quadrature_weights_values_j[j][0] *
                       this->quadrature_weights_values_i[i][0] *
                       h1_prime.getNorm()                      *
                       (u_u - u_l)                             *
                       (xi_k_vector[0] + xi_k_vector[1] + xi_k_vector[2]) *
                       std::exp(std::complex<double>(-1,0) * this->j * this->k * R_ij);                

                Ipq_xi += 
                          this->quadrature_weights_values_j[j][0] *
                          this->quadrature_weights_values_i[i][0] *
                          h1_prime.getNorm()                      *
                          (u_u - u_l)                             *
                          xi_k_vector[0]                          *
                          // xi1_ij                                  *
                          std::exp(std::complex<double>(-1,0) * this->j * this->k * R_ij);

                Ipq_eta += 
                           this->quadrature_weights_values_j[j][0] *
                           this->quadrature_weights_values_i[i][0] *
                           h1_prime.getNorm()                      *
                           (u_u - u_l)                             *
                           xi_k_vector[1]                          *
                           // xi2_ij                                  *
                           std::exp(std::complex<double>(-1,0) * this->j * this->k * R_ij);
                
                Ipq_zt += 
                           this->quadrature_weights_values_j[j][0] *
                           this->quadrature_weights_values_i[i][0] *
                           h1_prime.getNorm()                      *
                           (u_u - u_l)                             *
                           xi_k_vector[2]                          *
                           // xi3_ij                          *
                           std::exp(std::complex<double>(-1,0) * this->j * this->k * R_ij);
            }
        }
    }
    // Lets get the final I
    // Ipq_zeta = Ipq - Ipq_xi - Ipq_zeta;
    Ipq_zeta = Ipq_zt;

    // std::cout << Ipq << std::endl;
    // std::cout << Ipq_xi << std::endl;
    // std::cout << Ipq_eta << std::endl;
    // std::cout << Ipq_zeta << std::endl;
    // std::cout << Ipq_zt << std::endl;

    // Lets divide by the 2 * triangle area as noted in Equation 31 in RWG80
    Ipq = Ipq / (2.0 * this->triangles[p].getArea());
    Ipq_xi = Ipq_xi / (2.0 * this->triangles[p].getArea());
    Ipq_eta = Ipq_eta / (2.0 * this->triangles[p].getArea());
    Ipq_zeta = Ipq_zeta / (2.0 * this->triangles[p].getArea());

    // Now lets return the I values
    std::vector<std::complex<double>> i_vector;
    i_vector.push_back(Ipq);
    i_vector.push_back(Ipq_xi);
    i_vector.push_back(Ipq_eta);
    i_vector.push_back(Ipq_zeta);

    return i_vector;
}

std::vector<double> MoMSolverMPI::getMatrixXVector(std::vector<std::array<double, 3>> matrix,
                                                std::vector<double> vec_tor)
{
    // Lets get the multiplication of a 3x3 matrix by a 3x1 vector
    std::vector<double> return_vector;
    return_vector.resize(3);

    for(int i = 0; i < 3; i++)
    {
        return_vector[i] = matrix[i][0] * vec_tor[0] +
                           matrix[i][1] * vec_tor[1] +
                           matrix[i][2] * vec_tor[2];
    }

    return return_vector;
}


void MoMSolverMPI::printPartialZMN()
{
    std::cout << this->edges.size() << std::endl;
    for(int i = 0; i < this->edges.size(); i++)
    {
        for(int j = 0; j < this->edges.size(); j++)
        {
            if(std::isnan(this->zmn[i + this->edges.size() * j].real()))
            {
                std::cout << "NAN " << i << " <=> " << j <<std::endl;
            }   
        }
        
    }
}


































