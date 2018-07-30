#include "mom_solver.h"

/**
 *   \file mom_solver.cpp
 *   \brief Solve the method of moments(MoM) of an electromagnetic problem
 *
 *  Detailed description
 *  This class will first solve for Zmn using a face pair approach. Then Vm = ZmnIn will be solved
 * 	for In using LU decomposition on Zmn. 
 *
 *  Author:  Tameez Ebrahim
 *  Created: 27 July 2018
 *
 */

MoMSolver::MoMSolver(std::vector<Node> nodes,
				     std::vector<Triangle> triangles,
				     std::vector<Edge> edges,
			         std::vector<float> vrhs,
			         std::map<std::string, std::string> const_map)
{
	this->nodes = nodes;
	this->triangles = triangles;
	this->edges = edges;
	this->vrhs = vrhs;
	this->const_map = const_map;

	// Lets define some constants
	// Check for what the frequency is
	this->k = (2 * 3.141592653589793); 
	//		  (std::stof(this->const_map["C0"]) / std::stof(this->const_map["freqData"]));
	std::complex<float>complex_constant(0, 1);
	this->j = complex_constant;
}

void MoMSolver::calculateZmnByFace()
{
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

	// Declare Zmn TODO change
	std::complex<float> z_mn[this->edges.size()][this->edges.size()] = {std::complex<float>(0,0)}; 

	for(int p = 0; p < this->triangles.size(); p++)
	{
		for(int q = 0; q < this->triangles.size(); q++)
		{
			// Lets calculate A and Phi
			// To do so, the four I values need to be calculated
			// TODO comment this
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

				std::complex<float> phi = a_and_phi[3].getXComplexCoord() * this->edges[this->triangles[q].getEdges()[e]].getLength();
				
				if(this->edges[this->triangles[q].getEdges()[e]].getMinusTriangleIndex() == q)
				{
					a_pq = a_pq.getScalarMultiply(-1.0);
				}
				else
				{
					phi = phi * std::complex<float>(1.0, 0);
				}
				
				for(int r = 0; r < this->triangles[p].getEdges().size(); r++)
				{
					Node rho_c;
					float phi_sign;

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
					
					z_mn[this->triangles[p].getEdges()[r]][this->triangles[q].getEdges()[e]] =
					z_mn[this->triangles[p].getEdges()[r]][this->triangles[q].getEdges()[e]] +
						this->edges[this->triangles[p].getEdges()[r]].getLength() *
							(this->j * 
							std::complex<float>(2, 0) * 
							std::complex<float>(3.141592653589793, 0) *
							std::stof(this->const_map["C0"]) *
							a_pq.getDot(rho_c) / 
							std::complex<float>(0, 2) -
							phi * 
							phi_sign);	
					//std::cout << z_mn[this->triangles[p].getEdges()[r]][this->triangles[q].getEdges()[e]] << std::endl;
				}				 	

			} 

		}
	}

	for(int i = 0; i < this->edges.size(); i++)
	{
		for (int j = 0; j < this->edges.size(); j++)
		{
			std::cout << z_mn[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

std::vector<Node> MoMSolver::calculateAAndPhi(int p, int q)
{
	// Lets calculate the magnetic vector potential and the scalar potential
	// The formulae are equations 32 and 33 in RWG80
	// First, lets get the four I's calculated in calculateIpq 
	// Remember, the vector is ordered in [Ipq Ipq_xi Ipq_eta Ipq_zeta]
	std::vector<std::complex<float>> i_vector = this->calculateIpq(p, q);

	// There are three values for the magnetic vector potential (Apq) and one for the scalar
	// potential(Phipq)
	// Remember, from calculateIpq, r1 -> vertex_1, r2 -> vertex_2 and r3 -> vertex_3
	// Lets first create the return vector
	std::vector<Node> a_and_phi_vector;

	// Lets loop over the three triangle vertices
	for(int i = 0; i < 3; i++)
	{
		// TODO comment this
		Node a_pq_node_1 = this->nodes[this->triangles[q].getVertex1()].getScalarMultiply(i_vector[1]);
		Node a_pq_node_2 = this->nodes[this->triangles[q].getVertex2()].getScalarMultiply(i_vector[2]);
		Node a_pq_node_3 = this->nodes[this->triangles[q].getVertex3()].getScalarMultiply(i_vector[3]);
		Node a_pq_node_4 = this->nodes[this->triangles[q].getVertices()[i]].getScalarMultiply(i_vector[0]);

		Node a_pq_node = a_pq_node_1.getAddComplexNode(a_pq_node_2);
		a_pq_node = a_pq_node.getAddComplexNode(a_pq_node_3);
		a_pq_node = a_pq_node.getAddComplexNode(a_pq_node_4);
		

		float a_pq = (std::stof(this->const_map["MU_0"]) / (4 * 3.141592653589793));			 	
		a_pq_node = a_pq_node.getScalarMultiply(a_pq);
		a_and_phi_vector.push_back(a_pq_node);

	}

	std::complex<float> phi_pq_multiplier = std::complex<float>(1, 0) / 
											(this->j * std::complex<float>(2.0, 0) * 
											std::complex<float>(3.141592653589793, 0) *
											std::complex<float>(2.0, 0) * 
											std::complex<float>(3.141592653589793) * 
											std::stof(this->const_map["C0"]) *
											std::stof(this->const_map["EPS_0"])); 
	

	std::complex<float> phi_pq = i_vector[0] * phi_pq_multiplier;
	Node phi_pq_node(phi_pq, 0, 0);
	a_and_phi_vector.push_back(phi_pq_node);	

	// TODO DELETE
	//std::cout << phi_pq << std::endl;
	return  a_and_phi_vector;
	
}

std::vector<std::complex<float>> MoMSolver::calculateIpq(int p, int q)
{
	// Lets calculate the 4 Ipq integrals(equations 34a-d in RWG80)
	// The first three(a-c) will be calculated using numerical quadrature
	// The final one(d) is calculated from the first three
	// p and q, the triangle indices are the only arguments needed
	// All the other parameters are specified in const_map
	// If p == q, then singularity treatment may be applied
	// This is defined in const_map whether to use the singularity treatment or not
	std::vector<std::complex<float>> i_vector;


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
		std::complex<float> Ipq;		// RWG80 34a
		std::complex<float> Ipq_xi;		// RWG80 34b
		std::complex<float> Ipq_eta;	// RWG80 34c
		std::complex<float> Ipq_zeta;	// RWG80 34d


		// Now lets loop over the quadrature points
		for(int i = 0; i < this->quadrature_weights_values.size(); i++)
		{
			// Lets first get r' as noted in equation 30 in RWG80
			// So, it is [xi * vertex_1, eta * vertex_2, zeta * vertex_3]  
			// TODO Check if r' needs to be normalized
			float r_prime_x = this->quadrature_weights_values[i][1] * this->nodes[this->triangles[q].getVertex1()].getXCoord() +
					          this->quadrature_weights_values[i][2] * this->nodes[this->triangles[q].getVertex2()].getXCoord() +
				              this->quadrature_weights_values[i][3] * this->nodes[this->triangles[q].getVertex3()].getXCoord();
			
			float r_prime_y = this->quadrature_weights_values[i][1] * this->nodes[this->triangles[q].getVertex1()].getYCoord() +
					          this->quadrature_weights_values[i][2] * this->nodes[this->triangles[q].getVertex2()].getYCoord() +
				              this->quadrature_weights_values[i][3] * this->nodes[this->triangles[q].getVertex3()].getYCoord();
			
			float r_prime_z = this->quadrature_weights_values[i][1] * this->nodes[this->triangles[q].getVertex1()].getZCoord() +
					          this->quadrature_weights_values[i][2] * this->nodes[this->triangles[q].getVertex2()].getZCoord() +
				              this->quadrature_weights_values[i][3] * this->nodes[this->triangles[q].getVertex3()].getZCoord();
			
			Node r_prime(r_prime_x, r_prime_y, r_prime_z);


			// Now, lets get Rp as noted in equation 27 in RWG80
			// C++ doesn't have an easy way to take the norm of a vector
			// Lets do it manually from the node class
			// Here is Rp = norm(triangle_centre - r_prime)
			float R_p = this->triangles[p].getCentre().getDifferenceBetween(r_prime).getNorm();

			// Now lets get Greens function
			// That is, Green = (e^(-j*k*R_p)) / R_p
			// Remeber that the imaginary unit 'j' is defined as this->j	
			std::complex<float> greens_function = std::exp(std::complex<float>(-1.0, 0) * this->j * this->k * R_p) / R_p;

			// Finally, lets calculate three  of the four I's
			Ipq = Ipq + (this->quadrature_weights_values[i][0] * greens_function);  
			
			Ipq_xi = Ipq_xi + (this->quadrature_weights_values[i][0] * 
							   this->quadrature_weights_values[i][1] * 
							   greens_function);

			Ipq_eta = Ipq_eta + (this->quadrature_weights_values[i][0] * 
				                 this->quadrature_weights_values[i][2] * 
				                 greens_function);  
		}


		// Lets calculate the final I and push to the vector
		Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
		i_vector.push_back(Ipq);
		i_vector.push_back(Ipq_xi);
		i_vector.push_back(Ipq_eta);
		i_vector.push_back(Ipq_zeta);
		// TODO: DELETE
		//std::cout << Ipq <<"\t" << Ipq_xi <<"\t" << Ipq_eta <<"\t" << Ipq_zeta <<"\t" << std::endl;
	}

	return i_vector;
}




























