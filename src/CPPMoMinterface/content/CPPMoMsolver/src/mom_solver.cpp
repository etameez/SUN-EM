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