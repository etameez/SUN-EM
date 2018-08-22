#include "mom_file_writer.h"

/**
 *   \file mom_file_writer.cpp
 *   \brief Write results to a file that will be read by MATLAB
 *
 *  Detailed description
 *	A class to write Zmn, In and Vm to a file to be read by MATLAB 	
 *
 *  Author:  Tameez Ebrahim
 *  Created: 09 August 2018 
 *
 */

MoMFileWriter::MoMFileWriter()
{

}

void MoMFileWriter::writeIlhsToFile(std::string filename, std::vector<std::complex<double>> ilhs)
{
	// Lets write the current vector to a file
	// The form will be a+bi for MATLAB to read easily

	// Lets start by defining a temporary variable to store the imaginary piece(b)
	// and a variable for the sign
	double temp_imag;
	std::string sign;

	// Lets open the file
	std::ofstream file;	
	file.open(filename);

	for(int i = 0; i < ilhs.size(); i++)
	{
		temp_imag = ilhs[i].imag();

		// Check if the imaginary part is negative
		// If it is, set sign to - and multiply the imaginary part by -1
		// Else, set sign to +
		if(temp_imag < 0)
		{
			temp_imag *= -1.000;
			sign = "-";
		}
		else
		{
			sign = "+";
		}

		// Lets write a+bi to the file
		file << ilhs[i].real() << sign << temp_imag << "i" << std::endl;
	}

}	