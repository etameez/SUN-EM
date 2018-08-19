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
	double temp_imag;
	std::string sign;
	std::ofstream file;	
	file.open(filename);

	for(int i = 0; i < ilhs.size(); i++)
	{
		temp_imag = ilhs[i].imag();

		if(temp_imag < 0)
		{
			temp_imag *= -1.000;
			sign = "-";
		}
		else
		{
			sign = "+";
		}

		file << ilhs[i].real() << sign << temp_imag << "i" << std::endl;
	}

}	