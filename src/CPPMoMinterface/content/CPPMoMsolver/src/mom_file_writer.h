#ifndef MOM_FILE_WRITER
#define MOM_FILE_WRITER

#include <string>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>

class MoMFileWriter
{
	public:
		MoMFileWriter();
		void writeIlhsToFile(std::string filename, std::vector<std::complex<double>> ilhs);
};

#endif