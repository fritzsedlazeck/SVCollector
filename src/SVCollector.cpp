/*
 ============================================================================
 Name        : SVCollector.cpp
 Author      : Fritz Sedlazeck
 Version     :
 Copyright   : MIT license
 Description : A greedy framework to select samples based on initial observations.
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <map>
#include <math.h>

#include "Select_samples.h"

int main(int argc, char **argv) {

	if (argc == 1) {
		std::cerr <<"./SVCollector my_svs_vcf_file output_ranked"<<std::endl;
		std::cerr <<"my_svs_vcf_file: A valid uncompressed multisample VCF file."<<std::endl;
		std::cerr <<"output_ranked: The file to write out the ranked list with additional information."<<std::endl;
		return -1;
	} else {
		select_greedy(std::string(argv[2]), std::string(argv[3]));

	}

	return 0;
}
