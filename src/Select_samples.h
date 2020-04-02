/*
 * Select_samples.h
 *
 *  Created on: Feb 27, 2018
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_SELECT_SAMPLES_H_
#define ANALYSIS_SV_SELECT_SAMPLES_H_
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
//#include"../vcfs/Merge_VCF.h"
using namespace std;


void select_greedy(std::string vcf_file, int min_allele_count, int num_samples, int alleles,  std::string output, std::string preselected_file);
void select_topN(std::string vcf_file, int num_samples, bool use_alleles, std::string output);
void select_random(std::string vcf_file, int num_samples, bool use_alleles, std::string output);
void select_greedyv2(std::string vcf_file, int num_samples, std::string output);
void generate_matrix(std::string vcf_file, std::string output,std::string chosen);

#endif /* ANALYSIS_SV_SELECT_SAMPLES_H_ */
