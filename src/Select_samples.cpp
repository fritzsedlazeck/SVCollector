/*
 * Select_samples.cpp
 *
 *  Created on: Feb 27, 2018
 *      Author: sedlazec
 */

#include "Select_samples.h"

bool genotype_parse(char * buffer) {
	if ((buffer[0] == '.' || (buffer[0] == '0' && buffer[2] == '0')) || (strncmp(buffer, "./.:0", 5) == 0 || strncmp(buffer, "./.:.", 4) == 0)) {
		return false;
	}
	return true;
}

void parse_tmp_file(std::string tmp_file, int id_to_ignore, std::vector<double> & matrix,bool use_alleles) {
	std::string buffer;
	std::ifstream myfile;
	myfile.open(tmp_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << tmp_file.c_str() << std::endl;
		exit(0);
	}
	std::string output = tmp_file;
	output += "v2";
	FILE *file;
	file = fopen(output.c_str(), "w");
	getline(myfile, buffer);

	for (size_t i = 0; i < matrix.size(); i++) {
		matrix[i] = 0;
	}

	while (!myfile.eof()) {
		bool store = true;
		double freq=1;
		if(use_alleles){
			freq=atof(&buffer[0]); //restore the original allele count.
		}
		for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (buffer[i-1]==':' || buffer[i - 1] == ',') {//just parse the sample IDs.
				int id = atoi(&buffer[i]);
				if (id == id_to_ignore) { //if not yet selected This takes time!
					store = false;
					break;
				}
			}
		}

		if (store) {
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (buffer[i - 1] == ':' || buffer[i - 1] == ',') {
					int id = atoi(&buffer[i]);
					matrix[id]+=freq; //AF
				}
			}
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		}
		getline(myfile, buffer);
	}
	myfile.close();
	fclose(file);

	std::stringstream ss;
	ss << "mv ";
	ss << output;
	ss << " ";
	ss << tmp_file;
	system(ss.str().c_str());
}

std::vector<double> prep_file(std::string vcf_file, int min_allele_count, std::vector<std::string> & names, double &num_snp, std::string output,bool use_alleles) {
	std::cout << "Initial assessment of VCF:" << endl;
	std::string buffer;
	std::ifstream myfile;
	std::vector<double> matrix;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");

	getline(myfile, buffer);
	int line = 0;
	num_snp=0;
	while (!myfile.eof()) {
		if (names.empty() && (buffer[0] == '#' && buffer[1] == 'C')) { //parse names
			int count = 0;
			std::string id = "";
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i] != '\t') {
					id += buffer[i];
				}
				if (buffer[i] == '\t') {
					if (!id.empty()) {
						names.push_back(id);
						id = "";
					}
					count++;
				}
			}
			if (!id.empty()) {
				names.push_back(id);
			}

		} else if (buffer[0] != '#') { //parse variants;

			matrix.resize(names.size(), 0);
			line++;
			int count = 0;
			int num = 0;
			double alleles = 0;
			double af=-1; //just used for adams selection.

			std::string entries;
			entries.resize(names.size(), '0');
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if(count==7 && strncmp(&buffer[i],";AF=",4)==0){
					af=atof(&buffer[i+4]);
				}
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						entries[num] = '1';
						alleles++;
					}
					num++;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (alleles > (double)min_allele_count) {

				// use global AF !
				std::stringstream ss;
				double freq=1;
				if(use_alleles&& af>-1){ // if AF is givenin the format field.
					freq=af;
				}else if(use_alleles){ //else we compute it:
					freq=(alleles/(double(names.size())));
				}
		//		cout<<freq<< " " <<alleles<<endl;
				ss<<freq;
				num_snp+=freq;
				ss<<':';
				for (size_t j = 0; j < entries.size(); j++) {
					if (entries[j] == '1') {
						//print out j
						ss << j;
						ss << ',';
						matrix[j]+=freq;
					}
				}
				fprintf(file, "%s", ss.str().c_str());
				fprintf(file, "%c", '\n');
				if (line % 10000 == 0) {
					std::cout << "\tentries: " << line << std::endl;
				}
			}
		}

		getline(myfile, buffer);
	}
	fclose(file);
	myfile.close();
	return matrix;
}

void print_mat(std::vector<int> svs_count_mat) {
	for (size_t i = 0; i < svs_count_mat.size(); i++) {
		std::cout << svs_count_mat[i] << "\t";
	}
	std::cout << std::endl;
	std::cout << std::endl;
}

void select_greedy(std::string vcf_file, int min_allele_count, int num_samples, int alleles, std::string output) {
	std::vector<std::string> sample_names;
	double total_svs = 0;

	//we can actually just use a vector instead!
	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<double> svs_count_mat = prep_file(vcf_file,min_allele_count, sample_names, total_svs, tmp_file,(bool)(alleles==1));
	FILE *file;
	file = fopen(output.c_str(), "w");
	cout << "Total SV: " << total_svs << endl;

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\tSVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	double captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		//select max on main diag
		double max = 0;
		int max_id = -1;

		//select based on greedy:
		for (size_t j = 0; j < sample_names.size(); j++) {
			if (max < svs_count_mat[j]) {
				max = svs_count_mat[j];
				max_id = j;
			}
		}

		captured_svs += max;
		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", max);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", captured_svs / total_svs);
		fprintf(file, "%c", '\n');

		if (max == 0) {
			std::cerr << "No more samples to select from" << std::endl;
			break;
		}
		//erase joined svs
		parse_tmp_file(tmp_file, max_id, svs_count_mat,(bool)(alleles==1)); //span a  NxN matrix and stores the shared SVs
	}
	fclose(file);
	std::stringstream ss;
	ss << "rm ";
	ss << tmp_file;
	system(ss.str().c_str());
}

void select_topN(std::string vcf_file, int num_samples, std::string output) {

	std::vector<std::string> sample_names;
	double total_svs = 0;

	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<double> svs_count_mat = prep_file(vcf_file, 0,sample_names, total_svs, tmp_file,false);
	cout << "Total SV: " << total_svs << endl;
	std::vector<double> initial_counts = svs_count_mat;

	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\t%_SVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	double captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		//select max on main diag
		double max = 0;
		int max_id = -1;

		//select based on max num SVs:
		for (size_t j = 0; j < sample_names.size(); j++) {
			if (max < initial_counts[j]) {
				max = initial_counts[j];
				max_id = j;
			}
		}
		captured_svs += svs_count_mat[max_id];
		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", svs_count_mat[max_id]);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');

		if (max == 0) {
			std::cerr << "No more samples to select from" << std::endl;
			break;
		}

		initial_counts[max_id] = 0;
		parse_tmp_file(tmp_file, max_id, svs_count_mat,false); //span a  NxN matrix and stores the shared SVs
	}
	fclose(file);
	std::stringstream ss;
	ss << "rm ";
	ss << tmp_file;
	system(ss.str().c_str());
}

void select_random(std::string vcf_file, int num_samples, std::string output) {
	std::vector<std::string> sample_names;
	double total_svs = 0;
	srand(time(NULL));

//we can actually just use a vector instead!
	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<double> svs_count_mat = prep_file(vcf_file,0, sample_names, total_svs, tmp_file,false);
//print_mat(svs_count_mat);

	std::vector<int> ids;	// just to avoid that the same ID is picked twice.
	ids.resize(sample_names.size());
	for (size_t i = 0; i < ids.size(); i++) {
		ids[i] = i;
	}

	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\tSVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	double captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		//select max on main diag
		double max = 0;
		int max_id = -1;
		while (max_id == -1 || ids[max_id] == -1) {
			max_id = rand() % sample_names.size();
		}

		max = svs_count_mat[max_id];
		captured_svs += max;
		//std::cout<<"entry: "<<ids[max_id]<<" "<<max<<std::endl;

		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", max);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');
		ids[max_id] = -1;

		svs_count_mat[max_id] = 0;
		parse_tmp_file(tmp_file, max_id, svs_count_mat,false); //span a  NxN matrix and stores the shared SVs

	}
	fclose(file);
	std::stringstream ss;
	ss << "rm ";
	ss << tmp_file;
	system(ss.str().c_str());
}

std::map<std::string, bool> read_names(std::string chosen) {
	std::string buffer;
	std::ifstream myfile;
	std::map<std::string, bool> names;
	myfile.open(chosen.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << chosen.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	while (!myfile.eof()) {
		std::string name = "";
		for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\t'; i++) {
			name += buffer[i];
		}
		names[name] = true;
		getline(myfile, buffer);
	}

	myfile.close();
	return names;
}

void generate_matrix(std::string vcf_file, std::string output, std::string chosen) {
	std::map<std::string, bool> chosen_names = read_names(chosen);
	std::string buffer;
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");

	getline(myfile, buffer);
	int line = 0;
	while (!myfile.eof()) {
		if ((buffer[0] == '#' && buffer[1] == 'C')) { //parse names
			int count = 0;
			std::string id = "";
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i] != '\t') {
					id += buffer[i];
				}
				if (buffer[i] == '\t') {
					if (!id.empty()) {
						if (chosen_names.find(id) == chosen_names.end()) {
							id = " ";
						}
						fprintf(file, "%s", id.c_str());
						fprintf(file, "%c", '\t');
						id = "";
					}
					count++;
				}
			}
			if (chosen_names.find(id) == chosen_names.end()) {
				id = " ";
			}
			if (!id.empty()) {
				fprintf(file, "%s", id.c_str());
			}
			fprintf(file, "%c", '\n');

		} else if (buffer[0] != '#') { //parse sv;
			line++;
			int count = 0;
			int alleles = 0;
			std::stringstream ss;
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						ss << "1\t";
						alleles++;
					} else {
						ss << "0\t";
					}
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (alleles > 2) {
				size_t size = ss.str().size();
				//	cout<<ss.str().substr(0,size-1).c_str()<<endl;
				fprintf(file, "%s", ss.str().substr(0, size - 1).c_str());
				fprintf(file, "%c", '\n');
				if (line % 10000 == 0) {
					std::cout << "\tentries: " << line << std::endl;
				}
			}
		}

		getline(myfile, buffer);
	}
	fclose(file);
	myfile.close();

}

