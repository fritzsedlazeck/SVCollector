/*
 * Select_samples.cpp
 *
 *  Created on: Feb 27, 2018
 *      Author: sedlazec
 */

#include "Select_samples.h"

bool genotype_parse(char * buffer) {
	//TODO do we need to catch other situations?

	if ((buffer[0] == '0' && buffer[2] == '0') || (buffer[0] == '.' && buffer[2] == '.')) {
		return false;
	}
	//0/0 ./.
	return true;
}

std::vector<int> parase_matrix(std::string vcf_file, std::vector<std::string> & names, std::map<int, bool> taken_ids, int &num) {
	std::cout << "Parsing:" << endl;
	std::string buffer;
	std::ifstream myfile;
	std::vector<int> matrix;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	int line = 0;
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

		} else if (buffer[0] != '#') { //parse svs;
			if (matrix.empty()) { //init pairwise matrix;
				std::vector<int> tmp;
				matrix.resize(names.size(), 0);
			}
			line++;
			num++;
			//bool discard = false;
			int count = 0;
			bool include = false;
			int num = 0;

			std::string entries;
			entries.resize(names.size(), '0');
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						if (taken_ids.find(num) != taken_ids.end()) {
							include = false;
							break;
						}
						include = true;
						entries[num] = '1';
					}
					num++;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (include) {
				for (size_t j = 0; j < entries.size(); j++) {
					if (entries[j] == '1') {
						matrix[j]++;
					}
				}
			}
			if (line % 10000 == 0) {
				std::cout << "\tentries: " << line << std::endl;
			}
		}

		getline(myfile, buffer);
	}

	cout << "Fin" << endl;
	myfile.close();
	return matrix;
}
void parse_tmp_file(std::string tmp_file, std::map<int, bool> taken_ids, std::vector<int> & matrix) {
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
	int num = 0;
	while (!myfile.eof()) {
		bool store = true;
		for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (i == 0 || buffer[i - 1] == ',') {
				int id = atoi(&buffer[i]);
				if (taken_ids.find(id) != taken_ids.end()) { //if not yet selected This takes time!
					store = false;
				}
			}
		}

		if (store) {
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (i == 0 || buffer[i - 1] == ',') {
					int id = atoi(&buffer[i]);
					matrix[id]++;
				}
			}
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		}

		/*num++;
		 if (num % 10000 == 0) {
		 std::cout << "\tparsed entries:" << num << std::endl;
		 }*/
		getline(myfile, buffer);
	}
	myfile.close();
	fclose(file);

	std::stringstream ss;
	ss << "mv ";
	ss << output;
	ss << " ";
	ss << tmp_file;
	int i = system(ss.str().c_str());

}

std::vector<int> prep_file(std::string vcf_file, std::vector<std::string> & names, int &num, std::string output) {
	std::cout << "Initial assessment of VCF:" << endl;
	std::string buffer;
	std::ifstream myfile;
	std::vector<int> matrix;
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

		} else if (buffer[0] != '#') { //parse svs;

			std::vector<int> tmp;
			matrix.resize(names.size(), 0);

			line++;
			num++;

			int count = 0;
			int num = 0;

			std::string entries;
			entries.resize(names.size(), '0');
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						entries[num] = '1';
					}
					num++;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			std::stringstream ss;

			for (size_t j = 0; j < entries.size(); j++) {
				if (entries[j] == '1') {
					//print out j
					ss << j;
					ss << ',';
					matrix[j]++;
				}
			}
			fprintf(file, "%s", ss.str().c_str());
			fprintf(file, "%c", '\n');
			if (line % 10000 == 0) {
				std::cout << "\tentries: " << line << std::endl;
			}
		}

		getline(myfile, buffer);
	}
	fclose(file);
	//exit(0);
	//cout << "Fin" << endl;
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
/*
 void select_greedy(std::string vcf_file, int num_samples, std::string output) {
 std::map<int, bool> taken_ids;
 std::vector<std::string> sample_names;
 int total_svs = 0;

 //we can actually just use a vector instead!
 std::vector<int> svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids, total_svs); //span a  NxN matrix and stores the shared SVs
 //print_mat(svs_count_mat);

 FILE *file;
 file = fopen(output.c_str(), "w");

 fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\t%_SVs_captured\n");
 std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
 int captured_svs = 0;
 for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
 //select max on main diag
 int max = 0;
 int max_id = -1;

 //select based on greedy:
 for (size_t j = 0; j < sample_names.size(); j++) {
 //	cout << svs_count_mat[j][j] << "\t";
 if (max < svs_count_mat[j]) {
 max = svs_count_mat[j];
 max_id = j;
 }
 }

 captured_svs += max;
 //	std::cout <<"RANK:\t"<<i<<"\t"<<sample_names[max_id]<<"\t"<<max<<"\t"<<captured_svs<<"\t"<< (double) captured_svs / (double) total_svs<< std::endl;
 fprintf(file, "%s", sample_names[max_id].c_str());
 fprintf(file, "%c", '\t');
 fprintf(file, "%i", max);
 fprintf(file, "%c", '\t');
 fprintf(file, "%i", captured_svs);
 fprintf(file, "%c", '\t');
 fprintf(file, "%f", (double) captured_svs / (double) total_svs);
 fprintf(file, "%c", '\n');

 if (max == 0) {
 std::cerr << "No more samples to select from" << std::endl;
 break;
 }
 //	print_mat(svs_count_mat);
 //print max_id and max
 //erase joined svs over matrix given the pairwise matrix.
 taken_ids[max_id] = true;
 total_svs = 0;
 svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids, total_svs); //span a  NxN matrix and stores the shared SVs
 }
 fclose(file);
 }
 */
void select_greedy(std::string vcf_file, int num_samples, std::string output) {
	std::map<int, bool> taken_ids;
	std::vector<std::string> sample_names;
	int total_svs = 0;

	//we can actually just use a vector instead!
	//print_mat(svs_count_mat);
	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<int> svs_count_mat = prep_file(vcf_file, sample_names, total_svs, tmp_file);
	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\tSVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	int captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		std::cout << "\tSelecting sample: " << i + 1 << std::endl;
		//select max on main diag
		int max = 0;
		int max_id = -1;

		//select based on greedy:
		for (size_t j = 0; j < sample_names.size(); j++) {
			//	cout << svs_count_mat[j][j] << "\t";
			if (max < svs_count_mat[j]) {
				max = svs_count_mat[j];
				max_id = j;
			}
		}

		captured_svs += max;
		//	std::cout <<"RANK:\t"<<i<<"\t"<<sample_names[max_id]<<"\t"<<max<<"\t"<<captured_svs<<"\t"<< (double) captured_svs / (double) total_svs<< std::endl;
		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", max);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');

		if (max == 0) {
			std::cerr << "No more samples to select from" << std::endl;
			break;
		}
		//	print_mat(svs_count_mat);
		//print max_id and max
		//erase joined svs over matrix given the pairwise matrix.
		taken_ids[max_id] = true;

		parse_tmp_file(tmp_file, taken_ids, svs_count_mat); //span a  NxN matrix and stores the shared SVs
	}
	fclose(file);

	//TODO: delet tmp_file!!
}

void select_topN(std::string vcf_file, int num_samples, std::string output) {
	std::map<int, bool> taken_ids;
	std::vector<std::string> sample_names;
	int total_svs = 0;

//we can actually just use a vector instead!
	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<int> svs_count_mat = prep_file(vcf_file, sample_names, total_svs, tmp_file);

//	std::vector<int> svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids, total_svs); //span a  NxN matrix and stores the shared SVs
//print_mat(svs_count_mat);

	std::vector<int> initial_counts = svs_count_mat;

	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\t%_SVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	int captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		//select max on main diag
		int max = 0;
		int max_id = -1;

		//select based on max num SVs:
		for (size_t j = 0; j < sample_names.size(); j++) {
			//	cout << svs_count_mat[j][j] << "\t";
			if (max < initial_counts[j]) {
				max = initial_counts[j];
				max_id = j;
			}
		}
		//std::cout<<"entry: "<<sample_names[max_id]<<" "<<svs_count_mat[max_id]<<std::endl;
		captured_svs += svs_count_mat[max_id];
		initial_counts[max_id] = 0;

		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", svs_count_mat[max_id]);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');

		if (max == 0) {
			std::cerr << "No more samples to select from" << std::endl;
			break;
		}

		taken_ids[max_id] = true;
		parse_tmp_file(tmp_file, taken_ids, svs_count_mat); //span a  NxN matrix and stores the shared SVs

	}
	fclose(file);
}

void select_random(std::string vcf_file, int num_samples, std::string output) {
	std::map<int, bool> taken_ids;
	std::vector<std::string> sample_names;
	int total_svs = 0;

	srand(time(NULL));

//we can actually just use a vector instead!
	std::string tmp_file = output;
	tmp_file += "_tmp";
	std::vector<int> svs_count_mat = prep_file(vcf_file, sample_names, total_svs, tmp_file);
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
	int captured_svs = 0;
	for (size_t i = 0; i < sample_names.size() && i < num_samples; i++) {
		//select max on main diag
		int max = 0;
		int max_id = -1;
		while (max_id == -1 || ids[max_id] == -1) {
			max_id = rand() % sample_names.size();
		}

		max = svs_count_mat[max_id];
		captured_svs += max;
		//std::cout<<"entry: "<<ids[max_id]<<" "<<max<<std::endl;

		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", max);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');
		ids[max_id] = -1;
		//if (max == 0) {
		//		std::cerr << "No more samples to select from" << std::endl;
//			break;
//		}
		svs_count_mat[max_id] = 0;
		taken_ids[max_id] = true;
		parse_tmp_file(tmp_file, taken_ids, svs_count_mat); //span a  NxN matrix and stores the shared SVs

	}
	fclose(file);
}

