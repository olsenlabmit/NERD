#include<iostream>
#include<cstdio>
#include<vector>
#include<utility>
#include <vector>
#include <ctime>
#include <string>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <math.h> 
#include <cmath>
#include <random>
#include <chrono>
#include <set>
#include <map>
#include <numeric>
#include <fstream>
#include "subnet.h"

using namespace std;

vector <vector <vector <double>>> N_template;
double b; //Angstrom
int number_of_molecules_0, num_KMC, total_A, total_B;
vector <double> prob(MAX + 1, 0.0), conc; //M
double conv, V, f; //L
vector <int> num_molecules;

class KMC
{
private:
	vector <vector <unsigned short>> DAA, DBB, DAB;
	vector <int> crosslink_index_A, crosslink_index_B, crosslink_type_A, crosslink_type_B, reacted_index_A = vector<int>(total_A, -1), reacted_index_B = vector<int>(total_B, -1);
	std::set <int> unreacted_A, unreacted_B;
	vector <double> sum_A = vector<double>(total_A, 0.0);
	double sum_tot;
public:
	KMC()
	{
		sum_tot = 0;
	}
	void initialize()
	{
		vector <unsigned short> proxy;
		for (int i = 1; i < total_A; i++)
		{
			for (int j = 0; j < i; j++)
			{
				proxy.push_back(MAX);
			}
			DAA.push_back(proxy);
			proxy.clear();
		}

		for (int i = 1; i < total_B; i++)
		{
			for (int j = 0; j < i; j++)
			{
				proxy.push_back(MAX);
			}
			DBB.push_back(proxy);
			proxy.clear();
		}

		for (int i = 0; i < total_A; i++)
		{
			for (int j = 0; j < total_B; j++)
			{
				proxy.push_back(MAX);
			}
			DAB.push_back(proxy);
			proxy.clear();
		}

		//Initialize DAA, DAB, DBB 
		vector <int> A_or_B = { 0, 1 };
		for (int start = 0; start < 2; start++)
		{
			int counter = 0;
			int AB = A_or_B[start];
			for (int i = 0; i < num_molecules.size(); i++) //types of reactants
			{
				for (int j = 0; j < num_molecules[i]; j++) //number of NRx reactants
				{
					for (int k = 0; k < N_template[i][AB].size(); k++) //only As or Bs
					{
						for (int l = k + 1; l < N_template[i][AB].size(); l++)
						{
							//cout << std::to_string(l + counter + 1) + " " + std::to_string(k + counter + 1) << endl;
							//Adjust indices by subtracting 2 from the row and 1 from the column
							if (AB == 0)
							{
								DAA[l + counter - 1][k + counter] = N_template[i][AB][k] + N_template[i][AB][l];
							}
							else
							{
								DBB[l + counter - 1][k + counter] = N_template[i][AB][k] + N_template[i][AB][l];
							}
						}
					}
					counter += N_template[i][AB].size();
				}
			}
		}

		//Initialize crosslink_index_A, crosslink_index_B
		int counterA = 0, counterB = 0, c_index_A = 1, c_index_B = 1, c_type_A = 1, c_type_B = 1;
		for (int i = 0; i < num_molecules.size(); i++) //reactant x
		{
			for (int j = 0; j < num_molecules[i]; j++) //number of NRx reactants
			{
				for (int k = 0; k < N_template[i][0].size(); k++)
				{
					for (int l = 0; l < N_template[i][1].size(); l++)
					{
						DAB[k + counterA][l + counterB] = N_template[i][0][k] + N_template[i][1][l];
					}
				}
				counterA += N_template[i][0].size();
				counterB += N_template[i][1].size();
				for (int k = 0; k < N_template[i][0].size(); k++)
				{
					crosslink_index_A.push_back(c_index_A);
					crosslink_type_A.push_back(c_type_A);
				}
				for (int k = 0; k < N_template[i][1].size(); k++)
				{
					crosslink_index_B.push_back(c_index_B);
					crosslink_type_B.push_back(c_type_B);
				}
				c_index_A += 1;
				c_index_B += 1;
			}
			c_type_A += 1;
			c_type_B += 1;
		}

		//Initialize sum_A, sum_tot, unreacted_A, unreacted_B
		for (int i = 0; i < total_A; i++)
		{
			unreacted_A.insert(i);
			double sum = 0.0;
			for (int j = 0; j < total_B; j++)
			{
				sum += prob[DAB[i][j]];
				if (i == 0)
				{
					unreacted_B.insert(j);
				}
			}
			sum_A[i] = sum;
			sum_tot += sum;
		}
		//cout << "indexA" << endl; for (int i = 0; i < crosslink_index_A.size(); i++) cout << crosslink_index_A[i] << endl;
		//cout << "indexB" << endl; for (int i = 0; i < crosslink_index_B.size(); i++) cout << crosslink_index_B[i] << endl;
		//cout << "typeA" << endl; for (int i = 0; i < crosslink_type_A.size(); i++) cout << crosslink_type_A[i] << endl;
		//cout << "typeB" << endl; for (int i = 0; i < crosslink_type_B.size(); i++) cout << crosslink_type_B[i] << endl;
		//cout << "DAB" << endl; for (int i = 0; i < DAB.size(); i++) { for (int j = 0; j < DAB[i].size(); j++) cout << DAB[i][j] << " "; cout << endl; }
		//cout << "DAA" << endl; for (int i = 0; i < DAA.size(); i++) { for (int j = 0; j < DAA[i].size(); j++) cout << DAA[i][j] << " "; cout << endl; }
		//cout << "DBB" << endl; for (int i = 0; i < DBB.size(); i++) { for (int j = 0; j < DBB[i].size(); j++) cout << DBB[i][j] << " "; cout << endl; }
	}
	vector <int> choose()
	{
		double r = rand32() / 4294967296.0 * sum_tot;
		double sum = 0.0;
		int indexA = 0;
		for (int i = 0; i < sum_A.size(); i++)
		{
			sum += sum_A[i];
			if (sum >= r)
			{
				sum -= sum_A[i];
				indexA = i;
				break;
			}
		}
		for (int j = 0; j < DAB[indexA].size(); j++)
		{
			if (reacted_index_B[j] == -1)
			{
				sum += prob[DAB[indexA][j]];
				if (sum >= r)
				{
					return { indexA, j };
				}
			}
		}
	}
	void run_KMC()
	{
		int num_steps;
		if (total_A < total_B)
		{
			num_steps = total_A * conv;
		}
		else
		{
			num_steps = total_B * conv;
		}
		// cout << "Number of KMC Steps: " << num_steps << endl;

		for (int i = 0; i < num_steps; i++)
		{
			vector <int> ind = choose();
			int chosen_A = ind[0], chosen_B = ind[1];
			// cout << "IndexA: " << chosen_A << endl << "IndexB: " << chosen_B << endl;
			// cout << "sum_total " << sum_tot << endl;
			// for (int j = 0; j < sum_A.size(); j++) { cout << sum_A[j] << endl; }
			for (int iA = 0; iA < total_A; iA++)
			{
				if (reacted_index_A[iA] == -1)
				{
					sum_A[iA] -= prob[DAB[iA][chosen_B]];
					sum_tot -= prob[DAB[iA][chosen_B]];
				}
			}
			sum_tot -= sum_A[chosen_A];
			sum_A[chosen_A] = 0.0;
			reacted_index_A[chosen_A] = chosen_B;
			reacted_index_B[chosen_B] = chosen_A;
			bool intramolecular = DAB[chosen_A][chosen_B] > 0 && DAB[chosen_A][chosen_B] < MAX;
			bool intermolecular = DAB[chosen_A][chosen_B] == MAX;
			DAB[chosen_A][chosen_B] = 0;
			unreacted_A.erase(chosen_A);
			unreacted_B.erase(chosen_B);
			if (intramolecular)
			{
				vector <int> cAcB;
				int totA = 0;
				for (std::set<int>::iterator it = unreacted_A.begin(); it != unreacted_A.end(); it++)
				{
					if (DAB[*it][chosen_B] < MAX)
					{
						cAcB.push_back(*it);
						totA += 1;
					}
				}
				for (std::set<int>::iterator it = unreacted_B.begin(); it != unreacted_B.end(); it++)
				{
					if (DAB[chosen_A][*it] < MAX)
					{
						cAcB.push_back(*it);
					}
				}
				// cout << "totA: " << totA << endl;
				// for (int j = 0; j < cAcB.size(); j++)
				// {
				// 	cout << cAcB[j] << endl;
				// }

				for (int one = 0; one < cAcB.size(); one++)
				{
					for (int two = one + 1; two < cAcB.size(); two++)
					{
						unsigned short	AA_one = 0, AA_two = 0, AB_one = 0, BA_two = 0, BA_one = 0, AB_two = 0, BB_one = 0, BB_two = 0;
						if (one < totA)
						{
							AA_one = DAA[std::max(cAcB[one], chosen_A) - 1][std::min(cAcB[one], chosen_A)];
							AB_one = DAB[cAcB[one]][chosen_B];
						}
						else
						{
							BA_one = DAB[chosen_A][cAcB[one]];
							BB_one = DBB[std::max(cAcB[one], chosen_B) - 1][std::min(cAcB[one], chosen_B)];
						}
						if (two < totA)
						{
							AA_two = DAA[std::max(chosen_A, cAcB[two]) - 1][std::min(chosen_A, cAcB[two])];
							BA_two = DAB[cAcB[two]][chosen_B];
						}
						else
						{
							AB_two = DAB[chosen_A][cAcB[two]];
							BB_two = DBB[std::max(chosen_B, cAcB[two]) - 1][std::min(chosen_B, cAcB[two])];
						}

						unsigned short before = 0, path_one = 0, path_two = 0, path_three = 0, path_four = 0, min;
						if (one < totA && two < totA) //A on cluster 1 and A on cluster 2
						{
							before = DAA[std::max(cAcB[one], cAcB[two]) - 1][std::min(cAcB[one], cAcB[two])];
							path_one = AA_one + AA_two;
							path_two = AA_one + BA_two;
							path_three = AB_one + AA_two;
							path_four = AB_one + BA_two;
							vector <unsigned short> shortest = { before, path_one, path_two, path_three, path_four };
							min = *min_element(shortest.begin(), shortest.end());
							DAA[std::max(cAcB[one], cAcB[two]) - 1][std::min(cAcB[one], cAcB[two])] = min;
						}
						else if (one < totA && two >= totA) //A on cluster 1 and B on cluster 2
						{
							before = DAB[cAcB[one]][cAcB[two]];
							path_one = AA_one + AB_two;
							path_two = AA_one + BB_two;
							path_three = AB_one + AB_two;
							path_four = AB_one + BB_two;
							vector <unsigned short> shortest = { before, path_one, path_two, path_three, path_four };
							min = *min_element(shortest.begin(), shortest.end());
							DAB[cAcB[one]][cAcB[two]] = min;

							sum_A[cAcB[one]] = sum_A[cAcB[one]] - prob[before] + prob[min];
							sum_tot = sum_tot - prob[before] + prob[min];
						}
						else if (one >= totA && two < totA) //B on cluster 1 and A on cluster 2
						{
							before = DAB[cAcB[two]][cAcB[one]];
							path_one = BA_one + AA_two;
							path_two = BA_one + BA_two;
							path_three = BB_one + AA_two;
							path_four = BB_one + BA_two;
							vector <unsigned short> shortest = { before, path_one, path_two, path_three, path_four };
							min = *min_element(shortest.begin(), shortest.end());
							DAB[cAcB[two]][cAcB[one]] = min;

							sum_A[cAcB[two]] = sum_A[cAcB[two]] - prob[before] + prob[min];
							sum_tot = sum_tot - prob[before] + prob[min];
						}
						else if (one >= totA && two >= totA) //B on cluster 1 and B on cluster 2
						{
							before = DBB[std::max(cAcB[one], cAcB[two]) - 1][std::min(cAcB[one], cAcB[two])];
							path_one = BA_one + AB_two;
							path_two = BA_one + BB_two;
							path_three = BB_one + AB_two;
							path_four = BB_one + BB_two;
							vector <unsigned short> shortest = { before, path_one, path_two, path_three, path_four };
							min = *min_element(shortest.begin(), shortest.end());
							DBB[std::max(cAcB[one], cAcB[two]) - 1][std::min(cAcB[one], cAcB[two])] = min;
						}
					}
				}
			}
			else if (intermolecular)
			{
				vector <int> cASide, cBSide;
				int totA_Aside = 0, totA_Bside = 0;
				for (std::set<int>::iterator it = unreacted_A.begin(); it != unreacted_A.end(); it++)
				{
					if (DAA[std::max(chosen_A, *it) - 1][std::min(chosen_A, *it)] < MAX)
					{
						cASide.push_back(*it);
						totA_Aside += 1;
					}
					else if (DAB[*it][chosen_B] < MAX)
					{
						cBSide.push_back(*it);
						totA_Bside += 1;
					}
				}
				for (std::set<int>::iterator it = unreacted_B.begin(); it != unreacted_B.end(); it++)
				{
					if (DAB[chosen_A][*it] < MAX) 
					{ 
						cASide.push_back(*it); 
					}
					else if (DBB[std::max(chosen_B, *it) - 1][std::min(chosen_B, *it)] < MAX) 
					{ 
						cBSide.push_back(*it); 
					}
				}
				// cout << "totA_A: " << totA_Aside << endl;
				// cout << "totA_B: " << totA_Bside << endl;
				// cout << "cASide" << endl;
				// for (int j = 0; j < cASide.size(); j++) { cout << cASide[j] << endl; }
				// cout << "cBSide" << endl;
				// for (int j = 0; j < cBSide.size(); j++) { cout << cBSide[j] << endl; }

				for (int one = 0; one < cASide.size(); one++)
				{
					for (int two = 0; two < cBSide.size(); two++)
					{
						unsigned short	AA_one = 0, BA_one = 0, 
										BA_two = 0, BB_two = 0;
						if (one < totA_Aside)
						{
							AA_one = DAA[std::max(chosen_A, cASide[one]) - 1][std::min(chosen_A, cASide[one])];
						}
						else
						{
							BA_one = DAB[chosen_A][cASide[one]];
						}
						if (two < totA_Bside)
						{
							BA_two = DAB[cBSide[two]][chosen_B];
						}
						else
						{
							BB_two = DBB[std::max(chosen_B, cBSide[two]) - 1][std::min(chosen_B, cBSide[two])];
						}

						if (one < totA_Aside && two < totA_Bside) //A on Cluster A and A on Cluster B
						{
							DAA[std::max(cASide[one], cBSide[two]) - 1][std::min(cASide[one], cBSide[two])] = AA_one + BA_two;
						}
						else if (one < totA_Aside && two >= totA_Bside) //A on cluster A and B on cluster B
						{
							DAB[cASide[one]][cBSide[two]] = AA_one + BB_two;

							sum_A[cASide[one]] = sum_A[cASide[one]] - prob[MAX] + prob[AA_one + BB_two];
							sum_tot = sum_tot - prob[MAX] + prob[AA_one + BB_two];
						}
						else if (one >= totA_Aside && two < totA_Bside) //B on cluster A and A on cluster B
						{
							DAB[cBSide[two]][cASide[one]] = BA_one + BA_two;

							sum_A[cBSide[two]] = sum_A[cBSide[two]] - prob[MAX] + prob[BA_one + BA_two];
							sum_tot = sum_tot - prob[MAX] + prob[BA_one + BA_two];
						}
						else if (one >= totA_Aside && two >= totA_Bside) //B on cluster A and B on cluster B
						{
							DBB[std::max(cASide[one], cBSide[two]) - 1][std::min(cASide[one], cBSide[two])] = BA_one + BB_two;
						}
					}
				}
			}
		}
	}
	vector <double> analyze_reaction_graph(ofstream& to_MF)
	{
		//to_MF << "How Many A have not reacted? " << std::count(reacted_index_A.begin(), reacted_index_A.end(), -1) << endl << "How Many B have not reacted? " << std::count(reacted_index_B.begin(), reacted_index_B.end(), -1) << endl;
		//cout << "DAB" << endl; for (int i = 0; i < DAB.size(); i++) { for (int j = 0; j < DAB[i].size(); j++) cout << DAB[i][j] << " "; cout << endl; }
		//cout << "DAA" << endl; for (int i = 0; i < DAA.size(); i++) { for (int j = 0; j < DAA[i].size(); j++) cout << DAA[i][j] << " "; cout << endl; }
		//cout << "DBB" << endl; for (int i = 0; i < DBB.size(); i++) { for (int j = 0; j < DBB[i].size(); j++) cout << DBB[i][j] << " "; cout << endl; }
		//cout << "reacted index A" << endl; for (int i = 0; i < reacted_index_A.size(); i++) { cout << reacted_index_A[i] << endl; }
		//cout << "crosslink_index_A" << endl; for (int i = 0; i < crosslink_index_A.size(); i++) cout << crosslink_index_A[i] << endl;
		//cout << "crosslink_index_B" << endl; for (int i = 0; i < crosslink_index_B.size(); i++) cout << crosslink_index_B[i] << endl;
		//cout << "crosslink_type_A" << endl; for (int i = 0; i < crosslink_type_A.size(); i++) cout << crosslink_type_A[i] << endl;
		//cout << "crosslink_type_B" << endl; for (int i = 0; i < crosslink_type_B.size(); i++) cout << crosslink_type_B[i] << endl;

		std::set <int> already_in_network;
		vector <network> networks;
		vector <int> all_networks_size;
		double strand_ct = 0, primary_junct = 0, secondary_junct = 0, num_junctions = 0, v0 = 0;
		int i = 0;
		while (i < crosslink_index_A.size())
		{
			std::set <int> current_network;
			int start = crosslink_index_A[i];
			if (already_in_network.find(start) != already_in_network.end())
			{
				while (i < crosslink_index_A.size() && crosslink_index_A[i] == start)
				{
					i += 1;
				}
			}
			else
			{
				already_in_network.insert(crosslink_index_A[i]);
				current_network.insert(crosslink_index_A[i]);
				i += 1;
				build_network(current_network, already_in_network);

				vector <crosslink> cross_index_candidate;
				for (std::set<int>::iterator it = current_network.begin(); it != current_network.end(); it++)
				{
					crosslink k = crosslink(*it);
					k.primary_external(reacted_index_A, reacted_index_B, crosslink_index_A, crosslink_index_B, crosslink_type_A, crosslink_type_B);
					cross_index_candidate.push_back(k);
				}
				network n = network(cross_index_candidate, N_template, { V, conv }, reacted_index_A, reacted_index_B, crosslink_index_A, crosslink_index_B, crosslink_type_A, crosslink_type_B);
				int nj = n.primary();
				n.secondary();
				primary_junct += n.get_loop_count()[0];
				secondary_junct += n.get_loop_count()[1];
				num_junctions += nj;
				v0 += n.get_vo();
				strand_ct = n.get_strand_ct();
				networks.push_back(n);
				all_networks_size.push_back(nj);
			}
		}
		return { v0, primary_junct / num_junctions, secondary_junct / num_junctions };
	}
	bool build_network(std::set <int>& current_network, std::set <int>& already_in_network)
	{
		set<int>::iterator itr;
		vector <int> c_new;
		for (itr = current_network.begin(); itr != current_network.end(); itr++)
		{
			// find crosslink index that current crosslink connects to and add to set 
			bool found = false;
			for (int j = 0; j < crosslink_index_A.size(); j++)
			{
				if (*itr == crosslink_index_A[j])
				{
					found = true;
					if (reacted_index_A[j] != -1)
					{
						int b = crosslink_index_B[reacted_index_A[j]];
						if (current_network.find(b) == current_network.end())
						{
							c_new.push_back(b);
							already_in_network.insert(b);
						}
					}
				}
				else if (found)
				{
					break;
				}
			}
			found = false;
			for (int j = 0; j < crosslink_index_B.size(); j++)
			{
				if (*itr == crosslink_index_B[j])
				{
					found = true;
					if (reacted_index_B[j] != -1)
					{
						int a = crosslink_index_A[reacted_index_B[j]];
						if (current_network.find(a) == current_network.end())
						{
							c_new.push_back(a);
							already_in_network.insert(a);
						}
					}
				}
				else if (found)
				{
					break;
				}
			}
		}
		if (c_new.size() == 0) { return true; }
		for (int i = 0; i < c_new.size(); i++)
		{
			current_network.insert(c_new[i]);
		}
		return build_network(current_network, already_in_network);
	}
};

int main(int argc, char* argv[])
{
	if (argc >= 3)
	{
		conc[0] = std::stod(argv[1]);
		conc[1] = conc[0] / 2; // mol/L
		number_of_molecules_0 = std::stoi(argv[2]);
	}

	string myText;
	ifstream MyReadFile("../Info-CG.txt");
	vector <string> inputs;
	while (getline (MyReadFile, myText)) 
	{
		inputs.push_back(myText);
	}
	MyReadFile.close();

	string* a = &inputs[0];
	int n = 1;
	num_KMC = std::stoi(a[n]);
	n += 2;
	conv = stod(a[n]) / 100; // between 0 and 1
	n += 2;
	b = stod(a[n]); // Angstrom
	n += 2;
	int num_precursors = std::stoi(a[n]);
	n += 2;
	vector <int> number;
	for (int i = 0; i < num_precursors; i++)
	{
		// Now, computer reads "NEW"
		n += 2;
		conc.push_back(stod(a[n])); // mM
		n += 2;
		number.push_back(std::stoi(a[n]));
		n += 2;
		int num_A = std::stoi(a[n]);
		n += 1;
		vector <double> A;
		for (int j = 0; j < num_A; j++)
		{
			A.push_back(stod(a[n]));
			n += 1;
		}

		//Now, the computer reads "B"
		n += 1;
		int num_B = std::stoi(a[n]);
		n += 1;
		vector <double> B;
		for (int j = 0; j < num_B; j++)
		{
			B.push_back(stod(a[n]));
			n += 1;
		}
		n += 1;
		// Now, computer reads NEW 

		vector <vector <double>> N = {A, B};
		N_template.push_back(N);
	}

	number_of_molecules_0 = number[0];
	V = number_of_molecules_0 / (6.022 * pow(10, 23)) / (conc[0] / 1000); // molecules * (molecules / mol * mol / L)^ -1 = L
	num_molecules.push_back(number_of_molecules_0);
	for (int i = 1; i < conc.size(); i++)
	{
		num_molecules.push_back(conc[i] / conc[0] * number_of_molecules_0);
	}
	total_A = 0;
	total_B = 0;
	for (int i = 0; i < N_template.size(); i++)
	{
		total_A += N_template[i][0].size() * num_molecules[i];
		total_B += N_template[i][1].size() * num_molecules[i];
	}
	prob[0] = 0.0;
	for (int i = 1; i < MAX; i++)
	{
		prob[i] = 1.0 + pow(3 / (2 * 3.14159 * i * pow(b, 2)), 1.5) * V * pow(10, 27); // V[L] * m^3/1000L * 10^30 A^3 / m^3 * 1 / A^3 = 10^27 prefactor
	}
	prob[MAX] = 1.0;
	seed = time(NULL) + clock();
	RNG_initialize(seed);

	ofstream data_file("Data-File-" + std::to_string(int(conc[0])) + "mM.txt");
	data_file << "Concentration of Molecule 1: " << conc[0] << " mM" << endl;
	for (int i = 0; i < conc.size(); i++)
	{
		data_file << "Molecule " << std::to_string(i + 1) << endl;
		data_file << conc[i] << " mM" << endl;
		data_file << num_molecules[i] << " molecules" << endl;
	}

	data_file << "Total A: " << total_A << endl;
	data_file << "Total B: " << total_B << endl;
	vector <double> v(num_KMC, 0.0), primary(num_KMC, 0.0), secondary(num_KMC, 0.0);
	for (int i = 0; i < num_KMC; i++)
	{
		data_file << "-----------------KMC Simulation #" << (i + 1) << "--------------------" << endl;
		KMC trial = KMC();
		auto t1 = Clock::now();
		trial.initialize();
		trial.run_KMC();
		auto t2 = Clock::now();
		data_file << "Kinetic Monte Carlo for network formation took " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << endl;
		vector <double> sim_info = trial.analyze_reaction_graph(data_file);
		v[i] = sim_info[0];
		primary[i] = sim_info[1];
		secondary[i] = sim_info[2];
		data_file << "Chain Density (mol/L): " << v[i] << endl;
		data_file << "Primary Loop Fraction (# primary loop junctions per junction): " << primary[i] << endl;
		data_file << "Secondary Loop Fraction (# secondary loop junctions per junction): " << secondary[i] << endl;
		auto t3 = Clock::now();
		data_file << "Property Prediction takes " << std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count() << " seconds" << endl;
	}
	data_file.close();

	return 0;
}