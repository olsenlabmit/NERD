#include<iostream>
#include<cstdio>
#include<vector>
#include<utility>
#include<vector>
#include<ctime>
#include <string>
#include<algorithm>
#include <fstream>
#include <numeric>
#include <math.h> 
#include <cmath>
#include <chrono>
#include <set>
#include <map>
#define MATRIX_A 0x9908b0dfUL //for RNG
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

using namespace std;
const int MAX = pow(2, sizeof(short) * 8) - 1;
typedef std::chrono::high_resolution_clock Clock;

// RNG related
unsigned long mt[624];
int mti = 625;
double fn[128], wn[128];
int kn[128];
void RNG_initialize(unsigned long);
unsigned long rand32();
size_t seed;

class crosslink
{
private:
	int cross, type;
	vector <int> A_ext, B_ext, match_A_original, match_B_original, match_A_hold, match_B_hold;
	double external, primary, external_original, external_subtract;
	bool is_subnet, store, is_bound, cprimary;
public:
	crosslink(int cross_in)
	{
		cross = cross_in;
		external = 0;
		primary = 0;

		is_subnet = false;
		store = false;
		is_bound = false;
		cprimary = false;
	}
	void primary_external(vector <int> reacted_index_A, vector <int> reacted_index_B, vector <int> crosslink_index_A, vector <int> crosslink_index_B, vector <int> crosslink_type_A, vector <int> crosslink_type_B)
	{
		bool found = false;
		for (int i = 0; i < crosslink_index_A.size(); i++)
		{
			if (cross == crosslink_index_A[i])
			{
				type = crosslink_type_A[i];
				found = true;
				if (reacted_index_A[i] != -1 && cross != crosslink_index_B[reacted_index_A[i]])
				{
					external += 1;
					A_ext.push_back(reacted_index_A[i]);
				}
				else if (reacted_index_A[i] != -1 && cross == crosslink_index_B[reacted_index_A[i]])
				{
					primary += 0.5;
					A_ext.push_back(-1);
				}
				else
				{
					A_ext.push_back(-1);
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
			if (cross == crosslink_index_B[j])
			{
				type = crosslink_type_B[j];
				found = true;
				if (reacted_index_B[j] != -1 && cross != crosslink_index_A[reacted_index_B[j]])
				{
					external += 1;
					B_ext.push_back(reacted_index_B[j]);
				}
				else if (reacted_index_B[j] != -1 && cross == crosslink_index_A[reacted_index_B[j]])
				{
					primary += 0.5;
					B_ext.push_back(-1);
				}
				else
				{
					B_ext.push_back(-1);
				}
			}
			else if (found)
			{
				break;
			}
		}
		external_original = external;
		match_A_original = A_ext;
		match_B_original = B_ext;
	}
	vector <int> bifunctional_linker(bool landsOnA, int pos_AB_other)
	{
		if (landsOnA)
		{
			for (int pos_AB = 0; pos_AB < A_ext.size(); pos_AB++) { if (A_ext[pos_AB] != -1 && pos_AB != pos_AB_other) { return { pos_AB, true }; } }
			for (int pos_AB = 0; pos_AB < B_ext.size(); pos_AB++) { if (B_ext[pos_AB] != -1) { return { pos_AB, false }; } }
		}
		else
		{
			for (int pos_AB = 0; pos_AB < A_ext.size(); pos_AB++) { if (A_ext[pos_AB] != -1) { return { pos_AB, true }; } }
			for (int pos_AB = 0; pos_AB < B_ext.size(); pos_AB++) { if (B_ext[pos_AB] != -1 && pos_AB != pos_AB_other) { return { pos_AB, false }; } }
		}
	}
	void original()
	{
		external = external_original;
		A_ext = match_A_original;
		B_ext = match_B_original;
	}
	void is_subnetwork(bool isit)
	{
		is_subnet = isit;
	}
	bool is_subnetwork()
	{
		return is_subnet;
	}
	void is_boundary(bool isit)
	{
		is_bound = isit;
	}
	bool is_boundary()
	{
		return is_bound;
	}
	void contains_primary(bool isit)
	{
		cprimary = isit;
	}
	bool contains_primary()
	{
		return cprimary;
	}
	double get_external()
	{
		return external;
	}
	double get_primary()
	{
		return primary;
	}
	int get_cross()
	{
		return cross;
	}
	int get_type()
	{
		return type;
	}
	vector <int> get_A_con()
	{
		return A_ext;
	}
	vector <int> get_B_con()
	{
		return B_ext;
	}
	void hold(double subtract, vector <int> A, vector <int> B)
	{
		external_subtract = subtract;
		match_A_hold = A;
		match_B_hold = B;
		store = true;
	}
	bool is_holding()
	{
		return store;
	}
	void update()
	{
		external = external - external_subtract;
		A_ext = match_A_hold;
		B_ext = match_B_hold;
		store = false;
	}
	bool equals(crosslink a)
	{
		return cross == a.get_cross();
	}
};

class network
{
private:
	vector <crosslink> all_molecules;
	vector <vector <vector <double>>> N_template;
	double conv, V, v0, strand_ct;
	vector <int> reactedIndex_A, reactedIndex_B, crosslink_index_A, crosslink_index_B, crosslink_type_A, crosslink_type_B;
	vector <double> loop_count = vector <double>(5, 0.0);
public:
	network() = default;
	network(vector <crosslink> all_molecules_in,
		vector <vector <vector <double>>> N_template_in, vector <double> properties_in,
		vector <int> rA, vector <int> rB, vector <int> ciA, vector <int> ciB, vector <int> ctA, vector <int> ctB)
	{
		all_molecules = all_molecules_in;
		N_template = N_template_in;
		V = properties_in[0];
		conv = properties_in[1];
		reactedIndex_A = rA;
		reactedIndex_B = rB;
		crosslink_index_A = ciA;
		crosslink_index_B = ciB;
		crosslink_type_A = ctA;
		crosslink_type_B = ctB;
	}

	double primary()
	{
		std::set <int> primary_junct;
		double num_junctions = 0;
		for (int m = 0; m < all_molecules.size(); m++)
		{
			crosslink origin = all_molecules[m];
			if (origin.get_A_con().size() + origin.get_B_con().size() >= 3)
			{
				num_junctions += 1;
			}
			if (origin.get_external() >= 3)
			{
				for (int i = 0; i < origin.get_A_con().size(); i++)
				{
					if (origin.get_A_con()[i] != -1)
					{
						bool isA = true;
						int pos_AB = i;
						unsigned short kuhn_mon = 0;
						crosslink next = build_strand(origin, isA, pos_AB, kuhn_mon);
						if (next.equals(origin))
						{
							strand_ct += 0.5;
							primary_junct.insert(origin.get_cross());
							all_molecules[m].contains_primary(true);
						}
						else if (next.get_external() == 1)
						{
							strand_ct += 1.0;
						}
						else
						{
							strand_ct += 0.5;
						}
					}
				}
				for (int i = 0; i < origin.get_B_con().size(); i++)
				{
					if (origin.get_B_con()[i] != -1)
					{
						bool isA = false;
						int index_AB = i;
						unsigned short kuhn_mon = 0;
						crosslink next = build_strand(origin, isA, index_AB, kuhn_mon);
						if (next.equals(origin))
						{
							strand_ct += 0.5;
							primary_junct.insert(origin.get_cross());
							all_molecules[m].contains_primary(true);
						}
						else if (next.get_external() == 1)
						{
							strand_ct += 1.0;
						}
						else
						{
							strand_ct += 0.5;
						}
					}
				}
			}
			strand_ct += all_molecules[m].get_primary();
		}
		v0 = strand_ct / V / (6.022 * pow(10, 23));
		loop_count[0] = primary_junct.size();
		return num_junctions;
	}
	void secondary()
	{
		std::set <int> secondary_junct;
		for (int m = 0; m < all_molecules.size(); m++)
		{
			crosslink origin = all_molecules[m];
			if (origin.get_external() >= 3)
			{
				vector <int> next_cross;
				for (int A = 0; A < origin.get_A_con().size(); A++)
				{
					if (origin.get_A_con()[A] != -1)
					{
						bool isA = true;
						int pos_AB = A;
						unsigned short kuhn_mon = 0;
						crosslink next = build_strand(origin, isA, pos_AB, kuhn_mon);
						if (next.get_external() >= 3)
						{
							next_cross.push_back(next.get_cross());
						}
					}
				}
				for (int B = 0; B < origin.get_B_con().size(); B++)
				{
					if (origin.get_B_con()[B] != -1)
					{
						bool isA = false;
						int pos_AB = B;
						unsigned short kuhn_mon = 0;
						crosslink next = build_strand(origin, isA, pos_AB, kuhn_mon);
						if (next.get_external() >= 3)
						{
							next_cross.push_back(next.get_cross());
						}
					}
				}
				std::map<int, int> cross_to_number;
				for (int i : next_cross)
				{
					cross_to_number[i] += 1;
				}
				for (auto& dup : cross_to_number)
				{
					bool secondary_complete = dup.first != origin.get_cross() && dup.second >= 2;
					if (secondary_complete)
					{
						secondary_junct.insert(origin.get_cross());
						secondary_junct.insert(dup.first);
					}
				}
			}
		}
		loop_count[1] = secondary_junct.size();
	}
	crosslink build_strand(crosslink c, bool& startsfromA, int& pos_AB, unsigned short& kuhn_mon)
	{
		add_kuhn_length(c, startsfromA, pos_AB, kuhn_mon);
		crosslink k = AB_index_to_crosslink_index(c, startsfromA, pos_AB);
		bool terminate = k.equals(c) || k.get_external() >= 3 || k.get_external() == 1;
		bool landsOnA = !startsfromA;
		add_kuhn_length(k, landsOnA, pos_AB, kuhn_mon);
		while (!terminate)
		{
			vector <int> bf = k.bifunctional_linker(landsOnA, pos_AB);
			pos_AB = bf[0];
			startsfromA = bf[1];
			add_kuhn_length(k, startsfromA, pos_AB, kuhn_mon);
			k = AB_index_to_crosslink_index(k, startsfromA, pos_AB);
			landsOnA = !startsfromA;
			add_kuhn_length(k, landsOnA, pos_AB, kuhn_mon);
			terminate = k.equals(c) || k.get_external() >= 3 || k.get_external() == 1;
		}
		startsfromA = landsOnA;
		return k;
	}
	crosslink AB_index_to_crosslink_index(crosslink c, bool startsfromA, int& pos_AB)
	{
		int index_DAB, cross;
		vector <int> crosslink_index;
		if (startsfromA)
		{
			index_DAB = c.get_A_con()[pos_AB];
			cross = crosslink_index_B[index_DAB];
			crosslink_index = crosslink_index_B;
		}
		else
		{
			index_DAB = c.get_B_con()[pos_AB];
			cross = crosslink_index_A[index_DAB];
			crosslink_index = crosslink_index_A;
		}
		crosslink k = crosslink(0);
		for (int i = 0; i < all_molecules.size(); i++)
		{
			if (all_molecules[i].get_cross() == cross)
			{
				k = all_molecules[i];
				break;
			}
			else if (i == all_molecules.size() - 1)
			{
				cout << "Error! Could not Find crosslinker in network." << endl;
				cout << k.get_type() << endl;
			}
		}
		int counter = 0;
		while (index_DAB >= 1)
		{
			index_DAB -= 1;
			if (crosslink_index[index_DAB] == cross)
			{
				counter += 1;
			}
			else
			{
				pos_AB = counter;
				break;
			}
		}
		if (index_DAB == 0)
		{
			pos_AB = counter;
		}
		return k;
	}
	void add_kuhn_length(crosslink c, bool A, int pos_AB, unsigned short& kuhn_mon)
	{
		vector <vector <double>> type = N_template[c.get_type() - 1];
		if (A)
		{
			kuhn_mon += (unsigned short)type[0][pos_AB];
		}
		else
		{
			kuhn_mon += (unsigned short)type[1][pos_AB];
		}
	}

	double get_vo()
	{
		return v0;
	}
	double get_strand_ct()
	{
		return strand_ct;
	}
	vector <double> get_loop_count()
	{
		return loop_count;
	}
};

//Initializes random number generator with seed
//RNG is Mersenne Twister MT19937 algorithm
void RNG_initialize(unsigned long seed)
{
	mt[0] = seed & 0xffffffffUL;
	for (mti = 1; mti < 624; mti++)
	{
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
	double dn = 3.442619855899;
	int i;
	const double m1 = 2147483648.0;
	double q;
	double tn = 3.442619855899;
	const double vn = 9.91256303526217E-03;

	q = vn / exp(-0.5 * dn * dn);

	kn[0] = (dn / q) * m1;
	kn[1] = 0;

	wn[0] = q / m1;
	wn[127] = dn / m1;

	fn[0] = 1.0;
	fn[127] = exp(-0.5 * dn * dn);

	for (i = 126; i >= 1; i--)
	{
		dn = sqrt(-2 * log(vn / dn + exp(-0.5 * dn * dn)));
		kn[i + 1] = (dn / tn) * m1;
		tn = dn;
		fn[i] = exp(-0.5 * dn * dn);
		wn[i] = dn / m1;
	}
}
//Returns a random long between 0 and 4294967295
unsigned long rand32()
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= 624)
	{
		int kk;

		for (kk = 0; kk < 227; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + 397] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk < 623; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk - 227] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[623] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}
#pragma once
