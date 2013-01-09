#ifndef __DATABASE_H__
#define __DATABASE_H__


#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

struct in {
	double id;
	double ra;
	double dec;
	double mag;
};
struct dat1 {
	double id;
	double x;
	double y;
	double z;
	double mag;
};
struct dat2 {
	double id1;
	double id2;
	double id3;
	double alpha1;
	double alpha2;
	double beta;
};
struct coord {
	double x;
	double y;
};

struct coord_3d {
	double x;
	double y;
	double z;
};
struct angle_triple {
	double alpha1;
	double alpha2;
	double beta;
};
struct coordint {
	int x;
	int y;
};
struct pixel {
	unsigned char r;
	unsigned char g;
	unsigned char b;
};
struct quaternion {
	double q0;
	double q1;
	double q2;
	double q3;
};
struct vec3d {
	double x;
	double y;
	double z;
};

class Database
{
public:
	Database(string inputfile, string outputfile_table1, string outputfile_table2, double foc_len, double pix_size);
	~Database();

	//sucht alpha1, alpha2, beta in tabelle 2
	bool find_triple (angle_triple key, dat2 *result);
	//gibt koordinaten zu gesuchten stern zurück
	vec3d get_coords (double id);
	
	
private:
	vector<in> input_data;
	vector<dat1> table1_data;
	vector<dat2> table2_data;
	double focal_length;
	double pixel_size;

	void read_inputfile(string filename);
	void generate_table1(string filename);
	void generate_table2(string filename);
};

#endif

