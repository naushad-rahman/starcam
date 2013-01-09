// starcam.cpp : Defines the entry point for the console application.
// 
// LUFT UND RAUMFAHRTDYNAMIK 
// WS 12/13
// Eric Reinthal, Pouyan Azari, Serkan Dikmen
// Dezember 2012
// Version 1.0

#include "database.h"
using namespace std;

#define IMAGEFILE "bild.bmp"
#define XDIM 864
#define YDIM 640
#define NOISE_THRESHOLD 20
#define FOCAL_LENGTH 2.5e-2
#define PIXEL_SIZE 5.8e-6
#define DATABASE_INPUTFILE "hip_red_1.txt"
#define DATABASE_OUT_TABLE1 "table1.txt"
#define DATABASE_OUT_TABLE2 "table2.txt"






//2d array im speicher anlegen
//---rgb-bild
pixel** allocrgb (int x, int y) {
    pixel **img = new pixel*[y];
    for (int i = 0; i < y; i++) {
		img[i] = new pixel[x];
    }
    return img;
}
//---s/w-bild
int** allocsw (int x, int y) {
    int **img = new int*[y];
    for (int i = 0; i < y; i++) {
		img[i] = new int[x];
    }
    return img;
}
//speicher freigeben
void deletergb(pixel **map, int y) {
    for (int i=0; i<y; i++)
		delete map[i];
    delete map;
}
void deletesw(int **map, int y) {
    for (int i=0; i<y; i++)
		delete map[i];
    delete map;
}
void deletematrix(double **m, int y) {
    for (int i=0; i<y; i++)
		delete m[i];
    delete m;
}

//l�nge eines vektors berechnen
double v_length(vec3d v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
//kreuzprodukt von 2 vektoren
vec3d cross_product(vec3d v1, vec3d v2) {
	vec3d ret;
	ret.x = v1.y*v2.z - v1.z*v2.y;
	ret.y = v1.z*v2.x - v1.x*v2.z;
	ret.z = v1.x*v2.y - v1.y*v2.x;
	return ret;
}
//vektor durch skalar teilen
vec3d div_vec (vec3d v, double scalar) {
	v.x = v.x/scalar;
	v.y = v.y/scalar;
	v.z = v.z/scalar;
	return v;
}
//skalarprodukt von 2 vektoren
double dot_product(vec3d v1, vec3d v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}


//bild einlesen und daten speichern
pixel** read_img(vector<unsigned char>& head) {
	cout << "Lese Bild ein..." << endl;
	//bild �ffnen
	ifstream img;
	img.open(IMAGEFILE, fstream::binary);
	//begnn der pixel
	int start = 0;
	for(int a=0; a<10; a++)
    {
        img.get();
    }
	start = img.get();
	//header speichern
	img.seekg(0, ios_base::beg);
	for (int i=0; i<start; i++)
		head.push_back(img.get());
	//pixel speichern
	pixel** pix = allocrgb(XDIM, YDIM);
	
	for (int y=0; y<YDIM; y++)
		for (int x=0; x<XDIM; x++) {
			pix[y][x].r = img.get();
			pix[y][x].g = img.get();
			pix[y][x].b = img.get();
		}
	
	img.close();
	return pix;
}

//bild in graustufen umwandeln
int** img_to_greyscale(pixel** a) {
	int** b = allocsw(XDIM, YDIM);
	int color = 0;
	for (int y=0; y<YDIM; y++)
		for (int x=0; x<XDIM; x++) {
			color = a[y][x].r + a[y][x].g + a[y][x].b;
			b[y][x] = color/3;
		}
	return b;
}
//bild erzeugen aus pixeldaten
void write_img(int** data, vector<unsigned char> head, string filename) {
	ofstream o;
	o.open(filename, fstream::binary);
	//head schreiben
	for (vector<unsigned char>::iterator it = head.begin(); it != head.end(); ++it)
		o.put(*it);
	//pixel schreiben
	for (int y=0; y<YDIM; y++)
		for (int x=0; x<XDIM; x++) {
			o.put(data[y][x]);
			o.put(data[y][x]);
			o.put(data[y][x]);
		}
	o.close();
}

//pr�ft, ob eine operation auf ein bild im Bereich der Pixel geschieht
bool valid_range(int x, int y) {
	if (x>XDIM || x<0 || y>YDIM || y<0)
		return false;
	return true;
}

//schwerpunkte der sterne mit hilfe der erweiterten ROI-Methode herausfinden
vector<coord> star_emphasises(int** data, vector<unsigned char> head) {
	cout << "Berechne Schwerpunkte der Sterne mit Advanced ROI..." << endl;
	//alle schwerpunkte
	vector<coord> e;
	//rauschen entfernen und maximale helligkeit herausfinden
	int bright_min = NOISE_THRESHOLD;
	int bright_max = 0;
	for (int y=0; y<YDIM; y++)
		for (int x=0; x<XDIM; x++) {
			if (data[y][x] < bright_min)
				data[y][x] = 0;
			if (data[y][x] > bright_max)
				bright_max = data[y][x];
		}

	//ROIs herausfinden

	//map mit den ROIs
	int** roi_map = allocsw(XDIM, YDIM);
	for (int y=0; y<YDIM; y++)
		for (int x=0; x<XDIM; x++)
			roi_map[y][x] = 0;

	int roundmax; //helligkeit des aktuellen pixels
	int count = 0; //counter f�r die sterne
	int roundmaxnew; //evtl gr��ere helligkeit eines Nachbarpixels
	double max_rel; //relative helligkeit des pixels (im bezug auf hellsten)
	coordint max; //Hauptpixel des aktuellen ROI-Feldes
	int roisize; //Gr��e des ROI (Abh von max_rel)
	coord ecoord; //berechnete Schwerpunktkoordinate
	double ecoord_num_x; //Hilfsvariablen zur Schwerpunktberechnung
	double ecoord_denom_x;
	double ecoord_num_y;
	double ecoord_denom_y;
	for (int y=0; y<YDIM; y++)
			for (int x=0; x<XDIM; x++) {
				//nur in liste der sterne aufnehmen, falls aktueller pixel noch nicht in einem roi liegt
				if (data[y][x] != 0 && roi_map[y][x] == 0) {
					//lokal hellsten pixel herausfinden
					// dazu  im kreis un akt. pixel nach helleren suchen
					
					roundmaxnew=data[y][x];
					max.x=x;
					max.y=y;
					do {
						roundmax = roundmaxnew;
						for (int k=-3; k<=3; k++)
							for (int l=-3; l<=3; l++)
								if (valid_range(max.x+l, max.y+k) && data[max.y+k][max.x+l] > roundmaxnew && (k!=0 || l!=0)) {
									roundmaxnew = data[max.y+k][max.x+l];
									max.x = max.x+l;
									max.y = max.y+k;
								}
					}while (roundmaxnew != roundmax);
					//relative helligkeit des sterns
					max_rel = (double)(roundmaxnew-bright_min)/(bright_max-bright_min);

					//punkte je nach helligkeit in der ROI-Map markieren und schwerpunkt berechnen mit hilfe der Formel aus der Vorlesung f�r Advanced ROI
					if (max_rel >= 3.0/4.0)
						roisize=9;
					else if (max_rel >= 1.1/2.0)
						roisize=6;
					else if (max_rel >= 1.0/4.0)
						roisize=4;
					else
						roisize=2;

					ecoord_num_x = 0;
					ecoord_denom_x = 0;
					ecoord_denom_y = 0;
					ecoord_num_y = 0;

					for (int k=-roisize; k<=roisize; k++)
						for (int m=-roisize; m<=roisize; m++)
							if (valid_range(x+m, y+k)) {
								roi_map[y+k][x+m] = 255;
								ecoord_num_x += (double) data[y+k][x+m] * (x+m);
								ecoord_num_y += (double) data[y+k][x+m] * (y+k);
								ecoord_denom_x += (double) data[y+k][x+m];
								ecoord_denom_y += (double) data[y+k][x+m];
							}

					ecoord.x = ecoord_num_x/ecoord_denom_x;
					ecoord.y = ecoord_num_y/ecoord_denom_y;

					//punkt in vector der Schwerpunktliste einf�gen
					e.push_back(ecoord);
					//z�hler der sterne
					count++;
					//cout << count << ": " << ecoord.x << " / " << ecoord.y << endl;
				}
			}

	write_img(roi_map,head,"roi.bmp");
	return e;
}

//zentralen stern herausfinden und diesen aus der starlist l�schen
coord central_star (vector<coord>* starlist) {
	cout << "Finde zentralen Stern heraus..." << endl;
	coord c;
	//abstand aller sterne zum Mittelpunkt finden
	coord m;
	double distance = 0.0;
	double mindistance = -1.0;
	m.x=XDIM/2;
	m.y=YDIM/2;
	//iterator zwischenspeichern, um stern sp�ter zu l�schen aus der liste
	vector<coord>::iterator itmin;
	for (vector<coord>::iterator it = (*starlist).begin(); it != (*starlist).end(); it++) {
		distance = sqrt(pow((*it).x-m.x,2) + pow((*it).y-m.y,2));
		if (distance < mindistance || mindistance < 0) {
			mindistance = distance;
			itmin=it;
			c.x = (*it).x;
			c.y = (*it).y;
		}
	}
	(*starlist).erase(itmin);
	cout << "Zentraler Stern bei (x/y): " << c.x << " / " << c.y << endl;
	return c;
}

//alle alpha_i zwischen Mittelpunktstern und den Sternen der Liste berechnen (in rad)
vector<double> calc_alpha_i (vector<coord> starlist, coord central_star) {
	cout << "Berechne Alpha_i ..." << endl;
	vector<double> alpha; //r�ckgabevektor
	double angle; //einzelne winkel
	double dist_st1_st2; //abstand stern zu Mittelpunktstern
	double dist_st1_cam; //abstand Mittelpunktstern zu Camera
	double dist_st2_cam; //abstand Stern zu Camera
	for (vector<coord>::iterator it = starlist.begin(); it != starlist.end(); it++) {
		dist_st1_st2 = sqrt(pow( (central_star.x-(*it).x)*PIXEL_SIZE,2) + pow((central_star.y-(*it).y)*PIXEL_SIZE,2) );
		dist_st1_cam = sqrt(pow( central_star.x*PIXEL_SIZE,2) + pow( central_star.y*PIXEL_SIZE,2) + pow( FOCAL_LENGTH,2) );
		dist_st2_cam = sqrt(pow( (*it).x*PIXEL_SIZE,2) + pow( (*it).y*PIXEL_SIZE,2) + pow( FOCAL_LENGTH,2) );
		//winkelberechnung mit Kosinussatz
		angle = acos( (pow(dist_st1_st2,2) - pow(dist_st1_cam,2) - pow(dist_st2_cam,2)) / (-2*dist_st1_cam*dist_st2_cam) );
		alpha.push_back(angle);
	}
	return alpha;
}

//finde 2 kleinste alpha_i und berechne beta
angle_triple get_angles (vector<coord> starlist, vector<double> alphalist, coord central_star, coord* id2_coord) {
	cout << "Berechne 3 Winkel aus dem Bild (in rad)..." << endl;
	angle_triple res;
	//alphalist nach 2 kleinsten winkeln durchsuchen und zugeh�rige koordinaten speichern
	coord star1;
	coord star2;
	double alpha1=1000;
	double alpha2=1000;
	if (starlist.size() != alphalist.size())
		cout << "FEHLER: starlist.size() != alphalist.size()" << endl;

	vector<double>::iterator italpha = alphalist.begin();
	vector<coord>::iterator itstar = starlist.begin();

	//finde 2 kleinste alpha 
	while(italpha != alphalist.end()) {
		if ((*italpha) < alpha1) {
			alpha1 = *italpha;
			star1.x = (*itstar).x;
			star1.y = (*itstar).y;
		} else if ((*italpha) < alpha2) {
			alpha2 = *italpha;
			star2.x = (*itstar).x;
			star2.y = (*itstar).y;
		}
		itstar++;
		italpha++;
	}
	res.alpha1 = alpha1;
	res.alpha2 = alpha2;
	
	//berechne beta
	double dist_c_1; //winkel zwischen zentralstern und stern unter alpha1
	double dist_c_2; //winkel zwischen zentralstern und stern unter alpha2
	double dist_1_2; //winkel zwischen stern unter alpha2 und stern unter alpha1
	dist_c_1 = sqrt( pow(central_star.x-star1.x,2) + pow(central_star.y-star1.y,2) );
	dist_c_2 = sqrt( pow(central_star.x-star2.x,2) + pow(central_star.y-star2.y,2) );
	dist_1_2 = sqrt( pow(star2.x-star1.x,2) + pow(star2.y-star1.y,2) );
	//betaberechnung mit kosinussatz
	res.beta = acos((pow(dist_1_2,2)-pow(dist_c_1,2)-pow(dist_c_2,2)) / (-2*dist_c_1*dist_c_2) );
	*id2_coord = star1;
	return res;
}

//gibt sternenvektoren zur�ck (2 aus Datenbank (IDs �bergeben) und die 2 entsprechenden aus dem Bild (2 coordinaten �bergeben))
//vektoren v:bild w:db
void generate_starvector (vec3d* v1, vec3d* v2, vec3d* w1, vec3d* w2, double id1, double id2, coord c1, coord c2, Database* dbhandle) {
	//vektoren aus der DB holen
	*w1 = dbhandle->get_coords(id1);
	*w2 = dbhandle->get_coords(id2);
	//vektoren aus bild berechnen
	vec3d b1;
	vec3d b2;
	//z=brennweite
	b1.z = FOCAL_LENGTH;
	b2.z = b1.z;
	//x und y aus pixelentfernung zum mittelpunkt berechnen;
	coord center;
	center.x = XDIM/2;
	center.y = YDIM/2;
	//abstand = pxielentfernung*pixelgr��e
	b1.x = abs(center.x-c1.x)*PIXEL_SIZE;
	b2.x = abs(center.x-c2.x)*PIXEL_SIZE;
	b1.y = abs(center.y-c1.y)*PIXEL_SIZE;
	b2.y = abs(center.y-c2.y)*PIXEL_SIZE;

	*v1 = b1;
	*v2 = b2;
}

//rotationsmatrix aus 2 koordinaten berechnen
double** rot_matrix (vec3d v1, vec3d v2, vec3d w1, vec3d w2) {
	//3x3 matrix erstellen
	double** ret = new double*[3];
	for (int i=0; i<3; i++)
		ret[i] = new double[3];
	//transformation
	//sensorframe
	vec3d r1 = div_vec(v1, v_length(v1));
	vec3d r2 = div_vec(cross_product(v1,v2), v_length(cross_product(v1,v2)));
	vec3d r3 = cross_product(r2, r1);
	//irf frame
	vec3d t1 = div_vec(w1, v_length(w1));
	vec3d t2 = div_vec(cross_product(w1,w2), v_length(cross_product(w1,w2)));
	vec3d t3 = cross_product(t2, t1);

	//transpornierte vektoren von t;
	vec3d s1;
	vec3d s2;
	vec3d s3;

	s1.x = t1.x;
	s1.y = t2.x;
	s1.z = t3.x;

	s2.x = t1.y;
	s2.y = t2.y;
	s2.z = t3.y;

	s3.x = t1.z;
	s3.y = t2.z;
	s3.z = t3.z;

	//transformationsmatrix aufstellen
	ret[0][0] = dot_product(r1,s1);
	ret[0][1] = dot_product(r1,s2);
	ret[0][2] = dot_product(r1,s3);

	ret[1][0] = dot_product(r2,s1);
	ret[1][1] = dot_product(r2,s2);
	ret[1][2] = dot_product(r2,s3);

	ret[2][0] = dot_product(r3,s1);
	ret[2][1] = dot_product(r3,s2);
	ret[2][2] = dot_product(r3,s3);

	return ret;
}

//quaternionen aus einer transformationsmatrix berechnen
quaternion matrix_to_quaternion(double** a) {
	quaternion q;

	//spur der matrix berechnen
	double trace = a[0][0] + a[1][1] + a[2][2];

	//quaternionen
	q.q0 = abs(0.5*sqrt(1+trace));
	q.q1 = (0.5/q.q0)*(a[1][2]-a[2][1]);
	q.q2 = (0.5/q.q0)*(a[2][0]-a[0][2]);
	q.q3 = (0.5/q.q0)*(a[0][1]-a[1][0]);

	return q;
}

int main(int argc, char* argv[])
{

	//bild einlesen und header und pixelinfo speichern
	vector<unsigned char> image_header;
	pixel** image_data = read_img(image_header);
	int** image_grey = img_to_greyscale(image_data);

	 
	//schwerpunkt der sterne herausfinden und in liste speichern
	vector<coord> star_emphasis_list = star_emphasises(image_grey, image_header);
	//zentralen stern heruasfinden und diesen stern aus star_emphasis l�schen
	coord central = central_star(&star_emphasis_list);

	//alle winkel alpha_i im bild berechnen
	vector<double> alphalist = calc_alpha_i(star_emphasis_list, central);

	//alpha1, alpha2 und beta aus dem Bild
	coord coord2;
	angle_triple img_angles = get_angles(star_emphasis_list, alphalist, central, &coord2);

	cout << "Alpha1: " << img_angles.alpha1 << endl;
	cout << "Alpha2: " << img_angles.alpha2 << endl;
	cout << "Beta: " << img_angles.beta << endl;

	//textfile einlesen und beide Tabellen berechnen
	Database *catalog;
	catalog = new Database(DATABASE_INPUTFILE, DATABASE_OUT_TABLE1, DATABASE_OUT_TABLE2, FOCAL_LENGTH, PIXEL_SIZE);

	//nach winkel in tabelle 2 suchen und IDs ausgeben
	dat2 result;
	bool found = catalog->find_triple(img_angles, &result);

	if (!found) {
		cout << "Keine passende Sternenkonstellation in der DB gefunden!" << endl;
	} else {
		cout << "Passende IDs wurden in der DB gefunden. Gefundene Winkel:" << endl;
		cout << "Alpha1 (gemessen / DB): " << img_angles.alpha1 << " / " << result.alpha1 << endl;
		cout << "Alpha2 (gemessen / DB): " << img_angles.alpha2 << " / " << result.alpha2 << endl;
		cout << "Beta (gemessen / DB): " << img_angles.beta << " / " << result.beta << endl << endl;
		cout << "Gefundene Sterne:" << endl;
		cout << "ID1: " << result.id1 << endl;
		cout << "ID2: " << result.id2 << endl;
		cout << "ID3: " << result.id3 << endl;

		//vektoren zu 2 sternen aus bild und aus datenbank generieren
		vec3d v1;
		vec3d v2;
		vec3d w1;
		vec3d w2;
		generate_starvector(&v1, &v2, &w1, &w2, result.id1, result.id2, central, coord2, catalog);
		cout << endl;
		cout << "Vektoren (v: Sensorframe | w: IRF-Frame)(x y z):" << endl;
		cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << endl;
		cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << endl;
		cout << "w1: " << w1.x << " " << w1.y << " " << w1.z << endl;
		cout << "w2: " << w2.x << " " << w2.y << " " << w2.z << endl << endl;

	
	
		//rotationsmatrix berechnen
		double** rot = rot_matrix(v1, v2, w1, w2);
		cout << "Rotationsmatrix A: " << endl;
		cout << rot[0][0] << "\t" << rot[0][1] << "\t" << rot[0][2] << endl;
		cout << rot[1][0] << "\t" << rot[1][1] << "\t" << rot[1][2] << endl;
		cout << rot[2][0] << "\t\t" << rot[2][1] << "\t" << rot[2][2] << endl << endl;

		//quaternionen aus der Rotationsmatrix berechnen
		quaternion q = matrix_to_quaternion(rot);
		cout << "Quaternionen: " << endl;
		cout << "Q0:  " << q.q0 << endl;
		cout << "Q1:  " << q.q1 << endl;
		cout << "Q2:  " << q.q2 << endl;
		cout << "Q3:  " << q.q3 << endl;


		//speicher freigeben
		deletematrix(rot, 3);
	}

	//bild mit schwellenfilter schreiben
	write_img(image_grey, image_header, "out.bmp");
	//speicher freigeben
	deletesw(image_grey, YDIM);
	deletergb(image_data, YDIM);

	//system("PAUSE");
	return 0;
}

