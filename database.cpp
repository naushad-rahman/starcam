#include "database.h"


//erzeugt die Tabellen
Database::Database(string inputfile, string outputfile_table1, string outputfile_table2, double foc_len, double pix_size)
{
	focal_length = foc_len;
	pixel_size = pix_size;
	
	//falls datenbanken noch nicht erstellt, diese erzeugen
	//andernfalls textdateien einlesen
	ifstream tab1;
	tab1.open(outputfile_table1);
	if (tab1.good()) {
		cout << "Lese Daten aus Tabelle 1 ein..." << endl;
		//bereits generierte daten der tabelle einlesen
		dat1 data;
		while (tab1.good()) {
			tab1 >> data.id >> data.x >> data.y >> data.z >> data.mag;
			table1_data.push_back(data);
		}
		tab1.close();
	} else {
		tab1.close();
		read_inputfile(inputfile);
		generate_table1 (outputfile_table1);
	}
	
	ifstream tab2;
	tab2.open(outputfile_table2);
	if (tab2.good()) {
		cout << "Lese Daten aus Tabelle 2 ein..." << endl;
		//bereits generierte daten der tabelle einlesen
		dat2 data2;
		while (tab2.good()) {
			tab2 >> data2.id1 >> data2.id2 >> data2.id3 >> data2.alpha1 >> data2.alpha2 >> data2.beta;
			table2_data.push_back(data2);
		}
		tab2.close();
	} else {
		tab2.close();
		generate_table2 (outputfile_table2);
	}
}

Database::~Database()
{
}

// liest die Daten aus der Datei und schreibt sie in eine Vector
void Database::read_inputfile(string filename) {
	cout << "lese Daten aus Datei: " << filename << " ..." << endl;
	ifstream file;
	file.open (filename);
	//inhalt der datei
	double ra;
	double dec;
	double mag;
	//ausgabe
	in data;
	int counter = 1;
	while(file) {
		file >> ra >> dec >> mag;
		//nur in DB aufnehmen, falls magnitude unter schwennelwert
		if (mag < 6.5) {
			data.id = counter;
			data.dec = dec;
			data.ra = ra;
			data.mag = mag;
			input_data.push_back(data);
			counter++;
		}
	}

	file.close();
}

//rohdaten der gegebenen datenbank in tabelle1 schreiebn (Id, x, y, z, mag)
void Database::generate_table1(string filename) {
	cout << "Erstelle Tabelle 1(ID, X, Y, Z, Mag): " << filename << endl;
	//variablen für daten der tabelle
	double x;
	double y;
	double z;
	dat1 data;
	
	//datei für tabelle 1
	ofstream out;
	out.open(filename);

	for (vector<in>::iterator it = input_data.begin(); it != input_data.end(); it++) {
		x = cos((*it).dec) * cos ((*it).ra);
		y = cos((*it).dec) * sin ((*it).ra);
		z = (*it).mag * sin((*it).dec);
		data.id = (*it).id;
		data.x = x;
		data.y = y;
		data.z = z;
		data.mag = (*it).mag;
		//zeile in vector aufnehmen
		table1_data.push_back(data);
		//zeile in datei aufnehmen
		out << data.id << " " << x << " " << y << " " << z << " " << data.mag << " ";
	}
	out.close();
}

//aus den Daten von Tabelle 1 wird hier Tabelle 2 (id1, id2, id3, a1, a2, beta) generiert
void Database::generate_table2(string filename) {
	cout << "Erstelle Tabelle 2(ID1, ID2, ID3, Alpha1, Alpha2, Beta): " << filename << endl;
	//variablen für daten der neuen tabelle
	double a1;
	double a2;
	double beta;
	dat2 data;
	//sternvektoren
	double dist_it; //vektor von cam zum iterierenden stern
	double dist_s1; //vektor von cam zu nachbarsternen des it. stern
	double dist_it_s1; //vektor zwischen it. stern und nachbarsternen

	//datei für tabelle 2
	ofstream out;
	out.open(filename);

	// variablen für zwei kleinste alphas zu jedem stern
	double alpha; //alpha zu jedem stern
	double dist_a1_it; //vektor vom iterierenden stern zu stern unter a1
	double dist_a2_it; //vektor vom iterierenden stern zu stern unter a2
	double dist_a1_a2_squared; //vektor vom stern unter a1 zum stern unter a2 zum quadrat
	coord_3d c1; //koordinaten zum stern unter a1
	coord_3d c2; //koordinaten zum stern unter a2
	dat1 data2; //daten zum stern unter a1
	dat1 data3; //daten zum stern unter a2

	for (vector<dat1>::iterator it = table1_data.begin(); it!=table1_data.end(); it++) {
		//die 2 kleinsten alpha_i zu diesem Stern speichern
		a1 = 1000;
		a2 = 1000;
		dist_it = sqrt (pow((*it).x,2) + pow((*it).y,2) + pow((*it).z,2));
		//alle nachbarsterne durchlaufen
		for (vector<dat1>::iterator it_n = it+1; it_n!=table1_data.end(); it_n++) {
			
			dist_s1 = sqrt (pow((*it_n).x,2) + pow((*it_n).y,2) + pow((*it_n).z,2));
			dist_it_s1 = sqrt (pow((*it_n).x-(*it).x,2) + pow((*it_n).y-(*it).y,2) + pow((*it_n).z-(*it).z,2));
			//alpha berechnen
			alpha = acos((dist_it_s1*dist_it_s1-dist_it*dist_it-dist_s1*dist_s1) / (-2*dist_it*dist_s1));
			//falls kleinstes alpha, dann in liste aufnehmen
			if (alpha < a1) {
				a1 = alpha;
				data2 = (*it_n);
				dist_a1_it = dist_it_s1;
				c1.x = (*it_n).x;
				c1.y = (*it_n).y;
				c1.z = (*it_n).z;
			} else if (alpha < a2) {
				a2 = alpha;
				data3 = (*it_n);
				dist_a2_it = dist_it_s1;
				c2.x = (*it_n).x;
				c2.y = (*it_n).y;
				c2.z = (*it_n).z;
			}
			
		}

		//beta berechnen mit kosinussatz
		dist_a1_a2_squared = (c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z);
		beta = acos((dist_a1_a2_squared-dist_a1_it*dist_a1_it-dist_a2_it*dist_a2_it) / (-2*dist_a1_it*dist_a2_it));

		//zwei kleinsten nachbarn zu diesem stern in datenbank aufnehmen
		data.id1 = (*it).id;
		data.id2 = data2.id;
		data.id3 = data3.id;
		data.alpha1 = a1;
		data.alpha2 = a2;
		data.beta = beta;
		//daten in vector übernehmen
		table2_data.push_back(data);
		//daten in datei schreiben
		out << data.id1 << " " << data.id2 << " " << data.id3 << " " << data.alpha1 << " " << data.alpha2 << " " << data.beta << " ";
	}

	out.close();
}

bool Database::find_triple (angle_triple key, dat2 *result) {
	//durchsuche table2_data nach key
	(*result).id1=0;
	(*result).id2=0;
	(*result).id3=0;
	(*result).alpha1=0;
	(*result).alpha2=0;
	(*result).beta=0;

	double range = 0.01; //ungenauigkeitsbereich
	double variance; //quadratische abweichung vom aktuellen wert
	double variance_min = 1000; //quadratische abweichung vom bestwert (zum gesuchten wert)

	//db durchsuchen
	for (vector<dat2>::iterator it = table2_data.begin(); it != table2_data.end(); it++) {
		if ((*it).alpha1 < key.alpha1+range && (*it).alpha1 > key.alpha1-range)
			if ((*it).alpha2 < key.alpha2+range && (*it).alpha2 > key.alpha2-range)
				if ((*it).beta < key.beta+range && (*it).beta > key.beta-range) {
					//abweichung vom gesuchten wert berechnen
					variance = ((*it).alpha1-key.alpha1)*((*it).alpha1-key.alpha1)+((*it).alpha2-key.alpha2)*((*it).alpha2-key.alpha2)+((*it).beta-key.beta)*((*it).beta-key.beta);
					if (variance < variance_min) {
						*result = *it;
						variance_min = variance;
					}
				}
					
	}
	//falls in der datenbank was gefunden wurde, ist beta auf jeden fall != 0
	if ((*result).beta != 0)
		return true;
	return false;
}

vec3d Database::get_coords(double id) {
	vec3d ret;
	//id in tabelle1 suchen
	for (vector<dat1>::iterator it = table1_data.begin(); it!=table1_data.end(); it++) {
		if ((*it).id == id) {
			ret.x = (*it).x;
			ret.y = (*it).y;
			ret.z = (*it).z;
			break;
		}
	}
	
	return ret;
}

