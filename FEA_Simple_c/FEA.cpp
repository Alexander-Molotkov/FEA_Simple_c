#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>

using namespace std;

vector<double> cholesky_solve(vector<vector<double>> &a, vector<double>& b);
vector<double> forward_sub(vector<vector<double>> &l, vector<double> &b);
vector<double> backward_sub(vector<vector<double>>& u, vector<double>& y);
vector< vector<double> > cholesky_decomp(vector< vector<double> > &m1);
vector< vector<double> > matrix_transpose(vector < vector<double> >& m1);
vector<vector<double>> matrix_multiply(vector<vector<double>>& m1, vector<vector<double>>& m2);
vector<double> matrix_multiply(vector< vector<double> >& m1, vector< double>& m2);
vector< vector<double> > matrix_add(vector< vector <double> >& m1, vector< vector <double> >& m2);
vector<double> matrix_add(vector <double> & m1, vector <double> &m2);

void build_nodes(int sizeX, int sizeY, vector< vector<double> > &NODES);
void build_elems(int sizeX, int sizeY, vector< vector<double> > &ELEMS);
void build_supports(int sizeX, int sizeY, vector< vector<double> > &SUPPORTS);
void build_nodal_loads(int sizeX, int sizeY, vector< vector<double> > &NODALLOADS);
void build_support_disps(int sizeX, int sizeY, vector< vector<double> > &SUPPORTDISPS);

void build_local_basic_transform(vector< vector<double> >& abl, double L);

void debug(vector< vector<double> >);
void debug(vector< vector<int> >);
void debug(vector<int> v);
void debug(vector<double> v);

int main() {

	//Start Timer
	auto start = chrono::high_resolution_clock::now();

	const int NODE_X = 2;
	const int NODE_Y = 7;
	const int ELEMS_X = 7;
	const int ELEMS_Y = 10;
	const int SUPPORTS_X = 3;
	const int SUPPORTS_Y = NODE_Y;
	const int NODALLOADS_X = 3;
	const int NODALLOADS_Y = NODE_Y;
	const int SUPPORTDISPS_X = 3;
	const int SUPPORTDISPS_Y = NODE_Y;

	/*********************************************************************
	* The following functions are used to build the hardcoded input model
	*********************************************************************/

	printf("Calculating. . .\n");

	vector< vector<double> > NODES;
	for(int i = 0; i < NODE_Y; i++) {
		vector<double> temp;
		NODES.push_back(temp);
	}
	vector< vector<double> > ELEMS;
	for (int i = 0; i < ELEMS_Y; i++) {
		vector<double> temp;
		ELEMS.push_back(temp);
	}
	vector< vector<double> > SUPPORTS;
	for (int i = 0; i < SUPPORTS_Y; i++) {
		vector<double> temp;
		SUPPORTS.push_back(temp);
	}
	vector< vector<double> > NODALLOADS;
	for (int i = 0; i < NODALLOADS_Y; i++) {
		vector<double> temp;
		NODALLOADS.push_back(temp);
	}
	vector< vector<double> > SUPPORTDISPS;
	for (int i = 0; i < SUPPORTDISPS_Y; i++) {
		vector<double> temp;
		SUPPORTDISPS.push_back(temp);
	}

	build_nodes(NODES);
	//debug(NODES);
	build_elems(ELEMS);
	//debug(ELEMS);
	build_supports(SUPPORTS);
	//debug(SUPPORTS);
	build_nodal_loads(NODALLOADS);
	//debug(NODALLOADS);
	build_support_disps(SUPPORTDISPS);
	//debug(SUPPORTDISPS);

	/**************************************************************************
	*Everything below this line is generic and should work for any input above
	**************************************************************************/

	auto checkpoint = chrono::high_resolution_clock::now();
	auto checkpoint_start = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(checkpoint - start);
	cout << "Input model built in " << duration.count() << " milliseconds. \n";

	/***************************
	*	Assign equation numbers
	***************************/
	double nfdof = 0; 
	double ncdof = 0;

	//Create equation array
	vector< vector<double> > EQUATIONS;
	for (int i = 0; i < NODE_Y; i++) {

		vector<double> temp(3, 0);
		EQUATIONS.push_back(temp);
	}

	//Positive equation numbers are assigned to unconstrained DOFs
	//Negative equation numbers are assigned to constrained DOFs, i.e., where reactions where be developed
	for (int i = 0; i < NODE_Y; i++) {
		for (int j = 0; j < 3; j++) {

			if (SUPPORTS[i][j] == 0) {
				nfdof++;
				EQUATIONS[i][j] = nfdof;
			}
			else {
				ncdof++;
				EQUATIONS[i][j] = (ncdof * -1);
			}
		}
	}

	//debug(EQUATIONS);

	/**********************************************************************
	*	Create the joint load vector Pf and support displacement vector Uc
	**********************************************************************/

	vector <double> pf(nfdof, 0.0000);
	vector <double> uc(ncdof, 0.0000);

	for (int i = 0; i < NODE_Y; i++) {
		for (int j = 0; j < 3; j++) {

			if (EQUATIONS[i][j] > 0) { //Unconstrained DOF
				pf[EQUATIONS[i][j] - 1] = NODALLOADS[i][j];
			}
			else { //Constrained DOF
				uc[-EQUATIONS[i][j] - 1] = SUPPORTDISPS[i][j];
			}
		}
	}

	//Declaring nele, uf, kff
	int nele = ELEMS_Y;

	vector <double> uf(nfdof, 0);

	vector < vector<double> > kff;
	for (int i = 0; i < nfdof; i++) {
		vector<double> temp(nfdof, 0);
		kff.push_back(temp);
	}

	//Calculation
	for (int e = 0; e < nele; e++) {

		//I and J nodes for this element
		int ndI = (int)ELEMS[e][0];
		int ndJ = (int)ELEMS[e][1];

		//vector1.insert(vector1.end(), vector2.begin(), vector2.end());

		//Lookup vector
		vector <double> l;
		l = EQUATIONS[ndI - 1];
		l.insert(l.end(), EQUATIONS[ndJ - 1].begin(), EQUATIONS[ndJ - 1].end());

		//Difference in X-coordinates
		double DX = NODES[ndJ - 1][0] - NODES[ndI - 1][0];

		//Difference in y-coordinates
		double DY = NODES[ndJ - 1][1] - NODES[ndI - 1][1];

		//Element length and direction cosines
		double L = pow((pow(DX, 2) + pow(DY, 2)), 0.5);

		double c = DX / L;
		double s = DY / L;

		//Local-basic transformation
		vector< vector<double> > abl;
		for (int i = 0; i < 3; i++) {
			vector<double> temp;
			abl.push_back(temp);
		}
		build_local_basic_transform(abl, L);

		//Global-lobal transformation
		vector< vector<double> > alg;
		for (int i = 0; i < 6; i++) {
			vector<double> temp(6, 0);
			alg.push_back(temp);
		}
		alg[0][0] = c; alg[0][1] = s;
		alg[1][0] = -s; alg[1][1] = c;
		alg[2][2] = 1;
		alg[3][3] = c; alg[3][4] = s;
		alg[4][3] = -s; alg[4][4] = c;
		alg[5][5] = 1;

		//Select global displacements from Uf and Uc
		vector<double> u(6, 0);
		for (int i = 0; i < 6; i++) {
			if (l[i] > 0) {
				u[i] = uf[l[i] - 1];
			}
			else {
				u[i] = uc[-l[i] - 1];
			}
		}

		//Local Displacements
		vector<double> ul(6, 0);
		ul = matrix_multiply(alg, u);

		//Basic deformations
		vector<double> ub(3, 0);
		ub = matrix_multiply(abl, ul);

		//Element properties and member loads
		double E = ELEMS[e][2];
		double A = ELEMS[e][3];
		double I = ELEMS[e][4];
		double WX = ELEMS[e][5];
		double WY = ELEMS[e][6];

		//Basic Stiffness
		vector< vector<double> > kb;
		for (int i = 0; i < 3; i++) {
			vector<double> temp;
			kb.push_back(temp);
		}
		kb[0].push_back(E * A / L);
		kb[0].push_back(0);
		kb[0].push_back(0);
		kb[1].push_back(0);
		kb[1].push_back(4 * E * I / L);
		kb[1].push_back(2 * E * I / L);
		kb[2].push_back(0);
		kb[2].push_back(2 * E * I / L);
		kb[2].push_back(4 * E * I / L);

		//Local stiffness
		vector< vector<double> > x = matrix_transpose(abl);
		vector< vector<double> > y = matrix_multiply(kb, abl);
		vector< vector<double> > kl = matrix_multiply(x, y);

		//Global stiffness
		x = matrix_transpose(alg);
		y = matrix_multiply(kl, alg);
		vector< vector<double> > k = matrix_multiply(x, y);

		//Fixed-end basic forces
		vector<double> pb0;
		pb0.push_back(-WX * L / 2);
		pb0.push_back(-WY * pow(L, 2) / 12);
		pb0.push_back(WY * pow(L, 2) / 12);

		//Basic force-deformation relationship
		vector<double> pb = matrix_multiply(kb, ub);
		pb = matrix_add(pb, pb0);

		//"Reactions" due to member loads
		vector<double> plw;
		plw.push_back(-WX * L);
		plw.push_back(-WY * L / 2);
		plw.push_back(0);
		plw.push_back(0);
		plw.push_back(-WY * L / 2);
		plw.push_back(0);

		//Local forces
		x = matrix_transpose(abl);
		vector<double> pl = matrix_multiply(x, pb);
		pl = matrix_add(pl, plw);

		//Global forces
		x = matrix_transpose(alg);
		vector<double> p = matrix_multiply(x, pl);

		//Assemble Kff and Pf
		for (int j = 0; j < 6; j++){

			if (l[j] >= 0) {
				
				pf[ l[j] -1 ] -= p[j];

				for (int i = 0; i < 6; i++) {

					if (l[i] >= 0) {
						kff[ l[i]-1 ][ l[j]-1 ] += k[i][j];
					}
				}
			}
		}
	}

	/**************************************
	* Solve for the nodal displacements Uf
	**************************************/

	if (nfdof > 0) {
		uf = cholesky_solve(kff, pf);
	}

	/****************************
	* Assemble the reactions Pc
	****************************/

	checkpoint_start = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::microseconds>(checkpoint_start - start);
	cout << "Kff assembled, Pf modified, and Uf solved for in " << duration.count() << " milliseconds. ";

	vector<double> pc(ncdof, 0);

	for (int e = 0; e < nele; e++) {

		//I and J nodes for this element
		double ndI = ELEMS[e][0];
		double ndJ = ELEMS[e][1];

		//Lookup Vector
		vector<double> l;	
		l = EQUATIONS[ndI - 1];
		l.insert(l.end(), EQUATIONS[ndJ - 1].begin(), EQUATIONS[ndJ - 1].end());

		//Difference in X-coordinates
		double XI = NODES[ndI - 1][0];
		double XJ = NODES[ndJ - 1][0];
		double DX = XJ - XI;

		//Difference in Y-coordinates
		double YI = NODES[ndI - 1][1];
		double YJ = NODES[ndJ - 1][1];
		double DY = YJ - YI;

		//Element length and direction cosines
		double L = pow(pow(DX, 2) + pow(DY, 2), 0.5);
		double c = DX / L;
		double s = DY / L;

		//Local-basic transformation
		vector< vector<double> > abl;
		for (int i = 0; i < 3; i++) {
			vector<double> temp;
			abl.push_back(temp);
		}
		build_local_basic_transform(abl, L);

		//Global-lobal transformation
		vector< vector<double> > alg;
		for (int i = 0; i < 6; i++) {
			vector<double> temp(6, 0);
			alg.push_back(temp);
		}
		alg[0][0] = c; alg[0][1] = s;
		alg[1][0] = -s; alg[1][1] = c;
		alg[2][2] = 1;
		alg[3][3] = c; alg[3][4] = s;
		alg[4][3] = -s; alg[4][4] = c;
		alg[5][5] = 1;

		//Select global displacements from Uf and Uc
		vector<double> u(6, 0);
		for (int i = 0; i < 6; i++) {
			if (l[i] > 0) {
				u[i] = uf[l[i] - 1];
			}
			else {
				u[i] = uc[-l[i] - 1];
			}
		}

		//Local Displacements
		vector<double> ul(6, 0);
		ul = matrix_multiply(alg, u);

		//Basic deformations
		vector<double> ub(3, 0);
		ub = matrix_multiply(abl, ul);

		//Element properties and member loads
		double E = ELEMS[e][2];
		double A = ELEMS[e][3];
		double I = ELEMS[e][4];
		double WX = ELEMS[e][5];
		double WY = ELEMS[e][6];

		//Basic Stiffness
		vector< vector<double> > kb;
		for (int i = 0; i < 3; i++) {
			vector<double> temp;
			kb.push_back(temp);
		}
		kb[0].push_back(E* A / L);
		kb[0].push_back(0);
		kb[0].push_back(0);
		kb[1].push_back(0);
		kb[1].push_back(4 * E* I / L);
		kb[1].push_back(2 * E* I / L);
		kb[2].push_back(0);
		kb[2].push_back(2 * E* I / L);
		kb[2].push_back(4 * E* I / L);

		//Fixed-end basic forces
		vector<double> pb0;
		pb0.push_back(-WX * L / 2);
		pb0.push_back(-WY * pow(L, 2) / 12);
		pb0.push_back(WY* pow(L, 2) / 12);

		//Basic force-deformation relationship
		vector<double> pb = matrix_multiply(kb, ub);
		pb = matrix_add(pb, pb0);

		//"Reactions" due to member loads
		vector<double> plw;
		plw.push_back(-WX * L);
		plw.push_back(-WY * L / 2);
		plw.push_back(0);
		plw.push_back(0);
		plw.push_back(-WY * L / 2);
		plw.push_back(0);

		//Local forces
		vector< vector<double> > x = matrix_transpose(abl);
		vector<double> pl = matrix_multiply(x, pb);
		pl = matrix_add(pl, plw);

		//Global forces
		x = matrix_transpose(alg);
		vector<double> p = matrix_multiply(x, pl);

		//Assemble Pc
		for (int j = 0; j < 6; j++) {
			if( l[j] < 0 ){

				pc[-l[j] - 1] += p[j];
			}
		}
	}

	duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - checkpoint_start);
	cout << "\nPC assembled in " << duration.count() << " milliseconds.\n";

	return 0;
}

vector<double> cholesky_solve(vector< vector<double> > &a, vector<double> &b) {

	//Solve for Ax = b
	vector< vector<double> > l = cholesky_decomp(a);
	vector< vector<double> > u = matrix_transpose(l);

	vector<double> y = forward_sub(l, b);
	vector<double> x = backward_sub(u, y);

	return x;
}

vector<double> forward_sub(vector< vector<double> > &l, vector<double> &b) {
	
	vector<double> y(b.size(), 0);
	double s = 0;

	//On a lower-triangular matrix
	for (size_t i = 0; i < l.size(); i++) {

		double s = b[i];

		for (int j = 0; j < i; j++) {

			s -= l[i][j] * y[j];
		}
		y[i] = s / l[i][i];
	}
	return y;
}

vector<double> backward_sub(vector< vector<double> > &u, vector<double> &y) {

	vector<double> x(u.size(), 0);

	//On an upper-triangular matrix
	for (int i = u.size() -1; i >= 0; i--) {

		double s = y[i];

		for (size_t j = i; j < u.size(); j++) {

			s -= u[i][j] * x[j];
		}
		x[i] = s / u[i][i];
	}

	return x;
}

vector< vector<double> > cholesky_decomp(vector< vector<double> > &matrix) {

	int rows = matrix.size();
	int cols = matrix[0].size();

	if (rows != cols) {
		fprintf(stderr, "ERROR: You can only find the lower-triangular matrix of a square matrix.");
		exit(1);
	}

	vector< vector<double> > lower;
	for (int i = 0; i < rows; i++) {
		vector<double> temp(cols, 0);
		lower.push_back(temp);
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j <= i; j++) {

			double sum = 0;

			if (j == i) // summation for diagonals 
			{
				for (int k = 0; k < j; k++) {
					sum += pow(lower[j][k], 2);
				}
				lower[j][j] = sqrt(matrix[j][j] - sum);
			}
			else {

				// Evaluating L(i, j) using L(j, j) 
				for (int k = 0; k < j; k++) {
					sum += (lower[i][k] * lower[j][k]);
				}
				lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
			}
		}
	}
	return lower;
}

vector< vector<double> > matrix_transpose(vector< vector<double> > &m1) {

	int rows = m1.size();
	int cols = m1[0].size();

	vector< vector<double> > m2;
	for (int i = 0; i < cols; i++) {
		vector<double> temp(rows, 0);
		m2.push_back(temp);
	}

	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {

			m2[i][j] = m1[j][i];
		}
	}

	return m2;
}

vector< vector<double> > matrix_multiply(vector< vector<double> >& m1, vector< vector<double> >& m2) {

	//n x k by k x m = n x m matrix
	int rows1 = m1.size();
	int cols1 = m1[0].size();
	int rows2 = m2.size();
	int cols2 = m2[0].size();

	if (cols1 != rows2) {

		fprintf(stderr, "ERROR: Invalid matrix dims to matrix_multiply");
		exit(1);
	}

	vector< vector<double> > m3;
	for (int i = 0; i < rows1; i++) {
		vector<double> temp(cols2, 0);
		m3.push_back(temp);
	}
	
	for (int i = 0; i < rows1; i++) {
		for (int j = 0; j < cols2; j++) {

			for (int k = 0; k < cols1; k++) {
				m3[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return m3;
}

vector<double> matrix_multiply(vector< vector<double> >& m1, vector<double> &m2) {

	
	//n x k by k x m = n x m matrix
	int rows1 = m1.size();
	int cols1 = m1[0].size();
	int rows2 = 1;
	int cols2 = m2.size();

	if (cols1 != cols2) {

		fprintf(stderr, "ERROR: Invalid matrix dims to matrix_multiply");
		exit(1);
	}

	vector<double> m3(rows1, 0);
	
	for (int i = 0; i < rows1; i++) {
		for (int j = 0; j < cols2; j++) {

				m3[i] += m1[i][j] * m2[j];
			}
	}
	return m3;
}

vector< vector<double> > matrix_add(vector< vector <double> >& m1, vector< vector <double> > &m2) {

	if (m2.size() != m1.size()) {
		fprintf(stderr, "ERROR: Invalid matrix dims to matrix_add");
		exit(1);
	}

	vector< vector<double> > m3;
	for (size_t i = 0; i < m1.size(); i++) {
		vector<double> temp(m1[0].size(), 0);
		m3.push_back(temp);
	}

	for (size_t i = 0; i < m1.size(); i++) {
		for (size_t j = 0; j < m1[0].size(); j++) {
			m3[i][j] = m1[i][j] + m2[i][j];
		}
	}

	return m3;
}

vector<double> matrix_add(vector <double> &m1, vector <double> &m2) {

	if (m2.size() != m1.size()) {
		fprintf(stderr, "ERROR: Invalid matrix dims to matrix_add");
		exit(1);
	}

	vector<double> m3(m1.size(), 0);

	for (size_t i = 0; i < m1.size(); i++) {
		m3[i] = m1[i] + m2[i];
	}

	return m3;
}

//Prints out a 2d vector for debugging purposes
void debug (vector< vector<double> > v) {

	printf("\n");
	for (size_t i = 0; i < v.size(); i++) {
		printf("[ ");
		for (size_t j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}
void debug (vector< vector<int> > v) {

	printf("\n");
	for (size_t i = 0; i < v.size(); i++) {
		printf("[ ");
		for (size_t j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}
void debug(vector<double> v) {

	printf("\n");
	for (size_t i = 0; i < v.size(); i++) {
		printf("[");
		printf(" ");
		cout << v[i] << " ]\n";
	}
	return;
}
void debug(vector<int> v) {

	printf("\n[");
	for (size_t i = 0; i < v.size(); i++) {
		printf(" ");
		cout << v[i] << ", ";
	}
	printf("]\n");
	return;
}

//Functions that build the input model - Change these functions to change the input model
void build_nodes(vector< vector<double> > &NODES) {

	NODES[0].push_back(0.0);
	NODES[0].push_back(0.0);

	NODES[1].push_back(0.0);
	NODES[1].push_back(3099.0);

	NODES[2].push_back(0.0);
	NODES[2].push_back(5892.0);

	NODES[3].push_back(3048.0);
	NODES[3].push_back(3099.0);

	NODES[4].push_back(6096.0);
	NODES[4].push_back(0.0);

	NODES[5].push_back(6096.0);
	NODES[5].push_back(3099.0);

	NODES[6].push_back(6096.0);
	NODES[6].push_back(5892.0);

	return;
}
void build_elems(vector< vector<double> > &ELEMS) {

	ELEMS[0].push_back(1.0);
	ELEMS[0].push_back(2.0);
	ELEMS[0].push_back(200.0);
	ELEMS[0].push_back(10193.528);
	ELEMS[0].push_back(126118121.95679997);
	ELEMS[0].push_back(0.0);
	ELEMS[0].push_back(0.0);

	ELEMS[1].push_back(2.0);
	ELEMS[1].push_back(3.0);
	ELEMS[1].push_back(200.0);
	ELEMS[1].push_back(10193.528);
	ELEMS[1].push_back(126118121.95679997);
	ELEMS[1].push_back(0.0);
	ELEMS[1].push_back(0.0);

	ELEMS[2].push_back(5.0);
	ELEMS[2].push_back(6.0);
	ELEMS[2].push_back(200.0);
	ELEMS[2].push_back(10193.528);
	ELEMS[2].push_back(42871836.83679999);
	ELEMS[2].push_back(0.0);
	ELEMS[2].push_back(0.0);

	ELEMS[3].push_back(6.0);
	ELEMS[3].push_back(7.0);
	ELEMS[3].push_back(200.0);
	ELEMS[3].push_back(10193.528);
	ELEMS[3].push_back(42871836.83679999);
	ELEMS[3].push_back(0.0);
	ELEMS[3].push_back(0.0);

	ELEMS[4].push_back(2.0);
	ELEMS[4].push_back(4.0);
	ELEMS[4].push_back(200.0);
	ELEMS[4].push_back(10064.496);
	ELEMS[4].push_back(225181201.24959993);
	ELEMS[4].push_back(0.0);
	ELEMS[4].push_back(-0.012);

	ELEMS[5].push_back(4.0);
	ELEMS[5].push_back(6.0);
	ELEMS[5].push_back(200.0);
	ELEMS[5].push_back(10064.496);
	ELEMS[5].push_back(225181201.24959993);
	ELEMS[5].push_back(0.0);
	ELEMS[5].push_back(-0.012);

	ELEMS[6].push_back(3.0);
	ELEMS[6].push_back(7.0);
	ELEMS[6].push_back(200.0);
	ELEMS[6].push_back(10064.496);
	ELEMS[6].push_back(225181201.24959993);
	ELEMS[6].push_back(0.0);
	ELEMS[6].push_back(-0.008);

	ELEMS[7].push_back(1.0);
	ELEMS[7].push_back(4.0);
	ELEMS[7].push_back(300.0);
	ELEMS[7].push_back(3200.0);
	ELEMS[7].push_back(0.0);
	ELEMS[7].push_back(0.0);
	ELEMS[7].push_back(0.0);

	ELEMS[8].push_back(4.0);
	ELEMS[8].push_back(5.0);
	ELEMS[8].push_back(200.0);
	ELEMS[8].push_back(6283.8584);
	ELEMS[8].push_back(0.0);
	ELEMS[8].push_back(0.0);
	ELEMS[8].push_back(0.0);

	ELEMS[9].push_back(4.0);
	ELEMS[9].push_back(7.0);
	ELEMS[9].push_back(200.0);
	ELEMS[9].push_back(10580.623999999998);
	ELEMS[9].push_back(0.0);
	ELEMS[9].push_back(0.0);
	ELEMS[9].push_back(0.0);

	return;
}
void build_supports(vector< vector<double> > &SUPPORTS) {

	SUPPORTS[0].push_back(1.0);
	SUPPORTS[0].push_back(1.0);
	SUPPORTS[0].push_back(1.0);

	SUPPORTS[1].push_back(0.0);
	SUPPORTS[1].push_back(0.0);
	SUPPORTS[1].push_back(0.0);

	SUPPORTS[2].push_back(0.0);
	SUPPORTS[2].push_back(0.0);
	SUPPORTS[2].push_back(0.0);

	SUPPORTS[3].push_back(0.0);
	SUPPORTS[3].push_back(0.0);
	SUPPORTS[3].push_back(0.0);

	SUPPORTS[4].push_back(1.0);
	SUPPORTS[4].push_back(1.0);
	SUPPORTS[4].push_back(0.0);

	SUPPORTS[5].push_back(0.0);
	SUPPORTS[5].push_back(0.0);
	SUPPORTS[5].push_back(0.0);

	SUPPORTS[6].push_back(0.0);
	SUPPORTS[6].push_back(0.0);
	SUPPORTS[6].push_back(0.0);

	return;
}
void build_nodal_loads(vector< vector<double> > &NODALLOADS) {

	NODALLOADS[0].push_back(0.0);
	NODALLOADS[0].push_back(0.0);
	NODALLOADS[0].push_back(0.0);

	NODALLOADS[1].push_back(0.0);
	NODALLOADS[1].push_back(0.0);
	NODALLOADS[1].push_back(0.0);

	NODALLOADS[2].push_back(0.0);
	NODALLOADS[2].push_back(0.0);
	NODALLOADS[2].push_back(0.0);

	NODALLOADS[3].push_back(0.0);
	NODALLOADS[3].push_back(0.0);
	NODALLOADS[3].push_back(0.0);

	NODALLOADS[4].push_back(0.0);
	NODALLOADS[4].push_back(0.0);
	NODALLOADS[4].push_back(0.0);

	NODALLOADS[5].push_back(200.0);
	NODALLOADS[5].push_back(0.0);
	NODALLOADS[5].push_back(0.0);

	NODALLOADS[6].push_back(200.0);
	NODALLOADS[6].push_back(0.0);
	NODALLOADS[6].push_back(0.0);

	return;
}
void build_support_disps(vector< vector<double> > &SUPPORTDISPS) {

	SUPPORTDISPS[0].push_back(0.0);
	SUPPORTDISPS[0].push_back(0.0);
	SUPPORTDISPS[0].push_back(0.0);

	SUPPORTDISPS[1].push_back(0.0);
	SUPPORTDISPS[1].push_back(0.0);
	SUPPORTDISPS[1].push_back(0.0);

	SUPPORTDISPS[2].push_back(0.0);
	SUPPORTDISPS[2].push_back(0.0);
	SUPPORTDISPS[2].push_back(0.0);

	SUPPORTDISPS[3].push_back(0.0);
	SUPPORTDISPS[3].push_back(0.0);
	SUPPORTDISPS[3].push_back(0.0);

	SUPPORTDISPS[4].push_back(0.0);
	SUPPORTDISPS[4].push_back(0.0);
	SUPPORTDISPS[4].push_back(0.0);

	SUPPORTDISPS[5].push_back(0.0);
	SUPPORTDISPS[5].push_back(0.0);
	SUPPORTDISPS[5].push_back(0.0);

	SUPPORTDISPS[6].push_back(0.0);
	SUPPORTDISPS[6].push_back(0.0);
	SUPPORTDISPS[6].push_back(0.0);

	return;
}

//build the Local-basic transformation vector
void build_local_basic_transform(vector< vector<double> > &abl, double L) {

	abl[0].push_back(-1);
	abl[0].push_back(0);
	abl[0].push_back(0);
	abl[0].push_back(1);
	abl[0].push_back(0);
	abl[0].push_back(0);

	abl[1].push_back(0);
	abl[1].push_back(1/L);
	abl[1].push_back(1);
	abl[1].push_back(0);
	abl[1].push_back(-1/L);
	abl[1].push_back(0);

	abl[2].push_back(0);
	abl[2].push_back(1 / L);
	abl[2].push_back(0);
	abl[2].push_back(0);
	abl[2].push_back(-1 / L);
	abl[2].push_back(1);
}
