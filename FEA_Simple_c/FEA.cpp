#include <iostream>
#include <vector>
#include <numeric>

//http://eigen.tuxfamily.org/index.php?title=Main_Page
//License: MPL2 (free, copyleft)
#include <Eigen/Dense>

using namespace std;

vector< vector<float> > upper_triangular(vector< vector<float> > &m1);
vector< vector<float> > matrix_transpose(vector < vector<float> >& m1);
vector<vector<float>> matrix_multiply(vector<vector<float>>& m1, vector<vector<float>>& m2);

void build_nodes(vector< vector<float> > &NODES);
void build_elems(vector< vector<float> > &ELEMS);
void build_supports(vector< vector<float> > &SUPPORTS);
void build_nodal_loads(vector< vector<float> > &NODALLOADS);
void build_support_disps(vector< vector<float> > &SUPPORTDISPS);

void build_local_basic_transform(vector< vector<float> >& abl, float L);

void debug(vector< vector<float> >);
void debug(vector< vector<int> >);
void debug(vector<int> v);
void debug(vector<float> v);

int main() {

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

	printf("Building input model. . .\n");

	vector< vector<float> > NODES;
	for(int i = 0; i < NODE_Y; i++) {
		vector<float> temp;
		NODES.push_back(temp);
	}
	vector< vector<float> > ELEMS;
	for (int i = 0; i < ELEMS_Y; i++) {
		vector<float> temp;
		ELEMS.push_back(temp);
	}
	vector< vector<float> > SUPPORTS;
	for (int i = 0; i < SUPPORTS_Y; i++) {
		vector<float> temp;
		SUPPORTS.push_back(temp);
	}
	vector< vector<float> > NODALLOADS;
	for (int i = 0; i < NODALLOADS_Y; i++) {
		vector<float> temp;
		NODALLOADS.push_back(temp);
	}
	vector< vector<float> > SUPPORTDISPS;
	for (int i = 0; i < SUPPORTDISPS_Y; i++) {
		vector<float> temp;
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

	printf("Assembling Kff and modifying Pf. . .\n");

	/***************************
	*	Assign equation numbers
	***************************/
	int nfdof = 0; 
	int ncdof = 0;
	
	//Create equation array
	vector< vector<float> > EQUATIONS;
	for (int i = 0; i < NODE_Y; i++) {

		vector<float> temp(3, 0);
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

	vector <int> pf(nfdof, 0);
	vector <int> uc(ncdof, 0);

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

	//debug(pf);
	//debug(uc);

	/*********************************************************************************************************
	* Assemble stiffness matrix Kff and modify the load vector Pf for member loads and support displacements
	*********************************************************************************************************/

	//Declaring nele, uf, kff
	int nele = ELEMS_Y;

	vector <int> uf(nfdof, 0);

	vector < vector<float> > kff;
	for (int i = 0; i < nfdof; i++) {
		vector<float> temp(nfdof, 0);
		kff.push_back(temp);
	}

	//FOR DEBUG 
	/*
	vector < vector<float> > TEST;
	for (int i = 0; i < 5; i++) {
		vector<float> temp;
		TEST.push_back(temp);
		
		for (int j = 0; j < 5; j++) {
			TEST[i].push_back(i + j);
		}
	}
	vector < vector<float> > TEST2;
	for (int i = 0; i < 5; i++) {
		vector<float> temp;
		TEST2.push_back(temp);
		
		for (int j = 0; j < 5; j++) {
			TEST2[i].push_back(i + j);
		}
	}


	debug(TEST);
	debug(upper_triangular(TEST));

	debug( matrix_multiply(TEST, TEST2) );
	*/
	//END DEBUG 

	//Calculation
	for (int e = 0; e < nele; e++) {

		//I and J nodes for this element
		int ndI = (int)ELEMS[e][0];
		int ndJ = (int)ELEMS[e][1];

		//vector1.insert(vector1.end(), vector2.begin(), vector2.end());

		//Lookup vector
		vector <float> l;
		l = EQUATIONS[ndI - 1];
		l.insert( l.end(), EQUATIONS[ndJ - 1].begin(), EQUATIONS[ndJ-1].end() );

		//Difference in X-coordinates
		float dx = NODES[ndJ - 1][0] - NODES[ndI - 1][0];

		//Difference in y-coordinates
		float dy = NODES[ndJ - 1][1] - NODES[ndI - 1][1];

		//Element length and direction cosines
		float L = pow( (pow(dx, 2) + pow(dy, 2)), 0.5 );

		float c = dx / L;
		float s = dy / L;

		//Local-basic transformation
		vector< vector<float> > abl;
		for (int i = 0; i < 3; i++) {
			vector<float> temp;
			abl.push_back(temp);
		}
		build_local_basic_transform(abl, L);

		//Global-lobal transformation
		vector< vector<float> > alg;
		for (int i = 0; i < 6; i++) {
			vector<float> temp(6, 0);
			alg.push_back(temp);
		}
		alg[0][0] = c; alg[0][1] = s;
		alg[1][0] = -s; alg[1][1] = c;
		alg[2][2] = 1;
		alg[3][3] = c; alg[3][4] = s;
		alg[4][3] = -s; alg[4][4] = c;
		alg[5][5] = 1;

		//Select global displacements from Uf and Uc
		vector<float> u(6, 0);
		for(int i = 0; i < 6; i++) {
			if (l[i] > 0) {
				u[i] = uf[l[i] - 1];
			}
			else {
				u[i] = uc[-l[i] - 1];
			}
		}


		//Local Displacements
		vector<float> ul(6, 0);
		ul = matrix_multiply(alg, u);


	}

	//debug(kff);

	return 0;
}

vector< vector<float> > upper_triangular(vector< vector<float> > &m1) {

	int rows = m1.size();
	int cols = m1[0].size();

	if (rows != cols) {
		fprintf(stderr, "ERROR: You can only find the upper-triangular matrix of a square matrix.");
		exit(1);
	}

	vector< vector<float> > m2;
	for (int i = 0; i < rows; i++) {
		vector<float> temp(cols, 0);
		m2.push_back(temp);
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols - i; j++) {

			m2[i][cols - j - 1] = m1[i][cols - j - 1];
		}
	}
	return m2;
}

vector< vector<float> > matrix_transpose(vector< vector<float> > &m1) {

	int rows = m1.size();
	int cols = m1[0].size();

	vector< vector<float> > m2;
	for (int i = 0; i < cols; i++) {
		vector<float> temp(rows, 0);
		m2.push_back(temp);
	}

	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {

			m2[i][j] = m1[j][i];
		}
	}

	return m2;
}

vector< vector<float> > matrix_multiply(vector< vector<float> >& m1, vector< vector<float> >& m2) {

	
	//n x k by k x m = n x m matrix
	int rows1 = m1.size();
	int cols1 = m1[0].size();
	int rows2 = m2.size();
	int cols2 = m2[0].size();

	if (cols1 != rows2) {

		fprintf(stderr, "ERROR: Invalid matrix dims to matrix_multiply");
		exit(1);
	}

	vector< vector<float> > m3;
	for (int i = 0; i < rows1; i++) {
		vector<float> temp(cols2, 0);
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


//Prints out a 2d vector for debugging purposes
void debug (vector< vector<float> > v) {

	printf("\n");
	for (int i = 0; i < v.size(); i++) {
		printf("[ ");
		for (int j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}
void debug (vector< vector<int> > v) {

	printf("\n");
	for (int i = 0; i < v.size(); i++) {
		printf("[ ");
		for (int j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}
void debug(vector<float> v) {

	printf("\n[");
	for (int i = 0; i < v.size(); i++) {
		printf(" ");
		cout << v[i] << ", ";
	}
		printf("]\n");
	return;
}
void debug(vector<int> v) {

	printf("\n[");
	for (int i = 0; i < v.size(); i++) {
		printf(" ");
		cout << v[i] << ", ";
	}
	printf("]\n");
	return;
}

//Functions that build the input model - Change these functions to change the input model
void build_nodes(vector< vector<float> > &NODES) {

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
void build_elems(vector< vector<float> > &ELEMS) {

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
void build_supports(vector< vector<float> > &SUPPORTS) {

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
void build_nodal_loads(vector< vector<float> > &NODALLOADS) {

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
void build_support_disps(vector< vector<float> > &SUPPORTDISPS) {

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
void build_local_basic_transform(vector< vector<float> > &abl, float L) {

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
	abl[1].push_back(1/L);
	abl[1].push_back(0);

	abl[2].push_back(0);
	abl[2].push_back(1 / L);
	abl[2].push_back(0);
	abl[2].push_back(0);
	abl[2].push_back(-1 / L);
	abl[2].push_back(1);
}


