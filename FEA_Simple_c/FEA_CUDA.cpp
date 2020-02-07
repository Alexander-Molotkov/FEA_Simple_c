#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <assert.h>

//In order for CUDA to work you to link to the CUDA library in VS settings and set the project as a CUDA project
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolver_common.h>
#include <cusolverDn.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

using namespace std;

vector<double> cholesky_solve(vector<vector<double>>& a, vector<double>& b);
vector<double> cusolver_solve(vector<vector<double>>& a, vector<double>& b);

vector<double> forward_sub(vector<vector<double>>& l, vector<double>& b);
vector<double> backward_sub(vector<vector<double>>& u, vector<double>& y);
vector< vector<double> > cholesky_decomp(vector< vector<double> >& m1);
vector< vector<double> > matrix_transpose(vector < vector<double> >& m1);
vector<vector<double>> matrix_multiply(vector<vector<double>>& m1, vector<vector<double>>& m2);
vector<double> matrix_multiply(vector< vector<double> >& m1, vector< double>& m2);
vector< vector<double> > matrix_add(vector< vector <double> >& m1, vector< vector <double> >& m2);
vector<double> matrix_add(vector <double>& m1, vector <double>& m2);

double* vector_unwrap(vector<vector<double>>& v);
double* vector_unwrap(vector<double>& v);

vector< vector<double> > build_nodes(int MODEL_REPETITONS);
vector< vector<double> > build_elems(int MODEL_REPETITION);
vector< vector<double> > build_supports(int MODEL_REPETITIONS);
vector< vector<double> > build_nodal_loads(int MODEL_REPETITIONS);
vector< vector<double> > build_support_disps(int MODEL_REPETITIONS);

void build_local_basic_transform(int MODEL_REPETITIONS, vector< vector<double> >& abl, double L);

void debug(vector< vector<double> >);
void debug(vector< vector<int> >);
void debug(double* v, int size);
void debug(vector<int> v);
void debug(vector<double> v);

int main() {



	//END CUSOLVER SETUP

	//Repeats the base model x times for benchmarking purposes
	//When repeated, the model may not be accurate
	int MODEL_REPETITIONS = 50;

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

	vector< vector<double> > NODES = build_nodes(MODEL_REPETITIONS);
	vector< vector<double> > ELEMS = build_elems(MODEL_REPETITIONS);
	vector< vector<double> > SUPPORTS = build_supports(MODEL_REPETITIONS);
	vector< vector<double> > NODALLOADS = build_nodal_loads(MODEL_REPETITIONS);
	vector< vector<double> > SUPPORTDISPS = build_support_disps(MODEL_REPETITIONS);

	/**************************************************************************
	*Everything below this line is generic and should work for any input above
	**************************************************************************/

	auto checkpoint = chrono::high_resolution_clock::now();
	auto checkpoint_start = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(checkpoint - start);
	cout << "Input model built in " << duration.count() << " milliseconds. \n";

	/***************************
	*	Assign equation numbers
	***************************/
	double nfdof = 0;
	double ncdof = 0;

	//Create equation array
	vector< vector<double> > EQUATIONS;
	for (int i = 0; i < NODE_Y * MODEL_REPETITIONS; i++) {

		vector<double> temp(3, 0);
		EQUATIONS.push_back(temp);
	}

	//Positive equation numbers are assigned to unconstrained DOFs
	//Negative equation numbers are assigned to constrained DOFs, i.e., where reactions where be developed
	for (int i = 0; i < NODE_Y * MODEL_REPETITIONS; i++) {
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

	/**********************************************************************
	*	Create the joint load vector Pf and support displacement vector Uc
	**********************************************************************/

	vector <double> pf(nfdof, 0.0000);
	vector <double> uc(ncdof, 0.0000);

	for (int i = 0; i < NODE_Y * MODEL_REPETITIONS; i++) {
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
	for (int e = 0; e < nele * MODEL_REPETITIONS; e++) {

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
		build_local_basic_transform(MODEL_REPETITIONS, abl, L);

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
		for (int j = 0; j < 6; j++) {

			for (double x = 0; x < MODEL_REPETITIONS; x++) {

				if (l[j] >= 0) {

					pf[(l[j] - 1) + (x * (double)16)] -= p[j];

					for (int i = 0; i < 6; i++) {

						if (l[i] >= 0) {
							kff[(l[i] - 1) + (x * (double)16)][(l[j] - 1) + (x * (double)16)] += k[i][j];
						}
					}
				}
			}
		}
	}

	/**************************************
	* Solve for the nodal displacements Uf
	**************************************/

	if (nfdof > 0) {
		//uf = cholesky_solve(kff, pf);
		uf = cusolver_solve(kff, pf);
		debug(uf);
	}

	/****************************
	* Assemble the reactions Pc
	****************************/

	checkpoint_start = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(checkpoint_start - start);
	cout << "Kff assembled, Pf modified, and Uf solved for in " << duration.count() << " milliseconds. ";

	vector<double> pc(ncdof, 0);

	for (int e = 0; e < nele * MODEL_REPETITIONS; e++) {

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
		build_local_basic_transform(MODEL_REPETITIONS, abl, L);

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
		vector< vector<double> > x = matrix_transpose(abl);
		vector<double> pl = matrix_multiply(x, pb);
		pl = matrix_add(pl, plw);

		//Global forces
		x = matrix_transpose(alg);
		vector<double> p = matrix_multiply(x, pl);

		//Assemble Pc
		for (int j = 0; j < 6; j++) {

			for (int x = 0; x < MODEL_REPETITIONS; x++) {

				if (l[j] < 0) {
					pc[(-l[j] - 1) + (x * (double)5)] += p[j];
				}
			}
		}
	}

	duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - checkpoint_start);
	cout << "\nPC assembled in " << duration.count() << " milliseconds.\n";

	return 0;
}

vector<double> cholesky_solve(vector< vector<double> > & a, vector<double> & b) {

	//Solve for Ax = b
	vector< vector<double> > l = cholesky_decomp(a);
	vector< vector<double> > u = matrix_transpose(l);

	vector<double> y = forward_sub(l, b);
	vector<double> x = backward_sub(u, y);

	return x;
}

vector<double> cusolver_solve(vector<vector<double>>& kff, vector<double>& pf){

	//CUSOLVER SETUP
	cusolverDnHandle_t handle = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cudaStat1 = cudaSuccess;
	cudaError_t cudaStat2 = cudaSuccess;
	cudaError_t cudaStat3 = cudaSuccess;
	cudaError_t cudaStat4 = cudaSuccess;

	const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
	const int batchSize = 2;
	const int nrhs = 1;
	const int m = 3;
	const int lda = m;
	const int ldb = m;

	double* A1 = vector_unwrap(kff);
	double* B0 = vector_unwrap(pf);
	double X0[m]; /* X0 = A0\B0 */
	int infoArray[batchSize]; /* host copy of error info */

	double L0[lda * m]; /* cholesky factor of A0 */

	double* Aarray[batchSize];
	double* Barray[batchSize];

	double** d_Aarray = NULL;
	double** d_Barray = NULL;
	int* d_infoArray = NULL;

	/* step 1: create cusolver handle, bind a stream */
	status = cusolverDnCreate(&handle);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	assert(cudaSuccess == cudaStat1);

	status = cusolverDnSetStream(handle, stream);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/* step 2: copy A to device */
	for (int j = 0; j < batchSize; j++) {
		cudaStat1 = cudaMalloc((void**)& Aarray[j], sizeof(double) * lda * m);
		assert(cudaSuccess == cudaStat1);
		cudaStat2 = cudaMalloc((void**)& Barray[j], sizeof(double) * ldb * nrhs);
		assert(cudaSuccess == cudaStat2);
	}
	cudaStat1 = cudaMalloc((void**)& d_infoArray, sizeof(int) * batchSize);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMalloc((void**)& d_Aarray, sizeof(double*) * batchSize);
	cudaStat2 = cudaMalloc((void**)& d_Barray, sizeof(double*) * batchSize);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	cudaStat2 = cudaMemcpy(Aarray[1], A1, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	cudaStat1 = cudaMemcpy(Barray[0], B0, sizeof(double) * m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(Barray[1], B0, sizeof(double) * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	cudaStat1 = cudaMemcpy(d_Aarray, Aarray, sizeof(double*) * batchSize, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_Barray, Barray, sizeof(double*) * batchSize, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	cudaDeviceSynchronize();

	/* step 3: Cholesky factorization */
	status = cusolverDnDpotrfBatched(
		handle,
		uplo,
		m,
		d_Aarray,
		lda,
		d_infoArray,
		batchSize);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(infoArray, d_infoArray, sizeof(int) * batchSize, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(L0, Aarray[0], sizeof(double) * lda * m, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	//assert(0 == infoArray[0]);
	/* A1 is singular */
	//assert(2 == infoArray[1]);

	// step 4: solve A0*X0 = B0
	status = cusolverDnDpotrsBatched(
		handle,
		uplo,
		m,
		nrhs, /* only support rhs = 1*/
		d_Aarray,
		lda,
		d_Barray,
		ldb,
		d_infoArray,
		batchSize);


	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(infoArray, d_infoArray, sizeof(int), cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(X0, Barray[0], sizeof(double) * m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	cudaDeviceSynchronize();

	assert(0 == infoArray[0]);

	/* free resources */
	if (d_Aarray) cudaFree(d_Aarray);
	if (d_Barray) cudaFree(d_Barray);
	if (d_infoArray) cudaFree(d_infoArray);

	if (handle) cusolverDnDestroy(handle);
	if (stream) cudaStreamDestroy(stream);

	cudaDeviceReset();

	debug(X0, kff.size());

	return vector<double>();
}

vector<double> forward_sub(vector< vector<double> > & l, vector<double> & b) {

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

vector<double> backward_sub(vector< vector<double> > & u, vector<double> & y) {

	vector<double> x(u.size(), 0);

	//On an upper-triangular matrix
	for (int i = u.size() - 1; i >= 0; i--) {

		double s = y[i];

		for (size_t j = i; j < u.size(); j++) {

			s -= u[i][j] * x[j];
		}
		x[i] = s / u[i][i];
	}

	return x;
}

vector< vector<double> > cholesky_decomp(vector< vector<double> > & matrix) {

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

vector< vector<double> > matrix_transpose(vector< vector<double> > & m1) {

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

vector< vector<double> > matrix_multiply(vector< vector<double> > & m1, vector< vector<double> > & m2) {

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

vector<double> matrix_multiply(vector< vector<double> > & m1, vector<double> & m2) {


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

vector< vector<double> > matrix_add(vector< vector <double> > & m1, vector< vector <double> > & m2) {

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

vector<double> matrix_add(vector <double> & m1, vector <double> & m2) {

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

double* vector_unwrap(vector< vector<double> >& v) {

	double* arr = new double[v.size() * v[0].size()];

	for (int i = 0; i < v.size(); i++) {
		copy(v[i].begin(), v[i].end(), (arr+v[0].size()*i));
	}

	return arr;
}

double* vector_unwrap(vector<double> & v) {

	double* arr = new double[v.size()];

	copy(v.begin(), v.end(), arr);

	return arr;
}


//Prints out a 2d vector for debugging purposes
void debug(vector< vector<double> > v) {

	cout << endl << v.size() << "x" << v[0].size() << ":" << endl;

	for (size_t i = 0; i < v.size(); i++) {
		printf("[ ");
		for (size_t j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}

void debug(vector< vector<int> > v) {

	cout << endl << v.size() << "x" << v[0].size() << ":" << endl;

	for (size_t i = 0; i < v.size(); i++) {
		printf("[ ");
		for (size_t j = 0; j < v[i].size(); j++) {

			cout << v[i][j] << ", ";
		}
		printf("]\n");
	}

	return;
}
void debug(double* v, int size) {

	for (size_t i = 0; i < size; i++) {
		printf("[");
		printf(" ");
		cout << v[i] << " ]\n";
	}
	return;
}
void debug(vector<double> v) {

	cout << endl << v.size() << "x0:" << endl;

	for (size_t i = 0; i < v.size(); i++) {
		printf("[");
		printf(" ");
		cout << v[i] << " ]\n";
	}
	return;
}
void debug(vector<int> v) {

	cout << endl << v.size() << "x0:" << endl;

	printf("\n[");
	for (size_t i = 0; i < v.size(); i++) {
		printf(" ");
		cout << v[i] << ", ";
	}
	printf("]\n");
	return;
}

//Functions that build the input model - Change these functions to change the input model
vector< vector<double> > build_nodes(int MODEL_REPETITIONS) {

	vector< vector<double> > NODES;
	for (int i = 0; i < 7 * MODEL_REPETITIONS; i++) {
		vector<double> temp;
		NODES.push_back(temp);
	}

	for (int i = 0; i < MODEL_REPETITIONS; i++) {
		NODES[0 + (i * 7)].push_back(0.0);
		NODES[0 + (i * 7)].push_back(0.0);

		NODES[1 + (i * 7)].push_back(0.0);
		NODES[1 + (i * 7)].push_back(3099.0);

		NODES[2 + (i * 7)].push_back(0.0);
		NODES[2 + (i * 7)].push_back(5892.0);

		NODES[3 + (i * 7)].push_back(3048.0);
		NODES[3 + (i * 7)].push_back(3099.0);

		NODES[4 + (i * 7)].push_back(6096.0);
		NODES[4 + (i * 7)].push_back(0.0);

		NODES[5 + (i * 7)].push_back(6096.0);
		NODES[5 + (i * 7)].push_back(3099.0);

		NODES[6 + (i * 7)].push_back(6096.0);
		NODES[6 + (i * 7)].push_back(5892.0);
	}

	return NODES;
}
vector< vector<double> > build_elems(int MODEL_REPETITIONS) {

	vector< vector<double> > ELEMS;
	for (int i = 0; i < 10 * MODEL_REPETITIONS; i++) {
		vector<double> temp;
		ELEMS.push_back(temp);
	}

	for (int i = 0; i < MODEL_REPETITIONS; i++) {
		ELEMS[0 + (i * 10)].push_back(1.0);
		ELEMS[0 + (i * 10)].push_back(2.0);
		ELEMS[0 + (i * 10)].push_back(200.0);
		ELEMS[0 + (i * 10)].push_back(10193.528);
		ELEMS[0 + (i * 10)].push_back(126118121.95679997);
		ELEMS[0 + (i * 10)].push_back(0.0);
		ELEMS[0 + (i * 10)].push_back(0.0);

		ELEMS[1 + (i * 10)].push_back(2.0);
		ELEMS[1 + (i * 10)].push_back(3.0);
		ELEMS[1 + (i * 10)].push_back(200.0);
		ELEMS[1 + (i * 10)].push_back(10193.528);
		ELEMS[1 + (i * 10)].push_back(126118121.95679997);
		ELEMS[1 + (i * 10)].push_back(0.0);
		ELEMS[1 + (i * 10)].push_back(0.0);

		ELEMS[2 + (i * 10)].push_back(5.0);
		ELEMS[2 + (i * 10)].push_back(6.0);
		ELEMS[2 + (i * 10)].push_back(200.0);
		ELEMS[2 + (i * 10)].push_back(10193.528);
		ELEMS[2 + (i * 10)].push_back(42871836.83679999);
		ELEMS[2 + (i * 10)].push_back(0.0);
		ELEMS[2 + (i * 10)].push_back(0.0);

		ELEMS[3 + (i * 10)].push_back(6.0);
		ELEMS[3 + (i * 10)].push_back(7.0);
		ELEMS[3 + (i * 10)].push_back(200.0);
		ELEMS[3 + (i * 10)].push_back(10193.528);
		ELEMS[3 + (i * 10)].push_back(42871836.83679999);
		ELEMS[3 + (i * 10)].push_back(0.0);
		ELEMS[3 + (i * 10)].push_back(0.0);

		ELEMS[4 + (i * 10)].push_back(2.0);
		ELEMS[4 + (i * 10)].push_back(4.0);
		ELEMS[4 + (i * 10)].push_back(200.0);
		ELEMS[4 + (i * 10)].push_back(10064.496);
		ELEMS[4 + (i * 10)].push_back(225181201.24959993);
		ELEMS[4 + (i * 10)].push_back(0.0);
		ELEMS[4 + (i * 10)].push_back(-0.012);

		ELEMS[5 + (i * 10)].push_back(4.0);
		ELEMS[5 + (i * 10)].push_back(6.0);
		ELEMS[5 + (i * 10)].push_back(200.0);
		ELEMS[5 + (i * 10)].push_back(10064.496);
		ELEMS[5 + (i * 10)].push_back(225181201.24959993);
		ELEMS[5 + (i * 10)].push_back(0.0);
		ELEMS[5 + (i * 10)].push_back(-0.012);

		ELEMS[6 + (i * 10)].push_back(3.0);
		ELEMS[6 + (i * 10)].push_back(7.0);
		ELEMS[6 + (i * 10)].push_back(200.0);
		ELEMS[6 + (i * 10)].push_back(10064.496);
		ELEMS[6 + (i * 10)].push_back(225181201.24959993);
		ELEMS[6 + (i * 10)].push_back(0.0);
		ELEMS[6 + (i * 10)].push_back(-0.008);

		ELEMS[7 + (i * 10)].push_back(1.0);
		ELEMS[7 + (i * 10)].push_back(4.0);
		ELEMS[7 + (i * 10)].push_back(300.0);
		ELEMS[7 + (i * 10)].push_back(3200.0);
		ELEMS[7 + (i * 10)].push_back(0.0);
		ELEMS[7 + (i * 10)].push_back(0.0);
		ELEMS[7 + (i * 10)].push_back(0.0);

		ELEMS[8 + (i * 10)].push_back(4.0);
		ELEMS[8 + (i * 10)].push_back(5.0);
		ELEMS[8 + (i * 10)].push_back(200.0);
		ELEMS[8 + (i * 10)].push_back(6283.8584);
		ELEMS[8 + (i * 10)].push_back(0.0);
		ELEMS[8 + (i * 10)].push_back(0.0);
		ELEMS[8 + (i * 10)].push_back(0.0);

		ELEMS[9 + (i * 10)].push_back(4.0);
		ELEMS[9 + (i * 10)].push_back(7.0);
		ELEMS[9 + (i * 10)].push_back(200.0);
		ELEMS[9 + (i * 10)].push_back(10580.623999999998);
		ELEMS[9 + (i * 10)].push_back(0.0);
		ELEMS[9 + (i * 10)].push_back(0.0);
		ELEMS[9 + (i * 10)].push_back(0.0);
	}

	return ELEMS;
}
vector< vector<double> > build_supports(int MODEL_REPETITIONS) {

	vector< vector<double> > SUPPORTS;
	for (int i = 0; i < 7 * MODEL_REPETITIONS; i++) {
		vector<double> temp;
		SUPPORTS.push_back(temp);
	}

	for (int i = 0; i < MODEL_REPETITIONS; i++) {
		SUPPORTS[0 + (i * 7)].push_back(1.0);
		SUPPORTS[0 + (i * 7)].push_back(1.0);
		SUPPORTS[0 + (i * 7)].push_back(1.0);

		SUPPORTS[1 + (i * 7)].push_back(0.0);
		SUPPORTS[1 + (i * 7)].push_back(0.0);
		SUPPORTS[1 + (i * 7)].push_back(0.0);

		SUPPORTS[2 + (i * 7)].push_back(0.0);
		SUPPORTS[2 + (i * 7)].push_back(0.0);
		SUPPORTS[2 + (i * 7)].push_back(0.0);

		SUPPORTS[3 + (i * 7)].push_back(0.0);
		SUPPORTS[3 + (i * 7)].push_back(0.0);
		SUPPORTS[3 + (i * 7)].push_back(0.0);

		SUPPORTS[4 + (i * 7)].push_back(1.0);
		SUPPORTS[4 + (i * 7)].push_back(1.0);
		SUPPORTS[4 + (i * 7)].push_back(0.0);

		SUPPORTS[5 + (i * 7)].push_back(0.0);
		SUPPORTS[5 + (i * 7)].push_back(0.0);
		SUPPORTS[5 + (i * 7)].push_back(0.0);

		SUPPORTS[6 + (i * 7)].push_back(0.0);
		SUPPORTS[6 + (i * 7)].push_back(0.0);
		SUPPORTS[6 + (i * 7)].push_back(0.0);
	}

	return SUPPORTS;
}
vector< vector<double> > build_nodal_loads(int MODEL_REPETITIONS) {

	vector< vector<double> > NODALLOADS;
	for (int i = 0; i < 7 * MODEL_REPETITIONS; i++) {
		vector<double> temp;
		NODALLOADS.push_back(temp);
	}

	for (int i = 0; i < MODEL_REPETITIONS; i++) {
		NODALLOADS[0 + (i * 7)].push_back(0.0);
		NODALLOADS[0 + (i * 7)].push_back(0.0);
		NODALLOADS[0 + (i * 7)].push_back(0.0);

		NODALLOADS[1 + (i * 7)].push_back(0.0);
		NODALLOADS[1 + (i * 7)].push_back(0.0);
		NODALLOADS[1 + (i * 7)].push_back(0.0);

		NODALLOADS[2 + (i * 7)].push_back(0.0);
		NODALLOADS[2 + (i * 7)].push_back(0.0);
		NODALLOADS[2 + (i * 7)].push_back(0.0);

		NODALLOADS[3 + (i * 7)].push_back(0.0);
		NODALLOADS[3 + (i * 7)].push_back(0.0);
		NODALLOADS[3 + (i * 7)].push_back(0.0);

		NODALLOADS[4 + (i * 7)].push_back(0.0);
		NODALLOADS[4 + (i * 7)].push_back(0.0);
		NODALLOADS[4 + (i * 7)].push_back(0.0);

		NODALLOADS[5 + (i * 7)].push_back(200.0);
		NODALLOADS[5 + (i * 7)].push_back(0.0);
		NODALLOADS[5 + (i * 7)].push_back(0.0);

		NODALLOADS[6 + (i * 7)].push_back(200.0);
		NODALLOADS[6 + (i * 7)].push_back(0.0);
		NODALLOADS[6 + (i * 7)].push_back(0.0);
	}

	return NODALLOADS;
}
vector< vector<double> > build_support_disps(int MODEL_REPETITIONS) {

	vector< vector<double> > SUPPORTDISPS;
	for (int i = 0; i < 7 * MODEL_REPETITIONS; i++) {
		vector<double> temp;
		SUPPORTDISPS.push_back(temp);
	}

	for (int i = 0; i < MODEL_REPETITIONS; i++) {
		SUPPORTDISPS[0 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[0 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[0 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[1 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[1 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[1 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[2 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[2 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[2 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[3 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[3 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[3 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[4 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[4 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[4 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[5 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[5 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[5 + (i * 7)].push_back(0.0);

		SUPPORTDISPS[6 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[6 + (i * 7)].push_back(0.0);
		SUPPORTDISPS[6 + (i * 7)].push_back(0.0);
	}

	return SUPPORTDISPS;
}

//build the Local-basic transformation vector
void build_local_basic_transform(int MODEL_REPETITIONS, vector< vector<double> > & abl, double L) {

	abl[0].push_back(-1);
	abl[0].push_back(0);
	abl[0].push_back(0);
	abl[0].push_back(1);
	abl[0].push_back(0);
	abl[0].push_back(0);

	abl[1].push_back(0);
	abl[1].push_back(1 / L);
	abl[1].push_back(1);
	abl[1].push_back(0);
	abl[1].push_back(-1 / L);
	abl[1].push_back(0);

	abl[2].push_back(0);
	abl[2].push_back(1 / L);
	abl[2].push_back(0);
	abl[2].push_back(0);
	abl[2].push_back(-1 / L);
	abl[2].push_back(1);
}