/*****************************************************************************************************

 This file is where you'll put the source information for the functions you write.  Make sure it's 
 included in your project (shift-alt-A) and that the functions you add are declared in the header
 file, myEigenFunctions.h.  You can add any extra comments or checks that could be relevant to these
 functions as you need to.

====================================================================================================*/

#include <string>

#include "myEigenFunctions.h"

using namespace std;

double DotProduct(double **A, double **B, int n, int m)
{
	//
	//	This is a function that takes two matrices A and B of identical dimensions (n*m) and 
	//  calculates and returns their dot product.
	//
	double dot = 0.0;

	for(int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			dot += A[i][j]*B[i][j];
		}
	}
	return dot;
}

double DotProduct(double *A, double *B, int n) 
{
	//
	//	This is a function that takes two vectors A and B of identical length (n) and 
	//  calculates and returns their dot product.
	//

	double dot = 0.0;

	for(int i = 0; i < n; i++) {
		dot += A[i]*B[i];
	}
	return dot;
}

double* DotProduct(double **A, double *v, int n)
    {
          //
          //  This is a function that takes a nxn-matrix A and an n-dimensional vector v stores
          //  the product A.v at the original location of v
          //

    double *result = new double[n]; // pointer to result vector
  
    for(int i = 0; i < n; i++) {
      result[i] = 0.0; // initialize ith element of result v
      for(int j = 0; j <n; j++){
        result[i] += A[i][j] * v[j]; 
      }
    }

    return result;
  }



double** ReadData(string inputFileName, int n, int m)
    {
          //
          //  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
          //  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
          // The code for parsing the CSV file has been adapted from Bernardo Trindade
          // https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c

      double **spectra = new double* [n];

     
      for(int sample = 0; sample < n; sample++)
	spectra[sample] = new double [m];

    vector<vector<double> > data;
    cout << inputFileName << endl;
    ifstream inputFile(inputFileName);
    int l = 0;
      
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
	// cout << s << endl;
        if (l!=0) // ignore first line, which contains wavelengths
	  {
            istringstream ss(s);
            vector<double> record;

	    int column = 0;
            while (ss) {
	      // cout << "Row " << l << " " << ss << endl;
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
		  // cout << "Row " << l << " Column " << column << " line " << line << endl; 
                    spectra[l-1][column] = stod(line);
		    column++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
		}

        }
	l++;
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        //__throw_invalid_argument("File not found.");
    }      
      
    return spectra;
  }

void print_matrix(double **A, int n, int m)
{
  //
  // This procedure prints an nxm matrix.
  //
  for (int i=0; i < n; i++)
    {
    for (int j=0; j < m; j++)
      std::cout << A[i][j] << '\t';
    std::cout << '\n';
    }
}

Eigenpair power_method(double **A, double *v, int n, double tol)
{
  //
  // This function computes the largest eigenvalue and the corresponding eigenvector
  //

  Eigenpair eigenpair(n);


  
      
  return eigenpair;  
}

void deflate(double **A, Eigenpair eigenpair)
{
  //
  // This procedure removes eigenpair.vector from transformation matrix A in place
  //

  // Your code comes here
}

double** CenterMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and subtracts the emperical mean of each column.
          //  The function returns the point to the centered matrix
          //
      
      double **result = new double* [n];
      
      for(int row = 0; row < n; row++)
	    result[row] = new double [m];

      // Your code comes here
      vector<double> means(n);

      for(int col = 0; col< m; col++){
        double sum = 0;
        for(int row= 0; row < n; row++){
          sum += A[row][col];
        }
        double mean = sum/n;
        means[col] = mean;
      }

      for(int col = 0; col < m; col++){
        for(int row =0; row < n; row++){
          A[row][col] = A[row][col] - means[col];
        }
      }
      
    return result;
  }


double** CovarianceMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and 
          //  
          //
      
      double **cov = new double* [m];
      double **cA = CenterMatrix(A, n, m);
      
      for(int i = 0; i < m; i++)

	      cov[i] = new double [m];

      // Your code comes here

      
    return cov;
  }
