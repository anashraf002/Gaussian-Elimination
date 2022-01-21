#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

class PROGRAM_BANK{

   public:
  PROGRAM_BANK(int); // Default constructor
	PROGRAM_BANK(int, char *);//PROGRAM_BANK mat(size, file_name);
	void INVERT_MATRIX(char *);
   private:
	int N;
	float M[20][20];
};
PROGRAM_BANK::PROGRAM_BANK(int x)
{
	N = x;

	for (int i=0; i<N; i++)
	{
  		for (int j=0; j<N; j++)
		{
    			M[i][j]=0;
			//cout << M[i][j]<< " ";
			//cout << endl;
		}
	}
}

PROGRAM_BANK::PROGRAM_BANK(int x, char *y)
{
	N = x;
	ifstream infile(y, ios::in);
	for (int i=0; i<N; i++)
	{
  		for (int j=0; j<N; j++)
		{
			infile >> M[i][j];
			//cout << M[i][j]<< " ";
			//cout << endl;
		}
	}
}

void PROGRAM_BANK::INVERT_MATRIX(char *z)
{
	//program for matrix inversion
	ofstream outfile(z, ios::out);

	int d, e, i, j, k;
	float f, q[20][20], w[20][40], b, c;


	// Some stuff needed to be taken care of
	d = N;
	e = 2 * N;

	// Formulating the identity matrix
	for ( i = 0; i < N; i++ )
	{
		for ( j = d; j < e; j++ )
		{
			f = j - N;
			if ( i == f )
			{
		       		q[i][j] = 1;
			}
			else
			{
			       	q[i][j] = 0; 
			}
		}
	}

	// Putting A and the identity matrix, Q, together to make W
	for ( i=0; i < N; i++ )
	{
		for ( j = 0; j < e; j++ )
		{
			if ( j < N )
			{
			       	w[i][j] = M[i][j];
		       	}
			else if ( j >= N )
			{
			       	w[i][j] = q[i][j];
		       	}
		}
	}

	// The actual row reduction part
	//&Yacute; Comparing column entries
	for ( i = 0; i < N; i++ )
	{
		for ( j = 0; j < N; j++ )
		{
			if ( i != j )
			{
				b = w[i][i]; c = w[j][i];
				if (( b > 0 && c > 0 ) || ( b < 0 && c < 0 ))
				{
					for ( k = 0; k < e; k++ )
					{
						w[j][k]= ( fabs(c) * w[i][k] ) - ( fabs(b) * w[j][k] );
					}
				}
				else
				{
					for ( k = 0; k < e; k++ )
					{
						w[j][k]= ( fabs(c) * w[i][k] ) + ( fabs(b) * w[j][k] );
					}
				}
			}
		}
	}

	// The last divide
	for ( i = 0; i < N; i++ )
	{
		for ( j = d; j < e; j++ )
		{
			w[i][j] = w[i][j] / w[i][i];
		}
	}

	for ( i = 0; i < N; i++ )
	{
	
		// outfile << ">>> FROM_PROGRAM_BANK: ";
		for ( j = d; j < e; j++ )
		{
		       	outfile << setiosflags(ios::fixed | ios::showpoint) << setprecision(2) << w[i][j] <<" ";
	       	}
		outfile << endl;
	}
	//outfile << ">>> FROM_PROGRAM_BANK: END OF INVERSION";
	//outfile << endl;
}
