#include<iostream.h>
#include<fstream.h>
#include<stdlib.h>
#include<string.h>
#include "/run/media/Anas/Storage/Documents/Projects/inversion.h"

ofstream output_file("out.5", ios::out);

class MATRIX: public PROGRAM_BANK{

	friend MATRIX operator / (int, MATRIX &); 
				// FRIEND FUNCTION; example: 5 / X;
				// Compute inverse of X using INHERITANCE;
				// Multiply each element of inverted matrix
				// by 5; 
				// store the results into a temp MATRIX;
				// return the temp MATRIX;
	
	friend MATRIX operator * (int, MATRIX &); 
				// FRIEND FUNCTION; example: 5 * X;
				// each element of the matrix X is 
				// multiplied by 5;
				// store the results into a temp MATRIX;
				// return the temp MATRIX;

	public:	// public interfaces for this class
	
	  MATRIX (int, int); 	// constructor;
				// example: MATRIX X(d1,d2);
			   	// creates a matrix object with
			   	// private variables dim1 = d1 and dim2 = d2
			  	// and initializes each element to 0;
  
	  MATRIX(int, int, char *); 
	  			// constructor;
				// example: MATRIX X(d1,d2,"in_name");
				// read the elements of the matrix
				// from a given file in_name; d1 and d2 are 
				// the two dimensions of the matrix;
				

	  void PRINT();		// example: X.PRINT();
				// print the elements of the matrix X
				// into out.5 file;
				// returns no values;

	  int  operator < (int); //example: if( X < b)
				// check if the minimum element of the matrix X
				// is less than b
				// returns 1 if so; returns 0 otherwise;
		
	  void operator = (MATRIX &);
	  			// example X = Y;
				// Copy Y to X;
				// returns no values;		
		
	  void operator += (MATRIX &);  
	  			// example: X += Y;
                        	// matrix addition;
                        	// store the results in X;
				// values of Y remain unchanged;
				// returns no values;
		
	  void operator *= (MATRIX &); 
	  			// example: X *= Y;
				// matrix multiplication;
				// store the results in X;
				// values of Y remain unchanged;
				// returns no values;

	  MATRIX operator + (MATRIX &); 
	  			// example: X + Y;
				// Create a temp matrix object; 
				// add each element of Y to each element of X;
				// return this temp matrix object;
				// values of X and Y remain unchanged;

	  MATRIX operator * (MATRIX &); 
	  			// example: X * Y;
				// matrix multiplication;
				// store the results into a temp matrix
				// and return this temp matrix;
				// values of Y remain unchanged;

	  MATRIX operator / (MATRIX &); 
	  			// example: X / Y;
				// get inverse of Y and multiply it with X;
				// put the results into a temp MATRIX
				// and return temp MATRIX;
				// X and Y remain unchanged;
				// USE INHERITENCE FOR MATRIX INVERSION;
	
	private:

  	  int dim1, dim2; 	// dimensions of the matrix;
				// if dim2 is 1, a vector (no need for special
				// handling of vectors; everything should work
				// for both matrices and vectors)

  	  double A[20][20]; 	// 2 dimensional array for the matrix elements;
	  char in_file[11]; 	// hold the input file name
 };


// constructor
MATRIX::MATRIX(int d1, int d2)
        :PROGRAM_BANK(d1) // use the inheritance
{
     int i, j;
     dim1 = d1;
     dim2 = d2;
     
     if(d1 < 0 || d1 > 20 || d2< 0 || d2 >20)
     {
          exit(1);
     }
     
     for (i= 0; i<dim1; i++)
     {
         for (j= 0; j<dim2; j++)
         {
             A[i][j] = 0;
         }
     }

}

// constructor
MATRIX::MATRIX(int d1, int d2, char * file_name)
        :PROGRAM_BANK(d1, file_name) // use the inheritance
{
     int i, j;
     output_file << "+++ START +++++++++++++++++++++++++++++++++++++++++" << endl;
           
     if(d1 < 0 || d1 > 20 || d2< 0 || d2 >20)
     {
           output_file <<"+++ OUTPUT >>> INPUT ERROR"
           	       << endl
           	       <<"+++ OUTPUT >>> UNABLE TO CREATE MATRIX"
           	       << endl
           	       <<"+++ OUTPUT >>> EXITING PROGRAM"
           	       << endl
          	       << "+++ END +++++++++++++++++++++++++++++++++++++++++++"
          	       << endl;
           exit(1);
     }    
     else
     {
          dim1 = d1;
          dim2 = d2;
          
          // copy the input file name to private variable in_file
	  strcpy(in_file, file_name);

	  //define the input_filein this current block
 	  ifstream input_file(in_file, ios::in);
         
          
          for (i =0; i<dim1; i++)
          {
              for(j = 0; j<dim2; j++)
              {
                 input_file >> A[i][j];
              }
          }
          
          output_file <<"+++  _OUTPUT >>> CREATED A " << dim1 <<" X " << dim2
                      <<" MATRIX FROM " << in_file << endl;

     } 
	
     output_file << "+++   END +++++++++++++++++++++++++++++++++++++++++++" << endl;

}

 MATRIX operator / (int c, MATRIX & a_matrix)
{
     int i ,j;
     MATRIX temp_matrix(a_matrix.dim1, a_matrix.dim2);
     
     output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++" << endl;
     
     a_matrix.PROGRAM_BANK::INVERT_MATRIX("OUT_5_TEMP");
     
     
     ifstream input_file_temp("OUT_5_TEMP", ios::in);
     
     for(i = 0; i<temp_matrix.dim1; i++)
     {
         for(j = 0; j<temp_matrix.dim2; j++)
         {
             input_file_temp >> temp_matrix.A[i][j];
             temp_matrix.A[i][j] *= c;
         }
     }
     
     output_file <<"+++  _OUTPUT >>> THE RESULT OF " << c  << "/"
                 <<" X OPERATION IS:" 
                 << endl;
    
     for(i = 0; i<temp_matrix.dim1; i++)
     {
        for(j = 0; j<temp_matrix.dim2; j++)
        {
           output_file <<setprecision (2) << setiosflags(ios::fixed | ios::showpoint)
                       <<temp_matrix.A[i][j] << " ";                       
        }
        output_file << endl;
    }
     
    output_file <<"+++   END +++++++++++++++++++++++++++++++++++++++++++" << endl;

    return temp_matrix;
}

MATRIX operator * (int c, MATRIX & a_matrix)
{
     int i,j;
     
     MATRIX temp_matrix(a_matrix.dim1, a_matrix.dim2);
     
     output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++" << endl;
     
     
     for(i = 0; i<temp_matrix.dim1; i++)
     {
         for(j = 0; j<temp_matrix.dim2; j++)
         {
             temp_matrix.A[i][j] =  a_matrix.A[i][j];
             temp_matrix.A[i][j] *= c;
         }
     }
     
     output_file <<"+++  _OUTPUT >>> THE RESULT OF "<< c <<"* X OPERATION IS:"<< endl;
     
     for(i = 0; i<temp_matrix.dim1; i++)
     {
         for(j = 0; j<temp_matrix.dim2; j++)
         {
            output_file <<setprecision (2) << setiosflags (ios::fixed | ios::showpoint)
                        <<temp_matrix.A[i][j] <<" ";
         }
         
         output_file << endl;
     }
     
     output_file <<"+++   END +++++++++++++++++++++++++++++++++++++++++++" << endl;
     
     return temp_matrix;
}

void
MATRIX::PRINT( )
{
     output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
     		<< endl
     		<<"+++  _OUTPUT >>> CONTENTS OF MATRIX IS:"
     		<< endl;
     
     for (int i = 0; i<dim1; i++)
     {
         for (int j= 0; j<dim2; j++)
         {
             output_file <<setprecision(2) << setiosflags(ios::fixed |
             ios::showpoint) << A[i][j] << " ";
         }
         output_file << endl;
     }
	
     output_file<< "+++   END +++++++++++++++++++++++++++++++++++++++++++" << endl;


}

int
MATRIX::operator < (int b)
{
     int i, j, min_found;
     
     min_found=0;
     for(i= 0; i<dim1; i++)
     {
         for(j = 0; j<dim2; j++)
         {
            if(A[i][j] < b)
            {
               min_found = 1;
            }
            else
            {}
         }
     }
     
     if (min_found == 1)
     {
          output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
          	      << endl
          	      <<"+++  _OUTPUT >>> OVERLOADED OPERATOR < RETURNING TRUE"
          	      << endl
          	      <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
          	      << endl;
          return 1;
     }
     else
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR < RETURNING FALSE"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
	 return 0;
     }
         
     
}


void
MATRIX::operator = (MATRIX & a_matrix)
{
    if(dim1 != a_matrix.dim1 || dim2 != a_matrix.dim2)
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> ERROR"
         	     << endl
         	     <<"+++  _OUTPUT>>> INCOMPATIBLE MATRICES"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
     }
     else if( (dim1 == a_matrix.dim1) && (dim2 == a_matrix.dim2))
     {
         for(int i= 0; i<dim1; i++)
         {
             for(int j = 0; j<dim2; j++)
             {
                 A[i][j] = a_matrix.A[i][j];
                           
             }
         }
             
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR = COMPLETED"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
     }
     else{}
}

void
MATRIX::operator += (MATRIX & a_matrix)
{
     if(dim1 != a_matrix.dim1 || dim2 != a_matrix.dim2)
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> ERROR"
         	     << endl
         	     <<"+++  _OUTPUT >>> INCOMPATIBLE MATRICES"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
    
     }
     
     else if( (dim1 == a_matrix.dim1) && (dim2 == a_matrix.dim2))
     {
         for(int i= 0; i<dim1; i++)
         {
             for(int j = 0; j<dim2; j++)
             {
                 A[i][j] += a_matrix.A[i][j];
                           
             }
         }
             
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR += COMPLETED"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
     }
     else{}
}

void
MATRIX::operator *= (MATRIX & a_matrix)
{
     if(dim2 !=  a_matrix.dim1) 
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> ERROR"
         	     << endl
         	     <<"+++  _OUTPUT >>> INCOMPATIBLE MATRICES"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
     }
     
     else if(dim2 == a_matrix.dim1)
     {
     	 int i, j, k;
         MATRIX temp_matrix(dim1, a_matrix.dim2);
         
        
         for( i= 0; i<dim1; i++)
         {
             for( j = 0; j<a_matrix.dim2; j++)
             {
                 for(k = 0; k<dim2; k++)
                 {
                    temp_matrix.A[i][j] += A[i][k] * a_matrix.A[k][j];
                 }       
             }
         }
         
         dim2 = a_matrix.dim2;
         
         for( i= 0; i<dim1; i++)
         {
             for(j = 0; j<dim2; j++)
             {
                 A[i][j] = temp_matrix.A[i][j];
                           
             }
         }
                      
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR *= COMPLETED"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
     }
     else{}
}

MATRIX
MATRIX::operator + (MATRIX & a_matrix)
{
     int i, j;
     MATRIX temp_matrix(dim1, dim2);
     
     if(dim1 != a_matrix.dim1 || dim2 != a_matrix.dim2)
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> ERROR"
         	     << endl
         	     <<"+++  _OUTPUT >>> INCOMPATIBLE MATRICES"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
         	     
       	 MATRIX temp1_matrix(20,20);
      	 return temp1_matrix;
     }
     
     else
     {
         for(i= 0; i<dim1; i++)
         {
             for(j = 0; j<dim2; j++)
             {
                 temp_matrix.A[i][j] = A[i][j] + a_matrix.A[i][j];
                           
             }
         }
             
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR + COMPLETED"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
         
         return temp_matrix; 
    }
}

MATRIX
MATRIX::operator * (MATRIX & a_matrix)
{
     int i, j, k;
     MATRIX temp_matrix(dim1, a_matrix.dim2);

     if(dim2 != a_matrix.dim1 )
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	    << endl
         	    <<"+++  _OUTPUT >>> ERROR"
         	    << endl
         	    <<"+++  _OUTPUT >>> INCOMPATIBLE MATRICES"
         	    << endl
         	    <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	    << endl;
         	    
         MATRIX temp1_matrix(20,20);
      	 return temp1_matrix;
     }
     
     else if( dim2 == a_matrix.dim1) 
     {
         for(i= 0; i<dim1; i++)
         {
             for(j = 0; j<a_matrix.dim2; j++)
             {
                 for(k = 0; k<dim2; k++)
                 {
                     temp_matrix.A[i][j] += A[i][k] * a_matrix.A[k][j];
                 }         
             }
         }
             
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR * COMPLETED"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
         	     
         return temp_matrix;
         
     }
     else
     {
     	 MATRIX temp2_matrix(20,20);
      	 return temp2_matrix;
     }
    
}

MATRIX
MATRIX::operator / (MATRIX & a_matrix)
{
     int i, j, k;
     
     a_matrix.PROGRAM_BANK::INVERT_MATRIX("OUT.5.TEMP");
     
     MATRIX inv_matrix(20,20);
     
     ifstream input_file("OUT.5.TEMP", ios::in);
     
     inv_matrix.dim1 = a_matrix.dim1;
     inv_matrix.dim2 = a_matrix.dim2;
     
     for(i= 0; i<inv_matrix.dim1; i++)
     {
         for(j = 0; j<inv_matrix.dim2; j++)
         {
             input_file >> inv_matrix.A[i][j];
         }
     }
          
     
     
     if(dim2 != inv_matrix.dim1)
     {
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> ERROR"
         	     << endl
         	     <<"+++  _OUTPUT >>> INCOMPATIBLE MATRICES"
         	     << endl
         	     <<"+++   END +++++++++++++++++++++++++++++++++++++++++++"
         	     << endl;
         	     
      	 MATRIX temp1_matrix(20,20);
      	 return temp1_matrix;
     }
     
     
     
     else if( dim2 == inv_matrix.dim1)
     {
         MATRIX temp_matrix(dim1, inv_matrix.dim2);
         
         for(i = 0; i<dim1; i++)
         {
             for(j  = 0; j<inv_matrix.dim2; j++)
             {
                 for(k =0; k<dim2; k++) 
                 {
                    temp_matrix.A[i][j] += A[i][k] * inv_matrix.A[k][j];
                 }         
             }
         }
             
         output_file <<"+++   START +++++++++++++++++++++++++++++++++++++++++"
         	     << endl
         	     <<"+++  _OUTPUT >>> OVERLOADED OPERATOR / COMPLETED"
         	     << endl
         	     <<"+++  _OUTPUT >>> THE RESULT OF X/Y OPERATION IS:"
         	     << endl;
	 
	 for ( i = 0; i<temp_matrix.dim1; i++)
         {
             for ( j= 0; j<temp_matrix.dim2; j++)
             {
                 output_file <<setprecision(2) << setiosflags(ios::fixed |
                 ios::showpoint) << temp_matrix.A[i][j] << " ";
             }
             output_file << endl;
         }
	 output_file <<"+++   END +++++++++++++++++++++++++++++++++++++++++++" << endl;
	
	
	 return temp_matrix;
     }
     else
     {
     	 MATRIX temp2_matrix(20,20);
      	 return temp2_matrix;
     }
     
}
