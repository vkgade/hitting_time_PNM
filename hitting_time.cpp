#include<iostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>

#define step  1000
using namespace std;

char *set_start_file (void) {
  char *pt;
  pt = (char*) malloc (25);
  cout << "enter here --> ";
  cin  >> pt ;
  return pt;
}

const char *file_name = set_start_file();

int count_nodes(void) {
  unsigned int number_of_lines = 0;
  FILE *infile = fopen(file_name, "r");
  int ch;
  while (EOF != (ch=getc(infile)))
    if ('\n' == ch)
      ++number_of_lines;
  
  return number_of_lines;  
}

int nodes= count_nodes();

/**
 * PART@1
 */ 

int main() {
  cout << nodes << endl ;
  cout << file_name << endl;
  
  double **m,**a,sum, *sum_ROW;
  
  m= (double**) malloc (nodes*sizeof(double*));
  a= (double**) malloc (nodes*sizeof(double*));
  for (int x=0;x<nodes;x++) {
    *(m+x)= (double*) malloc (nodes*sizeof(double)) ; 
    *(a+x)= (double*) malloc (nodes*sizeof(double)) ;
  }

  sum_ROW = (double*) malloc (nodes*sizeof(double));
  
  if (remove("normalized_affinity_MAT.txt")!=0 || 
      remove("diagnol_MAT.txt")!=0 || 
      remove("THE_MATRIX.txt")!=0 || 
      remove("HITTING_TIME.txt")!=0) {
    perror ("Old output files could not be deleted or do not exist !!") ;
    cout << endl ;
  }
  else 
    cout << "All old outputfiles are deleted." << endl;  

  int i,j;
  
  ifstream inputFILE_affinity;
  inputFILE_affinity.open(file_name);
  if (!inputFILE_affinity) {
    cerr << "ERROR: File could not be opened" << endl ;
    exit(1);
  }
  while (!inputFILE_affinity.eof()) {
    for (j=0;j<=nodes-1;j++) {
      for (i=0;i<=nodes-1;i++) {
	inputFILE_affinity >> a[j][i];
      } 
    }  
  } inputFILE_affinity.close();    

  for (j=0;j<=nodes-1;j++) {
    sum =a[j][0];
    for (i=1;i<=nodes-1;i++) {
      sum = sum + a[j][i];
    } sum_ROW[j]=sum ;
  }

  ofstream outputFILE_normalized;
  outputFILE_normalized.open("normalized_affinity_MAT.txt",ios::app);
  for (j=0; j<=nodes-1;j++) {
    for (i=0; i<=nodes-1;i++) {
      m[j][i]=a[j][i]/sum_ROW[j];
      outputFILE_normalized << m[j][i] << "\t";  
    } outputFILE_normalized << endl;  
  }  outputFILE_normalized.close();
  
  ofstream outputFILE_diagnol;
  outputFILE_diagnol.open("diagnol_MAT.txt",ios::app);
  for (j=0; j<=nodes-1;j++) {
    for (i=0 ;i<=nodes-1;i++) {
      if (j!=i)  
	outputFILE_diagnol << 0 << "\t";
      else  
	outputFILE_diagnol << 1.0/sum_ROW[j] << "\t";
    }  outputFILE_diagnol << endl; 
  } outputFILE_diagnol.close(); 
  
  for (i=0;i<nodes;i++)
    free(m[i]);
  free(m);
  
  for (i=0;i<nodes;i++)
    free(a[i]);
  free(a);
   
  cout << "\nEND OF PART_1 !!" << endl ;

/**
 * PART@2
 */ 

  double *a_data, *b_data;
  a_data = (double*) malloc(((nodes-1)*(nodes-1)+1)*sizeof(double));
  b_data = (double*) malloc ((nodes-1)*sizeof(double));
  
  int k;
  for (k=0;k<nodes-1;k++) {
    b_data[k] = 1;
  }

  if (remove("nodes_ALL_INV_MAT.txt")!=0)
    perror("OLD nodes_ALL_INV_MAT.txt could not be deleted or file does not exist !!");
  else
    cout << "file deleted \n" ;
  double **M, **A ; 
  int z,e;
  
  M= (double**) malloc((nodes)*sizeof(double*));
  for (z=0;z<nodes;z++) {
    *(M+z)= (double*) malloc((nodes)*sizeof(double));  
  }

  ifstream the_M_MATRIX;
  the_M_MATRIX.open("normalized_affinity_MAT.txt");
  
  if (!the_M_MATRIX) {
    cerr << " the NORMALIZED MATRIX could not be found! Run part1.cpp again !" << endl ;
    exit (1);
  }

  while (!the_M_MATRIX.eof()) {
    for (j=0;j<nodes;j++) {
      for (i=0;i<nodes;i++) {
	the_M_MATRIX >> *(*(M+j)+i); 
      } 
    }
  } the_M_MATRIX.close();

  int s,t;
  ofstream all_nodes_matrix;
  all_nodes_matrix.open("nodes_ALL_INV_MAT.txt",ios::app);
  
  for (j=0;j<nodes;j++) {
    A= (double**) malloc ((nodes-1)*sizeof(double*));
    for (e=0;e<nodes-1;e++) {
      *(A+e)= (double*) malloc((nodes-1)*sizeof(double));
    }

    for (t=0;t<nodes-1;t++) {
      if (t<j) {
	for (s=0;s<nodes-1;s++) {
	  if (s<j) {
	    if (s==t) 
	      *(*(A+s)+t) = 1; 
	    else  
	      *(*(A+s)+t) =  -(*(*(M+t)+s)) ; 
	  } else  {
	    *(*(A+s)+t) =  -(*(*(M+t)+(s+1))) ; 
	  }
	}
      }  
 
      else {
	for (s=0;s<nodes-1;s++) {
	  if (s<j) {
	    *(*(A+s)+t) =  -(*(*(M+(t+1))+s));   
	  } else {
	    if (s==t) 
	      *(*(A+s)+t)= 1; 
	    else  
	      *(*(A+s)+t) =  -(*(*(M+(t+1))+(s+1))) ;
	  }
        }
      }
    }

    for (s=0;s<nodes-1;s++) {
      for (t=0;t<nodes-1;t++) {
	a_data[t+((nodes-1)*s)]=A[t][s] ;
      }
    }

    ofstream outputFILE;
    outputFILE.open("THE_MATRIX.txt",ios::app);
    gsl_matrix_view m = gsl_matrix_view_array (a_data, nodes-1, nodes-1);
    gsl_vector_view b = gsl_vector_view_array (b_data, nodes-1);
    gsl_vector *x = gsl_vector_alloc (nodes-1);
    int s;
    gsl_permutation * p = gsl_permutation_alloc (nodes-1);     
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    
    for(i=0;i<nodes-1;i++) { 
      outputFILE << gsl_vector_get(x,i) << "\t" ;
    }
    gsl_vector_fprintf (stdout, x, "%g"); /* comment this line for a clean stdout */        
    gsl_permutation_free (p);
    gsl_vector_free (x); 
    outputFILE << endl;
    outputFILE.close(); 
    
    
    all_nodes_matrix << endl ;
    for (i=0;i<nodes-1;i++) {
      free (A[i]);
    } free(A);
  } all_nodes_matrix.close();
         
  for (i=0;i<nodes;i++)
    free(M[i]);
  free(M);

  cout << "\nEND OF PART_2 !!" << endl;

/**
 * PART@2
 */ 

  int u ;
  double **DBZ, **H;
  DBZ= (double**) malloc (nodes*sizeof(double*));
  for (u=0;u<nodes;u++) {
    *(DBZ+u)= (double*) malloc (nodes*sizeof(double)) ;
  }
  H= (double**) malloc (nodes*sizeof(double*));
  for (u=0;u<nodes;u++) {
    *(H+u)= (double*) malloc (nodes*sizeof(double)) ;
  }

  ifstream the_matrix;
  the_matrix.open("THE_MATRIX.txt");
  if (!the_matrix)
    cerr << "The Matrix is gone !!" << endl ;
  else 
    cout << "The Matrix is good !!" << endl ;
  
  while (!the_matrix.eof()) {
    for (j=0;j<nodes;j++) {
      for (i=0;i<nodes-1;i++) {
	the_matrix >> DBZ[j][i] ;
      }
    }
  } the_matrix.close();

  ofstream matrix_revolution_final_part;
  matrix_revolution_final_part.open("HITTING_TIME.txt",ios::app);
  for (j=0;j<nodes;j++) {
    for (i=0;i<nodes;i++) {
      if (j>=i){
	if (i==j)
	  H[j][i]=0;
	else 
	  H[j][i]=DBZ[j][i];
      } else 
	H[j][i]=DBZ[j][i-1];

      matrix_revolution_final_part << H[j][i] << "\t" ;
    } matrix_revolution_final_part << endl; 
  } matrix_revolution_final_part.close();
  
  for (i=0;i<nodes;i++)
    free(H[i]);
  free(H);
  
  for (i=0;i<nodes;i++)
    free(DBZ[i]);
  free(DBZ);
  
  cout << "\nTHE END OF PART_3" << endl ;
  return 0;
}
