//Gvozdev 26.11.2019
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sched.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sched.h>
#include <pthread.h>
#include <cfloat>
#include <cstring>

#define ERR_CANNOT_OPEN      -1
#define ERR_READ             -2
#define ERR_UNKNOWN_ERROR    -3
#define N_LESS_THAN_M        -4
#define BLOCK_DEGENERATE     -5
#define MATRIX_DEGENERATE    -6
#define NOT_ENOUGH_MEMORY    -7
#define BAD_PTHREAD          -8
#define MUTEX_INIT           -9
#define MUTEX_DESTROY        -10

#define EPS 1e-16

#define LOG(...) std::cerr<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"

class Args
{
	public:
	    int p = 0;//# threads
	    int n = 0;//hight and width of matrix
	    int m = 0;//size of block
	    int k = 0;//[n/m]
	    int l = 0;//n%m
	    int num = 0;//num of the thread
	    int error = 0;//for error in process of solve
//	    int * index = 0;//array of index
	    int save_int = 0;//for another operations
	    int local_thread = 0;//#thread belongs this thread
	    double save_double = 0;
	    double * a = 0;//matrix a
	    double * b = 0;//vector for solution
	    double * x = 0;//vector of 0 and 1
	    double cpu_time = 0;//time of the thread
	    double total_time = 0;//time on solve
	    double residual = 0;
	    double error_solution = 0;
	    double norma = 0;
	    char * name = 0;//name of file with matrix *a
	    pthread_barrier_t *barrier = 0;//one for all
	    pthread_mutex_t *mutex = 0;//one for all
};


int read_matrix( double *a, int n , char * name);
void init_matrix(double *a, int n );
void init_x (double * x , int n );
void init_rhs(double *a, double *b, double *x, int n);
double compute_error(double *b ,double *x ,int n);
void init_matrix_file(double *a, int n, char *name);
void print_matrix_file(double *a, int n);
double compute_resideal(double *a,double *b, double *x,int n);

void* solve( void *arg );
int initialize_matrix(Args * arg , int * fatal_error);
void formula(double * a, int n, Args * arg);
void matrix_mult_vector( double * a , double * b , double * for_count , int n ,Args * arg );
void get_block(int i, int j, int width, int hight, int real_size_c , int n, double *c, double *a);
void put_block(int i, int j, int width, int hight, int real_size_c , int n, double *c, double *a);
void make_E( int real_size_c , double * c);
void make_0( int hight , int width , int real_size_c , double * c );
void swap_line( int i , int num , int * perest , int n , int m );

void multiply_line_on_inverse( int i , int column , double * a , double * b , double * c0 ,
                              double * c1 , double * c2 , int n , int m , int k , int p ,
                              double norm , int * perest_block , Args *arg );

void substract_lines( int i , int column , double * a , double * b , double * c0 ,
                      double * c1 , double * c2 , int n , int m , int k , int p ,
                      Args *arg , double * for_count , int size );
double norm_matrix ( double * a , int n , Args * arg);
void matrix_mult_matrix (double *a, double *b,double *c, int n , int N);
int search_inverse_matrix ( double * c0 , double * c1 , int m , double norm_eps_matrix , int * perest_block );//write inverse matrix for c0 in c1 . Jordan

void search_number_min_inverse_norm_in_column( int j , double * a , double * c0 , double * c1 ,
                                              int n , int m , int k , double norm , int * perest_block,
                                              Args * arg , int * index ); // n == k * m + p

void block_mult_block( double * c0 , double * c1 , double * c2 , int n );//n%3==0 c=a*b
/*c0 = c1 * c2 ; c0 - n*m ; c1 - n*k ; c2 - k*m ;
1st number - #lines ; 2d number - #column*/
void substract_block( double * c0 , double * c1 , double * c2 , int width , int hight , int real_size_c );//c0=c1-c2 hight*width
void recovery_solution( double * b , double * x , int k , int m , int p , Args * arg , int * index);
void print_matrix( double * a , int n );
template< typename T >
void print_vector( T * b , int n );
double get_full_time();
double get_time();
void search_error( double * x , double * b /*, int n */, Args * arg);
void search_residual( double * a , double * x , double * b , int n , Args * arg );

void init_matrix_paral( Args * arg );
void read_matrix_paral( Args * arg );
void reduce_max(int p, int *num_line,  double norm, int *lines, int num_column);

using namespace std;
