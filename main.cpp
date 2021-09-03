//Gvozdev 26.11.2019
#include "jordan.h"

int main(int argc, char ** argv)
{
    int n,m,p,k/*,*index*/; //Ax=b
    double * a, * b, * x , * save/*, norm*/;//initialize *b *x
    char * name = 0;
    /*cpu_set_t cpu;*/
    pthread_t *threads;
    pthread_barrier_t barrier;
    pthread_mutex_t mutex;
    Args *args;

    if ( argc<4 || argc > 5 || ((n = atoi(argv[1])) <= 0) || ((m = atoi(argv[2])) <= 0) || ((p = atoi(argv[3])) <= 0) )
    {
        printf("usage: %s n m p [name]\n", argv[0]); //n division matrix a
        return 1;
    }
    else if(argc == 5) { name = argv[4]; }

    if( m%3 != 0 )
    {
        printf("Error: m/3 != 0\n"); //m division block
        return 2;
    }

    if( pthread_mutex_init(&mutex, 0) )
    {
        cout<<"Cannot init mutex\n";
        return MUTEX_INIT;
    }

    a = new double [n*n];
    b = new double [n];
    x = new double [n];
    save = new double [n];
    threads = new pthread_t [p];
    args = new Args [p];
    k = n/m;

    if ( !a || !b || !x || !save )
    {
        printf("Not enough memory\n");
        delete [] a; delete [] b; delete []x; delete []save; delete [] threads; delete [] args ;
//        delete [] index; // Need it? No all from 27 to 32 is not.
        return 3;
    }

    int res = pthread_barrier_init (&barrier, nullptr, p);
    if (res != 0)
    {
        perror ("pthread_barrier_init");
        delete [] a; delete [] b; delete []x; delete []save; delete [] threads; delete [] args ;
        return BAD_PTHREAD;
    }

//    double t = get_full_time();

    for(int i = 0 ; i < p ; i ++ )/* in future last line to give p-1 thread */
    {
        args[i].mutex=&mutex;
        args[i].barrier=&barrier;
        args[i].p=p;
        args[i].num=i;
        args[i].k=k;
        args[i].l=n%m;
        args[i].name=name;
        args[i].x=x;
        args[i].n=n;
        args[i].m=m;
        args[i].a=a;
        args[i].b=b;
        res = pthread_create(&threads[i], 0, solve, (void*) (args+i));

        if (res != 0)
        {
            cout<<"pthread_create"<<endl;
            abort ();
        }
    }

    for( int i = 0 ; i < p ; i ++ )
    {
        pthread_join(threads[i],0);
    }
//    t = get_full_time() - t;
    pthread_barrier_destroy (&barrier);

    if( pthread_mutex_destroy(&mutex) )
    {
        cout<<"Cannot destroy mutex\n";
        return MUTEX_DESTROY;
    }

    if(args->error < 0)
    {
        switch (args->error)
        {
            case N_LESS_THAN_M:
                printf("Error n < m \n");
                break;
            case MATRIX_DEGENERATE:
                printf("Matrix degenerate \n");
                break;
            case ERR_CANNOT_OPEN:
                printf("Cannot open file %s \n",name);
                break;
            case ERR_READ:
                printf("Error read file %s \n", name);
                break;
            default:
                printf("Unknown error\n");
        }
        delete [] a; delete [] b; delete []x; delete []save; delete [] threads; delete [] args ;
        return 0 ;
    }

    for( int i =0 ; i < min(n,8) ; i ++ )
    {
        for( int j =0 ; j < min(n,8) ; j ++)
        {
            cout<<a[i*n+j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    double error = 0; // ||b-x||
    double residual = 0 ;
    for(int i = 0 ; i < p ; i ++ )
    {
        if(i==0)
        {
            error = args[0].error_solution;
            residual = args[0].residual;
        }
        if(error < args[i].error_solution)
        {
            error = args[i].error_solution;
        }
        if(residual < args[i].residual)
        {
            residual = args[i].residual;
        }
    }

    printf("Residual = %e  Error = %e  Total time: %.2lf  n = %d  m = %d  p = %d\n", residual, error, args[p-1].total_time, n, m, p);
    for(int i = 0 ; i < p ; i ++ )
    {
        printf("Time of %d thread = %.2lf\n", i , args[i].cpu_time );
    }

    delete [] a; delete [] b; delete []x; delete []save; delete [] threads; delete [] args ;
//    delete [] index;

    return 0;
}

void init_x (double * x , int n )
{
    for ( int i=0 ; i < n ; i++ ) { x[i] = i%2 ; }
}

void init_rhs(double *a, double *b, double *x, int n) //b=A*x
{
    for ( int i=0 ; i < n ; i++ )
    {
        b[i]=0;
        for ( int j = 0 ; j < n ; j++ )
        {
            b[i] += a[i*n+j] * x[j];
        }
    }

}

double compute_error(double *b, double *x,int n)// ||b-x||
{
    double y1,y2=0;
    for( int i = 0 ; i < n ; i ++ )
    {
        y1=fabs(b[i]-x[i]);
        if( y1 > y2 ){ y2 = y1 ; }
    }
    return y2;
}

double compute_resideal( double * a , double * b , double * x , int n )// ||Ab-x||
{
    double y, h =0,h2;
    for( int i = 0; i < n ; i ++ )
    {
        y=0;
        for(int j = 0 ; j < n ; j ++ )
        {
            y += a[ i*n+j ] * b[ j ] ;
        }
        y = y - x[i] ;
//        h = h + fabs(y); //another norm
        h2 = fabs(y);
        if( h2 > h ){ h = h2 ; }
    }
    return h;
}

void init_matrix(double *a, int n)
{
    int k;
    //srand(time(NULL));
    for(int i=0;i<n;i++)
    {
        k=i*n;
        for(int j=0;j<n;j++)
        {
            a[k+j]=fabs(i-j);//for example
            //a[k+j]=(double)(rand()%1001)/1000;
//            a[k+j]=1./(1+i+j);
            //a[k+j] = n-max(i,j);
        }
    }
}
