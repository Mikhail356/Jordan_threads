//Gvozdev 26.11.2019
#include "jordan.h"
void * solve( void * args )
{
    Args *arg = (Args*) args;

//    arg->cpu_time = get_time();

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET( get_nprocs() - arg->num - 1,&cpu);
    pthread_setaffinity_np(pthread_self(),sizeof(cpu),&cpu);

    static int fatal_error = 0;
    int /*num,*/ k, p, *perest_block, n, m;
    double *c0 , *c1 , *c2 /*, *x*/;
    double norm , y;
    static double max = 0;
    static int number = 0;
    double * for_count = nullptr;
    int * index = nullptr;

    k = arg->k;
    p = arg->l;
    n=arg->n;
    m=arg->m;
    int size = (/*k*m+p*/n+1);

    arg->local_thread = (k-arg->num)/arg->p;
    if((k-arg->num)%arg->p!=0){ arg->local_thread++ ; }
    if( k/arg->p==0 )
    {
        if(k%arg->p<arg->num){arg->local_thread = 0 ;}
    }

    initialize_matrix(arg, &fatal_error);
    pthread_barrier_wait(arg->barrier);
    if(fatal_error < 0)
    {
        arg->error = fatal_error;
        return 0 ;
    }

    if( n < m ) { arg->error = N_LESS_THAN_M; return 0; }

    for_count = new double [size*m];
    perest_block = new int [m];
    c0 = new double [m*m];
    c1 = new double [m*m];
    c2 = new double [m*m];
    index = new int [k];
    for(int i =0 ;i < k ; i ++ ){   index[i] = -1;  }
    if( !perest_block || !c0 || !c1 || !c2 )
    {
        fatal_error = NOT_ENOUGH_MEMORY;
    }
    pthread_barrier_wait( arg->barrier );
    if( fatal_error == NOT_ENOUGH_MEMORY ) { arg->error = NOT_ENOUGH_MEMORY; return 0 ;}

    init_x(for_count,n);
    matrix_mult_vector( arg->a , arg->b , for_count , n , arg ); //right hand side t e vector b Ax==b

    norm = EPS*norm_matrix( arg->a , n , arg );

    for( int i = 0 ; i < m ; i ++ ) { perest_block[i] = i ; }

    make_0(0,0,m,c2);
    make_0(0,0,m,c1);
    make_0(0,0,m,c0);

    if(arg->num==arg->p-1){arg->total_time = get_full_time();}
    arg->cpu_time = get_time();
    for( int i = 0 ; i < k ; i ++ ) // n == k * m + p // loop through columns
    {
//        LOG(i);
//        print_matrix(arg->a,n);
//        print_vector(arg->b,n);

        search_number_min_inverse_norm_in_column( i , arg->a , c0 , c1 , n , m , k ,norm ,
                                                            perest_block , arg , index);/*from 0 to k block
                                                                            min inverse norm in arg->norma*/
//        LOG(arg->save_int);
//        print_matrix(arg->a,n);
//        print_vector(arg->b,n);
//        LOG(arg->num);
        if( arg->error == MATRIX_DEGENERATE )
        {
            delete [] for_count ;
            delete [] perest_block ;
            delete [] c0 ;
            delete [] c1 ;
            delete [] c2 ;
            delete [] index ;
            return 0 ;
        }
/*-------------------------------------------------------sinhronizer---------------------------------------------*/
        if( arg->save_int != -1 )
        {
            max = arg->norma;/*will be initialize in search max in column*/
            number = arg->save_int;/*number of line*/
        }

        pthread_barrier_wait(arg->barrier);
        pthread_mutex_lock(arg->mutex);
        if( arg->norma < max && arg->save_int != -1 )
        {
            max = arg->norma ;
            number = arg->save_int;
        }
        pthread_mutex_unlock(arg->mutex);
        pthread_barrier_wait(arg->barrier);
//    pthread_mutex_lock(arg->mutex);
        index[i] = number ;
//        LOG(index[i]);
//    pthread_mutex_unlock(arg->mutex);
//    pthread_barrier_wait(arg->barrier);
/*-------------------------------------------------------sinhronizer---------------------------------------------*/

        multiply_line_on_inverse( index[i] , i , arg->a , arg->b , c0 , c1 , c2 , n , m , k , p , norm , perest_block , arg );

        pthread_mutex_lock(arg->mutex);
        print_matrix(arg->a,n);
        print_vector(arg->b, n);
        pthread_mutex_unlock(arg->mutex);
//        pthread_mutex_lock(arg->mutex);
//        LOG(arg->index[i]);
//        pthread_mutex_unlock(arg->mutex);
        pthread_barrier_wait(arg->barrier);

        substract_lines( index[i] , i , arg->a , arg->b , c0 , c1 , c2 ,
                             n , m , k , p , arg , for_count , size );

        pthread_barrier_wait(arg->barrier);

        pthread_mutex_lock(arg->mutex);
        print_matrix(arg->a,n);
        print_vector(arg->b, n);
        pthread_mutex_unlock(arg->mutex);
    }
//    print_matrix(arg->a,n);
    for( int i = arg->num ; i < n ; i += arg->p )
    {
        y = 0;
        for( int j = 0 ; j < n ; j ++ )
        {
            y += fabs( arg->a[i*n+j] );
        }

        if(y <= norm)
        {
            fatal_error = MATRIX_DEGENERATE;
        }
    }

    pthread_barrier_wait(arg->barrier);
    if(fatal_error == MATRIX_DEGENERATE)
    {
        delete [] for_count ;
        delete [] perest_block ;
        delete [] c0 ;
        delete [] c1 ;
        delete [] c2 ;
        delete [] index ;
        arg->error = MATRIX_DEGENERATE ;
        return 0;
    }

    if( p != 0 )//for down right side matrix a and vector b
    {
        if(arg->num == (k%arg->p))
        {
            get_block(k*n*m , k*m , p , p , m , n , c0 , arg->a);

            int schetchik = 0;
            for ( int i = 0 ; i < m ; i++ )
            {
                for( int j =0 ; j < m ; j++)
                {
                    if( i < p && j < p)
                    {
                        c0[ schetchik ]=c0[ i*m+j ];
                        schetchik ++ ;
                    }
                }
            }

            if( search_inverse_matrix( c0 , c1 , p , norm , perest_block ) == BLOCK_DEGENERATE)
            {
                fatal_error = MATRIX_DEGENERATE;
                arg->error = MATRIX_DEGENERATE;
            }

            if( arg->error == 0)
            {
                schetchik=0;
                for ( int i = 0 ; i < m ; i++ )
                {
                    for( int j =0 ; j < m ; j++)
                    {
                        c1[ i*m+j ]=c0[ schetchik ];
                        if( i < p && j < p)
                        {
                            schetchik ++ ;
                        }
                    }
                }
                make_0( p , p , m , c1 );
                get_block(k*m , 0 , 1 , p , m , 1 , c0 , arg->b);
                make_0( p , 1 , m , c0 );

                block_mult_block( c2 , c1 , c0 , m );//c2=c1*c0
                put_block(k*m , 0 , 1 , p , m , 1 , c2 , arg->b);

                make_E( m , c2);
                put_block(k*n*m , k*m , p , p , m , n , c2 , arg->a);//t k A^(-1)*A==E
            }
        }
        pthread_barrier_wait(arg->barrier);

        if(fatal_error==MATRIX_DEGENERATE)
        {
            delete [] for_count ;
            delete [] perest_block ;
            delete [] c0 ;
            delete [] c1 ;
            delete [] c2 ;
            delete [] index ;
            arg->error = MATRIX_DEGENERATE ;
            return 0 ;
        }

        for( int i = arg->num ; i < k ; i += arg->p )
        {
            get_block(i*n*m , k*m , p , m , m , n , c0 , arg->a);
            get_block( k*m  ,  0  , 1 , p , m , 1 , c1 , arg->b);
            make_0(m,p,m,c0);
            make_0(p,1,m,c1);
            block_mult_block( c2 , c0 , c1 , m);

            get_block(i*m , 0 , 1 , m , m , 1 , c1 , arg->b);

            substract_block( c0 , c1 , c2 , m , 1 , m);//c0=c1-c2

            put_block(i*m , 0 , 1 , m , m , 1 , c0 , arg->b);
        }
    }

    pthread_barrier_wait(arg->barrier);
    recovery_solution(arg->b , arg->x , k , m , p , arg , index);

    if(arg->num==arg->p-1){arg->total_time = get_full_time() - arg->total_time;}
    arg->cpu_time = get_time() - arg->cpu_time;
//    pthread_barrier_wait(arg->barrier);

    initialize_matrix(arg,&fatal_error);
    init_x(for_count,n);
    pthread_barrier_wait(arg->barrier);
    search_error( for_count , arg->x /*, n*/ , arg);
    search_residual( arg->a , for_count , arg->x , n , arg);

    delete [] for_count ;
    delete [] perest_block ;
    delete [] c0 ;
    delete [] c1 ;
    delete [] c2 ;
    delete [] index ;

    return 0;
}

int initialize_matrix(Args * arg , int * fatal_error)
{
    int n = arg->n , m = arg->m , k = arg->k , kol_pthread = arg->p , L = arg->l , num = arg->num;
    for(int i = num ; i < k ; i += kol_pthread)
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                arg->a[(i*m+j)*n+t]=0;
            }
        }
    }
    if(num==0&&L!=0)
    {
        for( int j = 0 ; j < L ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                arg->a[(k*m+j)*n+t]=0;
            }
        }
    }

    if (arg->name) //from file
    {
        if(num==0)
        {
            int res = read_matrix(arg->a, n, arg->name);
            if (res < 0)
            {
                *fatal_error = res;
            }
        }
        pthread_barrier_wait(arg->barrier);
        if(*fatal_error<0)
        {
            return 0 ;
        }
    }
    else // by formula
    {
        formula(arg->a, n, arg);
    }
    pthread_barrier_wait(arg->barrier);

    return 0 ;
}

int read_matrix(double *a , int n, char *name)
{
    ifstream file;
    file.open(name);

    if (!(file.is_open())) return ERR_CANNOT_OPEN;

    n=n*n;
    for(int i=0;i<n;i++)
    {
        if(!(file>>a[i])) return ERR_READ;
    }
    file.close();
    return 0;
}

void formula(double * a, int n, Args * arg)
{
    int m = arg->m , k = arg->k , kol_pthread = arg->p , L = arg->l , num = arg->num;

//    srand(time(NULL));
    for(int i = num ; i < k ; i += kol_pthread)
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                /*arg->*/a[(i*m+j)*n+t]=fabs(i*m+j-t);
//                arg->a[(i*m+j)*n+t]=(double)(rand()%1001)/1000;
//                arg->a[(i*m+j)*n+t]=1./(1+i*m+j+t);
//                arg->a[(i*m+j)*n+t] = n-max(i*m+j,t);
            }
        }
    }
    if(num==0&&L!=0)
    {
        for( int j = 0 ; j < L ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                /*arg->*/a[(k*m+j)*n+t]=fabs(k*m+j-t);
//                arg->a[(i*m+j)*n+t]=(double)(rand()%1001)/1000;
//                arg->a[(k*m+j)*n+t]=1./(1+k*m+j+t);
//                arg->a[(i*m+j)*n+t] = n-max(k*m+j,t);
            }
        }
    }
}

void matrix_mult_vector( double * a , double * b , double * x , int n , Args* arg ) //right hand side t e vector b Ax==b
{
    int m = arg->m , k = arg->k , kol_pthread = arg->p , L = arg->l , num = arg->num;
    double save = 0 ;

    for(int i = num ; i < k ; i += kol_pthread)
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            save = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                save += /*arg->*/a[(i*m+j)*n+t] * x[t];
            }
            /*arg->*/b[i*m+j] = save ;
        }
    }
    if( num==0&&L!=0)
    {
        for( int j = 0 ; j < L ; j ++ )
        {
            save = 0 ;
            for( int t = 0 ; t < n ; t ++ )
            {
                save += /*arg->*/a[(k*m+j)*n+t] * x[t];
            }
            /*arg->*/b[k*m+j] = save ;
        }
    }
}

void search_error( double * x , double * b ,/* int n ,*/ Args * arg)
{
    double y1,y2=0;
    int schet = 0 ;
    for(int i = arg->num ; i < arg->k ; i += arg->p )
    {
        schet = i*arg->m;
        for( int j = 0 ; j < arg->m ; j ++ )
        {
            y1=fabs(b[schet+j]-x[schet+j]);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    if(arg->num==0&&arg->l!=0)
    {
        schet = arg->k*arg->m;
        for( int j = 0 ; j < arg->l ; j ++ )
        {
            y1=fabs(b[schet+j]-x[schet+j]);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    arg->error_solution = y2;
}

void search_residual( double * a , double * x , double * b , int n , Args * arg )
{
    double y1=0,y2=0,save1=0,save2=0;
    int schet = 0 ;
    for(int i = arg->num ; i < arg->k ; i += arg->p )
    {
        schet = i*arg->m;
        for( int j = 0 ; j < arg->m ; j ++ )
        {
            save1 = 0; save2 = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                save1 += /*arg->*/a[(schet+j)*n+t] * x[t];
                save2 += /*arg->*/a[(schet+j)*n+t] * b[t];
            }
            y1=fabs(save1-save2);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    if(arg->num==0&&arg->l!=0)
    {
        schet = arg->k*arg->m;
        for( int j = 0 ; j < arg->l ; j ++ )
        {
            save1 = 0; save2 = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                save1 += /*arg->*/a[(schet+j)*n+t] * x[t];
                save2 += /*arg->*/a[(schet+j)*n+t] * b[t];
            }
            y1=fabs(save1-save2);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    arg->residual = y2;
}
void get_block(int i, int j, int width, int hight, int real_size_c, int n, double *c, double *a)
/*i,j koordinaty verhnego levogo ugla bloka c v matrice a
width, hight shirina i vysota bloka c
n size of matrix a */
{
    int p=0;
    p=i+j;

    for(int t=0;t<hight;t++)
    {
        for( int r = 0 ; r < width ; r ++ )
        {
            c[ t * real_size_c + r ] = a[ p + r ];
        }

        p+=n;
    }
}

void put_block(int i, int j, int width, int hight, int real_size_c , int n, double *c, double *a)
{
    int p;
    p=i+j;

    for(int t=0;t<hight;t++)
    {
        for(int r=0;r<width;r++)
        {
            a[p+r]=c[t*real_size_c+r];
        }

        p+=n;
    }
}

void make_0( int hight , int width , int real_size_c , double * c )
{
    for( int i = 0 ; i < real_size_c ; i ++ )
    {
        for( int j = 0 ; j < real_size_c ; j ++ )
        {
            if( i >= hight || j >= width )
            {
                c[ i * real_size_c + j ] = 0 ;
            }
        }
    }
}

void make_E( int real_size_c , double * c )
{
    for( int i = 0 ; i < real_size_c ; i ++ )
    {
        for( int j = 0 ; j < real_size_c ; j ++ )
        {
            c[ i * real_size_c + j ] = ( i != j ? 0 : 1 ) ;
        }
    }
}

double norm_matrix ( double * a , int n , Args * arg ) // a matrix   n * n , norm by line
{
    static double norm_of_matrix = 0;
    double max1 = 0 , max2 = 0 ;
    int schet = 0;
    int p = arg->p;/*# threads*/
    int m = arg->m;/*size of block*/
    int num = arg->num ; /*number of thread*/

//    pthread_mutex_lock(arg->mutex);
    p = p*m;
//    LOG("norm matrix");
//    LOG(arg->num);
//    LOG(n);
//    LOG(p);
//    LOG(m);
    for( int i = num*m ; i < n ; i += p )
    {
        schet = i * n ;
//        LOG(schet);
//        LOG(i);
        for( int temp = 0 ; temp < min(m,n-i) ; temp ++ )
        {
            max1 = 0;
            for( int j = 0 ; j < n ; j ++ )
            {
                max1 += fabs(a[ schet + temp*n + j ]);
            }

            if ( max1 > max2 ) {max2 = max1 ;}
        }
    }
//    LOG("ura");
    arg->save_double = max2 ;
    norm_of_matrix = max2 ;
//    LOG(arg->save_double);
//    pthread_mutex_unlock(arg->mutex);
    pthread_barrier_wait( arg->barrier );
    pthread_mutex_lock(arg->mutex);
    if ( arg->save_double > norm_of_matrix )
    {
        norm_of_matrix = arg->save_double ;
    }
    pthread_mutex_unlock(arg->mutex);
    pthread_barrier_wait( arg->barrier );

//pthread_mutex_lock(arg->mutex);
//LOG(norm_of_matrix);
//pthread_mutex_unlock(arg->mutex);

    return norm_of_matrix;
}
double norm_block( double * a , int n )
{
    double max1 , max2 = 0 ;
    int schet;

    for( int i = 0 ; i < n ; i ++ )
    {
        max1 = 0;
        schet = i * n ;

        for( int j = 0 ; j < n ; j ++ )
        {
            max1 = max1 + fabs(a[ schet + j ]);
        }

        if ( max1 > max2 ) {max2 = max1 ;}
    }
    return max2;
}

int search_inverse_matrix ( double * c0 , double * c1 , int m , double norm , int * perest_block )//write inverse matrix for c0 in c0 . Jordan
{
    double y,y0;
    int www=0 , scheti , schetj ;
    for( int i = 0 ; i < m ; i ++ )
    {
        perest_block[i]=i;
        for ( int  j = 0 ; j < m ; j ++ )
        {
            c1[ i * m + j ] = ( i != j ? 0 : 1 );
        }
    }

    for ( int i = 0 ; i < m ; i ++ )
    {
        y=0;
        for( int j = i ; j < m ; j ++ )
        {
            y0 = fabs(  c0[ perest_block[j] * m + i ] );
            if ( y0 > y )
            {
                y = y0;
                www = j;
            }
        }
        if( fabs(y) <= norm ) { return BLOCK_DEGENERATE ; }
        swap(perest_block [ i ] , perest_block [ www ]) ;

        scheti = perest_block [ i ] * m;
        y = 1./c0 [ scheti + i];

        for ( int  j = 0 ; j < m ; j ++ )
        {
            c0 [ scheti + j] *= y ;
            c1 [ scheti + j] *= y ;
        }

        for( int j = 0 ; j < m ; j ++ )//vichest is stroki j stroku www umnojuyu na koef j stroki
        {
            schetj = perest_block[j] * m;
            y=c0[ schetj + i ];
            www=0;

            for ( int t = 0 ; t < m ; t ++ )
            {
                if( scheti != schetj )
                {
                    c0[ schetj + t ] = c0[ schetj + t ] - ( y * c0[ scheti + t ] );
                    c1[ schetj + t ] = c1[ schetj + t ] - ( y * c1[ scheti + t ] );
                }
                if( fabs(c0[ schetj + t ]) < norm )
                {
                    www ++ ;
                }
            }

            if( www == m )
            {
                return BLOCK_DEGENERATE;
            }
        }
    }
    for ( int i = 0 ; i < m ; i ++ )
    {

        for( int j = 0 ; j < m ; j ++ )
        {
            c0 [ i * m + j ] = c1 [ perest_block[ i ] * m + j ] ;
        }
    }
    return 0;
}

void search_number_min_inverse_norm_in_column( int j , double * a , double * c0 , double * c1 ,
                                              int n , int m , int k , double norm , int * perest_block,
                                              Args * arg , int * index) // n == k * m + p
{
    int h,f3=0,t=0/*,kk=arg->k*/;
    static int u = 0;
    double f1,f2 = 0;
//    int num_thread = arg->num;
    int kol_pthread = arg->p;
    int local_count = 0;
//LOG(j);
//print_matrix(a,n);
//print_vector(arg->index,k);
    u=0; t=0;
    for(int i = arg->num  ; i < k ; i += kol_pthread )
    {
            get_block( i*n*m , j*m , m , m , m , n , c0 , a);

            h = search_inverse_matrix( c0 , c1 , m , norm , perest_block );

            if( h == BLOCK_DEGENERATE )
            {
                pthread_mutex_lock(arg->mutex);
                u ++ ;
                pthread_mutex_unlock(arg->mutex);
                local_count++;
            }
            else
            {
                f1 = norm_block( c0 , m );/*rewrite this plase*/
                h=0;
//                LOG(i);
                for(int www = 0 ; www <= i ; www++ )
                {
//                    LOG(index[www]);
//                    LOG(www);
                    if(index[www]==i){h=1; local_count++;}
                }
                if(h==0)
                {
                    if(t==0)
                    {
                        f2 = f1 ; f3 = i; t=1;
                    }
                    else if( f2 > f1 )
                    {
                        f3 = i ;
                        f2 = f1 ;
                    }
                }
            }
    }
    pthread_barrier_wait( arg->barrier );
//    pthread_mutex_lock(arg->mutex);
//    LOG(u);
    if ( u == k )
    {
        arg->error = MATRIX_DEGENERATE ;
    }
//    pthread_mutex_lock(arg->mutex);
//    LOG(n);
//    LOG(m);
//    LOG(arg->num);
//    LOG((n-(arg->num*m))/m);

    if( local_count < arg->local_thread )
    {
        arg->norma = f2;/*min norm in column*/
        arg->save_int = f3;/*num line with min norm*/
    }
    else
    {
        arg->save_int = -1;
    }
//    pthread_mutex_lock(arg->mutex);
//    LOG(arg->num);
//    LOG(local_count);
//    LOG(arg->local_thread);
//    LOG(arg->save_int);
//    LOG(f2);
//    LOG(f3);
//    pthread_mutex_unlock(arg->mutex);
}

void multiply_line_on_inverse( int i , int column , double * a , double * b , double * c0 ,
                               double * c1 , double * c2 , int n , int m , int k , int p ,
                               double norm , int * perest_block , Args * arg)// n == k * m + p
//i and j number lines of block
{
//    /*if(arg->num==(i%arg->p))
//    {*/
        //int step = arg->p;//need it?
//        int num_thread = arg->num;//need it?
//        LOG("hello");
//    pthread_mutex_lock(arg->mutex);
        get_block( i*n*m , column*m , m , m , m , n , c0 , a );
//        print_matrix(c0 , m);
        search_inverse_matrix( c0 , c1 , m , norm , perest_block );//inverse matrix in c0
//        print_matrix(c0 , m);
//        print_matrix(c1 , m);
        for( int j = arg->num ; j < k ; j += arg->p )
        {
            if(j>column)
            {
                get_block( i * n * m , j * m , m , m , m , n , c2 , a );
                block_mult_block( c1 , c0 , c2 , m ); // c1 = c0 * c2
                put_block( i * n * m, j * m , m , m , m , n , c1 , a );
            }
        }

//        if(p!=0&&arg->num==(i%arg->p))
//        {
//            get_block( i * n * m, k * m , p , m , m , n , c2 , a );
//            block_mult_block( c1 , c0 , c2 , m );
//            put_block( i * n * m, k * m , p , m , m , n , c1 , a );
//        }
//        print_matrix(a , n );
//        print_vector(b,n);
        if(arg->num==(i%arg->p))
        {
            if(p!=0)
            {
                get_block( i * n * m, k * m , p , m , m , n , c2 , a );
                block_mult_block( c1 , c0 , c2 , m );
                put_block( i * n * m, k * m , p , m , m , n , c1 , a );
            }
            get_block(i*m,0,1,m,m,1, c2, b);
    //        print_matrix(c2 , m );
            block_mult_block( c1 , c0 , c2 , m );
    //        print_matrix( c1 , m );
            put_block(i*m,0,1,m,m,1, c1, b);
        }
//        get_block(i*m,0,1,m,m,1, c2, b);
//        /*print_matrix(c2 , m );*/
//        block_mult_block( c1 , c0 , c2 , m );
//        /*print_matrix( c1 , m );*/
//        put_block(i*m,0,1,m,m,1, c1, b);
//        print_matrix(a , n );
//        print_vector(b , n );
//        pthread_mutex_unlock(arg->mutex);

//        make_E( m , c1 );/*uberi nahuj kogda vse budet ok*/
//        put_block( i*n*m , column*m , m , m , m , n , c1 , a );
//    /*}*/
}

void substract_lines( int i , int column , double * a , double * b , double * c0 ,
                      double * c1 , double * c2 , int n , int m , int k , int p,
                      Args * arg , double * for_count , int size )//n == k*m + p
// a[t][j] = a[t][j] - a[t][i]*a[i][j]
{
//    pthread_mutex_lock(arg->mutex);
//    print_matrix(a,n);
//    LOG(arg->num);
//    LOG(k);
//    LOG(arg->p);
//    pthread_mutex_unlock(arg->mutex);

//    /*for( int t = 0 ; t < k ; t ++ )
//    {
//        get_block( i*n*m , t*m , m , m , m ,   n  , c2 ,     a     );
//        put_block(   0   , t*m , m , m , m , size , c2 , for_count ); /* may be a problems because matrix for_count is not square */
//    }*/
//pthread_barrier_wait(arg->barrier);
//if(arg->num==0)
//{
//    for ( int count = 0 ; count < m ; count++ )
//    {
//        for( int j =0 ; j < size ; j++)
//        {
//            cout<<for_count[count*size+j]<<" ";
//        }
//        cout<<"\n";
//    }
//    cout<<"\n\n";
//}
//    /*get_block( i*n*m , k*m , p , m , m , n , c2 , a );
//    put_block( 0 , k*m , p , m , m , size , c2 , for_count );*/

    for(int t = 0 ; t < m ; t ++ )
    {
        memcpy( for_count+size*t , a+i*n*m+t*n , n*sizeof(double) );
    }

    get_block( i*m , 0 , 1 , m , m , 1 , c2 , b );
    put_block( 0 , k*m+p , 1 , m , m , size , c2 , for_count );

//pthread_barrier_wait(arg->barrier);
//if(arg->num==0)
//{
//    for ( int count = 0 ; count < m ; count++ )
//    {
//        for( int j =0 ; j < size ; j++)
//        {
//            cout<<for_count[count*size+j]<<" ";
//        }
//        cout<<"\n";
//    }
//    cout<<"\n\n";
//}
//pthread_barrier_wait(arg->barrier);
    for( int t = arg->num ; t < k ; t += arg->p )
    {
        if(t!=i)
        {

            for( int j = 0/*arg->num*/ ; j < k ; j ++/*= arg->p*/ ) /*to do better in future, need to begin ~ j==column+1*/
            {
                if( j>=column+1)
                {
                    get_block( t*n*m , column*m , m , m , m , n , c2 , a );

                    //for a
                    get_block( 0 , j*m , m , m , m , size , c1 , for_count );

                    block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1

                    get_block( t*n*m , j*m , m , m , m , n , c1 , a );
                    substract_block( c2 , c1 , c0 , m , m , m );//c2=c1-c0

                    put_block( t*n*m , j*m , m , m , m , n , c2 , a );
                }
            }
//            if( arg->num == 0 )
//            {
//                LOG(t);
                //for b
                get_block( t*n*m , column*m , m , m , m , n , c2 , a );
                get_block( 0 , k*m+p , 1 , m , m , size , c1 , for_count );
//                LOG("c2");
//                print_matrix(c2,m);
//                LOG("c1");
//                print_matrix(c1,m);
                block_mult_block( c0 , c2 , c1 , m );
//                LOG("c0 == c2*c1");
//                print_matrix(c0,m);
                get_block(t*m,0,1,m,m,1, c1, b);
//                LOG("c1");
//                print_matrix(c1,m);
                substract_block( c2 , c1 , c0 , m , 1 , m);
//                LOG(" c2 == c1 - c0 ");
//                print_matrix(c2,m);
                put_block(t*m,0,1,m,m,1, c2, b);
//            }

            if( p != 0 && (arg->num == (t%arg->p)) ) // for right side a
            {
//        pthread_mutex_lock(arg->mutex);
//                LOG(t);
//        pthread_mutex_unlock(arg->mutex);

                get_block( t*n*m , column*m , m , m , m , n , c2 , a );//p x m

                //for a

                get_block( 0 , k*m , p , m , m , size , c1 , for_count );//p x p
//                LOG("c2");
//                print_matrix(c2,m);
//                LOG("c1");
//                print_matrix(c1,m);

                block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1
//                LOG("c0 == c2*c1");
//                print_matrix(c0,m);
                get_block( t*n*m , k*m , p , m , m , n , c1 , a );//m x p
//                LOG("c1");
//                print_matrix(c1,m);
                substract_block( c2 , c1 , c0 , m , p , m);//c2=c1-c0
//                LOG(" c2 == c1 - c0 ");
//                print_matrix(c2,m);
                put_block( t*n*m , k*m , p , m , m , n , c2 , a );

                //for b already make
            }
        }
    }

    if( p != 0 && i != k && ( arg->num == (k%arg->p) )) // for down-right side block
    {
//        pthread_mutex_lock(arg->mutex);
//        LOG(" for down-right side block in substract_lines");
//        LOG(i);
//        LOG(arg->num);
//        pthread_mutex_unlock(arg->mutex);

        get_block( k*n*m , column*m , m , p , m , n , c2 , a );//p x m

        //for a

        get_block( 0 , k*m , p , m , m , size , c1 , for_count );//p x p

        block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1

        get_block( k*n*m , k*m , p , p , m , n , c1 , a );//m x p

        substract_block( c2 , c1 , c0 , p , p , m);//c2=c1-c0
        put_block( k*n*m , k*m , p , p , m , n , c2 , a );

        //for b
        get_block( k*n*m , column*m , m , p , m , n , c2 , a );

        get_block( 0 , k*m+p , 1 , m , m , size , c1 , for_count );//get_block( i*m , 0 , 1 , m , m , 1 , c1 , b );


        block_mult_block( c0 , c2 , c1 , m );

        get_block( k*m , 0 , 1 , p , m , 1 , c1 , b );
        substract_block( c2 , c1 , c0 , p , 1 , m);

        put_block( k*m , 0 , 1 , p , m , 1 , c2 , b );

    }

    if(p != 0 && ( arg->num == (k%arg->p) )) //for bottom except left side column
    {
        for(int j = column+1 ; j < k ; j ++ )
        {
            get_block( k*n*m , column*m , m , p , m , n , c2 , a);// p x m

            //for a
            get_block( 0 , j*m , m , m , m , size , c1 , for_count );

            block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1

            get_block( k*n*m , j*m , m , p , m , n , c1 , a );
            substract_block( c2 , c1 , c0 , p , m , m);//c1=c1-c0 p x m

            put_block( k*n*m , j*m , m , p , m , n , c2 , a );
        }
    }

//    for( int t = arg->num ; t < k ; t += arg->p )//for column left side //ubrat' eto pozje
//    {
//        if(t!=i)
//        {
//                //for a
//                get_block( t*n*m , column*m , m , m , m , n , c2 , a );
//                make_0( 0 , 0 , m , c2 );
//                put_block( t*n*m , column*m , m , m , m , n , c2 , a );
//        }
//    }
//
//    if(k!=i&&(arg->num==(k%arg->p)))//for left column from down //ubrat' eto pozje
//    {
//        get_block( k*n*m , column*m , m , p , m , n , c2 , a );
//        make_0( 0 , 0 , m , c2 );
//        put_block( k*n*m , column*m , m , p , m , n , c2 , a );
//    }
//    pthread_mutex_lock(arg->mutex);
//    print_matrix(a,n);
//    print_vector(b,n);
//    pthread_mutex_unlock(arg->mutex);
}

void substract_block( double * c0 , double * c1 , double * c2 , int hight , int width , int real_size_c)//c0=c1-c2 hight*width
{
    for(int i = 0 ; i < hight ; i ++ )
    {
        for( int j = 0 ; j < width ; j ++ )
        {
            c0[ i*real_size_c+j ] = c1[ i*real_size_c+j ] - c2[ i*real_size_c+j ];
        }
    }
}

void block_mult_block(double *c, double *a, double *b, int n) // n%3==0  c=a*b
{
	int i,j,k;
	double *pc = c , *pa = a , *pb = b , sum[9] ;

	for(i=0;i<n;i++)
    {
		for(j=0;j<n;j++)
        {
            c[i*n+j]=0.;
        }
    }

	for(i=0;i<n;i+=3)
    {
		for(j=0;j<n;j+=3)
        {
			sum[0]=sum[1]=sum[2]=sum[3]=sum[4]=sum[5]=sum[6]=sum[7]=sum[8]=0.;

			for(k=0;k<n;k++)
			{
				pa=a+i*n+k;
				pb=b+k*n+j;

				sum[0]+=pa[0]*pb[0];
				sum[1]+=pa[0]*pb[1];
				sum[2]+=pa[0]*pb[2];
				sum[3]+=pa[n]*pb[0];
				sum[4]+=pa[n]*pb[1];
				sum[5]+=pa[n]*pb[2];
				sum[6]+=pa[2*n]*pb[0];
				sum[7]+=pa[2*n]*pb[1];
				sum[8]+=pa[2*n]*pb[2];
			}

			pc=c+i*n+j;

			pc[0]			+=sum[0];
			pc[1]			+=sum[1];
			pc[2]			+=sum[2];
			pc[n]			+=sum[3];
			pc[n+1]			+=sum[4];
			pc[n+2]			+=sum[5];
			pc[2*n]			+=sum[6];
			pc[2*n+1]		+=sum[7];
			pc[2*n+2]		+=sum[8];
		}
    }
}

void recovery_solution( double * b , double * x , int k , int m , int p , Args * arg , int * index)
{
//    pthread_barrier_wait(arg->barrier);
    for( int i = arg->num ; i < k ; i += arg->p )
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            x[i*m+j]=b[index[i]*m+j];
        }
    }
    if( p!= 0 && arg->num == arg->p-1 )
    {
        for( int i = 0 ; i < p ; i ++ )
        {
            x[k*m+i]=b[k*m+i];
        }
    }
//    pthread_barrier_wait(arg->barrier);
}

void print_matrix( double * a , int n )
{
    for ( int i = 0 ; i < n ; i++ )
    {
        for( int j =0 ; j < n ; j++)
        {
            cout<<a[i*n+j]<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n\n";
}
template< typename T >
void print_vector( T *b, int n)
{
    for ( int i = 0 ; i < n ; i++ )
    {
        cout<<b[i]<<" ";
    }
    cout<<"\n\n";
}

void init_rhs_paral(double *a, double *b, double *x, int n) //b=A*x
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

