#include <stdio.h>
#include <time.h>
#include <omp.h>


int main (void)
{
    clock_t start = clock();

    size_t N = 100000000;
    #ifdef _OPENMP
    double wt_start = omp_get_wtime();    
    double sum = 0.0;
    /* printf("Hello from thread %d/%d\n",omp_get_thread_num(),omp_get_num_threads()); */
    #pragma omp parallel for reduction(+:sum)
    for (size_t ii = 0; ii < N; ii++){
        sum = sum + (ii * 2)/ 3;
    }
    printf("sum = %G\n",sum);
    double wt_end = omp_get_wtime();
    double wt = wt_end - wt_start;
    printf("Wall time elapsed (omp): %G\n",wt);
    #endif

    clock_t start_1 = clock();
    sum = 0.0;
    for (size_t ii = 0; ii < N; ii++){
        sum = sum + (ii * 2)/ 3;
    }
    printf("sum = %G\n",sum);
    clock_t end_1 = clock();
    printf("Wall time elapsed (serial): %G\n",(double)(end_1-start_1)/CLOCKS_PER_SEC);
    /* gettimeofday(&end,NULL); */

    clock_t end = clock();
        
    /* double tout = ((end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec)) / 1.0e6; */
    double tout = (double)(end-start) / CLOCKS_PER_SEC;
    printf("Wall time elapsed: %G\n",tout);
    return 0;
}

