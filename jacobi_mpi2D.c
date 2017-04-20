#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq)//lN is no of rows and columns
{
    int i;
    double  gres = 0.0, lres = 0.0;
    
    for (i = 0; i <(lN + 2)*(lN + 2); ++i) {
        int r=i/(lN+2);
        int c=i-r*(lN+2);
        int N=(lN+1);
        if(r==0||r==(lN+1)||c==0||c==(lN+1))
        {}
        else
            lres=lres+pow(invhsq*(4*lu[i]-(lu[(r-1)*(N+1)+c]+lu[(r)*(N+1)+c-1]+lu[(r+1)*(N+1)+c]+lu[(r)*(N+1)+c+1]))-1,2);
        
    }
 
    /* use allreduce for convenience; a reduce would also be sufficient */
    MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(gres);
}


int main(int argc, char * argv[])
{
    int mpirank, i, p, N, lN, iter, max_iters;//N is total no of rows and colums
    MPI_Status status, status1;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* get name of host running MPI process */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);
    
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &max_iters);
    
    /* compute number of unknowns handled by each process */
    lN = N / p;
    int flag=0;
    int j=p;
    while(j%4==0) {
        j=j/4;
        if(j==1)
            flag=1;
    }
    if(flag==0&&mpirank==0)
    {
        printf("Exiting. p must of the form 4^j");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if (((N*N) % p != 0) && mpirank == 0 ) {
        printf("N: %d, local N: %d\n", N, lN);
        printf("Exiting. N must be a multiple of p\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    /* Allocation of vectors, including left/upper and right/lower ghost points */
    double * lu    = (double *) calloc(sizeof(double), (lN + 2)*(lN + 2));
    double * lunew = (double *) calloc(sizeof(double), (lN + 2)*(lN + 2));
    double * lutemp;
    
    double h = 1.0 / (N + 1);
    double hsq = h * h;
    double invhsq = 1./hsq;
    double gres, gres0, tol = 1e-5;
    
    /* initial residual */
    gres0 = compute_residual(lu, lN, invhsq);
    gres = gres0;
    
    for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) {
        
        /* Jacobi step for local points */
        for (i = 0; i < (lN + 2)*(lN + 2); ++i) {
            int r=i/(lN+2);
            int c=i-r*(lN+2);
            if(r==0||r==(lN+1)||c==0||c==(lN+1))
            {}
            else{
                lunew[i]=.25*(hsq+lu[(r-1)*(lN+2)+c]+lu[(r)*(lN+2)+c-1]+lu[(r+1)*(lN+2)+c]+lu[(r)*(lN+2)+c+1]);
                //printf("%d\n",lunew[i]);
            }
            
            
        }
        
        /* communicate ghost values */ //change for 2D
        int np=sqrt(p);
        int mr=mpirank/(np);
        int mc=mpirank-mr*(np);
        
        if (mr <(np-1)) {
            int srrank=(mr+1)*np+mc;
            /* If not a top row process, send/recv bdry values to the top */
            for(i=0;i<(lN+2);i++)
            {
            MPI_Send(&(lunew[lN*(lN+2)+i]), 1, MPI_DOUBLE, srrank, 124, MPI_COMM_WORLD);
            MPI_Recv(&(lunew[(lN+1)*(lN+2)+i]), 1, MPI_DOUBLE, srrank, 123, MPI_COMM_WORLD, &status);
            }
        }
        if (mr > 0) {
            int srrank=(mr-1)*np+mc;
            /* If not a bottom row process, send/recv bdry values to the bottom */
            for(i=0;i<(lN+2);i++)
            {
                MPI_Send(&(lunew[(lN+2)+i]), 1, MPI_DOUBLE, srrank, 123, MPI_COMM_WORLD);
                MPI_Recv(&(lunew[i]), 1, MPI_DOUBLE, srrank, 124, MPI_COMM_WORLD, &status);
                        
            }
        }
       if (mc <(np-1)) {
          int srrank=mr*np+mc+1;
          /* If not a right-most column process, send/recv bdry values to the right */
          for(i=0;i<(lN+2);i++)
          {
            MPI_Send(&(lunew[i*(lN+2)+lN]), 1, MPI_DOUBLE, srrank, 128, MPI_COMM_WORLD);
            MPI_Recv(&(lunew[i*(lN+2)+lN+1]), 1, MPI_DOUBLE, srrank, 129, MPI_COMM_WORLD, &status);
                                          
           }
        }
       if (mc >0) {
         int srrank=mr*np+mc-1;
         /* If not a left-most column process, send/recv bdry values to the left */
         for(i=0;i<(lN+2);i++)
         {
            MPI_Send(&(lunew[i*(lN+2)+1]), 1, MPI_DOUBLE, srrank, 129, MPI_COMM_WORLD);
            MPI_Recv(&(lunew[i*(lN+2)]), 1, MPI_DOUBLE, srrank, 128, MPI_COMM_WORLD, &status);
                                      
         }
       }
                         
                         
        /* copy newu to u using pointer flipping */
        lutemp = lu; lu = lunew; lunew = lutemp;
        if (0 == (iter % 10)) {
            gres = compute_residual(lu, lN, invhsq);
            if (0 == mpirank) {
                printf("Iter %d: Residual: %g\n", iter, gres);
            }
        }
    }
    
    /* Clean up */
    free(lu);
    free(lunew);
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    if (0 == mpirank) {
        printf("Time elapsed is %f seconds.\n", elapsed);
    }
    MPI_Finalize();
    return 0;
}
