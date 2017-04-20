/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
    int *da = (int *)a;
    int *db = (int *)b;
    
    if (*da > *db)
        return 1;
    else if (*da < *db)
        return -1;
    else
        return 0;
}

int main( int argc, char *argv[])
{
    int rank;
    int i, N,c=0,p,l,j=0;
    int *vec,*lsplit,*splitter,*split;
    MPI_Request req[1];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Number of random numbers per processor (this should be increased
     * for actual tests or could be passed in through the command line */
    N = 20;
    
    vec = calloc(N, sizeof(int));
    /* seed random number generator differently on every core */
    srand((unsigned int) (rank + 393919));
    
    /* fill vector with random integers */
    for (i = 0; i < N; ++i) {
        vec[i] = rand();
    }

    printf("rank: %d, first entry: %d\n", rank, vec[0]);
    
    /* sort locally */
    qsort(vec, N, sizeof(int), compare);

    lsplit=(int *)malloc((p)*sizeof(int));
    for(i=0;i<N;i=i+(N/p))
    {
        lsplit[c]=vec[i];
        c=c+1;
    }
    
    /* randomly sample s entries from vector or select local splitters,
     * i.e., every N/P-th entry of the sorted vector */
    
    if(0==rank)
    {
        split=(int *)malloc(c*p*sizeof(int));
    }

    MPI_Gather(lsplit,c,MPI_INT,split,c,MPI_INT,0,MPI_COMM_WORLD);
    

    /* every processor communicates the selected entries
     * to the root processor; use for instance an MPI_Gather */
    splitter=(int *) calloc(p-1,sizeof(int));

    if(0==rank)
    {
        l=sizeof(split)/sizeof(int);
        qsort(split, l, sizeof(int), compare);
        
        
        for(i=0;i<(p-1);i++)
        {
            splitter[i]=split[i*c];
        }
        
        
    }

    /* root processor does a sort, determinates splitters that
     * split the data into P buckets of approximately the same size */
    
    MPI_Bcast (splitter, p-1, MPI_INT, 0, MPI_COMM_WORLD);

    /* root process broadcasts splitters */

    int *Bucket=(int *)calloc(N,sizeof(int));

    
    int *tracker=(int *) calloc (p,sizeof (int) * (p));
    int *sendcou=(int *) calloc (p,sizeof (int) * (p));
    int *reccou=(int *) calloc (p,sizeof (int) * (p));
    int *recdis=(int *) calloc (p,sizeof (int) * (p));
    int cou=0,k=0;
    for(i=0;i<N;++i)
    {
        if(j<(p-1))
        {
            if(splitter[j]>=vec[i])
            {
                Bucket[i]=vec[i];
                cou+=1;
            }
            else
            {
                j++;
                --i;
                tracker[j]=k;
                sendcou[j-1]=cou;
                cou=0;
            }
        }
        else
        {
            Bucket[i]=vec[i];
            cou+=1;
        }
        k++;
        
    }
    sendcou[p-1]=cou;
    //MPI_Allreduce(sendcou,reccou,p,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    int *temp=(int *) calloc (p,sizeof (int) * (p));
    MPI_Allgather(sendcou,p,MPI_INT,temp,p,MPI_INT,MPI_COMM_WORLD);
    /*
    for(i=0;i<p;i++)
    {
        MPI_Send(&sendcou[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Recv(&temp[i],1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
     */
    for(i=0;i<p;i++){
        recdis[i]=recdis[i-1]+temp[i-1];
        printf("%d\n",sendcou[i]);

    }

    int *recbuf=(int *) malloc (sizeof (int) * (reccou[rank]));

    /* every processor uses the obtained splitters to decide
     * which integers need to be sent to which other processor (local bins) */
    for(i=0;i<p;i++)
    {
        int *sendtemp=(int *) calloc(sendcou[i],sizeof(int));
        int *rectemp=(int *) calloc(temp[i],sizeof(int));

        for(j=0;j<sendcou[rank];j++)
            sendtemp[j]=Bucket[tracker[i]+j];
        MPI_Send(sendtemp,sendcou[i],MPI_INT,i,0,MPI_COMM_WORLD);
        MPI_Recv(rectemp,temp[i],MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(j=0;j<temp[i];j++){
            int d=recdis[i];
            recbuf[d+j]=rectemp[j];
        }
        free(sendtemp);
        free(rectemp);
    }
    //MPI_Alltoallv(Bucket,sendcou,tracker,MPI_INT,recbuf,reccou,recdis,MPI_INT,MPI_COMM_WORLD);
    /* send and receive: either you use MPI_AlltoallV, or
     * (and that might be easier), use an MPI_Alltoall to share
     * with every processor how many integers it should expect,
     * and then use MPI_Send and MPI_Recv to exchange the data */
    /*
    for(i=0;i<2*N;i++)
    {
        if(&recbuf[i]==NULL)
            recbuf[i]=-1;
    }
     */
    qsort(recbuf, reccou[rank], sizeof(int), compare);
    /* do a local sort */

    /* every processor writes its result to a file */
    free(vec);
    free(tracker);
    free(lsplit);
    free(splitter);
    if(rank==0)
        free(split);
    free(Bucket);
    free(sendcou);
    free(reccou);
    free(recbuf);
    free(temp);
    
    MPI_Finalize();
    return 0;
}
