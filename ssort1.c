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
    int i, N,p,l,j=0;
    int *vec,*lsplit,*splitter,*split;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Request request[p];

    /* Number of random numbers per processor (this should be increased
     * for actual tests or could be passed in through the command line */
    N = 100;
    
    vec = calloc(N, sizeof(int));
    /* seed random number generator differently on every core */
    srand((unsigned int) (rank + 393919));
    
    /* fill vector with random integers */
    for (i = 0; i < N; ++i) {
        vec[i] = rand();
        //printf("%d\n",vec[i]);
    }
    
    // printf("rank: %d, first entry: %d\n", rank, vec[0]);
    
    /* sort locally */
    qsort(vec, N, sizeof(int), compare);
    
    lsplit=(int *)malloc((p-1)*sizeof(int));
    for(i=0;i<(p-1);i++)
    {
        lsplit[i]=vec[(N/p)*(i+1)];
    }
    
    /* randomly sample s entries from vector or select local splitters,
     * i.e., every N/P-th entry of the sorted vector */
    
    if(0==rank)
    {
        split=(int *)malloc(p*(p-1)*sizeof(int));
    }
    
    MPI_Gather(lsplit,p-1,MPI_INT,split,p-1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){
        for(i=0;i<N;i++)
        {
        }
    }
    
    /* every processor communicates the selected entries
     * to the root processor; use for instance an MPI_Gather */
    splitter=(int *) calloc(p-1,sizeof(int));
    
    if(0==rank)
    {
        l=sizeof(split)/sizeof(int);
        qsort(split, l, sizeof(int), compare);
        
        
        for(i=0;i<(p-1);i++)
        {
            splitter[i]=split[p*(i+1)-1];
            // printf("%d\n",splitter[i]);
        }
        
    }
    
    /* root processor does a sort, determinates splitters that
     * split the data into P buckets of approximately the same size */
    
    MPI_Bcast (splitter, p-1, MPI_INT, 0, MPI_COMM_WORLD);
    /* root process broadcasts splitters */
    
    
    int *Bucket=(int *)calloc(N,sizeof(int));
    
    
    int *sendcou=(int *) calloc (p,sizeof (int) * (p));
    int *recdis=(int *) calloc (p,sizeof (int) * (p));
    int *reccou=(int *) calloc (p,sizeof (int) * (p));
    int *tracker=(int *) calloc (p,sizeof (int) * (p));

    int cou=0;
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
                sendcou[j-1]=cou;
                cou=0;
            }
        }
        else
        {
            Bucket[i]=vec[i];
            cou+=1;
        }
        
    }
    sendcou[p-1]=cou;
    for(j=1;j<p;j++){
        tracker[j]=tracker[j-1]+sendcou[j-1];
        //printf("%d\n",tracker[j]);
    }
    MPI_Allreduce(sendcou,reccou,p,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int *temp=(int *) calloc (p,sizeof (int) * (p));
    for(i=0;i<p;i++){
        MPI_Isend(&sendcou[i],1,MPI_INT,i,0,MPI_COMM_WORLD,&request[i]);
        MPI_Recv(&temp[i],1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
   // MPI_Allgather(sendcou,p,MPI_INT,temp,p,MPI_INT,MPI_COMM_WORLD);
    /*
   for(i=0;i<p;i++)
    {
        printf("rank:%d sending %d\n",rank,sendcou[i]);
        printf("rank:%d receiving %d\n",rank,temp[i]);

    }
     */
    for(i=1;i<p;i++)
        recdis[i]=recdis[i-1]+temp[i-1];
    int *recbuf=(int *) malloc (sizeof (int) * (reccou[rank]));
    for(i=0;i<p;i++)
    {
        int *sendtemp=(int *) calloc(sendcou[i],sizeof(int));
        int *rectemp=(int *) calloc(temp[i],sizeof(int));
        
        for(j=0;j<sendcou[i];j++)
        {
            sendtemp[j]=Bucket[tracker[i]+j];

        }
        MPI_Isend(sendtemp,sendcou[i],MPI_INT,i,0,MPI_COMM_WORLD,&request[i]);
        MPI_Recv(rectemp,temp[i],MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(j=0;j<temp[i];j++){
            int d=recdis[i];
            recbuf[d+j]=rectemp[j];
        }
        free(sendtemp);
        free(rectemp);
    }
    qsort(recbuf, reccou[rank], sizeof(int), compare);
    //for(i=0;i<reccou[rank];i++)
        //printf("rank:%d %d\n",rank, recbuf[i]);
    char str[100];
    sprintf(str, "ssp%d.txt", rank);
    FILE *f = fopen(str, "w");
    if(f== NULL)
    {
        printf("\nCan't open file or file doesn't exist.");
        exit(0);
    }
    fwrite(recbuf, sizeof(int),reccou[rank],f);
    fclose(f);
    free(tracker);
    free(recbuf);
    free(temp);
    free(recdis);
    free(vec);
    free(lsplit);
    free(splitter);
    if(rank==0)
        free(split);
    free(Bucket);
    free(sendcou);
    MPI_Finalize();
    return 0;
}
