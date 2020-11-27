#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

struct {
    pthread_mutex_t data_lock[5];
    int             int_val[5];
    float           float_val[5];
} Data;

void *Add_to_Value(void *arg);

int main()
{
    pthread_t pid[5];
    int values[5];

    /* initialize init struct */
    for (int i=0; i<5; i++)
    {
        pthread_mutex_init(&Data.data_lock[i], NULL);
        Data.int_val[i] = 0;
        Data.float_val[i] = 0;

        values[i] = 2*i+1;
    }

    /* set concurrency and create the thread */
    for (int i=0; i<5; i++)
    {
        pthread_create(&pid[i], NULL, Add_to_Value, &values[i]);
    }

    /* join all thread */
    for (int i=0; i<5; i++){
        pthread_join(pid[i], NULL);
    }

    for (int i=0; i<5; i++)
    {
        printf("%d\n", Data.int_val[i]);
    }

    printf("done\n");
    return 1;

}

void *Add_to_Value(void *arg)
{
    int *inval = (int*) arg;
    
    for (int i=0; i<10000; i++)
    {
        pthread_mutex_lock(&Data.data_lock[i%5]);
        Data.int_val[i%5] += *inval;
        Data.float_val[i%5] += (float) 1.5*(*inval);
        pthread_mutex_unlock(&Data.data_lock[i%5]);
    }

    return ((void*) 0);
}