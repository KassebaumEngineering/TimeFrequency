
#include <stdio.h>
#include <math.h>

main()
{
    int i,j;

    for(i=0; i< 32; i++ ){
        for( j=0; j< 128; j++){
            printf("%f\n", (float) sin((double)(i*128 + j) * (M_PI * i / 50.0) + M_PI/2.0) );
        }
    }
}
