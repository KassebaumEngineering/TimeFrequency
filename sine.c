
#include <stdio.h>
#include <math.h>

main()
{
    int i,j;

    for(i=0; i< 4096; i++ ){
            printf("%f\n", (float) sin((double)(M_PI  *  i * i / 40000.0) ));
    }
}
