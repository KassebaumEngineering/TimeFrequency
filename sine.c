
#include <stdio.h>
#include <math.h>

main()
{
    int i,j;
    double value;

    for(i=0; i< 4096; i++ ){
        value = sin((double)(4 * M_PI  *  i * i / 40000.0) );
        value += sin( 2*M_PI*i/10.0 );
        printf("%f\n", (float) value );
    }
}
