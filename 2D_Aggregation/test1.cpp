#include <math.h>
#include <stdio.h>
 
int main(void)
{
    double x, y, z;
    x = 4.27;
 
    y = j0(x);       /* y = -0.3660 is the order 0 bessel */
                     /* function of the first kind for x  */
    z = yn(3,x);     /* z = -0.0875 is the order 3 bessel */
                     /* function of the second kind for x */
 
    printf("y = %lf\n", y);
    printf("z = %lf\n", z);
}