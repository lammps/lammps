/* file for testing loading of shared objects and libraries */

int some_int_val = 12345;
double some_double_val = 6.78e-9;

int some_int_function(int arg)
{
    return arg*arg;
}

double some_double_function(double arg1, int arg2)
{
    double sum = 0;
    for (int i = 0; i < arg2; ++i)
        sum += arg1;
    return sum;
}
