// dummy file that gets compiled to a shared object but it not a plugin

#include <cstdio>

extern "C" {

void hello()
{
    puts("Hello, World!");
}
}
