#!/bin/bash
gcc -shared -fPIC -o ./lib/libmthread.so mthread.c
clang -O3 -c -mcpu=hip11 ldfslp_include.c -o lib/ldfslp.o
ar -crv lib/libldfslp.a lib/ldfslp.o
