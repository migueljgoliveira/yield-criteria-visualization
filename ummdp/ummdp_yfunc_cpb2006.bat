ECHO OFF
cls

del ummdp_yfunc_cpb2006.pyd >nul

python -m numpy.f2py -c -m ummdp_yfunc_cpb2006 ummdp_yfunc_cpb2006.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_cpb2006.log 2>&1

rename ummdp_yfunc_cpb2006.*.pyd ummdp_yfunc_cpb2006.pyd >nul