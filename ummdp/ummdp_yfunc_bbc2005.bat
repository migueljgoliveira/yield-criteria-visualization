ECHO OFF
cls

del ummdp_yfunc_bbc2005.pyd >nul

python -m numpy.f2py -c -m ummdp_yfunc_bbc2005 ummdp_yfunc_bbc2005.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_bbc2005.log 2>&1

rename ummdp_yfunc_bbc2005.*.pyd ummdp_yfunc_bbc2005.pyd >nul