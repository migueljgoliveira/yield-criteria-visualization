ECHO OFF
cls

del ummdp_yfunc_yld2004_18p.pyd >nul

python -m numpy.f2py -c -m ummdp_yfunc_yld2004_18p ummdp_yfunc_yld2004_18p.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_yld2004_18p.log 2>&1

rename ummdp_yfunc_yld2004_18p.*.pyd ummdp_yfunc_yld2004_18p.pyd >nul