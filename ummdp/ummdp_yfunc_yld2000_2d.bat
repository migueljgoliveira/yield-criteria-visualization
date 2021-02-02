ECHO OFF
cls

del ummdp_yfunc_yld2000_2d.pyd >nul

python -m numpy.f2py -c -m ummdp_yfunc_yld2000_2d ummdp_yfunc_yld2000_2d.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_yld2000_2d.log 2>&1

rename ummdp_yfunc_yld2000_2d.*.pyd ummdp_yfunc_yld2000_2d.pyd >nul