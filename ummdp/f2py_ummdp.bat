ECHO OFF
cls

del *.pyd >nul

python -m numpy.f2py -c -m ummdp_yfunc_mises ummdp_yfunc_mises.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_mises.log 2>&1
python -m numpy.f2py -c -m ummdp_yfunc_hill48 ummdp_yfunc_hill48.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_hill48.log 2>&1
python -m numpy.f2py -c -m ummdp_yfunc_yld2004_18p ummdp_yfunc_yld2004_18p.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_yld2004_18p.log 2>&1
python -m numpy.f2py -c -m ummdp_yfunc_cpb2006 ummdp_yfunc_cpb2006.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_cpb2006.log 2>&1
python -m numpy.f2py -c -m ummdp_yfunc_yld2000_2d ummdp_yfunc_yld2000_2d.f --fcompiler=intelvem --opt=/heap-arrays:0 1> ummdp_yfunc_yld2000_2d.log 2>&1

rename ummdp_yfunc_mises.*.pyd ummdp_yfunc_mises.pyd >nul
rename ummdp_yfunc_hill48.*.pyd ummdp_yfunc_hill48.pyd >nul
rename ummdp_yfunc_yld2004_18p.*.pyd ummdp_yfunc_yld2004_18p.pyd >nul
rename ummdp_yfunc_cpb2006.*.pyd ummdp_yfunc_cpb2006.pyd >nul
rename ummdp_yfunc_yld2000_2d.*.pyd ummdp_yfunc_yld2000_2d.pyd >nul
