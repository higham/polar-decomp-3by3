`polar-decomp-3by3` - MATLAB code for polar decomposition of a real 3-by-3 matrix
==========

About
-----

`polar-decomp-3by3` contains a MATLAB function `polar_quaternion` for
computing the polar decomposition of a real 3-by-3 matrix.
This decomposition has the form A = Q*H, where Q is orthogonal and H is
symmetric positive semidefinite.
The underlying algorithm is developed in the paper

N. J. Higham and V. Noferini.
"[An algorithm to compute the polar decomposition of a $3 \times 3$
matrix](http://link.springer.com/article/10.1007%2Fs11075-016-0098-7)". Numer. Algorithms, 73(2):349-369, 2016.

The function `test.m` tests that the code is working correctly.

Requirements
-------------

The codes have been developed under MATLAB 2016a
and have been tested with MATLAB 2016b.

License
-------

See `license.txt` for licensing information.

