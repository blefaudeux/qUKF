qUKF
=====

What is it ?
----------
An Eigen-based implementation of an UKF filter, allowing the use of quaternions to filter angles. This removes singularities bound to appear if you feed your filter with ]-pi, pi] values, contrary to most of the UKF filters out there.

* Standard implementation from the Julier and Uhlman initial UKF proposition (Julier, S. J., & Uhlmann, J. K. (1997). A New Extension of the Kalman Filter to Nonlinear Systems, 182–193.), the use of quaternions being initially proposed by Kraft (Kraft, E. (2003). A quaternion-based unscented Kalman filter for orientation tracking. Proceedings of the Sixth International Conference …, 1, 47–54. Retrieved from http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf)


What does it depends on ?
-------------------------
Eigen headers. The small demo requires OpenCV though, but really dispensable

How to build ?
--------------
after cloning, it's as simple as (unix, for windows you could probably do something with the CMake GUI):

1. `mkdir build`

2. `cd build`

3. `cmake ..` (or `cmake .. -DDEMO=1` if you want to build the demo executable)

4. `make`

General observations
--------------------
Code quality is not top notch, leaves a lot to be desired. Feel free to contribute, just check that 
Travis CI is still happy from time to time..

License
-------

The MIT License (MIT)
Copyright (c) 2013 Benjamin Lefaudeux

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
