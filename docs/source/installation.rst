*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{HODLR2lib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/sivaramambikasaran/HODLR2.git --recursive

The ``--recursive`` flag is argument to ensure that all submodules are also checked out.

Dependencies
-------------

- Eigen Linear Algebra Library (get it `here <https://bitbucket.org/eigen/eigen/>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.

**NOTE**: On MacOS, the default compiler is `clang` which doesn't have OpenMP support. You will have to use g++ to make use of the speedups from OpenMP::

    user@computer HODLR2$ brew install g++-8
    user@computer HODLR2$ export CXX=g++

Installation
-------------

Manually Installing
^^^^^^^^^^^^^^^^^^^

Then, set the environment variable ``EIGEN_PATH`` to the location of your Eigen installation. This is needed by the CMake script.::

    user@computer HODLR2$ export EIGEN_PATH=path/to/eigen/

Testing
-------

Now, we need to ensure that all the functions of the libraries function as intended. For this purpose, we will be running the script ``examples/testHODLR2D.cpp``.
Key in the file ```examples/testHODLR2D.cpp`` as input under ``INPUT_FILE`` in ``HODLR2lib/CMakeLists.txt``. Here you also set the name of the output executable, say ``test_HODLR2``, under ``OUTPUT_EXECUTABLE_NAME``.
Then create a directory called ``build`` and navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable.
To check this on your computer, run the following lines::

    user@computer HODLR2$ mkdir build && cd build
    user@computer build$ cmake ..
    user@computer build$ make
    user@computer build$ ./test_HODLR2

For a succesful test, the final line of output for this run would read:"Reached End of Test File Successfully! All functions work as intended!".

Building and Executing
----------------------

Key in the required ``.cpp`` to be used as input under ``INPUT_FILE`` in ``HODLR2lib/CMakeLists.txt``. Here you also set the name of the output executable under ``OUTPUT_EXECUTABLE_NAME``. Then navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable::

    user@computer build$ cmake path/to/HODLR2/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable
