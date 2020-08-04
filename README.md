# MCMCv4
## Markov Chain Monte Carlo for random noncommutative geometries
Code for simulating random fuzzy spaces. This is the version used to generate the data for the papers: [1612.00713](https://arxiv.org/abs/1612.00713) and [1902.03590](https://arxiv.org/abs/1902.03590) which is a somewhat updated version of the code used in [1510.01377](https://arxiv.org/abs/1510.01377)

This code relies on two external libaries.
- The Eigen::MatrixXcd class is implemented in the Eigen Matrix library https://eigen.tuxfamily.org/ and my code uses version 3.2.4 of it.

- The code also requires a random number generator which I used from here class https://www.agner.org/random/


# Pauls Guide on how to get this to work.

First steps, what version of the compilers/packages do we need?
  - We need to install Eigen. The version **I think** this code was written using is  version 3.2.4. To install this visit the website: https://eigen.tuxfamily.org/ and following the instruction there. At the time of writing this (August 2020) you can get version 3.2.4 from https://gitlab.com/libeigen/eigen/-/releases/3.2.4  
    - You download this is a .zip file or a .tar file (whichever your familiar with). You can follow these instructions: https://eigen.tuxfamily.org/dox/GettingStarted.html
    - Copy the contents of the .zip/.tar file to either your project directory or to somewhere your compiler has access to. For now, copy it to your project directory (unless you already know how to symlink it to /usr/local/include).
    - You need to edit a path variable in the makefile to be wherever you've copied your eigen library to (either the project directory or /usr/local/include). In the makefile, we have the variable: LDFLAGS. This needs to be set to the location of the eigen package.
  - We need to install some random number generators. This code uses some programs from: https://www.agner.org/random/
    - I have copied the relevant files into Randomc.
    - Note: If you download the files for yourself from the web, you will need to edit the "rancombi.cpp" file to remove the example. Compare the file you download to the one in this repo.


# Testing the setup
To check that everything works correctly, you should be able to compile and run the files "test.cpp" using the following commands:
```bash
g++ test.cpp -I eigen-3.2.4 -o test
./test
```

If all is well you should see the following:
```bash
m =
10.0008  55.865 14.7045
23.1538 63.2767 77.8865
85.5605 31.8959 77.9296
m * v =
165.844
383.367
383.141
test=
(0,0) (0,0) (0,0)
(0,0) (0,0) (0,0)
(0,0) (0,0) (0,0)
adjoint=
(0,-0) (0,-0) (0,-0)
(0,-0) (0,-0) (0,-0)
(0,-0) (0,-0) (0,-0)
```
# Running the program
Now that we have successfully installed the necessary packages, we can compile the code.
We do this easily by using the Makefile. Just type ```make``` in your terminal.

If all goes well, you should have the following output:
```bash
g++ -c -IEigen/Eigen3.24 -g3 -std=c++11 -Wno-ignored-attributes -Wno-deprecated-declarations -IEigen/Eigen3.24
 -I eigen-3.2.4 MCMCv4.cpp
g++ -c -IEigen/Eigen3.24 -g3 -std=c++11 -Wno-ignored-attributes -Wno-deprecated-declarations -IEigen/Eigen3.24
 progParams.cpp
g++ -c -IEigen/Eigen3.24 -g3 -std=c++11 -Wno-ignored-attributes -Wno-deprecated-declarations -IEigen/Eigen3.24
 -I eigen-3.2.4 Dirac.cpp
g++ MCMCv4.o progParams.o Dirac.o -o MCMCv4
```

Once we have compiled the code, we should check everything is working correctly.


## Input File Settings
To run the program after compilation, you need to choose what type of geometry you want to investigate. This is done by setting the various parameters in an input file.

An example of such an input file can be found at ```example_input.txt```. To set the variables you need to write the variable's name, followed by a space, followed by the input for that variable. Specifically, no equals sign or anything.

Here is a brief summary of the variables:
- ```matrixsize``` - This one is fairly obvious, it defines what size random matrices we are considering. This requires a positive integer, i.e. 1, 2, 3, 4,...
- ```Type``` - This variable dictates which Clifford Type we are considering for our random geometries. Note that not all types have been implemented yet. The only types available are: (p,q) =? (0,2). This variable requires the form ```pq``` for a type (p,q) Clifford type, i.e. for type (0,2), you need to set the ```Type``` variable to ```02```.
- ```steps``` - How long to run the sims
- ```outputfile``` - File name where is the data saved
- ```initialconfig``` -  There are options to start from the operator D being an identity matrix etc. If this variable is set to 5 then this starts with a random Dirac, check in Dirac.cpp to see what the others are
- ```measurement```  if this is set to 2 then it measures the eigenvalues of D, see progParams.cpp for explanation for other options
- ```couplingD4``` set it to non-zero if you want quartic terms in the action.
- ```couplingD2``` set to non-zero if you want quadratic terms in the action.
- ```couplingD22``` Set to zero so we have a simple model with only quadratic data (Paul has no idea what this means ??)
- ```finalmatrixfile``` - the code spits out the last matrix so that we can use it as an initial matrix in the future. This variable is the filename of where to save that, needs to be txt
