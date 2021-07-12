# Tomofast-x  v.1.1

Geophysical 3D potential field joint and constrained parallel inversion code.  
Vitaliy Ogarko, Jeremie Giraud, Roland Martin.

Tomofast-x is 3D parallel inversion platform to run single domain or joint inversion (gravity and magnetic data).
It can use local weighting of gradient regularization function, global and local petrophysical constraints (Gaussian mixture model and multiple disjoint bound constraints).
Tomofast-x can run in parallel on laptop and supercomputers, using distributed memory systems.

The source code is a companion to the publication detailing Tomofast-x geophysical aspects and examples of utilisation and realistic dataset:
J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell, 2021: 
Structural, petrophysical and geological constraints in potential field inversion using the Tomofast-x open-source code, 
Geoscientific Model Development Discussions, https://doi.org/10.5194/gmd-2021-14


If you are using the code please cite the following papers:

- V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021), 
"Disjoint interval bound constraints using the alternating direction method of multipliers for geologically constrained inversion: Application to gravity data," GEOPHYSICS 86: G1-G11.
https://doi.org/10.1190/geo2019-0633.1

- J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021), 
"Structural, petrophysical and geological constraints in potential field inversion using the Tomofast-x open-source code",
Geoscientific Model Development Discussions, https://doi.org/10.5194/gmd-2021-14


### Compiling and running

To compile the code you need: gcc v.4.9 or more recent, and the MPI library (such as OpenMPI).

The makefile is contained in the root folder and should be used to compile Tomofast-x. Compiling the code is a necessary step to be able to run inversions.  
To compile the code run the make comand in the code directory as:  
```shell
make
```

To run the code with your parameter file:
```shell
./tomofast3D -j <Parfile path>
```

For the parallel run on you machine execute:
```shell
mpirun -np <number CPUs> -j <Parfile path>
```

To run unit tests (serial and parallel):
```shell
./runtests.sh
mpirun -np 3 ./runtests.sh
```

### Running the examples

Input data is contained in folder ``data``.  
The input parameter file (Parfile) which contains all the paramters of the inversion, is stored in folder: ``parfiles``.  
Output data is stored in folder ``output``. A detailed path is specified in the Parfile parameter ``global.outputFolderPath``.


To run the test example:
```shell
./tomofast3D -j ./parfiles/Parfile_mansf_slice.txt
```

### Publications using the code

- J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021), 
"Structural, petrophysical and geological constraints in potential field inversion using the Tomofast-x open-source code",
Geoscientific Model Development Discussions, https://doi.org/10.5194/gmd-2021-14

- V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021), 
"Disjoint interval bound constraints using the alternating direction method of multipliers for geologically constrained inversion: Application to gravity data," GEOPHYSICS 86: G1-G11.
https://doi.org/10.1190/geo2019-0633.1



### Authors and contacts 

Vitaliy Ogarko, Jeremie Giraud, Roland Martin.  
For questions, contact Vitaliy Ogarko via vogarko@gmail.com


### License

Tomofast-x code is licensed under the Creative Commons license, 
the kind that allows commercial use but forces users to acknowledge usage of Tomofast-x and to cite the relevant work.


### Acknowledgments

The authors acknowledge Mark Jessell, Mark Lindsay, and Dimitri Komatitsch.

