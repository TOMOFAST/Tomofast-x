# Tomofast-x  v.1.1

Geophysical 3D potential field joint and constrained parallel inversion code.  
Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.

Tomofast-x is 3D parallel inversion platform to run single domain or joint inversion (gravity and magnetic data).
It can use local weighting of gradient regularization function, global and local petrophysical constraints (Gaussian mixture model and multiple disjoint bound constraints).
Tomofast-x can run in parallel on laptops and supercomputers, using distributed memory systems.

The source code is a companion to the publication detailing Tomofast-x geophysical calculations and examples of utilisation and realistic dataset:
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
To compile the code run the make command in the code directory as:  
```shell
make
```

To run the code with your parameter file:
```shell
./tomofastx -j <Parfile path>
```

To run the code in parallel execute:
```shell
mpirun -np <number CPUs> -j <Parfile path>
```

To run unit tests (serial and parallel):
```shell
./runtests.sh
mpirun -np 3 ./runtests.sh
```

### Running the examples

The input data is contained in the folder ``data``.  
The input parameter file (Parfile) which contains all the paramters of the inversion, is stored in the folder: ``parfiles``.  
The output data is stored in the folder ``output``. The full output folder path is specified in the Parfile parameter ``global.outputFolderPath``.


To run the test example:
```shell
./tomofastx -j ./parfiles/Parfile_mansf_slice.txt
```

If the code runs successfully you will see in the end of the screen log the messages "*Writing the full model...*", and "*THE END*".
To visualize the final model, open in Paraview the file ``Paraview/grav_final_model3D_full.vtk`` (or ``Paraview/mag_final_model3D_full.vtk``), located in the output folder.
For details on other output files, see the User Manual in the ``docs`` folder. 

### Publications using the code

- J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021), 
"Structural, petrophysical and geological constraints in potential field inversion using the Tomofast-x open-source code",
Geoscientific Model Development Discussions, https://doi.org/10.5194/gmd-2021-14

- V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021), 
"Disjoint interval bound constraints using the alternating direction method of multipliers for geologically constrained inversion: Application to gravity data," GEOPHYSICS 86: G1-G11, https://doi.org/10.1190/geo2019-0633.1

- J. Giraud, M. Lindsay, M. Jessell, and V. Ogarko (2020), "Towards plausible lithological classification from geophysical inversion: honouring geological principles in subsurface imaging", 
Solid Earth, 11: 419–436, https://doi.org/10.5194/se-11-419-2020

- J. Giraud, M. Lindsay, V. Ogarko, M. Jessell, R. Martin, and E. Pakyuz-Charrier (2019), "Integration of geoscientific uncertainty into geophysical inversion by means of local gradient regularization", 
Solid Earth, 10: 193–210, https://doi.org/10.5194/se-10-193-2019

- R. Martin, J. Giraud, V. Ogarko, S. Chevrot, S. Beller, P. Gégout, M. Jessell (2021), "Three-dimensional gravity anomaly data inversion in the Pyrenees using compressional seismic velocity model as structural similarity constraints",
Geophysical Journal International 225(2): 1063–1085, https://doi.org/10.1093/gji/ggaa414

- J. Giraud, V. Ogarko, M. Lindsay, E. Pakyuz-Charrier, M. Jessell, R. Martin (2019), "Sensitivity of constrained joint inversions to geological and petrophysical input data uncertainties with posterior geological analysis", 
Geophysical Journal International, 218(1): 666–688, https://doi.org/10.1093/gji/ggz152

- R. Martin, V. Ogarko, D. Komatitsch, M. Jessell (2018), "Parallel three-dimensional electrical capacitance data imaging using a nonlinear inversion algorithm and Lp norm-based model regularization", 
Measurement, 128: 428-445, https://doi.org/10.1016/j.measurement.2018.05.099

### Authors and contacts 

Vitaliy Ogarko, Jeremie Giraud, Roland Martin.  
For questions, contact Vitaliy Ogarko via vogarko@gmail.com


### License

Tomofast-x code is licensed under the MIT license, 
which allows commercial use but forces users to acknowledge usage of Tomofast-x and to cite the relevant work.


### Acknowledgments

The authors acknowledge Mark Jessell, Mark Lindsay, and Dimitri Komatitsch.

