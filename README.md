# Tomofast-x v.2.0

Geophysical 3D potential field joint and constrained parallel inversion code.  
Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin, Kim Frankcombe, Taige Liu

Tomofast-x is a 3D parallel inversion platform to run single domain or joint inversion of gravity and magnetic data.
It can invert multiple component magnetic data, and supports inversion for magnetisaion vector, including remanence.
It can use petrophysical constraints such as Gaussian mixture model and disjoint interval bound constraints.
Tomofast-x can run in parallel on desktops and supercomputers, and it has been optimised for both – shared abd distributed memory systems.
Tomofast-x supports the parallel wavelet compression of the sensitivity kernel for significant memory reduction and faster performance.
Tomofast-x allows to invert models with arbitrary surface topography.

### Citing Tomofast-x

If you are using the code (or some of its parts) please cite the following two papers which are detailing Tomofast-x's geophysical calculations and examples of utilisation:

- V. Ogarko, K. Frankcombe, T. Liu, J. Giraud, R. Martin, and M. Jessell (2024),
"Tomofast-x 2.0: an open-source parallel code for inversion of potential field data with topography using wavelet compression", Geosci. Model Dev., 17, 2325–2345, 
https://doi.org/10.5194/gmd-17-2325-2024

- J. Giraud, V. Ogarko, R. Martin, M. Jessell, and M. Lindsay (2021),
"Structural, petrophysical, and geological constraints in potential field inversion using the Tomofast-x v1.0 open-source code", 
Geosci. Model Dev., 14, 6681–6709, https://doi.org/10.5194/gmd-14-6681-2021

### Compiling and running

To compile the code you need: GCC or Intel compiler, and the MPI library (such as OpenMPI).

The makefile is contained in the root folder and should be used to compile Tomofast-x. Compiling the code is a necessary step to be able to run inversions.  
To compile the code run the make command in the code directory as:  
```shell
make
```

To run the code with your parameter file:
```shell
mpirun -np <Number-of-cores> ./tomofastx -j <Parfile path>
```

To run unit tests (serial and parallel):
```shell
./runtests.sh
mpirun -np 3 ./runtests.sh
```

### Running the examples

The input data is contained in the folder ``data``.  
The input parameter file (Parfile) which contains all the parameters of the inversion, is stored in the folder: ``parfiles``.  
The output data is stored in the folder ``output``. The full output folder path is specified in the Parfile parameter ``global.outputFolderPath``.


To run the test example:
```shell
./tomofastx -j ./parfiles/Parfile_mansf_slice.txt
```

If the code runs successfully you will see in the end of the screen log the messages "*Writing the full model...*", and "*THE END*".
To visualize the final model, open in Paraview the file ``Paraview/grav_final_model3D_full.vtk`` (or ``Paraview/mag_final_model3D_full.vtk``), located in the output folder.
For details on other output files, see the User Manual in the ``docs`` folder. 

### Publications using the code

- V. Ogarko, K. Frankcombe, T. Liu, J. Giraud, R. Martin, and M. Jessell (2024),
"Tomofast-x 2.0: an open-source parallel code for inversion of potential field data with topography using wavelet compression", Geosci. Model Dev., 17, 2325–2345, 
https://doi.org/10.5194/gmd-17-2325-2024

- J. Giraud, H. Seillé, M. Lindsay, G. Visser, V. Ogarko, and M. Jessell (2023),
"Utilisation of probabilistic magnetotelluric modelling to constrain magnetic data inversion: proof-of-concept and field application",
Solid Earth, 14, 43–68, https://doi.org/10.5194/se-14-43-2023

- J. Giraud, V. Ogarko, R. Martin, M. Jessell, and M. Lindsay (2021),
"Structural, petrophysical, and geological constraints in potential field inversion using the Tomofast-x v1.0 open-source code", 
Geosci. Model Dev., 14, 6681–6709, https://doi.org/10.5194/gmd-14-6681-2021

- V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021), 
"Disjoint interval bound constraints using the alternating direction method of multipliers for geologically constrained inversion: Application to gravity data", 
GEOPHYSICS 86: G1-G11, https://doi.org/10.1190/geo2019-0633.1

- R. Martin, J. Giraud, V. Ogarko, S. Chevrot, S. Beller, P. Gégout, M. Jessell (2021), "Three-dimensional gravity anomaly data inversion in the Pyrenees using compressional seismic velocity model as structural similarity constraints",
Geophysical Journal International 225(2): 1063–1085, https://doi.org/10.1093/gji/ggaa414

- M. Rashidifard, J. Giraud, M. Lindsay, M. Jessell, and V. Ogarko (2021), "Constraining 3D geometric gravity inversion with a 2D reflection seismic profile using a generalized level set approach: application to the eastern Yilgarn Craton", 
Solid Earth, 12, 2387–2406, https://doi.org/10.5194/se-12-2387-2021

- J. Giraud, M. Lindsay, M. Jessell, and V. Ogarko (2020), "Towards plausible lithological classification from geophysical inversion: honouring geological principles in subsurface imaging", 
Solid Earth, 11: 419–436, https://doi.org/10.5194/se-11-419-2020

- J. Giraud, M. Lindsay, V. Ogarko, M. Jessell, R. Martin, and E. Pakyuz-Charrier (2019), "Integration of geoscientific uncertainty into geophysical inversion by means of local gradient regularization", 
Solid Earth, 10: 193–210, https://doi.org/10.5194/se-10-193-2019

- J. Giraud, V. Ogarko, M. Lindsay, E. Pakyuz-Charrier, M. Jessell, R. Martin (2019), "Sensitivity of constrained joint inversions to geological and petrophysical input data uncertainties with posterior geological analysis", 
Geophysical Journal International, 218(1): 666–688, https://doi.org/10.1093/gji/ggz152

- R. Martin, V. Ogarko, D. Komatitsch, M. Jessell (2018), "Parallel three-dimensional electrical capacitance data imaging using a nonlinear inversion algorithm and Lp norm-based model regularization", 
Measurement, 128: 428-445, https://doi.org/10.1016/j.measurement.2018.05.099

### Authors and contacts 

Vitaliy Ogarko, Jeremie Giraud, Roland Martin, Kim Frankcombe, Taige Liu
For questions, contact Vitaliy Ogarko via vogarko@gmail.com


### License

Tomofast-x code is licensed under the MIT license. 
We request users to acknowledge the usage of Tomofast-x and to cite the relevant work.


### Acknowledgments

The authors acknowledge Mark Jessell, Mark Lindsay, Dimitri Komatitsch, Kim Frankcombe, and Taige Liu.

