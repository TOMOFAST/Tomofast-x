![Windows]([https://img.shields.io/Windows](https://img.shields.io/badge/any_text-you_like-blue))![WSL](https://img.shields.io)![macOS](https://img.shields.io)![Linux](https://img.shields.io)
# Tomofast-x v.2.0

A parallel inversion code for 3D geophysical potential field data, supporting joint and constrained inversion.  
Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin, Kim Frankcombe, Taige Liu

**Tomofast-x** is a powerful 3D parallel inversion platform designed for single-domain or joint inversion of gravity and magnetic data. It also supports the inversion of gravity gradiometry data (FTG) and can handle multiple-component magnetic data. The platform is capable of inverting for the magnetization vector, including remanence, while incorporating petrophysical constraints.

Optimized for both shared and distributed memory systems, Tomofast-x can run in parallel on desktops and supercomputers, ensuring scalability across different computing environments. It also features parallel wavelet compression of the sensitivity kernel, reducing memory usage and enhancing performance. Furthermore, Tomofast-x supports the inversion of models with arbitrary surface topography, making it highly versatile for various geological scenarios.

### Documentation

For detailed usage instructions and comprehensive documentation, please refer to the **[Tomofast-x User Manual](https://github.com/TOMOFAST/Tomofast-manual)**.

### Citing Tomofast-x

If you use the code (or any of its components), please cite the following two papers, which provide detailed explanations of Tomofast-x's geophysical calculations and examples of its application:

- V. Ogarko, K. Frankcombe, T. Liu, J. Giraud, R. Martin, and M. Jessell (2024),
"Tomofast-x 2.0: an open-source parallel code for inversion of potential field data with topography using wavelet compression", Geosci. Model Dev., 17, 2325–2345, 
https://doi.org/10.5194/gmd-17-2325-2024

- J. Giraud, V. Ogarko, R. Martin, M. Jessell, and M. Lindsay (2021),
"Structural, petrophysical, and geological constraints in potential field inversion using the Tomofast-x v1.0 open-source code", 
Geosci. Model Dev., 14, 6681–6709, https://doi.org/10.5194/gmd-14-6681-2021

### Installing and running

See the [Installation instructions](https://github.com/TOMOFAST/Tomofast-manual/blob/main/install.md) for detailed installation instructions.

To run the code with your parameter file:
```shell
mpirun -np <Number-of-cores> ./tomofastx -p <Parfile path>
```

To run unit tests (serial and parallel):
```shell
mpirun -np 1 ./runtests.sh
mpirun -np 3 ./runtests.sh
```

### Running the examples

#### Browser-based examples

You can run examples in your browser without local installation using our [Google Colab notebook](https://colab.research.google.com/github/TOMOFAST/Tomofast-examples/blob/main/Tomofast-x_examples.ipynb), which includes visualization of output models.

#### Local execution

To run the test example:
```bash
mpirun -np 1 ./tomofastx -p ./parfiles/Parfile_mansf_slice.txt
```

If the code runs successfully you will see at the end of the screen log the messages "Writing the full model...", and "THE END". To visualize the final model, open in Paraview the file `Paraview/grav_final_model3D_full.vtk`, located in the output folder. For details on other output files, see the [Tomofast-x User Manual](https://github.com/TOMOFAST/Tomofast-manual).


### Publications using the code

See the [Publication list](https://github.com/TOMOFAST/Tomofast-manual/blob/main/publications.md) for publications using the Tomofast-x code.

### Authors and contacts 

Vitaliy Ogarko, Jeremie Giraud, Roland Martin, Kim Frankcombe, Taige Liu

For questions, contact Vitaliy Ogarko via vogarko@gmail.com


### License

Tomofast-x code is licensed under the MIT license. 
We request users to acknowledge the usage of Tomofast-x and to cite the relevant work.


### Acknowledgments

The authors acknowledge Mark Jessell, Mark Lindsay, Dimitri Komatitsch, Kim Frankcombe, and Taige Liu.

