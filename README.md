# Tomofast-x  v.1.0.1

Geophysical 3D potential field joint and constrained inversion code.
Vitaliy Ogarko, Jeremie Giraud, Roland Martin.

TOMOFAST-x is inversion platform to run single domain or joint inversion (gravity and magnetic data). It can use local weighting of gradient regularization function, global and local petrophysical constraints (Gaussian mixture model and multiple disjoint bound constraints).
TOMOFAST-x is can run in parallel on laptop and supercomputers. It is licensed under the do-what-the-fuck-license-TODO: Creative Commons, the kind that allows commercial use but forces users to acknowledge usage of Tomofast-x and to cite the relevant work. 

The source code is a companion to the publication detailing Tomofast-x geophysical aspects and examples of utilisation and realistic dataset:
J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell, 2020: Structural, petrophysical and geological constraints in potential field inversion using the Tomofast-x open-source code, Geoscientific Model Development Discussions [to be submitted soon...]. 

Companion datasets to previous publications using Tomofast-x are provided in: 
<Zenodo Link on Smart Gradient Paper - done, only to copy and paste> 
<Zenodo Link on Mansfield example - done, only to copy and paste>  
<Zenodo Link on Mansfield example - > 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 
The source code can be downloaded from a repository at: https://github.com/TOMOFAST/
Details about the utilisation of Tomofast-x are provided in the User Manual. 

### Prerequisites

To be able to compile the code you need to have installed gcc v.4.9 or more recent. 
The example contained in the archive requires less than 1 Gb of memory to install, store and run. 


### Compiling

The makefile is contained in the root folder and should be used to compile Tomofast-x. Compiling the code is a necessary step to be able to run inversions. 


## Running the test

The data contained in the archive comprises all the necessary information and files to run successfully. 
The following is true if all the parameters are kept unchanged. 
Input data is contained in folder : /.../TODO
Output data is stored in folder: /.../TODO
The parameter file (called "parfile") which contains all the paramters of the inversion, is stored in folder: /.../TODO
Inversion can be run in the folder the archive was extracted in by executing the following command: [...]

Detailed information to run Tomofast-x can be found in the User Manual. 


## To run in parallel on local computer (example)
mpirun -np 4 ./tomofast3D -j ./parfiles/Parfile_MASTER.txt


### And coding style tests [remove?]

Explain what these tests test and why

```
Give an example. TODO OR TO REMOVE SECTION?
```

## Deployment

Add additional notes about how to deploy this on a live system


## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.


## Versioning [REMOVE?]

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 


## Authors and contacts 

Vitaliy Ogarko, Jeremie Giraud, Roland Martin. 
TODO: do we add a personal page or leave it like that? do we give email adresses? 

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project 
[?REMOVE THE ABOVE LINE?]


## License

This project is licensed under the do-what-the-fuck-license - see the [LICENSE.md](LICENSE.md) file for details
[TODO: LICENSE.md FILE]

## Acknowledgments

The authors acknowledge Mark Jessell, Mark Lindsay, and Dimitri Komatitsch.

