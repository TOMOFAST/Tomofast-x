<h1> TOMOFAST-x </h1> 

<h2>Tomofast User manual </h2>

<h3> Integrated inversion platform </h3> 

Contributing Writers of the User Manual

``` Ashwani Prabhakar, Jérémie Giraud ```

Main Developers of TOMOFAST-x

``` Vitaliy Ogarko, Jérémie Giraud, Roland Martin ```

Project supervision
```Mark Jessell```

Release date
30/07/2019

![01_first_page_image_logo](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/001-.JPG)

|Version | Date      | Revision Description and comments |
|:-----:|  :---:    |-----------|
|0.1	   |01/06/2019 | User’s manual template and checklist and initial work.|
|1.0	   |29/07/2019 | Release version of Tomofast used of predating publications.| 
|1.1     |08/08/2019 | Added Table with summary of parameters | 


<h3> Copyright </h3>

The information contained in this documentation is the exclusive property of CET, UWA. This work is protected under the copyright law and other conventions. No part of this work may be reproduced or transmitted in any form or by any means, whether it is mechanical, electronic, photocopying and recording, or by any information distribution system without the written permission from the concerned authorities. All concerned requests should be sent to **Jérémie Giraud** , Research Fellow, Centre for Exploration Targeting, SES, The University of Western Australia, Perth, WA 6009, Australia. Email Address as follows - jeremie.giraud@uwa.edu.au 

<h3> Disclaimer </h3>

The writers of this manual and of Tomofast-x have taken care to ensure the accuracy and quality of this User manual and of the software. However, all material is provided without any warranty whatsoever. While we make every effort to ensure that the manual is accurate and up to date, there remain the possibility of human error, and utilisation Tomofast-x may by no means provide exact, certain answers to any modelling problem. We do not guarantee, and accepts no legal liability whatsoever arising from or connected to, the accuracy or completeness of any material contained in this manual or obtained through usage of Tomofast-x. The information contained in this documentation is subject to be changed without any notice/ legal procedure. You are reading this manual and using Tomofast-x at your own risks. 

 <h3> PREFACE </h3>
 
 This User manual is intended for new Users with little or no experience using TOMOFAST-x. The goal of this documentation is to give a broad overview of the main functions of TOMOFAST-x and some basic instructions on how to set up and administer its functionality. This User manual includes a description of TOMOFAST-x functions and capabilities, contingencies, and step-by-step procedures for its access and use. This documentation will concentrate on demonstrating interaction with TOMOFAST-x.
Every effort has been made to ensure that this documentation remains an accurate representation of the functionality of TOMOFAST-x. As we all are familiar about the development of the software which continues even after its release, same goes with TOMOFAST-x. Our team is working on its development and will be doing its best to update this documentation in time.  
We would highly appreciate any feedback on this User manual. For any suggestion, please contact Jérémie Giraud, Research Fellow, Centre for Exploration Targeting, SES, The University of Western Australia, Perth, WA 6009, Australia. Email Address as follows - jeremie.giraud@uwa.edu.au 
 
 <h3> Acknowledgements </h3>

We would like to thank our sponsors and our funding bodies who participated directly or indirectly to the development of TOMOFAST-x and the redaction of this manual. Without them, it would not be possible to achieve this. Appreciation is expressed to the CALMIP supercomputing centre (Toulouse, France), for their support through Roland Martin’s computing projects no. P1138_2017 and no. P1138_2018 and for the computing time provided on the EOS and Olympe machines when testing Tomofast on large 3D models. The work has been supported by the Mineral Exploration Cooperative Research Centre whose activities are funded by the Australian Government's Cooperative Research Centre Programme. This is MinEx CRC Document 20**/ . The authors of this document and the developers of TOMOFAST-x are thankful to the Australian Federal Government for granting an International Postgraduate Research Scholarship to Jeremie Giraud. Vitaliy Ogarko acknowledges the Australian Research Council Centre of Excellence for All Sky Astrophysics in 3-D (ASTRO 3-D) for supporting some of his research efforts. They acknowledge the state government of Western Australia for supporting Mark Jessell through the Geological Survey of Western Australia, Royalties for Regions and the Exploration Incentive Scheme.

 <h3> How You Can Contribute to TOMOFAST-x </h3>

There are a number of ways that you can contribute to help in making TOMOFAST-x a better system. You can share ideas and suggestions with the authors of the code and the documentation. You are also welcome to contribute  through the github project: <link here>. For your information,

TOMOFAST-x will be available on the public platform from the first quarter of 2020. By 2022, Tomofast-x will be integrated in the [Loop platform]( https://pages.github.com/) for more information.

Potential contributions which includes the integration of spatial trends, geochemical information in the geological conditioning process, tests involving lithology- and location-dependent petrophysical uncertainty to better account for spatial variability of rock properties will be highly appreciated. You can contact the authors for your respective queries, suggestions and for further improvement of TOMOFAST-x.

If you find TOMOFAST-x useful or have questions regarding its potential usage or extensions, please do not hesitate while contacting the authors. 

 <h3> Distribution </h3>
 Tomofast is licensed under (discription here). To be completed upon public release. TODO. 
 
** Attribution-NonCommercial-ShareAlike **
** CC BY-NC-SA ** 

<h2> 1. A BRIEF INTRODUCTION TO TOMOFAST-x </h2>
<ol>
<li> Welcome to TOMOFAST-x. This manual is intended to help the User in getting started using TOMOFAST-x software and to illustrate the methods and procedures involved in conducting a Gravity, Magnetic and Joint Inversion projects. If the User is new to TOMOFAST-x, this manual is a great place to start— User can learn how to use TOMOFAST-x in order to solve problems related to geophysical inversion. It integrates both statistical petrophysical constraints, probabilistic geological models, cross-gradient constraints and local gradient constraints. It can be used to investigate how uncertainty propagates from the geological and petrophysical input measurements to the recovered lithological model.</li>

<li> For general information, the source code of TOMOFAST-x follows the object-oriented FORTRAN 2008 standard. The design of TOMOFAST-x utilizes classes derived to account for the mathematics of the problem. This permits to reduce software complexity, thereby facilitating the addition of new functionalities (Giraud et al. (2019a)). More information about the implementation and the code scalability on supercomputers is provided in [Giraud et al. 2020, ToBeWritten].</li>

<li> TOMOFAST-x operates on the least-square geophysical inverse problem equation as mentioned below</li>



  ![equation 01](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/002-eqn_01.JPG)	 (1)
  
  ***m = \[***m<sub>1</sub>, m<sub>2</sub>]*** i.e. ***m<sub>1</sub>*** stands for gravity and ***m<sub>2</sub>*** stands for magnetic 
 
 Here \,
   -***d*** refers to geophysical data
	
  -***m<sub>p</sub>*** refers to prior model
	  
  - **g(m)** refers to the respective model
	 
  - **‖W_d (d-g(m))‖<sub>2</sub><sup>2</sup>**  represents data term. For more information, please refer Lines & Treitel (1984).
	
  - **‖W_m (m-m_p )‖<sub>2</sub><sup>2</sup>** represents model term. For more information, please refer Hoerl & Kennard (1970).
	
  - **α‖W_H ∇m‖<sub>2</sub><sup>2</sup>**  represents structure term which stands for local gradient regularization. ***∇m** represents model gradient. ***α*** represents damping- gradient constraints which has been explained in the section ‘Damping- gradient constraints’ For more information, please refer Li and Oldenburg (1996).
   
 - **‖W_P P(m)‖<sub>2</sub><sup>2</sup>** represents petrophysics term i.e. clustering constraint. ***P<sub>(m)</sub>*** represents petrophysical distribution which stands as follows \, 
 
![equation 02](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/003-eqn_02.JPG)  (2)


Where,

**ω<sub>k</sub> = 1/n<sub>f</sub>**  always if no geol available
**ω<sub>k</sub>= ψ<sub>k, i</sub>**   in the i<sup>th</sup> cell otherwise	

Please refer Giraud et al. (2017) for further information.

![equation 03](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/004-eqn_03.JPG)  (3)

- ***W<sub>m</sub>=α<sub>m</sub>d<sub>m</sub> , W<sub>m</sub>*** represents model term weighting. d<sub>m</sub> is a diagonal matrix (covariance value) which is the last column of the Input model file.User can have a look in the section … **(add cross reference of Input Model Voxet File)** where the authors have explained about input model voxet file. **α<sub>m</sub>** is the model damping parameter. For more information about model damping parameter, please refer section.

- ***W<sub>H</sub>=d<sub>H</sub>, W<sub>H</sub>***  represents the weighting for the local gradient regularization. ***d<sub>H</sub>*** is a diagonal matrix which is the last column of input model file calculated from geological information. If no geological information is present, User can implement ***I*** and can come back to regular smoothness constraints. Please refer to Brown et al. (2012), Yan et al. (2017),Giraud et al. (2019b) for further information regarding spatially varying weights which affects local conditioning 

- ***W<sub>p</sub>=α<sub>p</sub> I, W<sub>P</sub>*** represents the weighting for the petrophysics term. This term is responsible for prior petrophysical information. α<sub>p</sub> is the clustering weight which can be visualized in the section of Clustering Constraints in the explanation of Parameter file.
***I*** stands for Identity Matrix

Note: items stricken out will be part of a future release. 

<li> More information about applications of TOMOFAST-x and case studies, please refer to Giraud et al. (2017), (2019a), (2019b), Martin et al. (2018).</li>

<li> The source code and the parameter files will be available on Github, project TOMOFAST-x (CET Geophy et al. 2019). In the testing presented here, we used the publicly available structural geological model of the Mansfield area (Victoria, Australia) of Pakyuz-Charrier (2018) as the reference geological structural framework.   The 2D section shown in this documentation corresponds to an extended version of the cross-section extracted from the Mansfield model using in Giraud et al. (2017).</li> 

<li> Additional reference property models, synthetic geophysical data, inversion model and recovered lithological models shown or discussed in this document are made available by Giraud et al. (2018) in an ASCII format usable by Tomofast-x using doi: 10.5281/zenodo.1003105.</li>
</ol> 


<h2> 2. Towards first TOMOFAST-x run </h2>

<h3>  Basic Requirements for Running TOMOFAST-x 3</h3>
<ol>
<li> Software/ Operating System 
<ol>
<li> Ubuntu  <b>or</b> </li>
<li> Windows ( need to install/setup Windows Linux subsystem)</li>
</ol>	
</li>
<li> Environment – GNU/ Linux </li>
<li> TOMOFAST-x need to be compiled using <b>gcc 4.9</b> or <b>higher</b> and the appropriate MPI libraries.</li>
</ol> 

**Linux (Ubuntu) users can jump to the _section 3.2_**
 
Windows users should install (activate) [Windows Subsystem](https://docs.microsoft.com/en-us/windows/wsl/install-win10)  for Linux. Then use bash on Ubuntu to run applications designed for Linux systems.
 
  The Windows Subsystem for Linux lets developers run GNU/Linux environment including most command-line tools, utilities, and applications directly on Windows, unmodified, without the overhead of a virtual machine. Its installation for TOMOFAST-x is detailed below

 <h2> 3. Installation Guide Windows users Start here </h2>

<h3> 3.1 Install Bash on Windows. </h3>
 With Windows Systems, the installation of the Windows Linux Subsystem is necessary (see description below). After it is installed, the procedure described is the same as with Linux systems. 
 Note: Windows Linux Subsystem can also be installed using Windows PowerShell.	
 
<ol>
<li> Open Settings.
![setting_icon](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/005-settings.png)</li> 

<li> Click on Update & security.</li>
<li> Click on For Developers as shown in Figure 1.</li>
<li> Under "Use developer features", select the Developer mode option to setup the environment to install Bash.</li>

![fig01](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/006_Figure_01.PNG)

<li> On the message box, click yes to turn on developer mode as shown in Figure 2.</li>

![fig02](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/007-Figure_02.PNG)

<li> After the necessary components install, User needs to restart their computer.</li>
<li> Once the computer reboots, open Control Panel.</li>
<li> Click on Programs.</li>
<li> Click on Turn Windows features on or off as shown in Figure 3.</li>

![fig03](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/008-Figure_03.PNG)

<li> Check the Windows Subsystem for Linux (beta) option as shown in Figure 4.</li>
<li> Click OK.</li>

![fig04](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/009-Figure_04.PNG)

<li> Once the components installed in the computer, click the restart now button to complete the task.</li>
<li> After the computer restarts, User will notice that Bash will not appear in the "Recently added" list of apps, this is because Bash is not   installed yet. Now that User has setup the necessary components, use the following steps to complete the installation of Bash.</li>
<li> Open Start, do a search for bash.exe, and press Enter.</li>
<li> On the command prompt, type Y and press Enter to download and install Bash from the Windows Store.</li>

![fig05](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/010-Figure_05.PNG)

<li> Then User will need to create a default UNIX user account as shown in Figure 5. This account does not have to be the same as their windows account. It is setup independently and will be used only for your Windows Linux Subsystem. Enter the Username in the required field and press Enter (User cannot use the Username "admin").</li>
<li> Close the "bash.exe" command prompt.</li>

</ol>
Now the User have completed the installation and setup, User can open the Bash tool from the Start menu like User would do with any other application.

<h3> 3.2  Prerequisites to run TOMOFAST-x  </h3>
 

After this point prerequisites installation is the same for both Windows and Linux users. Windows users should follow the installation commands below in bash on windows Linux subsystem. Ubuntu (or Linux) users should follow the same instructions in the Ubuntu (Linux) terminal.

Prior to the first TOMOFAST-x run, User needs to install the **[gcc](https://gcc.gnu.org/)** and **[gfortran](https://gcc.gnu.org/wiki/GFortran)** libraries in the Windows Linux Subsystem, and **[OPENMPI](https://www.open-mpi.org/)**  to run TOMOFAST-x in parallel on their computer.


 ```sudo apt-get install gcc ``` to install gcc.  
 
 ```sudo apt-get install gfortran ```to install gfortran.
 
 ``` sudo apt-get install g++ ```to get g++. 
 
 ``` sudo apt-install make ```to get make.

 ``` sudo apt-get install flex ```;to get flex (optional)

Please note that User needs to install all these dependencies in order to install **OPENMPI**, which will allow TOMOFAST-x to run in parallel. 

Get **[OPENMPI-2.1.1](https://www.open-mpi.org/software/ompi/v2.1/)**  using web browser. Then follow the [build](https://www.open-mpi.org/faq/?category=building) information.


``` gunzip -c openmpi-2.1.1.tar.gz | tar xf - ``` unzip using any unzip programme (eg gunzip)

``` cd openmpi-2.1.1 ```                          go to the openmpi unzipped folder 

``` ./configure --prefix=/usr/local ```           to configure the openmpi

<...lots of output...>

``` make all install ``` install opnmpi

Installation of OPENMPI-2.1.1 gets completed. 

``` sudo nano .bashrc ```

Add this line to your **.bascre file ** ```export LD_LIBRARY_PATH:=$PATH:/usr/lib/openmpi/lib ```

Now, User is ready to run the executable of TOMOFAST-x called **‘tomofast3D’** (provided in the archive containing this manual) or compiled using the source code downloaded from the github project **<addlink>** provided above. 
	
<h3> 3.3 Errors while/ after Installation</h3>

There may be some potential errors while/ after Installation. One of the example has been shown in this section. If User gets the error (libgfortran.so.4) as shown below in Figure 6.

![fig06](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/011-Figure_06.PNG)

User need to install libgfortran4 using the following commands:


```sudo apt-get install gcc-7``` 
```sudo apt-get install g++-7``` 
```sudo apt-get install libgfortran4```  to install libgfortran4 

After following the above steps, User will be able to solve the issue and will be able to run the executable of TOMOFAST-x.

<h2> 4. GETTING STARTED WITH TOMOFAST-x </h2>

In this section, we would like to describe about the starting of the TOMOFAST-x. We would also be discussing about the executable of TOMOFAST-x i.e. “tomofast3D” briefly. 

<h3> 4.1 Invoking TOMOFAST-x using a command line </h3>

<ol>
<li> Change your directory to the folder where the executable tomofast3D exists and run the same as shown below in Figure 7.</li>

![fig7](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/012-Figure_07.png)

<li> Command to run - mpirun ** -n 1 ./tomofast3D -j ./Parfile_mansf_slice.txt | tee out.txt ** The explanation of the command can be found in the Figure 8. </li>

![fig8](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/013-Figure_08.png)

<li> Now the executable is ready to run as shown below. If you are successful, a copy of the parameter file will be printed to the screen, followed by a log of the inversion as shown in Figure 9.</li>

![fig9](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/014-Figure_09.png)

</ol>

<h3> 4.2 About Executable/ compiling  </h3>


<h2> 5. PARAMETER FILE </h2> 

In this section, we would like to introduce the User to Parameter File for TOMOFAST-x which we also refer to as ‘Parfile’. This is the file which should be referred to visualize the type of inversions, different parameters, features, etc. which have been involved with TOMOFAST-x. Please refer to this section for further explanation of the Parameter File. 

<h3> 5.1 First part of parameter file </h3>

What we refer to here as the first part is shown below in Figure 10.

![fig 10](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/015-Figure_10.png)

The explanation of Figure 10 follows. 

<h3> 5.2 GLOBAL</h3> 

<i>Table 1. GLOBAL section of parfile.</i> 

| Parameter   | Value for example case |Range/remark  |
| :--- | :--- | :--- |
| \******* GLOBAL ******  	|N/A   |N/A |
|path to the output folder 	|output/mansf_slice/|N/A|

<h3>5.3 DIMENSIONS</h3>

This section concerns electrical capacitance tomography and is to be considered only for this kind of modelling. 

<ol>
<li>It contains the feature where user can edit the number of grid points. But for this instance, User need not to think about this as this has been set up at the beginning of the project TOMOFAST-x.</li>

<li>Also, User need not to think about ntheta and nel because these two haven’t got any relevance while running TOMOFAST-x. Just for information, these two are relatable while using ECT (Electrical Capacitance Tomography). For that, ntheta must be sufficiently larger multiple of nel (nel  represents number of electrodes). For further information, please refer to ECT, Martin et al. (2018).</li>
</ol>

<i>Table 2. DIMENSIONS section of parfile. </i>

|Parameter|Value for example case|Range/remark |
|:---  |---:|:---|
| \******* DIMENSIONS ******* 	|  	|N/A|
|nr (number of grid points)	|36 	|Survey dependant  |
|ntheta				|0	|Survey dependant  |
|nz				|0	|Survey dependant  |

<h3>5.4 GEOMETRY </h3>

This section concerns electrical capacitance tomography and is to be considered only for this kind of modelling. 

The features which are present in this section have been already set up in the starting of the project TOMOFAST-x.User need not to think about the parameters/ features in this section. These features should be remain unchanged while running TOMOFAST-x. Just for information, this section is related with the functionalities of the source code which is beyond the scope of this documentation.

<i>Table 3. GEOMETRY section of parfile.</i> 

|Parameter|	Value for example case|	Range/remark|
|:----|-----:|:----|
| ******* GEOMETRY ******** |  |	 	 N/A|
|nel (number of electrodes)|             	36|	Survey dependant | 
|nrings (elec rings 1 or 2) |            	3|	Survey dependant | 
|kguards (number of guards)|             	0|	Survey dependant | 
|fixed electrodes by geometry|           	0|	Survey dependant | 
|refinement (NO = 0, YES = 1) |          	0|	Survey dependant | 
|location R1  |                          	0.045|	Survey dependant | 
|location R2  |                          	0.06|	Survey dependant | 
|location R3 |                           	0.07|	Survey dependant |
|height sensor |                         	0.2|	Survey dependant | 
|space between guards|                   	0|	Survey dependant | 
|space between electrodes(deg)|          	0|	Survey dependant | 

<h3> 5.5 MODEL </h3>

This section concerns electrical capacitance tomography and is to be considered only for this kind of modelling. 

The features which are present in this section have been already set up in the starting of the project TOMOFAST-x. These parameters should remain unchanged while running TOMOFAST-x.

<i> Table 4. MODEL section of parfile. </i>

|Parameter|	Value for example case	|Range/remark| 
|:----|-----:|:----|
|******* MODEL *******     |                    | 	 N/A|
|num bubbles (0=no bubbles)|           	4	|Survey dependant|  
|location of the bubbles   |data/ECT/bubble_4vert.dat|	Survey dependant|  
|absolute permittivity     |           	1	|Survey dependant  |
|permittivity air          |          	1	|Survey dependant  |
|permittivity isol tube    |         	3.5	|Survey dependant  |
|permittivity oil          |        	2	|Survey dependant  |


<h3> 5.6 SOLVER parameters </h3>

This section concerns expert parameters that do not need to be adjusted for regular usage of Tomofast. 

<ol>

<li> This section contains the solver parameters which includes expert numerical tuning and have been applied in inversion process. User need not to change these parameters while running TOMOFAST-x. They have been set up in the beginning of the project TOMOFAST-x.</li>

<li> PCG stands for preconditioned conjugate gradients. Here, 0 stands for NO i.e. not applying PCG and values greater than 0 stand for YES i.e. apply PCG. In TOMOFAST-x, we have applied PCG and have kept it as 1.</li>

<li> For using L2 norm, apply 1. For using greater than L2 norm, apply 2. In TOMOFAST-x, we have kept it to L2 norm.</li>

<li> The remaining parameters in this section are expert parameters which are not needed to be changed by the User while running TOMOFAST-x.</li>

</ol>

<i> Table 5. SOLVER section of parfile </i>

|Parameter	                | Value for example case |Range/remark| 
|:----|-----:|:----|
| ******* SOLVER parameters ********* |                 |	N/A       |
| PCG precond (0=NO, YES>0)     |              	1	|Expert parameter |  
| relaxation omega1 PCG precon  |           	 0.8d0	|Expert parameter | 
| type of norms: 1=L2, 2=max    |         	1	|Survey dependant | 
| max num. of lin. solv. iters  |         	1000	|Survey dependant |  
| output freq linsolver (iters) |         	20	|Expert parameter |  
| tolerance of linear solver    |         	 1.d-12	|Expert parameter |


<h3> 5.7 GRAVITY / MAGNETISM parameters </h3> 

<ol>
	
<li> This section contains some of the gravity and magnetic parameters which are as follows.</li>

<li> Grid size – It can be changed according to the model. For example, tested input models (github_link_) contains got 8192 no. of cells, organised as  nx*ny*nz = 2*128*32. </li> 
 
<li> Model Format – It stands for the format of the input model file. Here, 
<ul>
<li> 1 → stands for Voxet data format </li>
<li> 2 → stands for Noddy data format </li>
<li> 3 → stands for GeoModeller data format </li>
</ul>
</li>

<li> User need not to change the input model file format in order to run TOMOFAST-x. Voxet format has been set as a default for TOMOFAST-x.</li> 
	
<li> User can edit the coarse size of the data format, if Noddy data format has been selected. In order to run TOMOFAST-x, user need not to look after Noddy data format. </li>

<li> User can put the location of the gravity model and grid file/ magnetic model and grid file which can subsequently be used to generate the forward data set which will be stored in the file as shown below in the Figure 11. </li>

![figure 11](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/016-Figure_11.png)

<li> User can also chose Depth Weighting according to – 

<ul>
<li>1→ Power – Here, weights are proportional to the inverse of the distance (i.e. between the data point and the prediction location) raised to the power value p. As a result, as the distance increases, the weights decrease rapidly. For more information, please refer to Li and Oldenburg (1996) & (1998) respectively. </li>

<li>2→ Sensitivity - The depth weighting function is selected based on sensitivity analysis. For more info, please refer – http://dx.doi.org/10.1190/1.1512749 </li>
<li>3→ Integrated Sensitivity technique – For more information, please refer to Portniaguine and Zhdanov (2002)</li>
</ul>
</li>
</ol>

<i> Table 6. GRAVITY / MAGNETISM parameters section of parfile. </i>

|Parameter			|Value for example case |Range/remark | 
|:----				|:-----			|	:----|
|******* GRAVITY / MAGNETISM parameters | 	 	| N/A |
| grid size (nx, ny, nz)	      |2 128 32		| Survey dependant|  
|model format (1-vox,2-Noddy,3-GeoMod)|  	1	|1-3		|
|coarse size x y z (for Noddy models) | 1 1 1 	 	|N/A		|
|model section beginning (Noddy) x y z|  0.d0 0.d0 0.d0	|N/A		|
|grav model and grid file             | mansf_slice_input/true_model_grav.txt	      |Survey dependant |
|mag  model and grid file             | mansf_slice_input/true_model_grav_replaced.txt|Survey dependant|  
|Noddy model .g00 file                | NILL		|N/A		|
|Noddy model .g12 file                | NILL		|N/A		|
|Depth weighting (1-pow,2-sens,3-isens)| 3		|1-3  		|


<h3>5.8 Second part of parameter file </h3>
                                                                  
The second part of the parameter file is given in Figure 12.

![figure12](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/017-Figure_12.png)

The explanation of Figure 12 follows.

<h3> 5.9 GRAV/ MAG DATA parameters </h3>

<ol>
<li>This section includes some of the features where user can edit the gravity and magnetic data parameters.</li>
<li>User can edit the gravity/ magnetic number of data according to their respective model.</li>
<li>User can vary the type of the data format in the parameter file according to their preference. For the inversion scheme described in this manual, ‘points data format’ should be used. </li>
<li>GDA stands for Geocentric Datum of Australia (GDA). For the inversion scheme described in this manual, GDA data format is not required at this stage. For running </li>
	
<li>TOMOFAST-x,User need not to change the data format. Keep the corresponding parameter as 1.</li>
	
<li>User can change the grids based on the grid file or model file.
<ul>
<li>1→ Grid based on grid file</li>
<li>2→ Grid based on model file</li>
</ul>
</li>

<li>User can put the location of gravity and magnetic grid files respectively.</li>
<li>User can also put the location of gravity and magnetic data files respectively.</li>
<li>User need not to worry about clipping threshold and calculation without sensitivity in order to run TOMOFAST-x.</li>

</ol>

<i>Table 7.  GRAV / MAG DATA parameters section of parfile. </i>

|Parameter	|Value for example case	|Range/remark| 
|:----				|:-----			|:----	|
|******* GRAV / MAG DATA parameters ******  | N/A 	| N/A	|
|grav number of data                	|256		|Survey dependant|
|mag number of data                     |256		|Survey dependant|
|data format (1-points, 2-GDA)          |1	 	|N/A|
|inverse data density (for GDA)         |1	 	|N/A|
|data section beginning x y z (for GDA) |0.d0 0.d0 0.d0	|N/A|
|grid based on grid file(1) or model(2) |1		|Survey dependant|
|grav grid file                         |mansf_slice_input/data_grid.txt	|Survey dependant|
|mag  grid file                         |mansf_slice_input/data_grid.txt	|Survey dependant|
|grav data file                         |output35/mansf_slice/grav_calc_read_data.txt	|Survey dependant|
|mag  data file                         |output35/mansf_slice/mag_calc_read_data.txt	|Survey dependant|
|grav data clipping threshold (0-no)    |0.d0		|N/A |
|mag  data clipping threshold (0-no)    |0.d0		|N/A |
|calc. data without sensit (debug)      |0		|N/A |


<h3> 5.10 PRIOR MODEL </h3>

<ol>
<li>This section represents the features of a Prior Model .</li>
<li>User can select the type of prior model from the below options available to TOMOFAST-x.
<ul>
<li>1→ TOMOFAST-x will smooth the model provided as grid file using a Gaussian filter and use it as prior model. This option is recommend only for tests on synthetic data.</li> 
<li>2→ Setting up the gravity/ magnetic model according to the desire.User can put the reasonable values across the respective options i.e. “set prior model grav/ mag” present in this section.</li>
<li>3→User can use this in order to read the prior model from the model file.</li>
<li>4→User can also chose this option in order to fix the value of their respective top laye   r while performing their respective inversion. But for now, this option should be avoided.</li>
</ul>
</li>	
<li>User can put the path of the  gravity/ magnetic prior models file relative to the executable.</li>
</ol>

<i>Table 8. PRIOR MODEL parameters section of parfile</i>

|Parameter				|Value for example case	|Range/remark|
|:----					|:-----		|:----		     |
|******* PRIOR MODEL *************** 	| N/A		| N/A		 |
|type(1-smooth,2-set,3-file,4-toplayer) |2		|Survey dependant|
|smooth grav prior model (# times)      |0		|N/A		 |
|smooth mag prior model (# times)       |0		|N/A		 |
|set prior model grav (if no smoothing) | 0.d0		| 		 |
|set prior model mag (if no smoothing)  | 0.d-9		|Survey dependant|
|grav prior model file                  |mansf_slice_input/mod_start_geol.txt	|Survey dependant|
|mag prior model file                   |mansf_slice_input/Wh.txt		|Survey dependant|


<h3>5.11 STARTING MODEL</h3>

<ol>
<li>This section represents the parameters  of a Starting Model.</li>
<li>User can select the type of starting model.
<ul>
<li>1→User can put the starting model equal to prior model using this option.</li>
<li>2→ Setting up the gravity/ magnetic model according to the desire.User can put the reasonable values across the respective options i.e. “set prior model grav/ mag” present in this section.</li>
<li>3→User can use this in order to read the prior model from the model file.</li>
</ul>
</li>
<li>User can put the location of gravity/ magnetic starting model file relative to the executable.</li>
</ol>

<i>Table 9. STARTING MODEL constants parameters section of parfile.</i>

|Parameter				|Value for example case	|Range				| 
|:----					|:-----			|:----	     			|
|******* STARTING MODEL ***************	|			|				| 	 
|model type (1-eq to prior,2-set,3-file)|1			|Survey dependant		|
|set prior model grav (if no smoothing) |0.d0			|Survey dependant		|
|set prior model mag (if no smoothing)  |10.d-9			|Survey dependant		|
|grav prior model file                  |mansf_slice_input/mod_start_geol.txt 	|Survey dependant|
|mag prior model file                   |mansf_slice_input/Wh.txt 		|Survey dependant|

<h3>5.12 MAGNETIC constants</h3>

This section includes some of the magnetic constants which are self-explanatory. They depend on the studied area. Generally, they are set at the beginning of a project, need not be changed while running TOMOFAST-x.

<i>Table 10. MAGNETIC constants parameters section of parfile.</i>

|Parameter			  |Suggested Input value|	Range/remark| 
|:----				  |:----- 	|:----|
|******* MAGNETIC constants ******| N/A	  	|N/A	|
|mag field inclination            | 75.d0	|Survey dependant|
|mag field declination            | 25.d0	|Survey dependant|
|ambient field inclination        | 75.d0	|Survey dependant|
|ambient field declination        | 25.d0	|Survey dependant|
|ambient field intensity (nT)     | 50000.d0	|Survey dependant|

<h3>5.13 GRAVITY constants </h3>

This section includes some of the gravity constants which are self-explanatory and model dependent. Generally, they are set at the beginning of a project, need not be changed while running TOMOFAST-x. 

<i>Table 11. GRAVITY constants parameters section of parfile.</i>

|Parameter	|Value for example case	|Range/remark|
|:----				  |:----- 	|:----|
|******* GRAVITY constants ****** |N/A		|N/A	|
|elevation (m, for GDA format)    | -0.1d0	|N/A	|
|Depth weighting power, beta**    | 1.4d0	|Survey dependant|
|Depth weighting constant, Z0**   | 0.0d0	|Survey dependant|


<h3>5.14 Third part of parameter file</h3>

The second part of the parameter file is given in Figure 13.

![figure 13](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/018-Figure_13.png)

The explanation of Figure 13 follows. 

<h3>5.15 MATRIX COMPRESSION parameters </h3>

<ol>
<li>This section includes some of the matrix compression parameters which are self-explanatory and model dependent. Generally, they are not be changed while running TOMOFAST-x. </li>

<li>For general information, threshold distance stands for the distance between source to cell while calculating radius.  It includes the feature of accepting the effects of the surrounding bodies/ rocks/ strata while calculating the Jacobian Matrix. </li>  

<li>Compression rate reduces the size of Jacobian Matrix according to the value which User puts across this parameter. For example, 1 stands for accepting the full Jacobian Matrix during Inversion as shown in the Parfile (Parameter File).User need not to change this parameter while running TOMOFAST-x.</li> 

</ol>

<i>Table 12. MATRIX COMPRESSION parameters section of parfile</i>

|Parameter	|Value for example case	|Range/remark|
|:----				  	|:----- 	|:----|
|******* MATRIX COMPRESSION parameters * |N/A		|N/A|
|distance threshold (source to cell)     |1.d+10 	|Survey dependant|
|compression rate (1.0 = full matrix)    |1.0d0		|Survey dependant|
|Depth weighting constant, Z0**          |0.0d0		|Survey dependant|


<h3>5.16 INVERSION parameters</h3>

<ol>
<li>This section includes the information about the inversion parameters.</li>

<li>In this,User can change the number of iterations for inversion  which has been refer as “number of inversions” in the Parfile (Parameter File).</li>

<li>User can also change the number of LSQR solver iterations. For better and accurate solving of the system of equations in TOMOFAST-x, User should increase the number of solver iterations.</li>

<li>User need not to think about the stopping criterion while performing inversions through TOMOFAST-x.</li>

<li>As TOMOFAST-x operates on LSQR, it has been kept 1 in the respective feature option.</li>

<li>As TOMOFAST-x is not using L1 norm, it has been kept 0 as shown.</li>

</ol>

<i>Table 13. INVERSION parameters section of parfile.</i>

|Parameter		|Value for example case	|Range/remark| 
|:----				  	|:----- |:----|
|******* INVERSION parameters ***	|N/A	|N/A|
|number of inversions                   |50	|Survey dependant|
|number of solver iterations            |100	|Survey dependant|
|stopping criterion                     |1.d-13	|Survey dependant|
|method (LSQR=1)                        |1	|Survey dependant|
|soft threshold ("L1-norm", no=0.)      |0	|Survey dependant|

<h3>5.17 MODEL DAMPING </h3>

<ol>

<li> In this section, user can put model damping coefficient (αm) for both gravity and magnetic models. For more information about (αm), please refer to the section ‘A BRIEF INTRODUCTION TO TOMOFAST-x’ where different types of weighting have been described.</li>

<li>In TOMOFAST-x, we have used L2 norm, so we have kept the power 2 of Lp norm. User need not to change this feature.</li>

</ol>

<i>Table 14. MODEL DAMPING parameters section of parfile.</i>

|Parameter	|Value for example case	|Range/remark| 
|:----				  	|:----- 	|:----|
|******* MODEL DAMPING (m - m_prior) *** |	 	|N/A|
|damping for model1 (grav/ECT)          |2.d-08		|Survey dependant|
|damping for model2 (mag)               |0.d-11		|Survey dependant|
|power p of Lp norm (for LSQR)          |2.0d0		|Survey dependant|
|method (LSQR=1)                        |1		|Survey dependant|
|soft threshold ("L1-norm", no=0.)      |0		|Survey dependant|

<h3>5.18 JOINT INVERSION parameters</h3>

<ol>

<li> In this section, User can change the weights of the problem1 and problem2 which stand for Gravity Inversion and Magnetic Inversion respectively. </li> 

<li> User can put reasonable weights in the respective features (for example – put the weight as 1 across the feature problem1 weight  and 0 across problem2 weight  in order run gravity inversion. </li>

<li> Similarly, User can put 0 across problem1 weight and 1.d-8 (an example that can be used as default value) across  problem2  weight in order to run magnetic inversion.</li>

<li> In order to run joint inversion,User can put 1 across problem1 weight and 1.d-8 across problem2 weight  simultaneously.</li>

<li> Column weight multipliers need not to be changed while running TOMOFAST-x. </li>

<li> User can change the number of iterations for both gravity and magnetic inversions. The numbers provided correspond to the number of separate domain inversions before applying joint inversion. </li>

</ol>

<i>Table 15. JOINT INVERSION parameters section of parfile.</i>

|Parameter	|Value for example case		|Range/remark|
|:----				  	|:----- |:----		|
|******* JOINT INVERSION parameters ***	|N/A	|N/A		|
|problem1 weight                        | 1.d0	|Survey dependant|
|problem2 weight                        | 1.d-8	|acceptable value|
|column weight1 multiplier              | 4.d+3	|keep as is	|
|column weight2 multiplier              | 1.d+0	|keep as is	|
|niter single for model1 (grav)         | 0	|Survey dependant|
|niter single for model2 (mag)          | 0	|Survey dependant|


<h3> 5.19 Damping- gradient constraints </h3>
	
<ol>
	
<li>User can vary the type of weight by putting 1 and 2 for global and local respectively across weight type feature. The local gradient feature is introduced in detail in Giraud et al. (2019b).</li>
	
<li>User can vary damping gradient (α), which is the part of linear gradient regularization, for both gravity and magnetic models. For more information about (α), please refer to the section ‘A BRIEF INTRODUCTION TO TOMOFAST-x’.</li>

</ol>

<i>Table 16. Damping-gradient constraints parameters section of parfile.  </i>

|Parameter	|Value for example case		|Range/remark	| 
|:----				  	|:----- |:----		|
|******* Damping-gradient constraints ** |	|N/A|
|weight type (1-global, 2-local)       | 1	|Survey dependant|
|damping gradient for model1 (grav)    |  1.d-7	|Survey dependant|
|damping gradient for model2 (mag)     |  0.d+4	|Survey dependant|


<h3> 5.20 Cross- gradient constraints </h3>

<ol>
<li>User can put cross- gradient weight (Ws) in order to run the structural term. It can be visualized in the Geophysical Inverse Problem Equation which has been mentioned in the section ‘A BRIEF INTRODUCTION TO TOMOFAST-x’. For more information about the cross-gradient constraint, please refer Gallardo and Meju (2003).</li>

<li>User need not to change ‘number of iterations in methods of weight’ and ‘x-grad derivative’. They have already been set up according to the requirement of TOMOFAST-x.</li>

</ol>


<i>Table 17. Cross-gradient constraints parameters section of parfile.</i>

|Parameter	|Value for example case	|Range/remark|
|:----				  	|:----- |:----		|
|******* Cross-gradient constraints *** |N/A	|N/A|
|cross-gradient weight                  |1.d-4	|Survey dependant|
|num of iterations in method of weights |0	|Survey dependant|
|x-grad deriv (1-fwd, 2-cent, 3-mixed)  |1	|Survey dependant|



<h3>  5.21 Clustering constraints </h3>

<ol>
	
<li>	User can put the clustering weights in order to involve petrophysics information for both gravity and magnetic inversions across the respective features. This weight is referred to as (αp) which can be seen in the petrophysics Term in Geophysical Inverse Problem Equation. For more information, please refer Giraud et al. (2017), (2019c).</li>

<li> User can vary the number of clusters according to the distribution matching the petrophysical measurements in the studied area or prior petrophysical information. </li>

<li> User should put the path of the input cluster file relative to the executable. This file contains the clustering mixture (see description is section of Cluster File. If User is not applying clustering then it is not necessary to provide a file and its path. </li>

<li> User can also put geological clustering weights across the respective feature. These weights can correspond to the probability of observation Pakyuz-Charrier et al. (2018) of the different lithologies in the same fashion as in Giraud et al. (2017), (2019c).User should put the path of the input weight file relative to the executable.  The input weight file organisation is detailed in section ‘Input Geological Weights File’.</li>

</ol>


<i> Table 18. Clustering constraints parameters section of parfile. </i>
	
|Parameter	         |Value for example case |Range/remark| 
|:----				  	|:-----  |:----		|
|******* Clustering constraints ******  | N/A	 |N/A |
|clustering problem1 (grav) weight      | 0.d-7	 |Survey dependant|
|clustering problem2 (mag)  weight      | 0.d-9	 |Survey dependant|
|number of clusters                     |4	 |Survey dependant|
|clustering mixtures                    |mansf_slice_input/clusters.txt	|Survey dependant |
|clustering geol weights per cell       |mansf_slice_input/weights_geol.txt |Survey dependant |
|type of optimization (1-normal,2-log)  |2	 |Survey dependant |
|type of constraints (1-global,2-local) |2	 |Survey dependant |




<h3> 5.22 ADMM constraints  </h3>

Note: full details about this functionality will be provided in a future release. 

<ol>

<li> User can opt using the feature of ADMM (alternating direction method of multipliers). User can enable this feature by 1 putting 1 and 2 for global and local respectively. User can put 0 for not using this feature. Features this part of the parfile parameterises are not yet in their final version in the source-code. </li>

<li>	If User selects local/ global as options mentioned above, User needs to put the location of the gravity/ magnetic bound constraints file across the respective features in this section.</li>
 
</ol>

<i> Table 19. ADMM constraints </i>

|Parameter			|Value for example case	|Range/remark| 
|:----				  	|:----- 	|:----		|
|******* ADMM constraints ************  | N/A		|N/A |
|enable admm? (0-no, 1-glob, 2-local)   | 0		|Survey dependant |
|grav local bound constraints file      | /grav_bound_constraints.txt	|Survey dependant |
|mag local bound constraints file       | NILL	|Survey dependant |
|rho xmin xmax for model1 (grav)        | 11.d-8 -30.d0 330	|Survey dependant |
|rho xmin xmax for model2 (mag)         | 1.d+5 -3.d-3 1.d+10	|Survey dependant |





<h2> 6. INPUT FILES FOR TOMOFAST- x  </h2>


In this section Input for TOMOFAST-x, we would like to introduce the User regarding the inputs which TOMOFAST-x takes in order to perform. In order to run TOMOFAST- x, user needs to have the basic input files. These input files are the files which user needs to provide in the Parfile (parameter file). After entering their respective paths across the respective parameters in the Parfile (Parameter File), they will be passed through TOMOFAST-x. Some example of the input files are shown below in the Figure 15. These input files have been explained in the following sub – sections  and their respective examples have also been shown throughout this manual.

![figure 14](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/019-Figure_14.png)

<h3> 6.1 Types of Basic Input Files  </h3>

In this section, we have tried to explain some of the basic input files for TOMOFAST-x. Please refer to this section for further information about the input files.

<h4> 6.1.1  Data Grid File </h4>

It contains the information about the grids which are being used during gravity and magnetic inversion. The format has been shown in the Figure 15. 

![figure 15](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/020-Figure_15.png)

<ol>

	
<li> The first row represents the number of cells/ data points.  </li>

<li> 	From second row onwards, the four columns represent X-axis data, Y-axis data, Z-axis data, Value_data. The example of data grid file is shown below in the Figure 16. </li>


![figure 16](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/021-Figure_16.png)	


<li> Example input for Data Grid File in Parfile (Parameter File) is shown below in the Figure 17.</li>

![figure 17](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/022-Figure_17.png)


</ol>


<h4>6.1.2 Cluster File  </h4>

<ol>
	
<li> The format of the cluster file is shown below in the Figure18. </li>

![figure 18](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/023-Figure_18.png)

<li> The first row represents the number of clusters. </li>

<li> From second row onwards, the five columns represent cluster weight, mean density contrast, standard deviation, mean magnetic susceptibility, standard deviation and correlation respectively as shown in the Figure 20. </li>

![figure 19](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/024-Figure_19.png)
![figure 20](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/025-Figure_20.png)

</ol>


<h4> 6.1.3  Input Model Voxet File </h4>

<ol>

<li> The first row represents the number of cells </li>

<li> From second row onwards, the columns represent X1_cell, X2_cell, Y1_cell, Y2_cell, Z1_cell, Z2_cell, density    contrast/magnetic susceptibility, value index_in_matrix_1, value index_in_matrix_2, value index_in_matrix_3 and covariance value/ model term weighting (d_m), repectively as shown in figure 6.8 below. Note that TOMOFAST-x uses finite differences, and that X1_cell, X2_cell, Y1_cell, Y2_cell, Z1_cell, Z2_cell represent to coordinates of the 6 faces of a right rectangular prism.  </li>

![figure 21](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/026-Figure_21.png)


</ol>


<h4> 6.1.4 Input Geological Weights File</h4>

<ol>

<li> The first row represents the dimension of the model i.e. number of cell  and number of lithology respectively.</li>

<li> From second row onwards, the four columns represent geological weights for 1st lithology, 2nd lithology, 3rd lithology and 4th lithology respectively as shown in the Figure 22. </li>

![figure 22](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/027-Figure_22.png)

<li> Example for input file of Geological Weight in the Parfile (Parameter File) is shown below in the Figure23.< /li>

![figure 23](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/028-Figure_23.png)

</ol>


<h3> 6.2 Terminal Input </h3>

<ol>

<li> In order to run TOMOFAST-x, user needs a terminal input using command line under linux operating systems or using the Windows Linux Subsystem. After following the sections of Basic Requirements and Installation, user need to provide required inputs to the Parfile (Parameter Files). After that, user needs to change the directory to the folder where the executable tomofast3D exists and run the same as shown below in the Figure24 and Figure25. </li>

![figure 24](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/029-Figure_24.png)

<li> Input Command Line - After changing directory in the terminal, run the following Input Command Line as shown below in the Figure 25. </li>

![figure 25](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/030-Figure_25.png)


<li> Explanation of the Input Command Line is shown below in the Figure 26.
<b> mpirun -n 1 ./tomofast3D -j ./Parfile_mansf_slice.txt | tee out.txt</b> </li>

![figure 26](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/031-Figure_26.png)

</ol>

<h2>7. OUTPUT FOR TOMOFAST-x </h2>

In this section, we have tried to describe about the output for TOMOFAST-x.User can get the output of the respective inversion at the respective location present in the Parfile (Parameter File) as shown in the Figure 27.

![fig27](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/032-Figure_27.png)

<h3> 7.1 Types of Output Files </h3>

<ol>

<li> In this section, user can have a look on the Samples of Output Files of the TOMOFAST-x as shown below in the Figure 28Figure 29. </li>
<li> User needs to follow the Voxet folder which appears in the output folder in order to get these output voxet files.</li>

</ol>


![fig 28](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/033-Figure_28.png)

<h3> 7.2 Sample Output Files </h3>	

In this section, a group of sample output files have been explained.	
	

<h4> 7.2.1 Clustering Data Output File </h4>

<ol>

<li> The first row represents the dimension of the model i.e. number of cell.</li>

<li> There are 16 columns present from the second row onwards. The first 6 columns represent X1_data, X2_data, Y1_data, Y2_data, Z1_data, Z2_data. For a particular cell/ row, X1 and X2 represent the limits of the x-axis data in which the respective values for mixture model have been obtained. Similarly, Y1 and Y2 along with Z1 and Z2 represent the limits across y-axis and z-axis respectively. </li>

<li> The columns 7<sup>th</sup>, 8<sup>th</sup> and 9<sup>th</sup> represent indices across x, y and z axis respectively. </li>

<li> The column 10<sup>th</sup> represents the calculated values for mixture model which have been used for petrophysical constraints.</li>

<li> The column 11<sup>th</sup> represents the values of the first derivative of the mixture model with respect to the property inverted for. Each of the following columns correspond to the value of the separate Gaussians making up the Gaussian mixture model used in the petrophysical constraints. They are given in the same order as they are defined in the cluster file. </li>

<li> <b> User can use this file for visualization only after deleting the first row of the file so that the file can acquire proper matrix format.</b> </li>

<li> Example output file of cluster data is shown below in the Figure 29. </li>

</ol>

![fig 29](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/034-Figure_29.png)


<h4> 7.2.2 Clustering Voxet Output File	</h4>

<ol>

<li> The first row represents the number of cells. From second row onwards, first 6 columns represent X1_data, X2_data, Y1_data, Y2_data, Z1_data, Z2_data. For a particular cell/ row, X1 and X2 represent the limits of the x-axis data in which the respective values for mixture model have been obtained. Similarly, Y1 and Y2 along with Z1 and Z2 represent the limits across y-axis and z-axis respectively. </li>
	 
<li> The 7<sup>th</sup> column represents the calculated values for mixture model which have been used for petrophysical constraints same as the column 10th of the clustering output data file. </li>

<li> The columns <sup>8th</sup>, 9<sup>th</sup> and 10<sup>th</sup> represent indices of the x- axis, y- axis and z-axis respectively.</li>

<li> <b> User can use this file for visualization only after deleting the first row of the file so that the file can acquire proper matrix format. Data can then be loaded using simple MATLAB/ Python codes. </b> </li>

<li> Example for cluster voxet output file of the TOMOFAST-x is shown below in the Figure 30.</li>

</ol>

![fig 30](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/035-Figure_30.png)

<h4> 7.2.3 Model Voxet Output File </h4>

<ol>

<li> The first row represents the number of cells. From second row onwards, first 6 columns represent X1_data, X2_data, Y1_data, Y2_data, Z1_data, Z2_data. For a particular cell/ row, X1 and X2 represent the limits of the x-axis data in which the respective values for mixture model have been obtained. Similarly, Y1 and Y2 along with Z1 and Z2 represent the limits across y-axis and z-axis respectively.</li>

<li> The column 7<sup>th</sup> represents density contrast/ magnetic susceptibility contrast values according to the respective gravity/ magnetic voxet output file.</li>

<li> The columns 8<sup>th</sup>, 9<sup>th</sup> and 10<sup>th</sup> represent indices of the x- axis, y- axis and z-axis respectively.</li>

<li><b> User can use this file for visualization only after deleting the first row of the file so that the file can acquire proper matrix format. Data can then be loaded using simple Matlab of Python codes. </b></li>

<li> Example gravity voxet output file of TOMOFAST-x is shown below in Figure 31.</li>

</ol>

![fig 31](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/036-Figure_31.png)

<h4> 7.2.4 Cost Output File </h4>

<ol>

<li> User can find the Cost Output File in the Output folder as shown below in the Figure 32. </li>

![fig 32](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/037-Figure_32.png)

<li> In this file, there are 8 columns present which are representing iteration number, gravity data, magnetic data, gravity model value, magnetic model value, gravity structure, magnetic structure, cross gradient value and petrophysical data respectively. </li>

<li><b> User can use this file for visualization only after deleting the last row of the file so that the file can acquire proper matrix format. Data can then be loaded using simple Matlab of Python codes.</b></li>

<li> Example of cost output file of TOMOFAST-x is shown below in the Figure 33.</li>

</ol>

![fig 33](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/038-Figure_33.png)



<h3> 7.3 Command screen output of TOMOFAST-x </h3>

<p>Command screen output is the output of the working of the TOMOFAST-x in the form of log file which can be referred for the visualization of running TOMOFAST-x as shown in Figure 34. It will be stored in the file namely <b>out.txt</b> in the output folder as specified in the command line example provided above. For more information regarding this file, please refer to ‘Terminal output of TOMOFAST-x’. </p>

![fig 34](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/039-Figure_34.png)

<h3> 7.4 Terminal output of TOMOFAST-x </h3>

<p> After running the input command line,User receives a text file as out.txt in the directory folder. This out.txt is the command log file which can be referred while debugging the inversion procedure. While running inversions, if TOMOFAST-x stops in between, then this out.txt which is a terminal output file can be referred to find the errors and to monitor the inversion while it is running. This file is also referred as command screen output. Example of terminal output is shown below in the respective figures in this section. </p>

![fig 35](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/040-Figure_35.png)

![fig 36](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/041-Figure_36.png)

![fig 37](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/042-Figure_37.png)

![fig 38](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/043-Figure_38.png)

<h3> 7.5 Working of the TOMOFAST-x </h3>

![fig 39](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/044-Figure_39.png)

Explanation of the Figure 39.

<ul>
	
<li> <b>ROW 135</b> – Precision has been kept DOUBLE (A type of floating point number).</li>
<li> <b>ROW 136</b> – allocates the rank where myrank is the rank of the processor which is writing the data. <br>Nbproc represents the number of processors which have been used for inversion.</li>
<li> <b>ROW 137</b> – statement for solving problem for gravity, magnetic and joint inversion.</li>
<li> <b>ROW 138</b> – statement for running the cross- gradient and clustering procedure respectively.<br> Here, <b>T</b> and <b>F</b> refers that you are applying AND not applying respectively. </li>
<li> <b>ROW 139</b> and <b>ROW 140</b> – show the allocation of the rank, nunber of data and nunber of elements for 2 cases respectively.</li>
<li> <b>ROW 141</b> – shows the statement for the allocation of the model </li>
<li> <b>ROW 142, 143, 144</b> and <b>145</b> – show that the allocation of model and grid for both cases has been completed respectively.</li>
<li> <b>ROW 146</b> and <b>149</b> – reading model from the input files respectively.</li>
<li> <b>ROW 147, 148, 150</b> and <b>151</b> – allocating the limits of the x-axis and y-axis i.e. the minimum and the maximum values of the limits respectively.</li>
<li> <b>ROW 152</b> and <b>153</b> – show the models that has been written in the designated folder in the form of text files </li>
<li> <b>ROW 154</b> – statement for data allocation </li>
<li> <b>ROW 155</b> and <b>156</b> – number of data allocated for both gravity and magnetic respectively </li>
<li> <b>ROW 157</b> and <b>158</b> – reading data from the respective output data grid files </li>
<li> <b>ROW 159</b> – statement for the sensitivity matrix allocation </li>
<li> <b>ROWS (160–166)</b> - show the allocation of inversion arrays, residuals, column weights, damping weights, sensitivity. <br>This information is useful to monitor the inversion. It writes the information after these different steps have been completed. It can be useful when debugging or when there is a problem to know when the inversion stops/crashes.<br> For visualization, these options can be seen in the Parfile (Parameter file) i.e. in the sections of JOINT INVERSION parameter and Model Damping respectively.</li>
<li> <b>ROW 167</b> – shows the allocation of inversion arrays i.e. for the case of Gravity Inversion.</li>
<li> <b>ROWS (168-174)</b> - show the allocation of inversion arrays, residuals, column weights, damping weights, sensitivity. For visualization, these options can be seen in the Parfile (Parameter file) i.e. in the sections of JOINT INVERSION parameter and Model Damping respectively. </li>

</ul>

![fig 40](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/045-Figure_40.png)

Explanation of Figure 40.

<ul>

<li> <b>ROW 175</b> – shows the allocation of inversion arrays i.e. for the case of Magnetic Inversion.</li>
<li> <b>ROWS 176 -179</b> – calculation of the sensitivity kernel  </li>
<li> <b>ROW 180 – 186</b> – statements showing that the respective data file have been read from the respective files and written to the respective files which are self - explanatory.</li>
<li> <b>ROW 187 – 203</b> – statements represents the reading of different models and their respective Xmin, Xmax, Ymin and Ymax. The remaining statements are self – explanatory.</li>
<li> <b>ROW 206</b> and <b>207</b> – represent the cost of the gravity and magnetic models respectively.</li> 

</ul>

![fig 41](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/046-Figure_41.png)

Explanation of Figure 41.

<ul>
	
<li> <b>ROW 209</b> and <b>258</b> – show the iteration number which is as same as inversion number </li>
<li> <b>ROWS 211</b> and <b>224</b> – represent the weights of the first and second problem i.e. for gravity inversion and magnetic inversion respectively </li>
<li> <b>ROWS 212</b> and <b>213</b> – represent misfit term cost and value of <b>nel</b> </li>
<li> <b>ROWS 214</b> and <b>227</b> – represent addition of model damping and ……… for the respective models </li>
<li> <b>ROWS 215</b> and <b>228</b> – represent the cost of the damping term</li>
<li> ROWS 217 and 230 – represent the value of damping gradient for the respective models. The type of weight i.e. whether it is global or local. As the type of weight is global over here, it has been mentioned 1. </li>
<li> <b>ROWS 218 – 221</b> and <b>231 – 234</b> – represent damping gradient cost term and their respective total cost.User need not to look after the same in order to run TOMOFAST-x.</li>
<li> <b>ROW 237</b> – represents the cross – gradient function and the type of derivative which has been set up in the parameter file. Here, 1 represents forward derivative. For more information, please refer the section of Cross gradients constraints in the Parfile (Parameter File).</li>
<li> <b>ROW 238</b> – represents cross – gradient cost which will be stored in the 8<sup>th</sup> column of the cost file if no clustering constraints i.e. petrophysics are used.</li>
<li> <b>ROWS 242 – 251</b> – represent the number of iteration, residual and gradient values in LSQR respectively. Here, each row represents these features after every 10 iterations. Here, it has been set up 100 number of solver iteration per inversion as the User can see in the section of INVERSION parameter in the Parfile (Parameter File). </li>
<li> <b>ROW 252</b> – represents the residual value, <b>r</b> which has come at the end of 100 iterations as the User can compare the same from ROW 251.</li>
<li> <b>ROW 253</b> – represents the value of gravity data cost which can be seen across the parameter ‘cost’ in the row. This value will be stored in the 2<sup>nd</sup> column of the cost file.</li>
<li> <b>ROW 254</b> – represents the value of magnetic data cost which can be seen across the parameter ‘cost’ in the row. This value will be stored in the 3<sup>rd</sup> column of the cost file.</li>
<li> <b>ROW 255</b> – represents the value of gravity model cost which can be seen across the parameter ‘model cost’ in the row. This value will be stored in the 4<sup>th</sup> column of the cost file.</li>
<li> <b>ROW 256</b> – represents the value of magnetic model cost which can be seen across the parameter ‘model cost’ in the row. This value will be stored in the 5<sup>th</sup> column. </li> 
	
</ul>

\** This procedure continues for every iteration for Inversion. Here, just for information, User can change the number of inversions according to their requirement in the INVERSION parameter section of the Parfile (Parameter File)

![fig 42](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/047-Figure_42.png)

Explanation of Figure 42.
<ul>

<li> <b>ROW 2658</b> – represents minimum and maximum for the gravity model values.  </li>
<li> <b>ROW 2659</b> – represents minimum and maximum for the magnetic model values. </li>
<li> <b>ROW 2661</b> and <b>2662</b> – represent the path of the final voxets of gravity and magnetic models respectively and the statements confirm that both have been written. </li>
<li> <b>ROWS 2663</b> and <b>2664</b> – represent the paths of final calculated gravity and magnetic data respectively and confirm that both have been written. </li>
<li> <b>ROWS 2665 – 2667</b> – represent the paths of final voxets for cross gradient, gravity sensitivity and magnetic sensitivity respectively </li>
<li> <b>ROW 2669</b> – represents total CPU accumulated times </li>
<li> <b>ROW 2670</b> – represents TOTAL CPU minimum time. Here total represents the total number of CPU used while running TOMOFAST-x.User can set the number of CPU in the input command line while invoking TOMOFAST-x. </li>
<li> <b>ROW 2671</b> – represents TOTAL CPU maximum time. </li> 

</ul>

<h2> 8. HOW TO RUN GEOPHYSICAL INVERSIONS IN TOMOFAST-x </h2>

<p>This section will briefly introduce about how to run the respective inversions either it is gravity, magnetic or joint inversion.</p>

<h3> For running Gravity Inversion</h3>

<p>User can put reasonable weights in the respective features in the section of JOINT INVERSION parameters, (for example – put the weight as 1 across the feature <b>problem1 weight</b>  and 0 across <b>problem2 weight</b> in order run gravity inversion.) </p>

<h3>For running Magnetic Inversion</h3>

<p>User can put 0 across <b>problem1</b> weight and 1.d-8 (an example that can be used as default value) across <b>problem2 weight</b> in order to run magnetic inversion in the section of JOINT INVERSION parameters in the Parameter File.</p>

<h3>For running Joint Inversion</h3>

<p>In order to run joint inversion, User can put 1 across <b>problem1</b> weight and 1.d-8 across <b>problem2</b> weight simultaneously across the respective features in the section of JOINT INVERSION parameters.</p>

<h3>Types of combinations during Inversion</h3>

<p>This section involves the introduction of some of the options available with TOMOFAST-x which can be utilized in order to test the respective inversion.</p>

<ul>

<li> Gravity </li> 
<li> Magnetic </li>
<li> Cross gradient </li>
<li> Petrophysical constraint</li>
<li> Local gradient regularization </li>
<li> Geological uncertainty (Wh)</li>
<li> Geological Weights </li>
<li> With/ without starting and prior model And their respective combinations</li>
	
</ul>



<h3>Simple Examples</h3> 

![fig 43](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/048-Figure_43.png)

![fig 44](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/049-Figure_44.png)

![fig 45](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/050-Figure_45.png)

![fig 46](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/051-Figure_46.png)

![fig 47](https://github.com/TOMOFAST/TOMOFASTx/blob/master/Docs/figures/052-Figure_47.png)



