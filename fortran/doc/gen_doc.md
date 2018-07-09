# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Stochastic Fortran implementation #

## About ##

(c) 2013-2018 Lesley De Cruz and Jonathan Demaeyer

See [LICENSE.txt](../../LICENSE.txt) for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
[doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

for the MAOOAM original code, and with

* Demaeyer, J. and Vannitsem, S.: Comparison of stochastic parameterizations in
the framework of a coupled ocean-atmosphere model, Nonlin. Processes Geophys.
Discuss., 2018.

for the stochastic part.

**Please cite both articles if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <lesley.decruz@meteo.be>, <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the MAOOAM [code repository](http://www.github.com/Climdyn/MAOOAM)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

A pdf version of this manual is available [here](../latex/Reference_manual.pdf).

------------------------------------------------------------------------

## Installation ##

The program can be installed with Makefile. We provide configuration files for 
two compilers : gfortran and ifort.

By default, gfortran is selected. To select one or the other, simply modify the 
Makefile accordingly or pass the COMPILER flag to `make`. If gfortran is
selected, the code should be compiled with gfortran 4.7+ (allows for
allocatable arrays in namelists).  If ifort is selected, the code has been
tested with the version 14.0.2 and we do not guarantee compatibility with older
compiler version.

To install, unpack the archive in a folder or clone with git:

```bash     
git clone https://github.com/Climdyn/MAOOAM.git
cd MAOOAM
```
     
and run:

```bash     
make
```     
By default, the inner products of the basis functions, used to compute the
coefficients of the ODEs, are not stored in memory. If you want to enable the
storage in memory of these inner products, run make with the following flag:

```bash
make RES=store
```

Depending on the chosen resolution, storing the inner products may result in a
huge memory usage and is not recommended unless you need them for a specific
purpose.
 
 Remark: The command "make clean" removes the compiled files.

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called aotensor which
includes all the coefficients. This tensor is computed once at the program
initialization.

The following files are part of the MAOOAM model alone:

* maooam.f90 : Main program.
* aotensor_def.f90 : Tensor aotensor computation module.
* IC_def.f90 : A module which loads the user specified initial condition.
* inprod_analytic.f90 : Inner products computation module.
* rk2_integrator.f90 : A module which contains the Heun integrator for the model equations.
* rk4_integrator.f90 : A module which contains the RK4 integrator for the model equations.
* Makefile : The Makefile.
* params.f90 : The model parameters module.
* tl_ad_tensor.f90 : Tangent Linear (TL) and Adjoint (AD) model tensors definition module
* rk2_tl_ad_integrator.f90 : Heun Tangent Linear (TL) and Adjoint (AD) model integrators module
* rk4_tl_ad_integrator.f90 : RK4 Tangent Linear (TL) and Adjoint (AD) model integrators module
* test_tl_ad.f90 : Tests for the Tangent Linear (TL) and Adjoint (AD) model versions
* README.md : A read me file.
* LICENSE.txt : The license text of the program.
* util.f90 : A module with various useful functions.
* tensor.f90 : Tensor utility module.
* stat.f90 : A module for statistic accumulation.
* params.nml : A namelist to specify the model parameters.
* int_params.nml : A namelist to specify the integration parameters.
* modeselection.nml : A namelist to specify which spectral decomposition will be used.

with the addition of the files:

* maooam_stoch.f90 : Stochastic implementation of MAOOAM.
* maooam_MTV.f90 : Main program - MTV implementation for MAOOAM. 
* maooam_WL.f90 : Main program - WL implementation for MAOOAM. 
* corrmod.f90 : Unresolved variables correlation matrix initialization module.
* corr_tensor.f90 : Correlations and derivatives for the memory term of the WL parameterization.
* dec_tensor.f90 : Tensor resolved-unresolved components decomposition module.
* int_comp.f90 : Utility module containing the routines to perform the integration of functions.
* int_corr.f90 : Module to compute or load the integrals of the correlation matrices.
* MAR.f90 : Multidimensional AutoRegressive (MAR) module to generate the correlation for the WL parameterization.
* memory.f90 : WL parameterization memory term \f$M_3\f$ computation module. 
* MTV_int_tensor.f90 : MTV tensors computation module.
* MTV_sigma_tensor.f90 : MTV noise sigma matrices computation module.
* WL_tensor.f90 : WL tensors computation module.
* rk2_stoch_integrator.f90 : Stochastic RK2 integration routines module.
* rk2_ss_integrator.f90 : Stochastic uncoupled resolved nonlinear and tangent linear RK2 dynamics integration module.
* rk2_MTV_integrator.f90 :  MTV RK2 integration routines module.
* rk2_WL_integrator.f90 :  WL RK2 integration routines module.
* sf_def.f90 : Module to select the resolved-unresolved components.
* SF.nml : A namelist to select the resolved-unresolved components.
* sqrt_mod.f90 : Utility module with various routine to compute matrix square root.
* stoch_mod.f90 : Utility module containing the stochastic related routines.
* stoch_params.f90 : Stochastic models parameters module.
* stoch_params.nml : A namelist to specify the stochastic models parameters.

which belong specifically to the stochastic implementation.

  
------------------------------------------------------------------------

## MAOOAM Usage ##

The user first has to fill the params.nml and int_params.nml namelist files according to their needs.
Indeed, model and integration parameters can be specified respectively in the params.nml and int_params.nml namelist files. Some examples related to already published article are available in the params folder.

The modeselection.nml namelist can then be filled : 
* NBOC and NBATM specify the number of blocks that will be used in respectively the ocean and
  the atmosphere. Each block corresponds to a given x and y wavenumber.
* The OMS and AMS arrays are integer arrays which specify which wavenumbers of
  the spectral decomposition will be used in respectively the ocean and the
  atmosphere. Their shapes are OMS(NBOC,2) and AMS(NBATM,2).
* The first dimension specifies the number attributed by the user to the block and the second
  dimension specifies the x and the y wavenumbers.
* The VDDG model, described in Vannitsem et al. (2015) is given as an example
  in the archive.
* Note that the variables of the model are numbered according to the chosen
  order of the blocks.

The Makefile allows to change the integrator being used for the time evolution.
The user should modify it according to its need.
By default a RK2 scheme is selected.

Finally, the IC.nml file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current IC.nml file (if it exists)
and run the program :

    ./maooam

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart :

    ./maooam

It will generate two files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)

The tangent linear and adjoint models of MAOOAM are provided in the
tl_ad_tensor, rk2_tl_ad_integrator and rk4_tl_ad_integrator modules. It is
documented [here](tl_ad_doc.md).

------------------------------------------------------------------------

## Stochastic code usage ##

The user first has to fill the MAOOAM model namelist files according to their needs (see the previous section).
Additional namelist files for the fine tuning of the parameterization must then be filled,
and some "definition" files (with the extension .def) must be provided. An example is provided with the code.

Full details over the parameterization options and definition files are available [here](sto_doc.md).

The program "maooam_stoch" will generate the evolution of the full stochastic dynamics with the command:

    ./maooam_stoch

or any other dynamics if specified as an argument (see the header of [maooam_stoch.f90](./maooam_stoch.f90)).
It will generate two files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)

The program "maooam_MTV" will generate the evolution of the MTV parameterization evolution, with the command:

    ./maooam_MTV

It will generate three files :
 * evol_MTV.dat : the recorded time evolution of the variables.
 * ptend_MTV.dat : the recorded time evolution of the tendencies (used for debugging).
 * mean_field_MTV.dat : the mean field (the climatology)

The program "maooam_WL" will generate the evolution of the MTV parameterization evolution, with the command:

    ./maooam_WL

It will generate three files :
 * evol_WL.dat : the recorded time evolution of the variables.
 * ptend_WL.dat : the recorded time evolution of the tendencies (used for debugging).
 * mean_field_WL.dat : the mean field (the climatology)

------------------------------------------------------------------------

## MAOOAM Implementation notes ##

As the system of differential equations is at most bilinear in \f$z_j\f$ (\f$j=1..n\f$), 
\f$\boldsymbol{z}\f$ being the array of variables, it can be expressed as a tensor
 contraction :

\f[  \frac{d z_i}{dt} = \sum_{j,k=0}^{ndim} \, \mathcal{T}_{i,j,k} \, z_k \; z_j \f]

with \f$z_0 = 1\f$.


The tensor aotensor_def::aotensor is the tensor \f$\mathcal{T}\f$ that encodes the differential equations is composed so that:

* \f$\mathcal{T}_{i,j,k}\f$ contains the contribution of \f$dz_i/dt\f$ proportional to \f$ z_j \, z_k\f$.
* Furthermore, \f$z_0\f$ is always equal to 1, so that \f$\mathcal{T}_{i,0,0}\f$ is the constant
contribution to \f$dz_i/dt\f$
* \f$\mathcal{T}_{i,j,0} + \mathcal{T}_{i,0,j}\f$ is the contribution to  \f$dz_i/dt\f$ which is linear in
\f$z_j\f$.

Ideally, the tensor aotensor_def::aotensor is composed as an upper triangular matrix 
(in the last two coordinates).

The tensor for this model is composed in the aotensor_def module and uses the
inner products defined in the inprod_analytic module.

------------------------------------------------------------------------

## Stochastic code implementation notes ##

A stochastic version of MAOOAM and two stochastic parameterization methods
(MTV and WL) are provided with this code.

The stochastic version of MAOOAM is given by

\f[  \frac{d \boldsymbol{z}}{dt} = f(\boldsymbol{z})+ \boldsymbol{q} \cdot \boldsymbol{dW} (t)\f]

where \f$\boldsymbol{dW}\f$ is a vector of standard Gaussian White noise and where several choice for \f$f(\boldsymbol{z})\f$ are available.
For instance, the default choice is to use the full dynamics: 
\f[ f(\boldsymbol{z}) = \sum_{j,k=0}^{ndim} \, \mathcal{T}_{i,j,k} \, z_k \; z_j . \f]
The implementation uses the tensorial framework described above and add some noise to it.
This stochastic version is further detailed [here](sto_doc.md).

The MTV parameterization for MAOOAM is given by

\f[ \frac{d\boldsymbol{x}}{dt} = F_x(\boldsymbol{x}) + \frac{1}{\delta} R(\boldsymbol{x}) + G(\boldsymbol{x}) + \sqrt{2} \,\, \boldsymbol{\sigma}(\boldsymbol{x}) \cdot \boldsymbol{dW} \f]

where \f$\boldsymbol{x}\f$ is the set of resolved variables and \f$\boldsymbol{dW}\f$ is a vector of standard Gaussian White noise.
\f$F_x\f$ is the set of tendencies of resolved system alone and \f$\delta\f$ is the timescale separation parameter.

The WL parameterizations for MAOOAM is given by

\f[\frac{d\boldsymbol{x}}{dt} = F_x(\boldsymbol{x}) + \varepsilon \, M_1(\boldsymbol{x}) + \varepsilon^2 \, M_2(\boldsymbol{x},t) + \varepsilon^2 \, M_3 (\boldsymbol{x},t)\f]

where \f$\varepsilon\f$ is the resolved-unresolved components coupling strenght and where the different terms \f$M_i\f$ account for different effect. 

The implementation for these two approaches uses the tensorial framework described above, with the addition of new tensors to account for the
terms \f$R,G,\sigma,M_1,M_2,M_3\f$. They are detailed more completely [here](sto_doc.md).



------------------------------------------------------------------------

## Final Remarks ##

The authors would like to thank Kris for help with the lua2fortran project. It
has greatly reduced the amount of (error-prone) work.

  No animals were harmed during the coding process.
