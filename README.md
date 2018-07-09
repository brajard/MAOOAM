# Modular arbitrary-order ocean-atmosphere model: MAOOAM stochastic parameterization branch

--------------------------------------------------------------------------------------------

## About ##

(c) 2013-2017 Lesley De Cruz and Jonathan Demaeyer

See [LICENSE.txt](LICENSE.txt) for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
[doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

for the MAOOAM original code, and with

* Demaeyer, J. and Vannitsem, S.: Comparison of stochastic parameterizations
in the framework of a coupled ocean-atmosphere model, Nonlin. Processes Geophys.
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

------------------------------------------------------------------------

## Model description ##

The atmospheric component of the model is based on the papers of Charney and
Straus (1980), Reinhold and Pierrehumbert (1982) and  Cehelsky and Tung (1987),
all published in the Journal of Atmospheric Sciences. The ocean component is
based on the papers of Pierini (2012), Barsugli and Battisti (1998). The
coupling between the two components includes wind forcings, radiative and heat
exchanges.

![](./misc/figs/movie_10x10_GMD.gif)

------------------------------------------------------------------------

## Branch description ##

This branch is a version of the MAOOAM model focused on stochastic methods
and more particularly stochastic parameterizations. 

------------------------------------------------------------------------

## Parameterizations description ##

For the moment, this code provide the implementation of two different
stochastic parameterizations for the MAOOAM ocean-atmosphere coupled
model.

The first one is based on homogenization and follow closely the 
approach proposed in Franzke et al. 2005. It is called hereafter "MTV", from the
names of the authors of the paper that first proposed to apply the method 
to climate systems (Majda et al., 2001).

The second one is based on the Ruelle response theory, and is called
hereafter "WL", from the names of the authors of the paper that proposed
to apply response theory to parameterization problems (Wouters and Lucarini, 2012).

Complete details about these parameterizations can be found in the article
to whom this code is attached, and in the references proposed below.

A short account is also proposed in the documentation.

------------------------------------------------------------------------

## Implementation notes ##

At the present moment, the parameterizations have been solely implemented
for the Fortran version of MAOOAM. Please read the README.md file in 
the fortran folder for the implementation notes. 

------------------------------------------------------------------------

## References ##

* Vannitsem, S., Demaeyer, J., De Cruz, L., and Ghil, M.: Low-frequency
variability and heat transport in a loworder nonlinear coupled ocean-atmosphere
model, Physica D: Nonlinear Phenomena, 309, 71-85, 2015. 

* De Cruz, L., Demaeyer, J., & Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0,
Geoscientific Model Development, 9(8), 2793-2808, 2016. 

* Majda, A. J., Timofeyev, I., & Vanden Eijnden, E.: A mathematical framework
for stochastic climate models,
Communications on Pure and Applied Mathematics, 54(8), 891-974, 2001.

* Franzke, C., Majda, A. J., & Vanden-Eijnden, E.: Low-order stochastic mode
reduction for a realistic barotropic model climate,
Journal of the atmospheric sciences, 62(6), 1722-1745, 2005.  

* Wouters, J., & Lucarini, V.: Disentangling multi-level systems: averaging,
correlations and memory.
Journal of Statistical Mechanics: Theory and Experiment, 2012(03), P03003, 2012.

* Demaeyer, J., & Vannitsem, S.: Stochastic parametrization of subgrid‐scale
processes in coupled ocean–atmosphere systems: benefits and limitations of response theory,
Quarterly Journal of the Royal Meteorological Society, 143(703), 881-896, 2017. 

Please see the main articles for the full list of references.
