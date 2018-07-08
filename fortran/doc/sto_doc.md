#  Modular arbitrary-order ocean-atmosphere model: The MTV and WL parameterizations #

## The stochastic version of MAOOAM ##

The stochastic version of MAOOAM is given by

\f[  \frac{d \boldsymbol{z}}{dt} = f(\boldsymbol{z})+ \boldsymbol{q} \cdot \boldsymbol{dW} (t)\f]

where \f$\boldsymbol{dW}\f$ is a vector of standard Gaussian White noise and where several choice for \f$f(\boldsymbol{z})\f$ are available.
For instance, the default choice is to use the full dynamics: 
\f[ f(\boldsymbol{z}) = \sum_{j,k=0}^{ndim} \, \mathcal{T}_{i,j,k} \, z_k \; z_j . \f]
The implementation uses thus the tensorial framework of MAOOAM and add some noise to it. To study parameterization methods in MAOOAM, the models variables \f$\boldsymbol z\f$ is divised in two components: the resolved component \f$\boldsymbol x\f$ and the unresolved component \f$\boldsymbol y\f$ (see below for more details).

 Since MAOOAM is a ocean-atmosphere model, it can be decomposed further into oceanic and atmospheric components:
\f[ \boldsymbol z = \{ \boldsymbol{x}_{\text{a}}, \boldsymbol{x}_{\text{o}}, \boldsymbol{y}_{\text{a}}, \boldsymbol{y}_{\text{o}}\} \f]
and in the present implementation, the noise amplitude can be set in each component:
\f[  \frac{d \boldsymbol{x}_{\text{a}}}{dt} = f_{x,\text{a}}(\boldsymbol{z})+ \boldsymbol{q}_{x,\text{a}} \cdot \boldsymbol{dW}_{x,\text{a}}  (t)\f]
\f[  \frac{d \boldsymbol{x}_{\text{o}}}{dt} = f_{x,\text{o}}(\boldsymbol{z})+ \boldsymbol{q}_{x,\text{o}} \cdot \boldsymbol{dW}_{x,\text{o}}  (t)\f]
\f[  \frac{d \boldsymbol{y}_{\text{a}}}{dt} = f_{y,\text{a}}(\boldsymbol{z})+ \boldsymbol{q}_{y,\text{a}} \cdot \boldsymbol{dW}_{y,\text{a}}  (t)\f]
\f[  \frac{d \boldsymbol{y}_{\text{o}}}{dt} = f_{y,\text{o}}(\boldsymbol{z})+ \boldsymbol{q}_{y,\text{o}} \cdot \boldsymbol{dW}_{y,\text{o}}  (t)\f]
through the parameters stoch_params::q_ar, stoch_params::q_au, stoch_params::q_or and stoch_params::q_ou.

------------------------------------------------------------------------

## The resolved-unresolved components ##

Due to the decomposition into resoved variables \f$ \boldsymbol x\f$ and unresolved variables \f$ \boldsymbol y\f$, the equation of the MAOOAM model can be rewritten:

\f[ \frac{d \boldsymbol x}{dt} = \boldsymbol{H}^x + \boldsymbol{L}^{xx}\cdot\boldsymbol{x} + \boldsymbol{L}^{xy}\cdot\boldsymbol{y} + \boldsymbol{B}^{xxx} : \boldsymbol{x} \otimes \boldsymbol{x} + \boldsymbol{B}^{xxy} : \boldsymbol{x} \otimes \boldsymbol{y} + \boldsymbol{B}^{xyy} : \boldsymbol{y} \otimes \boldsymbol{y} + \boldsymbol{q}_{x} \cdot \boldsymbol{dW}_x \f]
\f[ \frac{d \boldsymbol y}{dt} = \boldsymbol{H}^y + \boldsymbol{L}^{yx}\cdot\boldsymbol{x} + \boldsymbol{L}^{yy}\cdot\boldsymbol{y} + \boldsymbol{B}^{yxx} : \boldsymbol{x} \otimes \boldsymbol{x} + \boldsymbol{B}^{yxy} : \boldsymbol{x} \otimes \boldsymbol{y} + \boldsymbol{B}^{yyy} : \boldsymbol{y} \otimes \boldsymbol{y} + \boldsymbol{q}_{y} \cdot \boldsymbol{dW}_y \f]

where \f$\boldsymbol{q}_x= \{\boldsymbol{q}_{x,\text{a}},\boldsymbol{q}_{x,\text{o}}\}\f$ and \f$\boldsymbol{q}_y= \{\boldsymbol{q}_{y,\text{a}},\boldsymbol{q}_{y,\text{o}}\}\f$. We have thus also \f$\boldsymbol{dW}_x=\{\boldsymbol{dW}_{x,\text{a}},\boldsymbol{dW}_{x,\text{o}}\}\f$ and \f$\boldsymbol{dW}_y= \{\boldsymbol{dW}_{y,\text{a}},\boldsymbol{dW}_{y,\text{o}}\}\f$. The various terms of the equations above are accessible in the dec_tensor module. To specify which variables belong to the resolved (unresolved) component, the user must fill the SF.nml namelist file by setting the component of the vector sf_def::sf to 0 (1). This file must be filled before starting any of the stochastic and parameterization codes. If this file is not present, launch one of the programs. It will generate a new SF.nml file and then abort.

The purpose of the parameterization is to reduce the \f$\boldsymbol x\f$ equation by closing it while keeping the statisical properies of the full system. To apply the parameterizations proposed in this implementation, we consider a modified version of the equation above:

\f[ \frac{d \boldsymbol x}{dt} = F_x(\boldsymbol{x}) + \boldsymbol{q}_{x} \cdot \boldsymbol{dW}_x + \frac{\varepsilon}{\delta} \, \Psi_x(\boldsymbol{x},\boldsymbol{y}) \f]
\f[ \frac{d \boldsymbol y}{dt} = \frac{1}{\delta^2}\, \Big( F_y(\boldsymbol{y}) + \delta \, \boldsymbol{q}_{y} \cdot \boldsymbol{dW}_y \Big) + \frac{\varepsilon}{\delta} \, \Psi_y(\boldsymbol{x},\boldsymbol{y}) \f]

where \f$\varepsilon\f$ is the resolved-unresolved components coupling strength given by the parameter stoch_params::eps_pert. \f$\delta\f$ is the timescale separation parameter given by the parameter stoch_params::tdelta. By setting those to 1, one recover the first equations above.

The function \f$\Psi_x\f$ includes all the \f$\boldsymbol x\f$ terms, and thus \f$F_x\f$ and \f$\Psi_x\f$ are unequivocally defined.
On the other hand, depending on the value of the parameter stoch_params::mode, the terms regrouped in the function \f$F_y\f$ can be different.
Indeed, if stoch_params::mode is set to
      - 'qfst', then:
	\f[   F_y(\boldsymbol{y})= \boldsymbol{B}^{yyy} : \boldsymbol{y} \otimes \boldsymbol{y} \f]
      - 'ures', then:
	\f[   F_y(\boldsymbol{y})= \boldsymbol{H}^y + \boldsymbol{L}^{yy}\cdot\boldsymbol{y} + \boldsymbol{B}^{yyy} : \boldsymbol{y} \otimes \boldsymbol{y} \f]
 However, for the WL parameterization, this parameter must be set to 'ures' by definition. See the article accompagnying this code for more details.

------------------------------------------------------------------------

## The MTV parameterization ##

This parameterization is also called homogenization. Its acronym comes from the names of
the authors that proposed this approach for climate modes: Majda, Timofeyev and Vanden Eijnden
(Majda et al., 2001). It is given by 

\f[ \frac{d\boldsymbol{x}}{dt} = F_X(\boldsymbol{x}) + \frac{1}{\delta} R(\boldsymbol{x}) + G(\boldsymbol{x}) + \sqrt{2} \,\, \boldsymbol{\sigma}(\boldsymbol{x}) \cdot \boldsymbol{dW} \f]

where \f$\boldsymbol{x}\f$ is the set of resolved variables and \f$\boldsymbol{dW}\f$ is a vector of standard Gaussian White noise.
\f$F_x\f$ is the set of tendencies of resolved system alone and \f$\delta\f$ is the timescale separation parameter.

### Correlations specification ###

The ingredients needed to compute the terms\f$R,G,\boldsymbol\sigma\f$ of this parametrization are the unresolved variables covariance matrix and the integrated correlation matrices.
The unresolved variables covariance matrix is given by
\f[\boldsymbol\sigma_y = \langle \boldsymbol y \otimes \boldsymbol y \rangle \f]
and is present in the implementation through the matrices corrmod::corr_i and corrmod::corr_i_full.
Their inverses are also available through corrmod::inv_corr_i and corrmod::inv_corr_i_full .
The integrated correlation matrices are given by
\f[ \boldsymbol\Sigma = \int_0^\infty \, ds \langle \, \boldsymbol y \otimes \boldsymbol y^s \rangle \f]
\f[ \boldsymbol\Sigma_2 = \int_0^\infty ds \, \left(\langle \boldsymbol y \otimes \boldsymbol y^s \rangle \otimes \langle \boldsymbol y \otimes \boldsymbol y^s \rangle\right) \f]
and is present in the implementation through the matrices int_corr::corrint and int_corr::corr2int .

These matrices are computed from the correlation matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$ which is accessible
through the function corrmod::corrcomp.
For instance, the covariance matrix \f$\boldsymbol\sigma_y\f$ is then simply the correlation matrix at the lagtime 0,
and \f$\boldsymbol\Sigma\f$ and \f$\boldsymbol\Sigma_2\f$ can be computed via integration over the lagtime.

There exists three different ways to load the correlation matrix, specified by the value
of the parameters stoch_params::load_mode and stoch_params::int_corr_mode .
The stoch_params::load_mode specify how the correlation matrix is loaded can take three different values:
     - 'defi': from an analytical definition encoded in the corrmod module function corrmod::corrcomp_from_def .
     - 'spli': from a spline definition file 'corrspline.def' .
     - 'expo': from a fit with exponentials definition file 'correxpo.def'

The stoch_params::int_corr_mode specify how the correlation are integrated and can take two different values:
     - 'file': Integration results provided by files 'corrint.def' and 'corr2int.def'
     - 'prog': Integration computed directly by the program with the correlation matrix. Write 'corrint.def' and 'corr2int.def' files to be reused later.

These parameters can be set up in the namelist file stoch_params.nml .
Examples of the ".def" files specifying the integrals are provided with the code.

### Other MTV setup parameters ###

Some additional parameters complete the options possible for the MTV parameters :
     - stoch_params::mnuti :  Multiplicative noise update time interval -- Time interval over which the matrix \f$\boldsymbol\sigma(\boldsymbol x)\f$ is updated.
     - stoch_params::t_trans_stoch : Transient period of the stochastic model.
     - stoch_params::maxint : Specify the upper limit of the numerical integration if stoch_params::int_corr_mode is set to 'prog'.

### Definition files ###

The following definition files are needed by the parameterization, depending on the value of the parameters described above.
Examples of those files are joined to the code. The files include:
     - 'mean.def' : Mean \f$\langle \boldsymbol y \rangle\f$ of the unresolved variables.  
     - 'correxpo.def': Coefficients \f$a_k\f$ of the fit of the elements of the correlations matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$ with the function \f[   a_4+a_0 \, \exp\left(-\frac{s}{a_1}\right) \, \cos(a_2 \, s + a_3) \f]
                       where \f$t\f$ is the lag-time and \f$\tau\f$ is the decorrelation time. Used if stoch_params::load_mode is set to 'expo'.
     - 'corrspline.def': Coefficients \f$b_k\f$ of the spline used to model the elements of the correlation matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$. Used if stoch_params::load_mode is set to 'spli'.
     - 'corrint.def': File holding the matrix \f$\boldsymbol\Sigma\f$. Used if stoch_params::int_corr_mode is set to 'file'.
     - 'corr2int.def': File holding the matrix \f$\boldsymbol\Sigma_2\f$.

The various terms are then constructed according to these definition files.
More details on the format of the definition files can be found [here](def_doc.md). 

------------------------------------------------------------------------

## The WL parameterization ##

This parameterization is based on the Ruelle response theory. Its acronym comes from the names of
the authors that proposed this approach: Wouters and Lucarini (Wouters and Lucarini, 2012). It is given by 

\f[\frac{d\boldsymbol{x}}{dt} = F_x(\boldsymbol{x}) + \varepsilon \, M_1(\boldsymbol{x}) + \varepsilon^2 \, M_2(\boldsymbol{x},t) + \varepsilon^2 \, M_3 (\boldsymbol{x},t)\f]

where \f$\varepsilon\f$ is the resolved-unresolved components coupling strength and where the different terms \f$M_i\f$ account for average, correlation and memory effects. 

### Correlations specification ###

The ingredients needed to compute the \f$M_i\f$ terms of this parametrization are the unresolved variable covariance matrix \f$\langle \boldsymbol y \otimes \boldsymbol y \rangle\f$ and correlation matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$.
The unresolved variables covariance matrix is given by
\f[\boldsymbol\sigma_y = \langle \boldsymbol y \otimes \boldsymbol y \rangle \f]
and is present in the implementation through the matrices corrmod::corr_i and corrmod::corr_i_full.
Their inverses are also available through corrmod::inv_corr_i and corrmod::inv_corr_i_full.

The correlation matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$ is accessible
through the function corrmod::corrcomp.

As for the MTV case, there exists three different ways to load the correlation matrix, specified by the value
of the parameters stoch_params::load_mode and stoch_params::int_corr_mode .
The stoch_params::load_mode specify how the correlation matrix is loaded can take three different values:
     - 'defi': from an analytical definition encoded in the corrmod module function corrmod::corrcomp_from_def .
     - 'spli': from a spline definition file 'corrspline.def' .
     - 'expo': from a fit with exponentials definition file 'correxpo.def'

The correlation term \f$M_2\f$ is emulated by an order \f$m\f$ multidimensional AutoRegressive (MAR) process:
\f[ \boldsymbol u_n = \sum_{i=1}^m \boldsymbol u_{n-i} \cdot \boldsymbol W_i + \boldsymbol Q \cdot \boldsymbol\xi_n \f]
of which the \f$\boldsymbol W_i\f$ and \f$\boldsymbol Q\f$ matrices are also needed (the \f$\boldsymbol\xi_n\f$ are vectors of standard Gaussian white noise). It is implemented in the MAR module.

### Other WL setup parameters ###

Some additional parameters complete the options possible for the WL parameters :
     - stoch_params::muti : Memory term \f$M_3\f$ update time interval.
     - stoch_params::t_trans_stoch : Transient period of the stochastic model.
     - stoch_params::meml : Time over which the memory kernel is numericaly integrated.
     - stoch_params::t_trans_mem : Transient period of the stochastic model to initialize the memory term.
     - stoch_params::dts : Intrisic resolved dynamics time step.
     - stoch_params::x_int_mode : Integration mode for the resolved component (not used for the moment -- must be set to 'reso').

Note that the stoch_params::mode must absolutely be set to 'ures', by definition.

### Definition files ###

The following definition files are needed by the parameterization, depending on the value of the parameters described above.
Examples of those files are joined to the code. The files include:
     - 'correxpo.def': Coefficients \f$a_k\f$ of the fit of the elements of the correlations matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$ with the function \f[   a_4+a_0 \, \exp\left(-\frac{s}{a_1}\right) \, \cos(a_2 \, s + a_3) \f]
                       where \f$t\f$ is the lag-time and \f$\tau\f$ is the decorrelation time. Used if stoch_params::load_mode is set to 'expo'.
     - 'corrspline.def': Coefficients \f$b_k\f$ of the spline used to model the elements of the correlation matrix \f$\langle \boldsymbol y \otimes \boldsymbol y^s \rangle\f$. Used if stoch_params::load_mode is set to 'spli'.
     - 'MAR_R_params.def': File specifying the \f$\boldsymbol R = \boldsymbol Q^2\f$ matrix for the MAR.
     - 'MAR_W_params.def': File specifying the \f$\boldsymbol W_i\f$ matrices for the MAR.

The various terms are then constructed according to these definition files.
More details on the format of the definition files can be found [here](def_doc.md). 

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

* Demaeyer, J., & Vannitsem, S.: Stochastic parametrization of subgrid-scale
processes in coupled ocean–atmosphere systems: benefits and limitations of response theory,
Quarterly Journal of the Royal Meteorological Society, 143(703), 881-896, 2017. 

Please see the main article for the full list of references.
