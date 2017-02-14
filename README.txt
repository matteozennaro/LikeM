/******************************************************************************/
                  This is LikeM - version Feb 17
/******************************************************************************/

Developed by M. Zennaro - matteo.zennaro@unimi.it
This code is not (yet) intended for public use, and therefore it is still
subject to major modifications and tests. In the future it is intended to become
more flexible, making adding different likelihoods and data sets even easier.

What's this code?
/******************************************************************************/
LikeM is an MCMC sampler that implements a Metropolis-Hastings algorithm.
It is written in C and its main features are
1. MPI parallelization, it can run many chains at the same time
2. Outputs in getdist format for easy plotting
3. Param file for setting number of parameters, initial values and steps
4. Adaptive refinement of the jump probability, to speed up convergence
5. Convergence statistics (Gelman-Rubin) calculated at the end of the run
6. Check-pointing allows to restart a chain from where it stopped

Likelihoods already included:
1. Tinker model for fitting a Halo Mass Function
2. Q-model (like Cole et al 2005) for fitting non-linear bias

Likelihoods to be included:
1. Clustering Ratio (written for an older version, not yet included)
2. Planck covariance matrix (written for an older version, not yet included)
3. Full Planck likelihood (requires interface with PLC)

Known bugs and to-do list
1. Check-pointing needs to be updated for this parallel version
2. Gelman-Rubin needs to be extended for chains of different lengths
3. Gelman-Rubin test should be performed at regular intervals, stopping the
   chains if a convergence criterion has been fulfilled

Installation
/******************************************************************************/
Check the contents of the Makefile, in particular you should choose your
favorite compiler (default mpicc). For some likelihoods a link to the gsl
libraries might be required. Check that you have them installed on your machine

Then hit make and it should compile

Running chains
/******************************************************************************/
Some explanatory parameter files are provided (with extension .ini), you shoud
have a look at them.

/******************************************************************************/
