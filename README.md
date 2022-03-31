### Phonon SED
The *PhononSED* Fortran-90 code calculates phonon projected Spectral Energy Densities (SEDs) from molecular dynamics atomic velocity data and phonon eigenvectors calculated with GULP.

Compile with
`make`

Run as:
`./PhononSED.x < input/PhononSED.inp`
`python fitter.py`

To compile a parallel version with MPI use
`make parallel`

To run the parallel version, use
`mpirun -stdin all -np *number_of_processors* ./PhononSED.x < input/PhononSED.inp`

*fitter.py* is an example python code calculates phonon lifetimes by fitting Lorentzians to the SED data. A newer version has been developed in Matlab.


### References
* G. Kumar, F. G. VanGessel, D. C. Elton, P. W. Chung, *MRS Adv.*, **4(40)**, 2191-2199 (2019)
* G. Kumar, F. G. VanGessel, P. W. Chung, *Propellants, Explos. Pyrotech.*, **45(2)**, 169-176 (2020)
* J. M. Larkin, Ph.D. thesis, Carnegie Mellon University, 2013
* Larkin, et al., *Phys. Rev. B* **81**, 081411(R) (2010)
* A. J. H. McGaughey and M. Kaviany, *Phys. Rev. B* **69**, 094303 (2004).
* J. E. Turney, E. S. Landry, A. J. H. McGaughey, and C. H. Amon, *Phys. Rev. B* **79**, 064301 (2009).
