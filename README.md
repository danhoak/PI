# PI
Parametric Instabilities

Collection of scripts to model parametric instability in aLIGO and AdVirgo.

PI_tools.py contains functions that implement the matrix calculation of parametric gain in Evans et al (arXiv:0910.2716).

PI_example.py uses the PI_tools.py functions to reproduce Figs 4 and 7 from Evans et al. (Note, there is a factor of two error in Eq. 9 from Evans et al. which was corrected in later publications, e.g. arXiv:1502.06058.)

PI_montecarlo.py calculates the parametric gain for mechanical modes in the Advanced Virgo interferometer, using as inputs the mechanical mode frequencies and the overlaps between the mechanical modes and optical modes.  It marginalizes over uncertainties in Q, Bmn, and cavity Gouy phase to generate an estimate of the probability to encounter PI at various circulating powers in the arms.