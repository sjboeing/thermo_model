# A simple model to explore whether different thermodynamic variables are similar

**Click  on one of the sliders to start exploring**

![example](https://github.com/sjboeing/thermo_model/blob/master/example_deltatv.png)

Inspired by George Bryan's work on this, but including entrainment and precipitation. For precipitation, a fall-out 'length-scale' is used for simplicity (as there is no prognostic-vertical velocity).

When entrainment is set proportional to  1/z, the entrainment constant chosen is typical for 1000m height. Entrainment only occurs one condensation takes place, and only above 200m for model stability.

Other options:
- MSE correction for the buoyancy production term in the prognostic equation for MSE
- No ice phase at all
- Freezing at zero degrees or over a temperature range
- Saturation with respect to liquid water or ice

Variables:
- s: entropy
- mse: moist static energy
- ths: entropy based potential temperature
- Various approximations (the main ones use constant values for Lv and cp)

Inputs:
- RH: uniform throughout atmsophere
- Virtual temperature lapse rate: idem
- Parcel temperature and moisture excess at surface with respect to this reference profile

There is a python visual interface which uses wx, but it should also be possible to compile the code in stand-alone mode and plot using e.g. gnuplot.

gfortran: f95 -fdefault-double precision-8 simpletest.f90 -o simpletest
./simpletest > gpfile

