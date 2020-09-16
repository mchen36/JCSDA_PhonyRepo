The Community Surface Emissivity Model (CSEM) system is a feature-rich highly modularized software package for accurately modeling the surface emissivity and reflectivity of various Earth’s surfaces in the wavelength range from visible to microwave. Enclosed in CSEM are not only the physical models based on sound radiative transfer equations, but also a variety of empirical and semi-empirical models, surface-type-based lookup tables, and global emissivity atlases, which were generally developed from specific satellite or ground observations. CSEM evolved from the surface modules of the Community Radiative Transfer Model (CRTM), but with completely reorganized data structures and new module interfaces so that it can be used as a subsystem of a general radiate transfer model, e.g., CRTM, or as an independent surface radiative transfer modelling system. More importantly, the object-oriented software structure provides a very flexible platform for implementing and testing new model components.  Multiple options of the same model kind (e.g., microwave land model) may be also easily implemented on the CSEM platform, which enables very easy comparisons among the models alike. 


How to build CSEM:

1) Go to Build/env.setup and source the specific compiler config file

   e.g., source gfortran.setup

   note: make sure the following three environment variables are
   already defined or included in the setup file.
 
   B shell (sh, bash)


         export NETCDF_HOME=path to the netcdf

         export HDF5_HOME=path to the hdf5

    C shell (csh)


         setenv NETCDF_HOME  path to the netcdf

         setenv HDF5_HOME    path to the hdf5

    run step 1) for the first time fresh installation or as long as
    these three ENV variables are changed.

2) Generate the file "configure"
    ./autogen.sh
3) Generate the file "Makefile", you may specify where the CSEM library will be
   installed. The default is the current directory
    ./configure --prefix=path for the CSEM library to be installed
4) make
5) make install

