# XLab-FFTBarotropic
# Author: Hsu, Tien-Yiao
# Purpose:
This program is to run the 2-dimensional bartropic model with additional output of
  (1) filamentation time (Rozoff et al. 2006) computed and
  (2) effective eddy diffusivity (Hendricks and Schubert 2009)
  (3) deformation factor (developed by Hsu, Tien-Yiao)
# Method
  (1) Pseudo-spectral method (2/3 rule applied)
  (2) Runge-Kutto of order 4. (RK4)
  
# Compatible Compiler
  (1) CygwinGcc 4.8.3 with C++ 11 standard (lambda expression involved)
  
# Dependency
  (1) FFTW3 3.3.4 or above (we use float version -> fftw3f library)
    # For yum users:
    yum install -y fftw-devel fftw-libs-single

# Initial test
    make
    ./setup_env.sh
    ./example.sh

# Build
    make

# Run Program
    # Generate example initial field
    ./bin/makefield.out   # output to "output/initial_vorticity.bin"
    
    # Execute main program
    ./bin/main.out
    
    # Customize input/output folder and initial file name
    ./bin/main.out -imyinput -omyoutput -Imyinitialfile

# Environment Variable
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:[path_to_workspace]/lib
    
    # or simply
    cd [path_to_workspace]
    . setup_env.sh
