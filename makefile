
%.o : %.f90
	pgf90 -O3 -C -c $<

OBJECTS = kinds.o common_variables.o forces.o calc_bond.o calc_angle.o calc_vdw.o coulomb.o\
      coulomb_sf.o coulomb_dsf.o fcheck.o dists.o velocity_verlet.o read_data.o thermostat.o read_input.o\
      initvel.o boltz_vel.o dump.o thermo_dump.o main.o leapfroginte.o pseudo_e_tb.o run_engines.o\
	qm_forces.o qm_dump.o write_restart.o planczos.o rank3d.o indexx.o qmsetup.o grid3d_new.o kinetic.o


main: $(OBJECTS)
	pgf90 -O3 -C $(OBJECTS) -o main -llapack -lblas -lm 

