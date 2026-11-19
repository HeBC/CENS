FC = gfortran
FFLAGS = -O3 -fPIC -ffree-line-length-none
OBJS = constants.o ang_mom_module.o spin_isospin_operators.o chiral_nnlo_opt_simple.o source2_api.o projection_math.o partialwaves.o

.PHONY: all clean
all: libn2loopt.a libn2loopt.so

constants.o: vendor/constants.f90
	$(FC) $(FFLAGS) -c $< -o $@
ang_mom_module.o: vendor/ang_mom_module.f90
	$(FC) $(FFLAGS) -c $< -o $@
spin_isospin_operators.o: vendor/spin_isospin_operators.f90
	$(FC) $(FFLAGS) -c $< -o $@
chiral_nnlo_opt_simple.o: vendor/chiral_nnlo_opt_simple.f90 constants.o ang_mom_module.o spin_isospin_operators.o
	$(FC) $(FFLAGS) -c $< -o $@
source2_api.o: source2_api.f90 chiral_nnlo_opt_simple.o
	$(FC) $(FFLAGS) -c $< -o $@
projection_math.o: projection_math.f90
	$(FC) $(FFLAGS) -c $< -o $@
partialwaves.o: partialwaves.f90 source2_api.o projection_math.o
	$(FC) $(FFLAGS) -c $< -o $@
libn2loopt.a: $(OBJS)
	$(AR) rcs $@ $^
libn2loopt.so: $(OBJS)
	$(FC) -shared -o $@ $^
clean:
	rm -f $(OBJS) *.mod libn2loopt.a libn2loopt.so
