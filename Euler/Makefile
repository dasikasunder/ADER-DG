FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8
SRC=ader_dg.f90 \
    analyze_soln.f90 \
    basis.f90 \
    calc_time_step.f90 \
    element_update.f90 \
    initialize.f90 \
    la.f90 \
    limiter.f90 \
    main.f90 \
    pde.f90 \
    quadrature.f90 \
    riemann_problem.f90 \
    space_time_predictor.f90 \
    surface_integral.f90 \
    volume_integral.f90 \
    write_data.f90
    
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<
aderdg: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ}
clean: 
	@rm -f *.mod *.o aderdg
distclean: 
	@rm -f *.mod *.o *.vtk aderdg
