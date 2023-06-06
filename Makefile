#	File name		: Makefile
#	Date			: December 2022
#	Version			: 3.0.0
#	Author			: Adhika
#	Co-Author		: Angga

DEST			= .

LIBS			= -larmadillo

INPS			= 

COMPILER		= g++

OPT_STD 		= -std=c++11

OPT_O 			= -O3

OPT_OMP			= -fopenmp

#OPT_MODEL		= -mcmodel=medium

OPTFLAG			= $(OPT_STD) $(OPT_O) $(OPT_DEPEND) $(OPT_OMP) $(OPT_MODEL)
 
MAKEFILE		= Makefile

PROGRAM			= main

# List of .cpp source file path
SRCS			= main.cpp \
				global.cpp \
				src/geometry/geometry.cpp	\
				src/geometry/2D_objects_generator.cpp\
				src/geometry/3D_objects_generator.cpp\
				src/geometry/u_var.cpp \
				src/geometry/moving_body.cpp \
				src/initialization/initialization.cpp \
				src/initialization/init2D.cpp \
				src/initialization/init3D.cpp \
				src/save_data/state_data.cpp \
				src/save_data/summary_log.cpp \
				src/save_data/force_output.cpp \
				src/save_data/force_linear_impulse.cpp \
				src/neighbor/neighbor.cpp \
				src/neighbor/cell_list_utils.cpp \
				src/neighbor/cell_list_init.cpp \
				src/neighbor/cell_list_ngh.cpp \
				src/neighbor/cell_list_adaptive.cpp \
				src/neighbor/direct_find.cpp \
				src/neighbor/link_list.cpp \
				src/neighbor/link_list_utils.cpp \
				src/neighbor/hashing_grid.cpp \
				src/neighbor/spatial_hashing.cpp \
				src/remeshing/kernel.cpp \
				src/remeshing/inter_search.cpp \
				src/remeshing/redistribute_particles.cpp \
				src/remeshing/dc_remeshing.cpp \
				src/remeshing/main_remesh.cpp \
				src/penalization/penalization-kai_def.cpp \
				src/penalization/penalization-no_slip_bc.cpp	\
				src/penalization/penalization.cpp \
				src/velocity_poisson/velocity_poisson.cpp \
				src/velocity_poisson/treeCell.cpp \
				src/velocity_poisson/fmm2D.cpp \
				src/velocity_poisson/regularization_function_2d.cpp \
				src/velocity_poisson/biotsavart_direct_2d.cpp \
				src/velocity_poisson/biotsavart_fmm_2d.cpp	\
				src/velocity_poisson/subroutine_all_fmm_2d.cpp	\
				src/pressure_poisson/pressure_poisson.cpp	\
				src/advection/main_advection.cpp \
				src/diffusion/main_diffusion.cpp \
				src/LSMPS/LSMPSa.cpp \
				src/LSMPS/LSMPSb.cpp \
				# src/solver_poisson/solver.cpp \
				# src/solver_poisson/solver_poisson.cpp \
				# src/diffusion/pse.cpp \
				# src/DC_operator/vandermonde.cpp \
				# src/DC_operator/dc_pse.cpp \
				# src/DC_operator/dc_gradient.cpp \
				# src/remeshing/sign_particles.cpp \
				# src/remeshing/basic_remeshing.cpp \
				# src/remeshing/particle_adaption/multiblock_adaption.cpp \
				# src/remeshing/particle_adaption/particle_adaption.cpp \
				# src/remeshing/particle_adaption/imr_removal_insertion.cpp \
				# src/remeshing/particle_adaption/imr_resolution_field.cpp \
				# src/remeshing/particle_adaption/imr_resolution_interp.cpp \
				# src/remeshing/particle_adaption/imr_searching.cpp \
				# src/remeshing/particle_adaption/imr_set_properties.cpp \
				# src/remeshing/particle_adaption/imr_steepest_descent.cpp \

# List of target .o object file path
OBJS			= $(SRCS:.cpp=.o)

# List of target .o object file path
.cpp.o:
			$(COMPILER) $(OPTFLAG) -c $*.cpp -o $*.o

# COMPILERF 		= gfortran

# F90FLAGS 		= -L/usr/include -L/usr/lib/x86_64-linux-gnu -I/usr/include -I/usr/local/include 

# SRCSF 			= src/Fortran/Utils

# OBJSF 			= $(SRCSF)/memory_fmm_2d.o \
# 				$(SRCSF)/subroutine_all_fmm_2d.o \
# 				$(SRCSF)/biotsavart_fmm_2d.o \

# %.o: %.f90
# 				# $(COMPILERF) -c $< -o $@
# 				$(COMPILERF) $(F90FLAGS) -c $< -o $@

all:			$(PROGRAM)

$(PROGRAM):	$(OBJS)
			$(COMPILER) $(OPT_OMP) -o $(PROGRAM) $(OBJS)
#			@echo -n "Loading Program $(PROGRAM) ... "
#			@$(COMPILER) $(OPTFLAG) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
			@echo "\n--------------- DONE COMPILE --------------"

run:
			@echo "--------------- RUN PROGRAM ---------------"
			./$(PROGRAM)

clean:;		@rm -f $(OBJS) # $(SRCS:.cpp=.il) $(SRCS:.cpp=.d) 
#			@rm -f $(SRCSF)/*.o $(DEST)/*.mod
			@rm -f $(PROGRAM)
			@echo "----------------- CLEANED -----------------"

delete_data:
			@while [ -z "$$CONTINUE" ]; do \
			read -r -p "Do you really want to delete data files [Y/N] ? " CONTINUE; \
			done ; \
			if [ $$CONTINUE != "y" ] && [ $$CONTINUE != "Y" ]; then \
			echo "Exiting." ; exit 1 ; \
			fi
			@rm -f output/*.dat output/*.jpg output/*.png output/*.csv
			@echo "-------------- DATA DELETED ---------------"
