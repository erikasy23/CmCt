# Makefile for CMCT
#
# I'm not sure this is currently using the best way to build .mod files,
# especially since gfortran doesn't create new ones if the interfaces
# haven't changed.  Maybe they should have their own rules, perhaps using
# rm foo.mod, gfortran -fsyntax-only foo.f08?
#
# AUTHOR: Jeff Guerber, GSFC 615/SigmaSpace, Nov. 2015
# Last SVN commit on atlas: $Id: Makefile 110 2016-07-12 07:28:00Z jguerber $

OBJS =  cmct_main.o                \
	cmct_runcfg_mod.o          \
	cism_polarstereo_mod.o     \
	datasets_mod.o             \
	histogram_class_mod.o      \
	icesat_cism_class_mod.o    \
	icesat_cism_rec_mod.o      \
	icesat_cmctnc_class_mod.o  \
	icesat_cmct_rec_mod.o      \
	icesat_config_mod.o        \
	kinds_mod.o                \
	model_cism_txt_mod.o       \
	model_netcdf_mod.o         \
	moments2d_class_mod.o      \
	polar_stereographic_mod.o

#FSON_INC = /home/cmct/JSON/include  # fson*.mod are here: /usr/local/include/
#FSON_LIB = /home/cmct/JSON/lib      # libfson.so is here: /usr/local/lib
#NCDF_INC = /opt/csw/include             # netcdf.mod is here: /usr/local/include

FSON_INC = /usr/local/include  # fson*.mod are here
FSON_LIB = /usr/local/lib      # libfson.so is here
NCDF_INC = /usr/local/include  # netcdf.mod is here
NCDF_LIB = /usr/local/lib      # netCDF libs, e.g. libnetcdff.so

FC = gfortran
FFLAGS = -g -O0 -Wall
PROGRAM = cmct_main

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FFLAGS) -o cmct_main $(OBJS) -L$(FSON_LIB) -lfson -L$(NCDF_LIB) -lnetcdff
	chmod a+rx cmct_launch.ksh cmct_launch_config.ksh

cmct_main.o: cmct_main.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod \
		model_cism_txt_mod.mod    moments2d_class_mod.mod     \
		datasets_mod.mod          histogram_class_mod.mod     \
		cmct_runcfg_mod.mod       model_netcdf_mod.mod
	$(FC) -c $(FFLAGS) cmct_main.f08 -I$(FSON_INC) -I$(NCDF_INC)


cmct_runcfg_mod.o  cmct_runcfg_mod.mod: \
		cmct_runcfg_mod.f08
	$(FC) -c $(FFLAGS) cmct_runcfg_mod.f08 -I$(FSON_INC)


cism_polarstereo_mod.o  cism_polarstereo_mod.mod: \
		cism_polarstereo_mod.f03 \
		kinds_mod.mod  polar_stereographic_mod.mod
	$(FC) -c $(FFLAGS) cism_polarstereo_mod.f03


datasets_mod.o  datasets_mod.mod: \
		datasets_mod.f08 \
		icesat_cism_rec_mod.mod   icesat_cism_class_mod.mod  \
		icesat_cmct_rec_mod.mod   icesat_cmctnc_class_mod.mod
	$(FC) -c $(FFLAGS) datasets_mod.f08 -I$(FSON_INC)


histogram_class_mod.o  histogram_class_mod.mod:  histogram_class_mod.f08
	$(FC) -c $(FFLAGS) histogram_class_mod.f08


icesat_cism_class_mod.o  icesat_cism_class_mod.mod: \
		icesat_cism_class_mod.f08 \
		icesat_cism_rec_mod.mod  icesat_config_mod.mod
	$(FC) -c $(FFLAGS) icesat_cism_class_mod.f08


icesat_cism_rec_mod.o  icesat_cism_rec_mod.mod:  \
		icesat_cism_rec_mod.f08
	$(FC) -c $(FFLAGS) icesat_cism_rec_mod.f08


icesat_cmct_class_mod.o  icesat_cmct_class_mod.mod: \
		icesat_cmct_class_mod.f08 \
		icesat_cmct_rec_mod.mod  icesat_config_mod.mod
	$(FC) -c $(FFLAGS) icesat_cmct_class_mod.f08 -I$(FSON_INC)


icesat_cmctnc_class_mod.o  icesat_cmctnc_class_mod.mod: \
		icesat_cmctnc_class_mod.f08 \
		icesat_cmct_rec_mod.mod  icesat_config_mod.mod
	$(FC) -c $(FFLAGS) icesat_cmctnc_class_mod.f08 -I$(NCDF_INC)


icesat_cmct_rec_mod.o  icesat_cmct_rec_mod.mod:  \
		icesat_cmct_rec_mod.f08
	$(FC) -c $(FFLAGS) icesat_cmct_rec_mod.f08


icesat_config_mod.o  icesat_config_mod.mod: \
		icesat_config_mod.f08
	$(FC) -c $(FFLAGS) icesat_config_mod.f08 -I$(FSON_INC)


kinds_mod.o  kinds_mod.mod:  kinds_mod.f90
	$(FC) -c $(FFLAGS)  kinds_mod.f90


model_cism_txt_mod.o  model_cism_txt_mod.mod: \
		model_cism_txt_mod.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod
	$(FC) -c $(FFLAGS) model_cism_txt_mod.f08


model_netcdf_mod.o  model_netcdf_mod.mod: \
		model_netcdf_mod.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod
	$(FC) -c $(FFLAGS) model_netcdf_mod.f08 -I$(NCDF_INC)


moments2d_class_mod.o  moments2d_class_mod.mod: \
		moments2d_class_mod.f08
	$(FC) -c $(FFLAGS) moments2d_class_mod.f08


polar_stereographic_mod.o  polar_stereographic_mod.mod: \
		polar_stereographic_mod.f90  \
		kinds_mod.mod
	$(FC) -c $(FFLAGS) polar_stereographic_mod.f90

clean: 
	rm -f $(PROGRAM) *.o *.mod
