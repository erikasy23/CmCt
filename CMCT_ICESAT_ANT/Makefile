# Makefile for CMCT
# to build on ggsghpcc
#
# I'm not sure this is currently using the best way to build .mod files,
# especially since gfortran doesn't create new ones if the interfaces
# haven't changed.  Maybe they should have their own rules, perhaps using
# rm foo.mod, gfortran -fsyntax-only foo.f08?
#
#
# AUTHOR: Jeff Guerber, GSFC 615/SigmaSpace, Nov. 2015
# HISTORY:
# 2016-03-21 JRG: FSON_INC, FSON_LIB, NCDF_INC macros.
# 2016-03-24 JRG: icesat_cmct_rec, icesat_cmct_class (but not yet in build)
# 2016-04-06 JRG: cmct_runcfg_mod
# 2016-04-14 JRG: model_netcdf_mod
# 2016-05-12 JRG: model_cism_txt_mod no longer requires NetCDF.
# 2016-05-25 JRG: icesat_config_mod
# 2016-05-27 JRG: icesat_cmctnc_class_mod.  icesat_cism_class_mod no longer
#    needs -I$FSON_INC.
# 2016-06-09 JRG: New datasets_mod now depends on icesat_*_mod.mod,
#    cmct_main now doesn't.
# 2016-07-11 JRG: Set permissions of launch scripts.
#
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


FSON_INC = ${HOME}/CmCt/externalpackages/fson/build     # fson*.mod are here
FSON_LIB = ${HOME}/CmCt/externalpackages/fson/dist      # libfson.so is here
NCDF_INC = ${HOME}/.conda_envs/CmCt/include             # netcdf.mod is here
NCDF_LIB = ${HOME}/.conda_envs/CmCt/lib                 # libnetcdff.so is here
 
##all:
##	@echo "Please specify exact target (for now)"


cmct_main: $(OBJS)
	gfortran -g -O0 -o cmct_main $(OBJS) -L$(FSON_LIB) -lfson -L$(NCDF_LIB) -lnetcdff
	chmod a+rx cmct_launch.ksh cmct_launch_config.ksh

cmct_main.o: cmct_main.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod \
		model_cism_txt_mod.mod    moments2d_class_mod.mod     \
		datasets_mod.mod          histogram_class_mod.mod     \
		cmct_runcfg_mod.mod       model_netcdf_mod.mod
	gfortran -c -g -O0 cmct_main.f08 -I$(FSON_INC) -I$(NCDF_INC)


cmct_runcfg_mod.o  cmct_runcfg_mod.mod: \
		cmct_runcfg_mod.f08
	gfortran -c -g -O0 cmct_runcfg_mod.f08 -I$(FSON_INC)


cism_polarstereo_mod.o  cism_polarstereo_mod.mod: \
		cism_polarstereo_mod.f03 \
		kinds_mod.mod  polar_stereographic_mod.mod
	gfortran -c -g -O0 cism_polarstereo_mod.f03


datasets_mod.o  datasets_mod.mod: \
		datasets_mod.f08 \
		icesat_cism_rec_mod.mod   icesat_cism_class_mod.mod  \
		icesat_cmct_rec_mod.mod   icesat_cmctnc_class_mod.mod
	gfortran -c -g -O0 datasets_mod.f08 -I$(FSON_INC)


histogram_class_mod.o  histogram_class_mod.mod:  histogram_class_mod.f08
	gfortran -c -g -O0 histogram_class_mod.f08


icesat_cism_class_mod.o  icesat_cism_class_mod.mod: \
		icesat_cism_class_mod.f08 \
		icesat_cism_rec_mod.mod  icesat_config_mod.mod
	gfortran -c -g -O0 icesat_cism_class_mod.f08


icesat_cism_rec_mod.o  icesat_cism_rec_mod.mod:  \
		icesat_cism_rec_mod.f08
	gfortran -c -g -O0 icesat_cism_rec_mod.f08


icesat_cmct_class_mod.o  icesat_cmct_class_mod.mod: \
		icesat_cmct_class_mod.f08 \
		icesat_cmct_rec_mod.mod  icesat_config_mod.mod
	gfortran -c -g -O0 icesat_cmct_class_mod.f08 -I$(FSON_INC)


icesat_cmctnc_class_mod.o  icesat_cmctnc_class_mod.mod: \
		icesat_cmctnc_class_mod.f08 \
		icesat_cmct_rec_mod.mod  icesat_config_mod.mod
	gfortran -c -g -O0 icesat_cmctnc_class_mod.f08 -I$(NCDF_INC)


icesat_cmct_rec_mod.o  icesat_cmct_rec_mod.mod:  \
		icesat_cmct_rec_mod.f08
	gfortran -c -g -O0 icesat_cmct_rec_mod.f08


icesat_config_mod.o  icesat_config_mod.mod: \
		icesat_config_mod.f08
	gfortran -c -g -O0 icesat_config_mod.f08 -I$(FSON_INC)


kinds_mod.o  kinds_mod.mod:  kinds_mod.f90
	gfortran -c -g -O0  kinds_mod.f90


model_cism_txt_mod.o  model_cism_txt_mod.mod: \
		model_cism_txt_mod.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod
	gfortran -c -g -O0 model_cism_txt_mod.f08


model_netcdf_mod.o  model_netcdf_mod.mod: \
		model_netcdf_mod.f08 \
		cism_polarstereo_mod.mod  polar_stereographic_mod.mod
	gfortran -c -g -O0 model_netcdf_mod.f08 -I$(NCDF_INC)


moments2d_class_mod.o  moments2d_class_mod.mod: \
		moments2d_class_mod.f08
	gfortran -c -g -O0 moments2d_class_mod.f08


polar_stereographic_mod.o  polar_stereographic_mod.mod: \
		polar_stereographic_mod.f90  \
		kinds_mod.mod
	gfortran -c -g -O0 polar_stereographic_mod.f90

clean:
	rm -rf *.o *.mod cmct_main

