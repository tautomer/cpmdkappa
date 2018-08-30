######  Fortran Complier  ##################################
#
#
#

objects = \
modules.o init.o readtraj.o cpmdkappa.o

wham.out : $(objects)
	${FC} ${FLAGS} -o kappa.out $(objects)

$(objects): %.o : %.F90
	${FC} ${FLAGS} -c $<

#
#
############################################################
