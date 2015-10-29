PROG = mpl
FC = gfortran
#FFLAGS = -O3 -fno-automatic 
FFLAGS = -O3 
LDFLAGS = 
SRCS = nrtype.f90 dvmlm.f parser.f90 units.f90 command_line.f90 data.f90 model.f90 dvmlm_wrapper.f90 scrs.f90 mpl.f90 
OBJS = nrtype.o dvmlm.o parser.o units.o command_line.o data.o model.o dvmlm_wrapper.o scrs.o mpl.o 

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90 .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<


nrtype.o: nrtype.f90
dvmlm.o: dvmlm.f
parser.o: parser.f90
units.o: units.f90 
command_line.o: command_line.f90 nrtype.o
data.o: data.f90 nrtype.o parser.o
model.o: model.f90 nrtype.o data.o 
dvmlm_wrapper.o: dvmlm_wrapper.f90 dvmlm.o model.o nrtype.o
scrs.o: scrs.f90 nrtype.o data.o model.o units.o
mpl.o: mpl.f90 nrtype.o units.o parser.o command_line.o data.o model.o scrs.o dvmlm_wrapper.o 

.PHONY: clean realclean debug showtargets

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod

# realclean: remove all regenerable files, including executables
realclean: clean
	$(RM) $(PROG) $(OBJS)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g -fbounds-check" "OSTYPE=$(OSTYPE:-gnu=)"

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean realclean debug fcheck fdepend 

