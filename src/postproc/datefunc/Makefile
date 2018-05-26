# The following three directory may need to set.
POSTPPATH  = .
FC      = ifort

OBJS = stringfunc.o time_manager.o datefunc.o

datefunc: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.f90
	$(FC) -c $<

clean:
	rm -f *.o *.mod datefunc
