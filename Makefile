PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CPPFLAGS=-I"/nfs/users/nfs_j/jl11/R-modules/RInside/include"
CXXFLAGS=-Wall -O3 -std=c++11
LDFLAGS=-L"/nfs/users/nfs_j/jl11/R-modules/RInside/lib" -lRInside

all: mlst_sa

clean:
	$(RM) *.o ~* mlst_sa

install: all
	install -d $(BINDIR)
	install mlst_sa $(BINDIR)

