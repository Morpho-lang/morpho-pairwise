src = $(wildcard *.c) $(wildcard */*.c)
obj = $(src:.c=.dylib)

LDFLAGS  = -lm -lcblas -llapack -lcxsparse -L/opt/homebrew/lib
CFLAGS   =  -undefined dynamic_lookup -dynamiclib $(EXTCFLAGS) -O3 -I/opt/homebrew/include

pairwise: 
	$(CC) -o $@.dylib $(src) $(LDFLAGS) $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) 

.PHONY: install
install: pairwise
	cp pairwise.dylib ../lib/pairwise.dylib
	make clean 
	