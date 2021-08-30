CC          = gcc
CFLAGS      = -Wall -g -O2
LDFLAGS     = -lz -lm
prefix      = /usr/local
exec_prefix = $(prefix)/bin

src = $(wildcard *.c)
obj = $(src:.c=.o)

trim: $(obj)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) trim

.PHONY: install
install:
	cp trim $(exec_prefix)







