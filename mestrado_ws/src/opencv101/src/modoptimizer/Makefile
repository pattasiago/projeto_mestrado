src = $(wildcard *.c)
obj = $(src:.c=.o)

LDFLAGS = -lGL -lglut -lpng -lz -lm

myprog: $(obj)
	$(CC) -o $@ $^ meschach.a $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) myprog
