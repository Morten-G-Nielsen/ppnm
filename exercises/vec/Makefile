CXXFlags = -Wall -Wextra -Werror -std=c++23 -O
LINK.o = $(CXX) $(LDFLAGS)
LDLLIBS = -lstdc++ -lm

all: Out.txt
	cat $<

Out.txt : main
	./$< > $@

main : main.o 
main.o : vec.h 

.PHONEY:clean
clean:
	$(RM) --force main Out.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
