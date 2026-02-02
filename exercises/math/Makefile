CXXFlags = -Wall -Wextra -Werror -std=c++23 -O
LINK.o = $(CXX) $(LDFLAGS)
LDLLIBS = -lstdc++ -lm
Out.txt. : main
	./main > Out.txt

main : main.o sfuns.o
	$(CXX) main.o sfuns.o -o main $(LDLLIBS)

main.o : sfuns.h main.cc
	$(CXX) $(CXXFlags) -c main.cc -o main.o

sfuns.o : sfuns.h sfuns.cc
	$(CXX) $(CXXFlags) -c sfuns.cc -o sfuns.o

.PHONEY:clean
clean:
	$(RM) --force main Out.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
