CXX = c++
CXXFlags = -Wall -Wextra -Werror -std=c++23 -O
LDLLIBS = -lstdc++ -lm
Out.txt. : main
	./main > Out.txt

main : main.o hello.o
	$(CXX) main.o hello.o -o main $(LDLLIBS)

main.o : hello.h main.cc
	$(CXX) $(CXXFlags) -c main.cc -o main.o

hello.o : hello.h hello.cc
	$(CXX) $(CXXFlags) -c hello.cc -o hello.o

.PHONEY:clean
clean:
	$(RM) --force main main.o Out.txt

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
