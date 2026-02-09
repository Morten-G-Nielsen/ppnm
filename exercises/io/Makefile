CXXFlags = -Wall -Wextra -Werror -std=c++23 -O
LDLLIBS = -lstdc++ -lm

all: Out.txt 
	cat Out.txt

Out.txt : main data.txt n.txt
	./main --input data.txt --output $@ $$(cat n.txt)  >> $@ 
	echo 1 2 3 | ./main >> $@

data.txt : Makefile
	echo 23 55 > $@
	echo 2.1 2.66 >> $@
	echo 3.1415 0.33333 >> $@

N=10
n.txt : Makefile
	>$@
	for i in $$(seq $N); do echo "-n " $$i >> $@; done

.PHONEY:clean
clean:
	$(RM) --force main *.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
