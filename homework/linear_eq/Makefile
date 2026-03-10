CXXFLAGS = -std=c++23 -Wall -I../lib/include -O3 -march=native
LDLIBS = -lstdc++ -lm
LIB_SRCS = $(wildcard ../lib/src/core/*.cc)
SRCS = $(wildcard src/*.cc)

all: times.svg Out.txt

times.svg: Out.times
	echo '\
		set terminal svg ;\
		set title "QR decomp for NxN matrix" ;\
		set output "$@" ;\
		set xlabel "N" ;\
		set ylabel "time [s]" ;\
		a = 1;\
		f(x) = a*(x**3);\
		fit f(x) "$<" using 1:2 via a ;\
		plot \
			"$<" using 1:2 with points title "times" \
			, f(x) title "fit with n³" ;\
		' | gnuplot

Out.times: main
	>$@
	for N in $$(seq 100 60 1000); do \
		time --format "$$N %U" --output $@ --append \
		./$< -size $$N > /dev/null ;\
		done

Out.txt: main
	./$< > $@ 

main: $(SRCS) $(LIB_SRCS)
	$(CXX) $^ $(CXXFLAGS) $(LDLIBS) -o $@


clean:
	$(RM) Out* *.svg *.log main
