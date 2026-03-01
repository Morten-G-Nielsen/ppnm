CXXFLAGS = -std=c++23 -Wall
LDLIBS = -lstdc++ -lm

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
	for N in $$(seq 100 20 400); do \
		time --format "$$N %e %U" --output $@ --append \
		./$< -size $$N ;\
		done

Out.txt: main
	./$< >$@ 

main: main.cc


clean:
	$(RM) Out* main
