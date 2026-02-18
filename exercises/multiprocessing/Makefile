CXXFLAGS = -std=c++23 -Wall -O
LDLIBS = -lstdc++ -lm -lpthread

times.svg: out.times
	echo '\
		set terminal svg ;\
		set output "$@" ;\
		set xlabel "number of threads" ;\
		set ylabel "real time" ;\
		a=1; b=10;\
		f(x)=a+b/x;\
		fit f(x) "$<" using 1:2 via a,b ;\
		plot \
			"$<" using 1:2 with linespoints title "times" \
			,f(x) title "fit with 1/n" ;\
		' | tee log.times | gnuplot

TIME = time --append --output $@ --format "$$nthreads %e %U" 
N=1e11
out.times : main
	>$@
	>log.main
	for nthreads in $$(seq 8); \
			do $(TIME) ./main -terms $(N) -threads $$nthreads >> log.main; \
		done

clean:
	$(RM) log* out* main

test:
	echo $$0
