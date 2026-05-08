CXXFLAGS = -std=c++23 -I../lib/include -O3 -march=native
SRCS = ../lib/src/core/*.cc ../lib/src/qr.cc ../lib/src/ode.cc

all: Wavefunction.svg conv_eps.svg conv_acc.svg conv_rmax.svg conv_rmin.svg out.txt
	mv *.svg Plots/
	mv conv* Data/
	mv data.txt Data/

Wavefunction.svg: data.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Numerical wavefunction comparison"; \
		set xlabel "r"; \
		set ylabel "f(r)"; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "Calculated points", x*exp(-x) t "Analytical solution"; \
			' | gnuplot

conv_eps.svg: conv_eps.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Convergence for varying eps"; \
		set xlabel "eps"; \
		set ylabel "Error"; \
		set logscale xy; \
		set xrange [*:*] reverse; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "|E_{calc} - E_{exact}|"; \
			' | gnuplot

conv_acc.svg: conv_acc.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Convergence for varying acc"; \
		set xlabel "acc"; \
		set ylabel "Error"; \
		set logscale xy; \
		set xrange [*:*] reverse; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "|E_{calc} - E_{exact}|"; \
			' | gnuplot

conv_rmax.svg: conv_rmax.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Convergence for varying rmax"; \
		set xlabel "rmax"; \
		set ylabel "Error"; \
		set logscale y; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "|E_{calc} - E_{exact}|"; \
			' | gnuplot

conv_rmin.svg: conv_rmin.txt
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Convergence for varying rmin"; \
		set xlabel "rmin"; \
		set ylabel "Error"; \
		set logscale xy; \
		set xrange [*:*] reverse; \
		set grid; \
		plot \
		 "$<" u 1:2 w p t "|E_{calc} - E_{exact}|"; \
		 ' | gnuplot


conv_eps.txt: main
	> $@
	for e in $$(seq 2 10); do \
		echo -n "1e-$$e " >> $@; \
		./$< -eps 1e-$$e >> $@; done

conv_acc.txt: main
	> $@
	for a in $$(seq 2 10); do \
		echo -n "1e-$$a " >> $@; \
		./$< -acc 1e-$$a >> $@; done

conv_rmax.txt: main
	> $@
	for r in $$(seq 4 25); do \
		echo -n "$$r " >> $@; \
		./$< -rmax $$r >> $@; done

conv_rmin.txt: main
	> $@
	for r in $$(seq 1 10); do \
		echo -n "1e-$$r " >> $@; \
		./$< -rmin 1e-$$r >> $@; done

data.txt: main
	./$< -id 2 -acc 1e-6 -eps 1e-6 -rmin 1e-3 -rmax 8 > $@

out.txt: test
	./$< > $@

test: test.cc root.cc
	$(CXX) $^ $(SRCS) $(CXXFLAGS) -o $@	

main: main.cc root.cc
	$(CXX) $^ $(SRCS) $(CXXFLAGS) -o $@

clean:
	$(RM) Plots/* Data/* main test *.txt
