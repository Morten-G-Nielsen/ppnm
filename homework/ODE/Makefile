CXXFLAGS = -std=c++23 -I../lib/include

PLOTS = Three-Body.svg Relativistic-precession.svg Elliptical.svg Circular.svg Lotka-Volterra.svg Oscillator.svg

all: $(PLOTS)
	mv *svg plots/
	$(RM) *.txt

Three-Body.svg: Three-Body.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		plot \
		"$<" u 2:3 w p t "Planet 1", \
		"$<" u 4:5 w p t "Planet 2", \
		"$<" u 6:7 w p t "Planet 3"; \
			' | gnuplot

Three-Body.txt: main
	./$< -id 4 -a 0 -b 10 -y -0.970 0.243 0.970 -0.243 0 0 0.466 0.432 0.466 0.432 -0.932 -0.864 > $@

Relativistic-precession.svg: RPPO3.txt interpolate
	./interpolate --input $< --output 1.txt
	echo "0 0" > temp.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set size ratio -1; \
		set title "Relativistic precession"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		plot "temp.txt" u 1:2 w p title "Star", \
		"$<" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w p t "Calculatet points", \
		"1.txt" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w l lw 2 t "Relativistic precission"; \
			' | gnuplot
	$(RM) 1.txt temp.txt

Elliptical.svg: RPPO2.txt interpolate
	./interpolate --input $< --output 1.txt
	echo "0 0" > temp.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Elliptical orbit"; \
		set size ratio -1; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		plot "temp.txt" u 1:2 w p title "Star", \
		"$<" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w p t "Calculatet points", \
		"1.txt" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w l lw 2 t "Planet"; \
			' | gnuplot
	$(RM) 1.txt temp.txt

Circular.svg: RPPO1.txt interpolate
	./interpolate --input $< --output 1.txt
	echo "0 0" > temp.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Circular orbit"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set xrange [-1.5:1.5]; set yrange [-1.5:1.5]; \
		set grid; \
		plot \
			"temp.txt" u 1:2 w p title "Star", \
			"$<" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w p t "Calculatet points", \
			"1.txt" u (1/$$2)*cos($$1):(1/$$2)*sin($$1) w l lw 2 t "Planet"; \
			' | gnuplot
	$(RM) 1.txt temp.txt

RPPO3.txt: main
	./$< -id 3 -a 0 -b 100 -epsilon 0.01 -y 1 -0.5 -h 1e-3 -acc 1e-6 -eps 1e-6 > $@
RPPO2.txt: main
	./$< -id 3 -a 0 -b 50 -y 1 -0.5 -h 1e-3 -acc 1e-6 -eps 1e-6 > $@
RPPO1.txt: main
	./$< -id 3 -a 0 -b 10 -y 1 0 -h 1e-3 -acc 1e-6 -eps 1e-6 > $@

Lotka-Volterra.svg: Lotka-Volterra.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Lotka Volterra"; \
		set xlabel "time"; \
		set ylabel "x"; \
		set grid; \
		plot \
			"$<" u 1:2 w l t "Prey", "$<" u 1:3 w l t "Predetor"; \
			' | gnuplot
Lotka-Volterra.txt: main
	./$< -id 2 -a 0 -b 10 -y 10 5 > $@

Oscillator.svg: Oscillator.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Oscillator with friction"; \
		set xlabel "time"; \
		set ylabel "x"; \
		set grid; \
		plot \
			"$<" u 1:2 w l t "Position over time", "$<" u 1:3 w l t "Velocity over time"; \
			' | gnuplot
Oscillator.txt: main
	./$< -id 1 -a 0 -b 20 -y 3 0 > $@

interpolate: interpolate.cc ../lib/src/interp.cc ../lib/src/core/vector.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

main: main.cc ode.cc ../lib/src/core/vector.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) *.txt plots/*.svg interpolate main
