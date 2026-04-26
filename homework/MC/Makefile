CXXFLAGS = -std=c++23 -O3 -march=native -Iinclude

all: Convergence.svg Comparison.svg out.txt
	mv Convergence.svg Comparison.svg Plots/
	$(RM) data*

Convergence.svg: data2.txt
	echo ' \
		set terminal svg background "white"; \
		set output "$@"; \
		set title "True error/Estimated error\nFor x^2 + y^2 <= 1"; \
		set xlabel "Points N"; \
		set ylabel "Error"; \
		set logscale xy; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "True error", \
			"$<" u 1:3 w p t "Estimated error", \
			2/sqrt(x) w l lc rgb "black" dt 2 t "2/sqrt(N)"; \
		' | gnuplot

Comparison.svg: data.txt
	echo ' \
		set terminal svg background "white"; \
		set output "$@"; \
		set title "Comparison in methods\nFor f(x, y)=x^2 + y^2;  x,y∈[0, 1]"; \
		set xlabel "Points N"; \
		set ylabel "Error from result"; \
		set logscale xy; \
		set grid; \
		plot \
			"$<" u 1:2 w p t "Plain", \
			"$<" u 1:3 w p t "Quasi", \
			"$<" u 1:4 w p t "Stratified"; \
	' | gnuplot

data2.txt: scaling
	for ((n=200; n<=2000000; n*=2)); do \
		./$< -N $$n >> $@; \
	done

data.txt: main
	for ((n=200; n<=2000000; n*=2)); do \
		./$< -id 2 -N $$n >> $@; \
	done

out.txt: main
	./$< -N 1e5 > $@

scaling: src/scaling.cc src/mc.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

main: src/main.cc src/mc.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) main scaling *.txt Plots/*.svg
