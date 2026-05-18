CXXFLAGS = -std=c++23 -I../lib/include -O3 -march=native

FILES = core/*.cc qr.cc minimize.cc quad.cc ode.cc
SRCS = $(addprefix ../lib/src/, $(FILES))

all: ODE.svg Interp.svg
	mv *.svg Plots/
	mv *.txt Data/

ODE.svg: data2.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "DE: y″ = -sin(y)"; \
		set xlabel "Time"; \
		set ylabel "y"; \
		set grid; \
		plot \
		"$<" u 1:2 w l t "NN", "$<" u 3:4 w p t "RK45"; \
		' | gnuplot

Interp.svg: data.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set title "Interpolation of cos(5x-1)*exp(-x*x)"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set grid; \
		plot \
		"$<" u 1:2 w l t "True function", "$<" u 1:3 w p t "NN interpolation"; \
		' | gnuplot

data2.txt: main
	./$< -id 2 > $@

data.txt: main
	./$< > $@

main: main.cc ann.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) main *.txt *.svg
