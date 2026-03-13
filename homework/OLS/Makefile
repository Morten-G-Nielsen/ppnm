CXXFLAGS = -std=c++23 -Wall -I../lib/include

SRCS_CORE = $(wildcard ../lib/src/core/*.cc)

all: fit.svg out.txt

fit.svg: data.txt main
	./main > param.gp
	( \
	echo "load 'param.gp'"; \
	echo 'f(x) = exp(a+b*x)'; \
	echo 'f1(x) = exp((a+da)+(b+db)*x)'; \
	echo 'f2(x) = exp((a-da)+(b-db)*x)'; \
	echo ' \
		set terminal svg; \
		set title "Decay"; \
		set output "$@"; \
		set xlabel "time [days]"; \
		set ylabel "Activity"; \
		set grid; \
		plot \
			"$<" u 1:2:3 w yerrorbars title "data" \
			, f(x) w l title "best fit" \
			, f1(x) w l title "c_i+δc_i" \
			, f2(x) w l title "c_i-δc_i"; \
	') | gnuplot
	$(RM) param.gp

out.txt: calc
	./$< > $@

calc: src/calc.cc $(SRCS_CORE) ../lib/src/qr.cc src/ols.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

main: src/main.cc $(SRCS_CORE) ../lib/src/qr.cc src/ols.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) *.svg out.txt calc main 
