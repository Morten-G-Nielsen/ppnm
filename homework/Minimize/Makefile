CXXFLAGS = -std=c++23 -I../lib/include -O3 -march=native
SRCS = ../lib/src/core/vector.cc ../lib/src/core/matrix.cc ../lib/src/qr.cc

all: out.txt Fit.svg

Fit.svg: higgs.txt main
	echo ' \
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "Energy [GeV]"; \
		set ylabel "Signal"; \
		set title "Best fit for Higgs data"; \
		set grid; \
		plot \
			"$<" u 1:2:3 w ye t "Data", \
			"tmp.txt" u 1:2 w l t "Fit"; \
	'	| gnuplot
	$(RM) tmp.txt

out.txt: test main
	./main < higgs.txt 1> tmp.txt 2> $@
	./$< >> $@

test: test.cc minimize.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o  $@

main: main.cc minimize.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) main test out.txt Fit.svg
