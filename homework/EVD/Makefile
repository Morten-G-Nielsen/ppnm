CXXFLAGS = -std=c++23 -I../lib/include/core -Iinclude -O3 -march=native

SRCS = $(wildcard ../lib/src/core/*.cc)

all: time.svg wf.svg rmax.svg dr.svg out.txt

time.svg: out.times.txt
	echo '\
				set terminal svg; \
				set output "plots/$@"; \
				set title "Time vs N"; \
				set xlabel "N";\
				set ylabel "time s"; \
				a=1;\
				f(x)=a*x**3;\
				fit f(x) "$<" u 1:2 via a; \
				plot "$<" u 1:2 w points title "Calc time" \
						,f(x) title "fit with N³"; \
				' | gnuplot 
	$(RM) $<


wf.svg: wf.data.txt
	echo '\
			set terminal svg; \
			set output "plots/$@"; \
			set title  "Numerical vs Analytical Wavefunctions (l=0)"; \
			set xlabel "r (Bohr)"; \
			set ylabel "u(n,r)=r*R(n,r)"; \
			set grid; \
			plot "$<" u 1:2 w lines title "Numerical n=1"\
			,2*x*exp(-x) w lines dt 2 title "Theory n=1" \
			, "$<" u 1:3 w lines title "Numerical n=2" \
			, (1/sqrt(2))*x*(1-x/2)*exp(-x/2) w lines dt 2 title "Theory n=2" \
			, "$<" u 1:4 w lines title "Numerical n=3" \
			, 2*x*(1-2*x/3+2*(x**2)/27)*exp(-x/3)/3/sqrt(2) w lines dt 2 title "Theory n=3"; \
	' | gnuplot
	$(RM) $<
	

dr.svg: dr.data.txt
	echo '\
			set terminal svg; \
			set output "plots/$@"; \
			set title  "Ground state energy vs dr with r_{max}=20"; \
			set xlabel "dr (Bohr)"; \
			set ylabel "E_0 (Hartree)"; \
			set grid; \
			plot "$<" with linespoints title "Numerical E_0"\
			,-0.5 linetype 2 title "Theory (-0.5)"; \
	' | gnuplot
	$(RM) $<

rmax.svg: rmax.data.txt
	echo '\
			set terminal svg; \
			set output "plots/$@"; \
			set title  "Ground state energy vs r_{max} with dr=0.1"; \
			set xlabel "r_{max} (Bohr)"; \
			set ylabel "E_0 (Hartree)"; \
			set grid; \
			plot "$<" with linespoints title "Numerical E_0"\
			,-0.5 linetype 2 title "Theory (-0.5)"; \
	' | gnuplot
	$(RM) $<

out.times.txt: main
	>$@
	for N in $$(seq 100 50 1000); do \
		time -ao $@ -f "$$N %e" ./$< -n $$N& \
	done; wait
	sort $@ -o $@

wf.data.txt: main
	./$< -rmax 40 -dr 0.04 -wf true > $@

rmax.data.txt: main
	>$@
	for N in $$(LC_NUMERIC=C seq 2 0.5 10); do \
		echo "Calculating for rmax=$$N..."; \
		echo -n "$$N " >> $@; \
		./$< -rmax $$N -dr 0.1 >> $@; \
	done

dr.data.txt: main
	>$@
	for N in $$(LC_NUMERIC=C seq 1 -0.05 0.05); do \
		echo "Calculating for dr=$$N..."; \
		echo -n "$$N " >> $@; \
		./$< -rmax 20 -dr $$N >> $@; \
	done

out.txt: test
	./$< > $@

test: src/test.cc src/evd.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o $@

main: src/main.cc src/evd.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) plots/*.svg main test *.log
