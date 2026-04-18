CXXFLAGS = -std=c++23 -Iinclude -O3 -march=native

# erf(1) from wolfram alpha
val = 0.8427007929497148693412206350826092592960669979663029084599378978

all: erfofone.svg erf.svg out.txt

erfofone.svg: erfofone.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "acc"; \
		set ylabel "err"; \
		set logscale xy; \
		set grid; \
		plot \
		"$<" u 1:(abs($$3 - $(val))) w p t "Error"; \
		' | gnuplot
	mv $@ plots
	$(RM) $<


erf.svg: erfdata.txt
	echo '\
		set terminal svg background rgb "white"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "erf(x)"; \
		set grid; \
		plot \
		"$<" u 1:2 w l dt 1 lw 2 t "My erf", erf(x) w l dt 2 lw 2 t "Gnuplot erf(x)"; \
		' | gnuplot
	mv $@ plots
	$(RM) $<

erfofone.txt: erf
	> $@
	for i in $$(LC_NUMERIC=C seq 1 15); do \
		echo -n "1e-$$i " >> $@; \
		./$< -z 1 -acc 1e-$$i -eps 0 >> $@; done

erfdata.txt: erf
	> $@
	for i in $$(LC_NUMERIC=C seq -3 0.01 3); do ./$< -z $$i >> $@; done	

out.txt: main
	./$< > $@
	python src/comp.py >> $@

erf: src/erf.cc src/quad.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

main: src/main.cc src/quad.cc
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) *.txt *.svg erf main
