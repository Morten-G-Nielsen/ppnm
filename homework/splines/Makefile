CXXFLAGS = -std=c++23 -I../lib/include

all: CubicSpline.svg QuadraticSpline.svg LinearSpline.svg TrueShape.svg
	mv $^ plots/

TrueShape.svg: True_shape.txt
	echo ' \
		set terminal svg; \
		set title "True shape"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "y=cos(x)+cos(2x+0.5)"; \
		set grid; \
		plot \
			"$<" u 1:2 w l title "Line" \
			, "$<" u 1:3 w l title "Deriv" \
			, "$<" u 1:4 w l title "Integ"; \
		' | gnuplot	
	$(RM) True_shape.txt
CubicSpline.svg: cubic_spline.txt data.txt
	echo ' \
		set terminal svg; \
		set title "Cubic spline"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "y=cos(x)+cos(2x+0.5)"; \
		set grid; \
		plot \
			"data.txt" u 1:2 w p title "data" \
			,"$<" u 1:2 w l title "Line" \
			, "$<" u 1:3 w l title "Deriv" \
			, "$<" u 1:4 w l title "Integ"; \
		' | gnuplot	
	$(RM) cubic_spline.txt
QuadraticSpline.svg: quad_spline.txt data.txt
	echo ' \
		set terminal svg; \
		set title "Quadratic spline"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "y=cos(x)+cos(2x+0.5)"; \
		set grid; \
		plot \
			"data.txt" u 1:2 w p title "data" \
			,"$<" u 1:2 w l title "Line" \
			, "$<" u 1:3 w l title "Deriv" \
			, "$<" u 1:4 w l title "Integ"; \
		' | gnuplot	
	$(RM) quad_spline.txt
LinearSpline.svg: line_spline.txt data.txt
	echo ' \
		set terminal svg; \
		set title "Linear spline"; \
		set output "$@"; \
		set xlabel "x"; \
		set ylabel "y=cos(x)+cos(2x+0.5)"; \
		set grid; \
		plot \
			"data.txt" u 1:2 w p title "data" \
			,"$<" u 1:2 w l title "Line" \
			, "$<" u 1:3 w l title "Deriv" \
			, "$<" u 1:4 w l title "Integ"; \
		' | gnuplot	
	$(RM) line_spline.txt

True_shape.txt: main
	./$< > $@
cubic_spline.txt: main
	./$< -c > $@
quad_spline.txt: main
	./$< -q > $@
line_spline.txt: main
	./$< -l > $@

main: main.cc interp.cc ../lib/src/core/vector.cc
	$(CXX) $^ $(CXXFLAGS) -o $@
