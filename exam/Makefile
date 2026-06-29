CXXFLAGS = -std=c++23 -I../homework/lib/include -O3 -march=native
FILES = core/*.cc qr.cc
SRCS = $(addprefix ../homework/lib/src/, $(FILES))

plots: main
	echo "Generating numerical dataset streams..."
	./$< -id 0 > data/data_scaling.txt
	./$< -id 1 > data/data_evals.txt
	./$< -id 2 > data/data_stress_BB.txt
	./$< -id 3 > data/data_stress_GB.txt
	./$< -id 4 > data/data_trajectory.txt
	./$< -id 5 > data/data_chain.txt
	echo ' \
		set terminal svg size 800,600 dynamic enhanced font "Verdana,10" background rgb "white"; \
		set output "plots/plot_scaling.svg"; \
		set title "Computational Runtime Complexity Scaling (O(N^2) vs O(N^3))" font "Verdana,12,Bold"; \
		set xlabel "Dimension (N)"; \
		set ylabel "Execution Time (ms)"; \
		set grid lc rgb "#d6dbdf" lt 1 lw 0.5; \
		plot "data/data_scaling.txt" u 1:2 w lp lw 2 t "Good Broyden", \
			"data/data_scaling.txt" u 1:3 w lp lw 2 t "Bad Broyden", \
			"data/data_scaling.txt" u 1:4 w lp lw 2 t "Newton Baseline"; \
		' | gnuplot
	echo ' \
		set terminal svg size 800,600 dynamic enhanced font "Verdana,10" background rgb "white"; \
		set output "plots/plot_evals.svg"; \
		set title "Total Target Function Evaluations" font "Verdana,12,Bold"; \
		set xlabel "Dimension (N)"; \
		set ylabel "Number of f(x) Calls"; \
		set grid lc rgb "#d6dbdf" lt 1 lw 0.5; \
		plot "data/data_evals.txt" u 1:2 w lp lw 2 t "Good Broyden", \
			"data/data_evals.txt" u 1:3 w lp lw 2 t "Bad Broyden", \
			"data/data_evals.txt" u 1:4 w lp lw 2 t "Newton Baseline"; \
		' | gnuplot
	echo 'set terminal svg size 800,600 dynamic enhanced background rgb "white"; \
	       set output "plots/plot_stability_BB.svg"; \
	       set title "Convergence Stability Basin Map using Bad Broyden\n{/*0.8 f_1(x,y)=x^2+y^2-4, f_2(x,y)=sin(x+y^2)-0.5}" font "Verdana,12,Bold"; \
	       set xlabel "Initial Guess X_0"; set ylabel "Initial Guess Y_0"; \
				 set cblabel "Iteration Count" font "Verdana,10"; \
	       set view map; \
	       set pm3d at b; \
	       set palette rgb 33,13,10; \
	       splot "data/data_stress_BB.txt" using 1:2:3 with pm3d notitle; \
		' | gnuplot
	echo 'set terminal svg size 800,600 dynamic enhanced background rgb "white"; \
	       set output "plots/plot_stability_GB.svg"; \
	       set title "Convergence Stability Basin Map using Good Broyden\n{/*0.8 f_1(x,y)=x^2+y^2-4, f_2(x,y)=sin(x+y^2)-0.5}" font "Verdana,12,Bold"; \
	       set xlabel "Initial Guess X_0"; set ylabel "Initial Guess Y_0"; \
				 set cblabel "Iteration Count" font "Verdana,10"; \
	       set view map; \
	       set pm3d at b; \
	       set palette rgb 33,13,10; \
	       splot "data/data_stress_GB.txt" using 1:2:3 with pm3d notitle; \
		' | gnuplot
	echo ' \
		set terminal svg size 800,600 enhanced background rgb "white"; \
		set output "plots/plot_trajectory.svg"; \
		\
		set title "Backtracking vs. Quadratic Interpolation"; \
		set xlabel "x"; \
		set ylabel "y"; \
		set xrange [-6:6]; \
		set yrange [-6:6]; \
		set size ratio -1; \
		set grid; \
		set key top left; \
		\
		set samples 400; \
		set isosamples 400; \
		\
		R(x,y)=100*(y-x**2)**2+(1-x)**2; \
		\
		set contour; \
		unset surface; \
		unset clabel; \
		set view map; \
		set cntrparam levels discrete 0,1,2,5,10,20,50,100,200,400,800,1000,2000,3000; \
		\
		set table $$contours; \
		splot R(x,y); \
		unset table; \
		\
		plot \
				$$contours using 1:2 with lines lc rgb "#b0b7bd" title "Rosenbrock contours", \
				"data/data_trajectory.txt" using 1:2 with linespoints lw 1 pt 7 ps 1 lc rgb "#e74c3c" title "Backtracking", \
				"data/data_trajectory.txt" using 3:4 with linespoints lw 1 pt 6 ps 1 lc rgb "#2ecc71" title "Quadratic interpolation"; \
		' | gnuplot
	echo 'set terminal svg size 800,600 dynamic enhanced background rgb "white"; \
      set output "plots/plot_chain.svg"; \
      set title "Static Equilibrium Catenary Arc of an N=100 Hanging Chain" font "Verdana,12,Bold"; \
      set xlabel "Horizontal Position X"; set ylabel "Vertical Position Y"; \
      set grid lc rgb "#d6dbdf" lt 1 lw 0.5; \
      set xrange [-5.5:5.5]; \
      set yrange [-140:10]; \
      set key top right; \
      plot "data/data_chain.txt" using 1:2 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#9b59b6" title "Solved Catenary Profile";' | gnuplot
	echo "SVG plots generated in plots directory!"

main: main.cc root.cc $(SRCS)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	$(RM) main data/*.txt plots/*.svg
