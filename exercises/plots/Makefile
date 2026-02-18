all: fig1.svg fig2.svg fig3.svg fig4.svg

fig4.svg: gpi/fig4.gpi
	gnuplot $<

fig3.svg: gpi/fig3.gpi 
	gnuplot $<

fig2.svg: gpi/fig2.gpi data/gamma.data
	gnuplot $<

fig1.svg: gpi/fig1.gpi data/erf.data
	gnuplot $<

clean:
	$(RM) *.svg
