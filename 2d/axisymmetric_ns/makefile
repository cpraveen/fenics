all: gmsh xml

gmsh: annulus.msh

annulus.msh: annulus.geo
	gmsh -2 annulus.geo

xml: annulus.xml

annulus.xml: annulus.msh
	dolfin-convert -i gmsh -o xml annulus.msh annulus.xml

clean:
	rm -f *.pvd *.vtu

allclean:
	rm -f *.pvd *.vtu *.msh *.xml *.bak
