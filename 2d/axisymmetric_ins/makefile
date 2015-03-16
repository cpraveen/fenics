all: gmsh xml mark

gmsh: annulus.msh

annulus.msh: annulus.geo
	gmsh -2 annulus.geo

xml: annulus.xml

annulus.xml: annulus.msh
	dolfin-convert -i gmsh -o xml annulus.msh annulus.xml

mark: subdomains.xml

subdomains.xml: annulus.xml
	python ./subdomains.py

clean:
	rm -f *.pvd *.vtu

allclean:
	make clean
	rm -f *.msh *.xml *.bak
