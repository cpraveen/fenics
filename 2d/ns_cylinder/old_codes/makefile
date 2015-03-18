P=cylinder_in_channel

all: gmsh xml subdomains

gmsh: $(P).geo
	gmsh -2 $<

xml: $(P).msh
	dolfin-convert -i gmsh -o xml $< $(P).xml

subdomains: $(P).xml
	python ./subdomains.py

clean:
	rm -rf *.msh *.xml *.pvd *.vtu *.bak *.pos
