P=cylinder

all: $(P).xml

$(P).msh: $(P).geo
	gmsh -2 $<

$(P).xml: $(P).msh
	dolfin-convert -i gmsh -o xml $< $(P).xml

clean:
	rm -rf *.msh *.xml *.pvd *.vtu *.bak *.pos
