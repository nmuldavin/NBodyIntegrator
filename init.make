init : makegalaxy.o helperfunctions.o Vector.o treenode.o maketree.o treeforce.o integrate.o dataio.o
	g++ -o init makegalaxy.o helperfunctions.o Vector.o treenode.o maketree.o treeforce.o integrate.o dataio.o

makegalaxy.o : makegalaxy.cpp library/helperfunctions.hpp library/dataio.hpp library/integrate.hpp library/maketree.hpp library/treeforce.hpp library/treenode.hpp library/Vector.hpp
	g++ -c makegalaxy.cpp

helperfunctions.o : library/helperfunctions.cpp library/helperfunctions.hpp library/dataio.hpp
	g++ -c library/helperfunctions.cpp

dataio.o : library/dataio.cpp library/dataio.hpp library/integrate.hpp
	g++ -c library/dataio.cpp

integrate.o : library/integrate.cpp library/integrate.hpp library/treeforce.hpp
	g++ -c library/integrate.cpp

maketree.o : library/maketree.cpp library/maketree.hpp library/treenode.hpp
	g++ -c library/maketree.cpp

treeforce.o : library/treeforce.cpp library/treeforce.hpp library/treenode.hpp
	g++ -c library/treeforce.cpp

treenode.o : library/treenode.cpp library/treenode.hpp library/Vector.hpp
	g++ -c library/treenode.cpp

Vector.o : library/Vector.cpp library/Vector.hpp
	g++ -c library/Vector.cpp