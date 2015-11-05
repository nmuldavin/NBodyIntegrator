nbody : main.o Vector.o treenode.o maketree.o treeforce.o integrate.o dataio.o
	g++ -o nbody main.o Vector.o treenode.o maketree.o treeforce.o dataio.o integrate.o
    
main.o : main.cpp library/treenode.hpp library/Vector.hpp library/maketree.hpp library/treeforce.hpp library/dataio.hpp library/integrate.hpp
	g++ -c main.cpp

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