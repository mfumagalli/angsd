CC=gcc
CXX=g++

FLAGS=-O3 


all: supersim thetaStat emOptim2 msToGlf

bfgs.o: bfgs.c bfgs.h
	$(CC) $(FLAGS) -c bfgs.c -Wno-vexing-parse
simnextgen: simnextgen10.c
	$(CC) $(FLAGS) simnextgen10.c -o simnextgen -lm -lz

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim
emOptim2: emOptim2.cpp
	$(CXX) $(FLAGS) emOptim2.cpp bfgs.cpp -o emOptim2 -lz -lpthread -Wno-write-strings
thetaStat:	thetaStat.cpp stats.cpp ../kstring.o ../bgzf.o
	$(CXX) $(FLAGS) thetaStat.cpp -I../  ../kstring.o ../bgzf.o ../knetfile.o -lz -o thetaStat
msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm

clean:
	rm  -f *.o optimSFS simnextgen supersim emOptim calcStat bgid emOptim2 msToGlf thetaStat *~
