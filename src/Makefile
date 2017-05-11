CC = g++

Relation: Relation.hpp Relation.cpp
	$(CC) Relation.cpp -o Relation

Distr_join.o: Relation.hpp Relation.cpp
	mpic++ -c Relation.cpp -o Distr_join.o

Distr_join: Distr_join.o Distr_join.cpp
	mpic++ Distr_join.o Distr_join.cpp -o Distr_join

test1: 
	mpirun -np 2 Distr_join "../data/facebook.dat" "../data/facebook.dat"

test2:
	mpirun -np 4 Distr_join "../data/data1.txt" "../data/data2.txt"

test3:
	mpirun -np 4 Distr_join "../data/dblp.dat" "../data/dblp.dat"