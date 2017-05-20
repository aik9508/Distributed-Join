CC = g++

Distr_join.o: Relation.hpp Relation.cpp
	mpic++ -c Relation.cpp -o Distr_join.o

Distr_join: Distr_join.o Distr_join.cpp
	mpic++ Distr_join.o Distr_join.cpp -o Distr_join

HyperCube.o: Relation.hpp Relation.cpp
	mpic++ -c Relation.cpp -o HyperCube.o

HyperCube: HyperCube.o HyperCube.cpp
	mpic++ HyperCube.o HyperCube.cpp -o HyperCube

Multi_join.o: Relation.hpp Relation.cpp
	mpic++ -c Relation.cpp -o Multi_join.o

Multi_join: Multi_join.o Multi_join.cpp
	mpic++ Multi_join.o Multi_join.cpp -o Multi_join

test_sequential.o: Relation.hpp Relation.cpp
	$(CC) -c Relation.cpp -o test_sequential.o

test_sequential: test_sequential.o test_sequential.cpp
	$(CC) test_sequential.o test_sequential.cpp -o test_sequential

test1: 
	mpirun -np 4 Distr_join "../data/facebook.dat" "../data/facebook.dat"

test2:
	mpirun -np 2 Distr_join "../data/data1.txt" "../data/data2.txt"

test3:
	mpirun -np 8 Distr_join "../data/dblp.dat" "../data/dblp.dat"

test4:
	mpirun -np 4 Distr_join "../data/data3.txt" "../data/data3.txt"

test5:
	mpirun -np 1 Multi_join 3 4 "../data/data1.txt" 2 0 1 "../data/data2.txt" 2 1 2 "../data/data2.txt" 2 2 3

test6:
	mpirun -np 4 Multi_join 3 3 "../data/facebook.dat" 2 0 1 "../data/facebook.dat" 2 1 2 "../data/facebook.dat" 2 2 0

hypercube_test:
	mpirun -np 8 HyperCube "../data/dblp.dat"

create_data: create_data.cpp
	$(CC) create_data.cpp -o create_data

clean:
	rm -f Distr_join create_data HyperCube Relation test_sequential *.o

