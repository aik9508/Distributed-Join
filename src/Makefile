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
	mpirun -np 2 Distr_join "../data/twitter.dat" "../data/twitter.dat"

test3:
	mpirun -np 4 Distr_join "../data/dblp.dat" "../data/dblp.dat"

test4:
	mpirun -np 4 Distr_join "../data/data3.txt" "../data/data3.txt"

multi_join_test1:
	mpirun -np 8 Multi_join 4 4 "../test_commande/test_fb2.txt"

multi_join_test2:
	mpirun -np 4 Multi_join 3 3 "../test_commande/test_fb1.txt"

multi_join_test3:
	mpirun -np 8 Multi_join 3 3 "../test_commande/test_tw.txt"

multi_join_test4:
	mpirun -np 8 Multi_join 3 3 "../test_commande/test_dblp.txt"

multi_join_test5:
	mpirun -np 4 Multi_join 6 4 "../test_commande/test_sm3.txt"

hypercube_test1:
	mpirun -np 8 HyperCube 3 3 3 2 4 1 "../test_commande/test_fb1.txt"

hypercube_test2:
	mpirun -np 8 HyperCube 3 3 3 2 2 2 "../test_commande/test_tw.txt"

hypercube_test3:
	mpirun -np 8 HyperCube 3 3 3 2 2 2 "../test_commande/test_dblp.txt"

hypercube_test4:
	mpirun -np 16 HyperCube 4 4 2 2 2 2 "../test_commande/test_fb2.txt"

hypercube_test5:
	mpirun -np 16 HyperCube 4 4 2 2 2 2 "../test_commande/test_sm1.txt"

hypercube_test6:
	mpirun -np 32 HyperCube 5 5 2 2 2 2 2 "../test_commande/test_sm2.txt"

hypercube_test7:
	mpirun -np 16 HyperCube 6 4 2 2 2 2 "../test_commande/test_sm3.txt"

create_data: create_data.cpp
	$(CC) create_data.cpp -o create_data

clean:
	rm -f Distr_join create_data HyperCube Relation Multi_join test_sequential *.o

