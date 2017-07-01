CC = g++

default: All

Relation.o: Relation.hpp Relation.cpp
	$(CC) -c Relation.cpp -o Relation.o
	
utils.o: utils.hpp utils.cpp Relation.o
	$(CC) -c utils.cpp -o utils.o

Seq_join: utils.o Relation.o Seq_join.cpp
	$(CC) Relation.o utils.o Seq_join.cpp -o Seq_join

Distr_join: utils.o Distr_join.cpp Relation.o
	mpic++ Relation.o utils.o Distr_join.cpp -o Distr_join

Hypercube: utils.o Relation.o HyperCube.cpp
	mpic++ Relation.o utils.o HyperCube.cpp -o HyperCube

Multi_join: utils.o Relation.o Multi_join.cpp
	mpic++ Relation.o utils.o Multi_join.cpp -o Multi_join

All: Relation.hpp Relation.cpp utils.hpp utils.cpp Seq_join.cpp Distr_join.cpp Multi_join.cpp HyperCube.cpp
	$(CC) -c Relation.cpp -o Relation.o
	$(CC) -c utils.cpp -o utils.o
	$(CC) Relation.o utils.o Seq_join.cpp -o Seq_join
	mpic++ Relation.o utils.o Multi_join.cpp -o Multi_join
	mpic++ Relation.o utils.o Distr_join.cpp -o Distr_join
	mpic++ Relation.o utils.o HyperCube.cpp -o HyperCube	

seq_join_test:
	./Seq_join "../data/data.txt" 2 0 1 "../data/data.txt" 2 1 2

seq_join_fb:
	./Seq_join "../data/facebook.dat" 2 0 1 "../data/facebook.dat" 2 1 2

distr_join_test_loc:
	mpirun -np 4 Distr_join "../data/data.txt" 2 0 1 "../data/data.txt" 2 1 2

distr_join_fb_loc: 
	mpirun -np 4 Distr_join "../data/facebook.dat" 2 0 1 "../data/facebook.dat" 2 1 2

distr_join_tw_loc:
	mpirun -np 4 Distr_join "../data/twitter.dat" 2 0 1 "../data/twitter.dat" 2 1 2

distr_join_db_loc:
	mpirun -np 4 Distr_join "../data/dblp.dat" 2 0 1 "../data/dblp.dat" 2 1 2

distr_join_fb: 
	salloc --ntasks=4 --ntasks-per-node=1 mpirun Distr_join "../data/facebook.dat" 2 0 1 "../data/facebook.dat" 2 1 2

distr_join_tw:
	salloc --ntasks=4 --ntasks-per-node=1 mpirun Distr_join "../data/twitter.dat" 2 0 1 "../data/twitter.dat" 2 1 2

distr_join_db:
	salloc --ntasks=4 --ntasks-per-node=1 mpirun Distr_join "../data/dblp.dat" 2 0 1 "../data/dblp.dat" 2 1 2

multi_join_fb_loc:
	mpirun -np 6 Multi_join 3 3 "../test_command/fb_triangle.txt"

multi_join_tw_loc:
	mpirun -np 6 Multi_join 3 3 "../test_command/tw_triangle.txt"

multi_join_db_loc:
	mpirun -np 6 Multi_join 3 3 "../test_command/db_triangle.txt"

multi_join_fb:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun Multi_join 3 3 "../test_command/fb_triangle.txt"

multi_join_tw:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun Multi_join 3 3 "../test_command/tw_triangle.txt"

multi_join_db:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun Multi_join 3 3 "../test_command/db_triangle.txt"

hypercube_fb:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun HyperCube 3 3 2 2 2 "../test_command/fb_triangle.txt"

hypercube_tw:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun HyperCube 3 3 2 2 2 "../test_command/tw_triangle.txt"

hypercube_db:
	salloc --ntasks=8 --ntasks-per-node=1 mpirun HyperCube 3 3 2 2 2 "../test_command/db_triangle.txt"

hypercube_fb4:
	salloc --ntasks=64 --ntasks-per-node=2 mpirun HyperCube 6 4 2 4 2 4 "../test_command/fb_4complet.txt"

hypercube_db4:
	salloc --ntasks=64 --ntasks-per-node=2 mpirun HyperCube 6 4 2 4 2 4 "../test_command/db_4complet.txt"

clean:
	rm -f Relation Seq_join Distr_join Multi_join HyperCube *.o

