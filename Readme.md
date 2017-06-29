# Distributed join processing on social network data

## Structure
This project consists of 3 folders:
/data: includes data files: facebook.data twitter.dat dblp.dat data.txt(made manually)
/src: includes all header and source code files
/test_command: includes arguments for executing programs

In the root directory, a Makefile is provided to compile and execute programs.

## Compiling
For compiling programs, run: make All.
Alternatively, for compiling each algorithm individually, run:
make Seq_join
make Distr_join
make Multi_join
make Hypercube

## Testing
For testing the performance, we provides Makefile items for each algorithm and each data file.

### Testing sequential algorithm for 4 data files:
make seq_join_test
make seq_join_fb

### Testing distributed join algorithm with single machine
make distr_join_test_loc
make distr_join_fb_loc
make distr_join_tw_loc
make distr_join_db_loc

### Testing distributed join algorithm with machine cluster
make distr_join_fb
make distr_join_tw
make distr_join_db

### Testing multi way join algorithm with single machine
make multi_join_fb_loc
make multi_join_tw_loc
make multi_join_db_loc

### Testing multi way join algorithm with machine cluster
make multi_join_fb
make multi_join_tw
make multi_join_db

### Testing hypercube join algorithm with machine cluster
make hypercube_fb
make hypercube_tw
make hypercube_db

### Testing searching complete subgraph of 4 vertices
make hypercube_fb4
make hypercube_db4

### Delete all compiled files
make clean