# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc

# Include paths...
Release_Include_Path= -I/usr/include -Iinclude/ANN -Iinclude/NPN -Iinclude -I/usr/include/c++/4.7.2

# Library paths...
Release_Library_Path= -L/usr/lib/boost 


# Additional libraries...
Release_Libraries=  

# Preprocessor definitions...
Release_Preprocessor_Definitions=-D GCC_BUILD -D NDEBUG -D _CONSOLE 

# Implictly linked object files...
Release_Implicitly_Linked_Objects=

# Compiler flags...
Release_Compiler_Flags= -O3 -Wall  
# extra optimization -O3 -march=native

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Release

# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/src/ANN/ANN.o gccRelease/src/ANN/brute.o gccRelease/src/ANN/checkActiveClass3.o gccRelease/src/ANN/kd_dump.o gccRelease/src/ANN/kd_fix_rad_search.o gccRelease/src/ANN/kd_pr_search.o gccRelease/src/ANN/kd_search.o gccRelease/src/ANN/kd_split.o gccRelease/src/ANN/kd_tree.o gccRelease/src/ANN/kd_util.o gccRelease/src/ANN/perf.o gccRelease/src/ANN/checkActiveClass.o gccRelease/src/ANN/checkActiveClass1.o gccRelease/src/ANN/checkActiveClass2.o gccRelease/src/NPN/DisjointSets.o gccRelease/src/NPN/NPN.o gccRelease/src/NPN/myTiming.o 
	g++ gccRelease/src/ANN/ANN.o gccRelease/src/ANN/brute.o gccRelease/src/ANN/checkActiveClass3.o gccRelease/src/ANN/kd_dump.o gccRelease/src/ANN/kd_fix_rad_search.o gccRelease/src/ANN/kd_pr_search.o gccRelease/src/ANN/kd_search.o gccRelease/src/ANN/kd_split.o gccRelease/src/ANN/kd_tree.o gccRelease/src/ANN/kd_util.o gccRelease/src/ANN/perf.o gccRelease/src/ANN/checkActiveClass.o gccRelease/src/ANN/checkActiveClass1.o gccRelease/src/ANN/checkActiveClass2.o gccRelease/src/NPN/DisjointSets.o gccRelease/src/NPN/NPN.o gccRelease/src/NPN/myTiming.o  $(Release_Library_Path) $(Release_Libraries) -Wl,-rpath,./ -o ../bin/Linux/NPN.o

# Compiles file src/ANN/ANN.cpp for the Release configuration...
-include gccRelease/src/ANN/ANN.d
gccRelease/src/ANN/ANN.o: src/ANN/ANN.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/ANN.cpp $(Release_Include_Path) -o gccRelease/src/ANN/ANN.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/ANN.cpp $(Release_Include_Path) > gccRelease/src/ANN/ANN.d


# Compiles file src/ANN/brute.cpp for the Release configuration...
-include gccRelease/src/ANN/brute.d
gccRelease/src/ANN/brute.o: src/ANN/brute.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/brute.cpp $(Release_Include_Path) -o gccRelease/src/ANN/brute.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/brute.cpp $(Release_Include_Path) > gccRelease/src/ANN/brute.d

# Compiles file src/ANN/checkActiveClass3.cpp for the Release configuration...
-include gccRelease/src/ANN/checkActiveClass3.d
gccRelease/src/ANN/checkActiveClass3.o: src/ANN/checkActiveClass3.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/checkActiveClass3.cpp $(Release_Include_Path) -o gccRelease/src/ANN/checkActiveClass3.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/checkActiveClass3.cpp $(Release_Include_Path) > gccRelease/src/ANN/checkActiveClass3.d

# Compiles file src/ANN/kd_dump.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_dump.d
gccRelease/src/ANN/kd_dump.o: src/ANN/kd_dump.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_dump.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_dump.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_dump.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_dump.d

# Compiles file src/ANN/kd_fix_rad_search.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_fix_rad_search.d
gccRelease/src/ANN/kd_fix_rad_search.o: src/ANN/kd_fix_rad_search.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_fix_rad_search.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_fix_rad_search.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_fix_rad_search.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_fix_rad_search.d

# Compiles file src/ANN/kd_pr_search.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_pr_search.d
gccRelease/src/ANN/kd_pr_search.o: src/ANN/kd_pr_search.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_pr_search.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_pr_search.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_pr_search.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_pr_search.d

# Compiles file src/ANN/kd_search.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_search.d
gccRelease/src/ANN/kd_search.o: src/ANN/kd_search.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_search.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_search.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_search.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_search.d

# Compiles file src/ANN/kd_split.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_split.d
gccRelease/src/ANN/kd_split.o: src/ANN/kd_split.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_split.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_split.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_split.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_split.d

# Compiles file src/ANN/kd_tree.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_tree.d
gccRelease/src/ANN/kd_tree.o: src/ANN/kd_tree.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_tree.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_tree.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_tree.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_tree.d

# Compiles file src/ANN/kd_util.cpp for the Release configuration...
-include gccRelease/src/ANN/kd_util.d
gccRelease/src/ANN/kd_util.o: src/ANN/kd_util.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/kd_util.cpp $(Release_Include_Path) -o gccRelease/src/ANN/kd_util.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/kd_util.cpp $(Release_Include_Path) > gccRelease/src/ANN/kd_util.d

# Compiles file src/ANN/perf.cpp for the Release configuration...
-include gccRelease/src/ANN/perf.d
gccRelease/src/ANN/perf.o: src/ANN/perf.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/perf.cpp $(Release_Include_Path) -o gccRelease/src/ANN/perf.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/perf.cpp $(Release_Include_Path) > gccRelease/src/ANN/perf.d

# Compiles file src/ANN/checkActiveClass.cpp for the Release configuration...
-include gccRelease/src/ANN/checkActiveClass.d
gccRelease/src/ANN/checkActiveClass.o: src/ANN/checkActiveClass.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/checkActiveClass.cpp $(Release_Include_Path) -o gccRelease/src/ANN/checkActiveClass.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/checkActiveClass.cpp $(Release_Include_Path) > gccRelease/src/ANN/checkActiveClass.d

# Compiles file src/ANN/checkActiveClass1.cpp for the Release configuration...
-include gccRelease/src/ANN/checkActiveClass1.d
gccRelease/src/ANN/checkActiveClass1.o: src/ANN/checkActiveClass1.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/checkActiveClass1.cpp $(Release_Include_Path) -o gccRelease/src/ANN/checkActiveClass1.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/checkActiveClass1.cpp $(Release_Include_Path) > gccRelease/src/ANN/checkActiveClass1.d

# Compiles file src/ANN/checkActiveClass2.cpp for the Release configuration...
-include gccRelease/src/ANN/checkActiveClass2.d
gccRelease/src/ANN/checkActiveClass2.o: src/ANN/checkActiveClass2.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ANN/checkActiveClass2.cpp $(Release_Include_Path) -o gccRelease/src/ANN/checkActiveClass2.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/ANN/checkActiveClass2.cpp $(Release_Include_Path) > gccRelease/src/ANN/checkActiveClass2.d

# Compiles file src/NPN/DisjointSets.cpp for the Release configuration...
-include gccRelease/src/NPN/DisjointSets.d
gccRelease/src/NPN/DisjointSets.o: src/NPN/DisjointSets.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/NPN/DisjointSets.cpp $(Release_Include_Path) -o gccRelease/src/NPN/DisjointSets.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/NPN/DisjointSets.cpp $(Release_Include_Path) > gccRelease/src/NPN/DisjointSets.d

# Compiles file src/NPN/myTiming.cpp for the Release configuration...
-include gccRelease/src/NPN/myTiming.d
gccRelease/src/NPN/myTiming.o: src/NPN/myTiming.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/NPN/myTiming.cpp $(Release_Include_Path) -o gccRelease/src/NPN/myTiming.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/NPN/myTiming.cpp $(Release_Include_Path) > gccRelease/src/NPN/myTiming.d

# Compiles file src/NPN/NPN.cpp for the Release configuration...
-include gccRelease/src/NPN/NPN.d
gccRelease/src/NPN/NPN.o: src/NPN/NPN.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/NPN/NPN.cpp $(Release_Include_Path) -o gccRelease/src/NPN/NPN.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/NPN/NPN.cpp $(Release_Include_Path) > gccRelease/src/NPN/NPN.d


# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccRelease/src
	mkdir -p gccRelease/src
	mkdir -p gccRelease/src/ANN
	mkdir -p gccRelease/src/NPN
	mkdir -p ../bin/Linux

# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f -r gccRelease/src/ANN/*.o
	rm -f -r gccRelease/src/ANN/*.d
	rm -f -r gccRelease/src/NPN/*.o
	rm -f -r gccRelease/src/NPN/*.d
	rm -f -r ../bin/Linux/*.a
	rm -f -r ../bin/Linux/*.so
	rm -f -r ../bin/Linux/*.dll
	rm -f -r ../bin/Linux/*.exe
