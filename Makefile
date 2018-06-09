default:
	nvcc *.cu *.cpp -lineinfo -g -maxrregcount 64 -arch=sm_35 --std=c++11 -O3 --compiler-options "-Wall -Wextra"
gcc:
	g++ *.cpp -g --std=c++11 -O3 -Wall -Wextra -DGCC 
