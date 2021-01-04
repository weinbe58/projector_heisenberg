

current_dir = $(shell pwd)
CXX=g++
FLAGS=--std=c++17 -O3 -I$(current_dir)/../../

chain : chain.o
	g++ chain.o -o chain.exe



%.o: %.cxx
	$(CXX) -g $(FLAGS) -c $<


