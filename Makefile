CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall

dac_sub_delaunay: dac_sub_delaunay.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

test: dac_sub_delaunay
	./dac_sub_delaunay --gen sphere --dim 4 --n 300 --s 3 --eps 0.9 --beta 1 --verify
	./dac_sub_delaunay --gen torus --n 600 --s 2 --eps 1.0 --beta 0.5 --verify
	./dac_sub_delaunay --gen random --dim 3 --n 250 --s 3 --eps 0.4 --beta 0 --verify

clean:
	rm -f dac_sub_delaunay

.PHONY: test clean
