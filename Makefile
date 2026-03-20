CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic

all: simplify

simplify: simplify.cpp
	$(CXX) $(CXXFLAGS) simplify.cpp -o simplify

clean:
	rm -f simplify