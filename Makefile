CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2
TARGET   := simplify
SRCS     := main.cpp ring.cpp priority_queue.cpp topology.cpp io.cpp
OBJS     := $(SRCS:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)

# Dependencies (headers)
main.o:          main.cpp ring.h priority_queue.h topology.h io.h
ring.o:          ring.cpp ring.h
priority_queue.o: priority_queue.cpp priority_queue.h ring.h
topology.o:      topology.cpp topology.h ring.h
io.o:            io.cpp io.h ring.h
