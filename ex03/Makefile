CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wshadow -Werror -fopenmp -O3 -DNDEBUG -ftree-vectorizer-verbose=1 -ffast-math -march=native

TARGET =rbgs
HXX=RBGS.h Timer.h

OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) $(HXX) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET)
