CXX = mpic++
CXXFLAGS = -std=c++11
CXXFLAGSP = -std=c++11 -Wall -Wextra -Wshadow -O3 -DNDEBUG  -ffast-math -march=native
CXXFLAGSVECTO = -ftree-vectorizer-verbose=1

TARGET =cg
HXX=Vector.h Timer.h

OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) $(HXX) 
	$(CXX) $(CXXFLAGS) $(CXXFLAGSP) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGSP) -c $<

clean:
	rm -rf *.o $(TARGET)
