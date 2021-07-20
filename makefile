CXX = g++
TARGET = mocca
SOURCES = mocca.cpp aperture.cpp ray.cpp
OBJS = $(SOURCES:.cpp=.o)
CXXFLAGS = -g -Wall
INCLUDES = -I "eigen/"

all: $(TARGET)

$(TARGET): $(OBJS)
		$(CXX) $(INCLUDES) $(OBJS) -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

clean:
		rm *.o mocca
