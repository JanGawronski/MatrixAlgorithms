CXX = g++
CXXFLAGS = -std=c++23 -O3 -Wall -Wextra
DEBUGFLAGS = -std=c++23 -g -O0 -Wall -Wextra
TARGET = compression
SOURCES = main.cpp SupportFunctions.cpp Compression.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean run debug batch

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

run: $(TARGET)
	./$(TARGET)

batch: $(TARGET)
	@mkdir -p output
	@parallel --colsep ' ' './$(TARGET) {1} {2} > output/svd_r={1}_e={2}.txt' :::: input.txt

debug: CXXFLAGS = $(DEBUGFLAGS)
debug: clean all