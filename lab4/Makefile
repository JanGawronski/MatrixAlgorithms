CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -I../eigen
LDFLAGS = 

TARGET = hmatrix
TEST_VECTOR = test_vector
SIMPLE_TEST = simple_test
VISUALIZE_EXAMPLE = visualize_example

# Main sources (including Compression.cpp from lab3)
SOURCES = main.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Test vector sources
TEST_VECTOR_SOURCES = test_vector.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
TEST_VECTOR_OBJECTS = $(TEST_VECTOR_SOURCES:.cpp=.o)

# Simple test sources
SIMPLE_TEST_SOURCES = simple_test.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
SIMPLE_TEST_OBJECTS = $(SIMPLE_TEST_SOURCES:.cpp=.o)

# Visualization example sources
VISUALIZE_EXAMPLE_SOURCES = visualize_example.cpp HMatrix.cpp SupportFunctions.cpp Compression.cpp
VISUALIZE_EXAMPLE_OBJECTS = $(VISUALIZE_EXAMPLE_SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

$(TEST_VECTOR): $(TEST_VECTOR_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TEST_VECTOR) $(TEST_VECTOR_OBJECTS) $(LDFLAGS)

$(SIMPLE_TEST): $(SIMPLE_TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(SIMPLE_TEST) $(SIMPLE_TEST_OBJECTS) $(LDFLAGS)

$(VISUALIZE_EXAMPLE): $(VISUALIZE_EXAMPLE_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(VISUALIZE_EXAMPLE) $(VISUALIZE_EXAMPLE_OBJECTS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TEST_VECTOR_OBJECTS) $(SIMPLE_TEST_OBJECTS) $(VISUALIZE_EXAMPLE_OBJECTS) \
	      $(TARGET) $(TEST_VECTOR) $(SIMPLE_TEST) $(VISUALIZE_EXAMPLE) *.txt *.exe *.png

run: $(TARGET)
	./$(TARGET)

test_vector: $(TEST_VECTOR)
	./$(TEST_VECTOR)

simple_test: $(SIMPLE_TEST)
	./$(SIMPLE_TEST)

visualize: $(VISUALIZE_EXAMPLE)
	./$(VISUALIZE_EXAMPLE)

.PHONY: all clean run test_vector simple_test visualize
