CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET = matmul.exe
SOURCES = main.cpp SupportFunctions.cpp Binet.cpp Strassen.cpp AI.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean run batch

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) wynik.txt

run: $(TARGET)
	./$(TARGET)

batch: $(TARGET)
	./$(TARGET) sizes.txt 1

# Dependencies (optional, helps with incremental builds)
main.o: main.cpp SupportFunctions.h Mnozenie.h Binet.h Strassen.h AI.h
SupportFunctions.o: SupportFunctions.cpp SupportFunctions.h
Binet.o: Binet.cpp Binet.h Mnozenie.h SupportFunctions.h
Strassen.o: Strassen.cpp Strassen.h Mnozenie.h SupportFunctions.h
AI.o: AI.cpp AI.h Mnozenie.h SupportFunctions.h