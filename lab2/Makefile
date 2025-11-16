CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
DEBUGFLAGS = -std=c++17 -g -O0 -Wall -Wextra
TARGET = main.exe
SOURCES = main.cpp SupportFunctions.cpp Binet.cpp Strassen.cpp Inverse.cpp LUfactorization.cpp GaussElimination.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean run batch debug parrallel

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
	./$(TARGET) 4 sizes.txt strassen.txt

parallel: $(TARGET)
	@seq 1 6 | parallel './$(TARGET) {} sizes.txt {}.txt'

debug: CXXFLAGS = $(DEBUGFLAGS)
debug: clean all

main.o: main.cpp SupportFunctions.h Mnozenie.h Binet.h Strassen.h
SupportFunctions.o: SupportFunctions.cpp SupportFunctions.h
Binet.o: Binet.cpp Binet.h Mnozenie.h SupportFunctions.h
Strassen.o: Strassen.cpp Strassen.h Mnozenie.h SupportFunctions.h
Inverse.o: Inverse.cpp Inverse.h Mnozenie.h SupportFunctions.h
LUfactorization.o: LUfactorization.cpp LUfactorization.h Mnozenie.h SupportFunctions.h
GaussElimination.o: GaussElimination.cpp GaussElimination.h Mnozenie.h SupportFunctions.h