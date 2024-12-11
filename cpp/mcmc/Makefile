CC = g++
CFLAGS = -O2 -Wall -std=gnu++14
LDFLAGS =
INCLUDES =
CLEANC=rm -rf *.o *.out

SOURCES = main.cpp LoadPrevalence.cpp McmcObject.cpp Period.cpp Episode.cpp utils.cpp  
HEADERS = LoadPrevalence.h McmcObject.cpp Period.h Episode.h utils.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = program.out

covid: $(EXECUTABLE) $(SOURCES) $(HEADERS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(INCLUDES) $(OBJECTS) -o $@ -lstdc++fs

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	$(CLEANC)
