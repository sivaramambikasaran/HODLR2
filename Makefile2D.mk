CC			=/usr/local/bin/g++-9
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I/usr/local/include #-I/usr/local/opt/icu4c/include #-ffast-math -Ofast -ffinite-math-only
LDFLAGS		=-fopenmp -std=c++17 -I/usr/local/include #-L/usr/local/opt/icu4c/lib #-Ofast
SOURCES		=./testHODLR2D.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testHODLR2D

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(KERNEL) $(HOMOG) $< -o $@

clean:
	rm a.out testHODLR2D *.o
