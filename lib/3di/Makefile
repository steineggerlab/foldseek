CPP = g++
# FLAGS = -03
FLAGS = -g -Wall

LIBS = -Iinclude -I.

test: main
	echo "./example/d1qd6.1.ent" | ./main

main: main.o transform.o load.o
	$(CPP) main.o load.o transform.o -o main

%.o: %.cpp
	$(CPP) $(FLAGS) $(LIBS) -c $<

clean:
	rm *.o
	rm main
