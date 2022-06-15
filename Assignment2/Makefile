CXX = g++
OBJECTS = main.o 
FLAGS = -O3 

exec:  $(OBJECTS)
	$(CXX) $(FLAGS) liblbfgs/lbfgs.c main.cpp -o main
	
clean :
	rm $(OBJECTS) main