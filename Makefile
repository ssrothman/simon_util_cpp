test1: test1.o combinatorics.o iterating.o
	g++ $^ -o $@

%.o: %.cc
	g++ -c -o $@ $^ -I/usr/local/include/Minuit2 

clean:
	find . -maxdepth 1 -type f -executable -exec rm {} +
	rm -f *.o
