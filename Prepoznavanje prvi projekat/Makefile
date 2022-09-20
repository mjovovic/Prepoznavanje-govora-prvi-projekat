
default: compile run

compile:
	#mora export DISPLAY=:0 u terminalu
	g++ main.c -o main -lm -lsndfile -I/usr/include/python2.7 -lpython2.7  -I kissfft kissfft/kiss_fft.c 
run:
	./main

clean:
	rm main