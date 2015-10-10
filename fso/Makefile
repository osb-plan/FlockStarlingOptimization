# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS= -Wall -std=c++11 -O3

all: 
	$(CC) $(CFLAGS) -o fso_exec main.cpp fso.cpp -I./
clean:
	rm  fso_exec global_best.txt position1.txt  velocity.txt 

