CC = g++
CFLAGS = -g -Wall -std=c++17

.PHONY: all

all: FiniteAbelianClassifier

FiniteAbelianClassifier: FiniteAbelianClassifier.cpp
	$(CC) $(CFLAGS) FiniteAbelianClassifier.cpp -o FiniteAbelianClassifier