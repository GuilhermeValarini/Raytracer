EXEC=raytracer
INPUT=$(sort $(wildcard *.in))
OUTPUT=$(wildcard *.ppm)

CXX=g++
FLAGS=-O3 -fopenmp -Wno-unused-result

RM=rm -f

GEN=gen_scene
PI=python

all: $(EXEC)

$(EXEC):
	$(CXX) $(FLAGS) $(EXEC).cpp -o $(EXEC)

clean:
	$(RM) $(EXEC)
	for file in $(OUTPUT); do $(RM) $$file; done

exec:clean all
	for file in $(INPUT); do ./$(EXEC) $$file "$${file%%.*}.ppm"; done

rm_test:
	for file in $(INPUT); do $(RM) $$file; done

gen_test:
	$(PI) $(GEN).py

new_test: rm_test gen_test
