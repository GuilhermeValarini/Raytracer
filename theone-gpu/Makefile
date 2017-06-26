EXEC=raytracer
INPUT=$(sort $(wildcard *.in))
OUTPUT=$(wildcard *.ppm)

CXX=nvcc
FLAGS=-O3 -g -G -gencode arch=compute_50,code=sm_50

RM=rm -f

GEN=gen_scene
PI=python

HEADER="Resolution, Width, Height, Number spheres, Maximum depth, Run time"
RESFILE=results.csv

all: $(EXEC)

$(EXEC):
	$(CXX) $(FLAGS) $(EXEC).cu -o $(EXEC)

clean:
	$(RM) $(EXEC) $(RESFILE)
	for file in $(OUTPUT); do $(RM) $$file; done

exec:clean all
	echo $(HEADER) >> $(RESFILE)
	for file in $(INPUT); do echo Testing $$file; ./$(EXEC) $$file "$${file%%.*}.ppm" >> $(RESFILE); done

rm_test:
	for file in $(INPUT); do $(RM) $$file; done

gen_test:
	$(PI) $(GEN).py

new_test: rm_test gen_test
