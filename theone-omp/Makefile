EXEC=raytracer
INPUT=$(abspath $(sort $(wildcard ../input/*.in)))
OUTPUT=$(basename $(notdir $(INPUT)))

CXX=g++
FLAGS=-fopenmp -Wno-unused-result

RM=rm -f

HEADER="Resolution,Width,Height,Number_spheres,Run_time"
RESFILE=results.csv

all:
	$(CXX) $(FLAGS) $(EXEC).cpp -o $(EXEC)

clean:
	$(RM) $(EXEC) *.ppm *.csv

exec:clean all run

run:
	$(RM) $(RESFILE)
	echo $(HEADER) >> $(RESFILE)
	for file in $(INPUT); do y=$${file%.in}; x=$${y##*/}; echo Testing $$x.in; ./$(EXEC) $$file $$x.ppm >> $(RESFILE); done
