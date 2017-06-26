$(MAKE)= make
$(RM)=rm

all:
	echo "Compiling on CPU with no parallel code"
	$(MAKE) -C theone-cpu all
	echo "Compiling on CPU with OpenMP anotations"
	$(MAKE) -C theone-omp all
	echo "Compiling on GPU"
	$(MAKE) -C theone-gpu all

run:
	echo "Running on CPU with no parallel code"
	$(MAKE) -C theone-cpu run
	echo "Running on CPU with OpenMP anotations"
	$(MAKE) -C theone-omp run
	echo "Running on GPU"
	$(MAKE) -C theone-gpu run


clean:
	echo "Cleanning all images and compiled files"
	$(RM) *.csv
	$(MAKE) -C theone-cpu clean
	$(MAKE) -C theone-omp clean
	$(MAKE) -C theone-gpu clean
