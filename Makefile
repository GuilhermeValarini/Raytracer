$(MAKE)= make

exec:
	echo "Running on CPU with no parallel code"
	$(MAKE) -C theone-cpu exec
	echo "Running on CPU with OpenMP anotations"
	$(MAKE) -C theone-omp exec
	echo "Running on GPU"
	$(MAKE) -C theone-gpu exec

clean:
	echo "Cleanning all images and compiled files"
	$(MAKE) -C theone-cpu clean
	$(MAKE) -C theone-omp clean
	$(MAKE) -C theone-gpu clean
