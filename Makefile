$(MAKE)= make

exec:
	echo "Running on CPU with no parallel code"
	$(MAKE) -C theone-cpu exec
	echo "Running on CPU with OpenMP anotations"
	$(MAKE) -C theone-omp exec
	echo "Running on GPU"
	$(MAKE) -C theone-gpu exec
