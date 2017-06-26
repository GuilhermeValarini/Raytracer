$(MAKE)= make
$(RM)=rm

exec:
	$(RM) *.csv
	echo "Running on CPU with no parallel code"
	$(MAKE) -C theone-cpu exec
	cp	./theone-cpu/results.csv ./result-cpu.csv
	echo "Running on CPU with OpenMP anotations"
	$(MAKE) -C theone-omp exec
	cp	./theone-omp/results.csv ./result-omp.csv
	echo "Running on GPU"
	$(MAKE) -C theone-gpu exec
	cp	./theone-gpu/results.csv ./result-gpu.csv

clean:
	echo "Cleanning all images and compiled files"
	$(MAKE) -C theone-cpu clean
	$(MAKE) -C theone-omp clean
	$(MAKE) -C theone-gpu clean
