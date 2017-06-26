import csv

fieldnamesOutput = ["Filename","Resolution","Width","Height","Number_spheres","Run_time","Speedup"]

cpuFilename = "./theone-cpu/results.csv"
ompFilename = "./theone-omp/results.csv"
gpuFilename = "./theone-gpu/results.csv"
cpuOutFilename = "results-cpu.csv"
ompOutFilename = "results-omp.csv"
gpuOutFilename = "results-gpu.csv"

with open(cpuFilename) as cpuCSV:
    cpuReader = csv.DictReader(cpuCSV)
    with open(ompFilename) as ompCSV:
        ompReader = csv.DictReader(ompCSV)
        with open(gpuFilename) as gpuCSV:
            gpuReader = csv.DictReader(gpuCSV)

            cpuTimes = []

            with open(cpuOutFilename, 'w') as cpuOut:
                cpuWriter = csv.DictWriter(cpuOut, fieldnames=fieldnamesOutput)
                cpuWriter.writeheader();
                print "<<<------------------------- Serial ------------------------->>>"
                for cpu in cpuReader:
                    cpuTimes.append(cpu["Run_time"])
                    speedup = 1.0
                    cpu["Speedup"] = "{:.2f}".format(speedup)
                    print "---->>> {}".format(cpu["Filename"])
                    print "Speedup: {}".format(cpu["Speedup"])
                    cpuWriter.writerow(cpu)

            with open(ompOutFilename, 'w') as ompOut:
                ompWriter = csv.DictWriter(ompOut, fieldnames=fieldnamesOutput)
                ompWriter.writeheader();
                print "<<<------------------------- OpenMP ------------------------->>>"
                i = 0
                for omp in ompReader:
                    speedup = float(cpuTimes[i]) / float(omp["Run_time"])
                    omp["Speedup"] = "{:.2f}".format(speedup)
                    print "---->>> {}".format(omp["Filename"])
                    print "Speedup: {}".format(omp["Speedup"])
                    ompWriter.writerow(omp)
                    i+=1

            with open(gpuOutFilename, 'w') as gpuOut:
                gpuWriter = csv.DictWriter(gpuOut, fieldnames=fieldnamesOutput)
                gpuWriter.writeheader();
                print "<<<------------------------- OpenMP ------------------------->>>"
                i = 0
                for gpu in gpuReader:
                    speedup = float(cpuTimes[i]) / float(gpu["Run_time"])
                    gpu["Speedup"] = "{:.2f}".format(speedup)
                    print "---->>> {}".format(gpu["Filename"])
                    print "Speedup: {}".format(gpu["Speedup"])
                    gpuWriter.writerow(gpu)
                    i+=1
