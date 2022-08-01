EXEC:=final_project
HOSTFILE:=hostfile.txt

all:
	mpicc -fopenmp project.c -o $(EXEC)
	
run:
	mpiexec -np 2 ./$(EXEC)
	
run_on2:
	mpiexec -np 2 -hostfile $(HOSTFILE) ./$(EXEC)
	
clean:
	rm $(EXEC)
