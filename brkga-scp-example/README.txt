To build the executable run "make". IMPORTANT: You must have CUDA installed!!

It will be created a folder build/ and inside it the executable cuda-scp.
To test the program you can execute 

	./build/cuda-scp ./instances/orlib/scp41.txt


To set up the BRKGA parameters, edit file config.txt

The file ./config.txt has runtime configuration parameters for cuda-tsp. By default it
was set decode type to be done in the CPU (HOST_DECODE,decode_type 1). 