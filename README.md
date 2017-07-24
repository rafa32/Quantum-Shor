To compile the C program, use the command:

 gcc  -o  <program_name>  shor.c  -lm

-lm simply links the math library (with trigonometry functions, exponents, etc.)

To run:
./<program_name>

The program will ask for a number N, to be factored. This number must be a product of two different primes (for reasons given in the report's Conclusions).
The last output of the program should be the two found factors of N. It will print some progress while running.
