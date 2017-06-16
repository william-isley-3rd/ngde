Authors: William C Isley III

Heavily modified from version written and published by: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, May 14, 2005.

** NGDE.exe is the executable file that one should use on a windows machine. Note that the code needs to be recompiled if one wants to make changes to it.

This C code is used to solve for nucleation, coagulation and surface growth problems. This example problem solves for characteristics of an Aluminum aerosol. However, it can be used
for any other material whose properties are known. There are four different sections in this code:

1) Coagulation, 
2) Nucleation (classical) + coagulation, 
3) Surface growth 
4) Unified GDE with all the three phenomena combined. 

The above phenomena constitute four different sections of the code. Each of the four sections are independent of each other. They have especially been written in such manner so 
that, one can identify the contribution of each phenomena to the GDE. Also, if one requires to use nucleation, coagulation or surface growth alone, it can be easily done.   

The solution algorithm involves a new approach to solve the Aerosol GDE, where the particle volume space is divided into nodes of zero width. Using nodes allows us to cover the 
entire size range of aerosol (1 nm to 10 micrometer) using just 40 nodes.A more detailed description of the theory and algorithm can be found in the reference:

The main theoretical constraints of this implementation are:
1. Use of the Self Consistent Classical nucleation Model
2. Free molecule collision kernel
3. Free molecule surface growth

To compile the code "ngde.c" the following command should be used to compile the code : "gcc ngde.c -lm"

To run the code one can modify the input files as directed and then use the command "a.out" to run the code. Using the command "a.out > outputfile" would direct the output to a file 
called "outputfile".


COMMENTS ON INPUT FILES:

The input file "main.inp" should contain all the property data, as shown in the example input file. If the user wants to run the code for a specific coagulation or a surface growth 
problem, the corresponding input files "coag.inp"/"grow.inp" must be modified accordingly.


IMPORTANT NOTE:

The input file format should not be modified unless the user is willing to make the corresponding changes in the code. Only the numerical entities after the colon (":") may be changed as desired.

Finally, the code has been well commented and making changes by the user should not be difficult.


