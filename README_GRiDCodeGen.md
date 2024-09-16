# Adding New Algorithms to GRiDCodeGenerator

GRiDCodeGenerator provides a framework to write GPU based algorithms for robots with different amounts of limbs/joints. 

## Using crba as an example algorithm, here are the steps to take: 

1. Create a python file for the algorithm, ```generateGRiD.py``` uses ```_crba.py``` and ```GRiDCodeGenerator.py``` to create ```grid.cuh``` 
2. To run ```generateGRiD.py``` use the command: 
```python
python3 generateGRiD.py iiwa_simple.urdf
``` 
```iiwa_simple.urdf``` can be replaced with any URDF parser you wish to test. This is helpful because ```grid.cuh``` will look different depending on the URDF because the python code in ```_crba.py``` will do different things based on the amount of bfs levels etc 
3. Create a CUDA C++ file ```test_crba.cu``` in order to test the crba algorithm from ```grid.cuh``` 
4. To run ```grid.cuh``` use the command: 
```python
nvcc test_crba.cu -o test
./test
```
5. ```_crba.p```y consists of different functions that are all used in GRiDCodeGenerator so be sure to add them to the algorithms section also be sure to add the ```CRBA_SHARED_MEM_COUNT``` param to ```gen_add_constant_helpers()``` and don’t forget to add the ```gen_crba``` function call to gen_all_code 
6. This ensures that ```GRiDCodeGenerator.py``` can call ```_crba.py``` and all of the functions in it
7. Import your algorithm to ```__init__.py``` in algorithms.  
8. Add the CPU/ GPU inputs and outputs to ```gen_init_gridData```  

##  Useful Functions in ```_crba.py```:

1. ```Gen_crba_inner_temp_mem_size``` sets the size of all of the temporary variables within the function 
2. ```Gen_crba_inner_function_call``` includes all of the input variables including gravity and temp and generates the function call
3. ```Gen_crba_inner``` contains everything that will be outputted to ```grid.cuh``` 
+  Function call, parameters and parameter notes, the actual algorithm
+  Use the helper functions to add code lines or to add loops 
+  The code outputted must be correct CUDA C++ code so be mindful of syntax and punctuation
4. ```Gen_crba_device_temp_mem_size``` keeps track of ```Gen_crba_inner_temp_mem_size``` + wrapper size of XI matrices 
5. ```Gen_crba_device``` loads the XI matrices and calls the ```Gen_crba_inner_function_call```
6. ```Gen_crba_kernel``` keeps track of temp memory and XI matrices and calls ```Gen_crba_inner_function_call```
7. ```Gen_crba_host``` calls ```Gen_crba_kernel``` and keeps track of allocating and freeing memory 
8. ```Gen_crba``` calls ```Gen_crba_inner```, ```Gen_crba_device```, ```Gen_crba_kernel``` and ```Gen_crba_host``` 
+ This is the function that is being called by ```GRiDCodeGenerator``` so it called all of the other functions and wrappers 
9. Host calls kernel which calls inner function and wrappers for inner function come from device 

## Useful Functions in ```helpers/_code_generation_helpers.py```
1. ```Gen_add_code_line``` adds the line of code, and indents for the next lines if necessary (when ```add_indent_after``` is true)
2. ```Gen_add_code_lines``` adds multiple lines of code
3. ```Gen_add_end_control_flow``` decreases the indent and adds closing bracket
4. ```Gen_add_end_function``` used the end of the function
5. ```Gen_add_func_doc``` adds function doc
6. ```Gen_add_parallel_loop``` this is how you add a parallel loop to the code. 
7. ```Var_name``` gives the name of the for loop variable, max_val gives the max value
8. ```Gen_add_serial_ops``` adds an if statement for only one thread (this is very helpful for print statements while debugging)
9. ```Gen_add_sync``` adds a sync for the threads
10. ```Gen_add_multi_threaded_select``` adds an if statement

## Switch out all of “crba” for the name of the algorithm you are working on. 
