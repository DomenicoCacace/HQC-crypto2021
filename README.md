# Cryptography and Architectures for Computer Security Project

## Introduction
This work describes a masking-based side channel resistent implementation of HQC, a third round alternative
candidate for the NIST Post-Quantum Cryptography competition. 

The implementation on this repository has been tested on
an ARM Cortex-M4 processor, showing a significant decrease in terms of leaked information, while still keeping a relatively small overhead.

## Contents
This repository contains:


- report: the documentation about the work done, the presentation and the benchmark results
- scripts: python and bash scripts to automate the code generation, testing and analysis processes
- src: C implementation of the cryptosystem
 
## Installation
In order to compile and run HQC, the following softwares are needed: <code>cmake, make, stm32programmer, gcc-arm-none-eabi</code>.

This implementation is tailored for the STM32F401RE board; to change the target, you need to generate the configuration files on STM32CodeMX and add them in the <code>src/stm32</code> folder.

### Compilation
We use CMake to manage all the executables we generate; assuming to be in the <code>build</code> folder:

<code>
cmake .. -DSECLVL=X -DMASKLVL=Y -DMODE="MODE" -DCROSSCOMPILE=CROSS -DVERBOSE=VERB
</code>


<list>
  <li>X: security level (128, 192, 256)
  <li>Y: number of shares of the masking scheme (1, 2, 3, 4)
  <li>MODE: the executable to be compiled (<code>CONST-KEM, CONST-PKE, TIMING-KEM, TIMING-PKE, FUNCTIONAL</code>)
    <li> CROSS: 1 to compile for the stm32 board, 0 for the native architecture
    <li> VERB: the verbosity level of the log messages (1, 2)
</list>
