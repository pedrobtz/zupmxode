# zupmxode: Population Pharmacokinetic and Pharmacodynamic Simulation Without a Compiler

Simulates population pharmacokinetic and pharmacodynamic models written
in a subset of the 'rxode2' model language over 'NONMEM'-style event
records, with no C compiler on the machine. Models are parsed in R and
compiled to native code through 'MIR', or run by its interpreter where
native code cannot be made. Closed-form compartment models are solved
exactly, linear systems by a matrix exponential, and the rest by 'LSODA'
or Dormand-Prince; subjects are solved in parallel with 'OpenMP'.

## See also

Useful links:

- <https://pedrobtz.github.io/zupmxode/>

- <https://github.com/pedrobtz/zupmxode>

- Report bugs at <https://github.com/pedrobtz/zupmxode/issues>

## Author

**Maintainer**: Pedro Baltazar <pedrobtz@gmail.com> \[copyright holder\]

Authors:

- Pedro Baltazar <pedrobtz@gmail.com> \[copyright holder\]
