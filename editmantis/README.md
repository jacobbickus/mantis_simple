# mantisworking
work in progress version of mantis

make using cmake
To Run:
  mantis >out.out

Should make a cherenkov.root file with 4 root trees
Very possible a seg fault occurs. If seg fault doesn't occur all of the G4cout will be in out.out file
mantis defaults to read from the mantis.in macro which can be adjusted to change the number of particles to run.
Should only take a few seconds to run default settings.

Changes Added:

Should be able to print G4cout to out file automatically now.
