#!/bin/bash

# Move scripts to lib and link to them from bin.
for script in *.py; do
  mv utils/$script $PREFIX/lib
  ln -s ../lib/$script $PREFIX/bin
done