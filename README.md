# Planetary Ephemeris Program (PEP)

To build:

```
$ git clone https://gitlab.com/jbattat/pep_core.git
$ cd pep_core
$ make
```

This compiles the code in the `pep/`, `peputil/`, and `verify/` subdirectories
and creates the binary files in the `bigtest/` subdirectory necessary for
running the validation script `bigtest`.

To run `bigtest` to verify (part of) the installation:
```
$ make bigtest
```

Caution: the output from this validation procedure is compared with the
standard output included in this package and distilled to ignore trivial
differences in time stamps, home directories, and processor speeds, but
discrepancies due to hardware and compiler differences may appear to be
distressingly profuse.  The condensed output is a set of files named
`t*.verout`, which will each consist of 1116 bytes if this package is
built on exactly the same platform as the original, but may be much
larger in some cases.  In such cases, expert help must be enlisted to
confirm the success of the build.

Then, to "install" the code (symlinking into `~/bin/`), run:
```
$ make install
```
You can undo this (delete the symlinks in `~/bin`) with `$ make uninstall`.  
If you prefer to symlink into a different directory, set the TARGET_DIR parameter:
```
$ make TARGET_DIR=~/my/other/dir install
```
note, that you also have to specify the directory by hand for `uninstall`:
```
$ make TARGET_DIR=~/my/other/dir uninstall
```

If you wish to remove the `.o` files and other detrius from the build:  
```
$ make clean
```
If you also wish to remove the executables, then run `make cleanall`.
