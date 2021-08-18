# LCModel compiled binaries for Windows & Mac

In this repository, I am collecting readily compiled and executable binaries for the [LCModel software](http://s-provencher.com/lcmodel.shtml), a popular method for linear-combination modeling of magnetic resonance spectroscopy data.

## Quick download links

[MacOS Catalina](LCModel_macos_catalina.7z)

[Windows 10](LCModel_win10.7z)

## Contents of this repository

- `binaries` contains the compiled LCModel executables, sorted by operating system (first level) and version (second level). As the compiled executables are just a little over 100 MB in file size (the maximum allowed size on GitHub if we don't want to use GitLFS), they are compressed into [.7z format, which you can open with the free, open-source 7-Zip tool](https://www.7-zip.org/).
- `source` contains the LCModel Fortran77 source code, as downloaded from the [LCModel website](http://s-provencher.com/pub/LCModel/source.zip).
- `test_lcm` contains a test dataset that Dr. Martin Wilson (https://github.com/martin3141) has shared along with the first compilation instructions in an [MRSHub forum thread](https://forum.mrshub.org/t/building-lcmodel/317).
- `README.md` is the file you are reading.
- `LICENSE.md` contains the licensing agreement.

## Compilation procedures

This is how the executables in the `binaries` folder were generated

### Windows 10

Windows 10 Enterprise Edition

GFortran 11.1.0 from http://www.equation.com/ftpdir/gcc/gcc-11.1.0-64.exe

Enter at Windows command prompt:

`gfortran.exe -ffpe-summary=none -std=legacy -O3 LCModel.f -o LCModel.exe`

### MacOS Catalina

Mac OS Catalina 10.15.7

GFortran 10.2 from https://github.com/fxcoudert/gfortran-for-macOS/releases/download/10.2/gfortran-10.2-Catalina.dmg

Enter at Terminal prompt:

```
gfortran -c -fno-backslash -fno-f2c -O3 -fall-intrinsics -std=legacy -Wuninitialized -ffpe-summary=none LCModel.f
gfortran LCModel.o -o lcmodel
```

## Testing procedures

Copy the compiled binary into the `test_lcm` folder and run

`./lcmodel < control.file` (Mac)

`lcmodel < control.file` (Windows)

Then, compare the output (`out.ps`) with the gold-standard output (`out_ref_build.ps`).
Apart from the date and the version (`6.3-N` vs. `6.3-R`), the output should be identical.

I will add a more sophisticated testing routine in the future to ensure full numerical agreement.

## Limitations

This repository *only* contains the core LCModel executable. None of the other LCModel tools (LCMgui, makebasis etc.) are available.

There are multiple open-source tools available that include functions that will help you with the generation of LCModel basis and control files.
- [Osprey](https://github.com/schorschinho/osprey)
- [FID-A](https://github.com/cic-methods/fid-a)
- [FSL-MRS](https://github.com/wexeee/fsl_mrs)
- [spant](https://github.com/martin3141/spant)

Please contact the developers of these software packages for help.

## Contact, feedback, suggestions

Please contact me for any sort of feedback about the binaries I have provided here.
If you have questions about how to use LCModel, please refer to the [official manual](http://s-provencher.com/pub/LCModel/manual/manual.pdf), or submit a question on the [MRSHub LCModel forum](https://forum.mrshub.org/c/mrs-software/lcmodel/8).

## Acknowledgements

The LCModel source code was made available free of charge by Dr. Stephen Provencher. Please see the licensing agreement for details.
