.PHONY: package

test_lcm/out.ps: binaries/linux/lcmodel
	cd test_lcm && \
	../binaries/linux/lcmodel < control.file

package: binaries/linux/lcmodel.xz
binaries/linux/lcmodel.xz: binaries/linux/lcmodel
	xz -k $^

binaries/linux/lcmodel: | binaries/linux
	gfortran -c -fno-backslash -fno-f2c -O3 -fall-intrinsics -std=legacy -Wuninitialized -ffpe-summary=none source/LCModel.f
	gfortran LCModel.o -o binaries/lcmodel

binaries/linux:
	mkdir -p $@
