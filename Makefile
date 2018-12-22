SHELL=/bin/bash

.PHONY: test unzip

.INTERMEDIATE: skdef1.test skdef2.test skdef3.test sktp1.test sktp2.test sktd1.test sktd2.test skodd1.test skodd2.test qcgett.test qcgetp.test

all: supcrt92 dprons92.dat

cprons92: cprons92.for
	gfortran $^ -o $@

dprons92.dat: sprons92.dat cprons92
	rm -f $@
	echo -e '$<\n$@' | ./cprons92 >/dev/null

supcrt92: sup92d.f rep92d.f reac92d.f h2o92d.for
	gfortran $^ -o $@

unzip: supcrt92.zip
	unzip $^
	ls | grep -v Makefile | while read -r f; do mv "$$f" "$$(echo "$$f" | tr '[:upper:]' '[:lower:]')"; done

supcrt92.zip:
	curl 'https://pages.uoregon.edu/palandri/data/Supcrt92.zip' > $@

test: all skdef1.test skdef2.test skdef3.test sktp1.test sktp2.test sktd1.test sktd2.test skodd1.test skodd2.test qcgett.test qcgetp.test
	# the header is different between the files, so skip it
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' skdef1.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' skdef1.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' skdef2.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' skdef2.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' skdef3.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' skdef3.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' sktp1.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' sktp1.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' sktp2.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' sktp2.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' sktd1.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' sktd1.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' sktd2.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' sktd2.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' skodd1.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' skodd1.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' skodd2.tab) \
		<(sed -n '/ REACTION  1 /,$$ p' skodd2.test)
	# the remaining two tests have slight differences in output (what was
	# simply 0.0 can become -0.0). We use some regex magic to convert negative
	# zeros to normal zeros.
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' qcgett.tab) \
		<(sed -n 's/ -0\.\(0*\) / 0.\1 /g; / REACTION  1 /,$$ p;' qcgett.test)
	diff --ignore-space-change \
		<(sed -n '/ REACTION  1 /,$$ p' qcgetp.tab) \
		<(sed -n 's/ -0\.\(0*\) / 0.\1 /g; / REACTION  1 /,$$ p;' qcgetp.test)

skdef1.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n1\n1\n1\nskarn.rxn\n$@\n1' | ./supcrt92 >/dev/null

skdef2.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n1\n2\n1\nskarn.rxn\n$@\n1' | ./supcrt92 >/dev/null

skdef3.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n1\n3\n1\nskarn.rxn\n$@\n1' | ./supcrt92 >/dev/null

sktp1.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nsktp1.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

sktp2.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nsktp2.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

sktd1.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nsktd1.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

sktd2.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nsktd2.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

skodd1.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nskodd1.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

skodd2.test: skarn.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nskodd2.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

qcgett.test: qc.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nqcgett.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

qcgetp.test: qc.rxn supcrt92 dprons92.dat
	echo -e 'y\n2\nqcgetp.con\n1\n$<\n$@\n1' | ./supcrt92 >/dev/null

clean:
	rm -f cprons92 dprons92.dat supcrt92 skdef1.test skdef2.test skdef3.test sktp1.test sktp2.test sktd1.test sktd2.test skodd1.test skodd2.test qcgett.test qcgetp.test
