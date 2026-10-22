.PHONY: all test clean

all: test

# Deletes the environment's Cargo lock first: rnaseq/FINDINGS.md R27.
test:
	rm -f $${MORLOC_STATE:-/opt/morloc-state}/cache/rust-env.lock
	morloc make -o test test.loc
	./test | tee test.out
	tail -n 1 test.out | grep -qx true

clean:
	rm -rf test test-build test.out
