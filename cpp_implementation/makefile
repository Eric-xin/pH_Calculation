all:
	g++ -std=c++11 -o exec/CBE src/CBE.cpp
	g++ -std=c++11 -o exec/PBE src/PBE.cpp
	ln -sf exec/CBE CBE
	ln -sf exec/PBE PBE

clean:
	rm exec/CBE exec/PBE
	rm CBE PBE
	rm io_test/CBE.sam.tmp.out io_test/PBE.sam.tmp.out

check:
	@echo "Checking files..."
	@ls src

test:
	@echo "Checking files..."
	@ls io_test
	@echo "Testing CBE..."
	./CBE < io_test/CBE.sam.in > io_test/CBE.sam.tmp.out
	@echo "Testing PBE..."
	./PBE < io_test/PBE.sam.in > io_test/PBE.sam.tmp.out