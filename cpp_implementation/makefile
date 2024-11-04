all:
	g++ -std=c++11 -o exec/CBE src/CBE.cpp
	g++ -std=c++11 -o exec/PBE src/PBE.cpp
	ln -sf exec/CBE CBE
	ln -sf exec/PBE PBE

clean:
	rm exec/CBE exec/PBE
	rm CBE PBE