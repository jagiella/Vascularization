all:
	g++ Source/main.cpp Source/Perfusion.cpp Source/VesselGraph.cpp Source/Tumor.cpp ../libnix/src/IO/EPS.cpp ../libnix/src/IO/Povray.cpp -ISource -I../libnix/src -I/opt/local/include/ -L/usr/local/lib/ -lnix