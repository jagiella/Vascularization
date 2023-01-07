//#define USE_GUI
#ifdef USE_GUI

#include "Perfusion.hpp"
#include <QApplication>
#include "window.hpp"

// SIMULATION THREAD

#include <QThread>
#include <QIcon>
class SimulationThread : public QThread
  {
  public:
     VesselGraph *vg;
     Window *window;
     int argc;
     char **argv;

     void run(){
        fprintf(stderr, "TEST\n");
        perfusion(argc, argv, vg);
     }
  };


int main(int argc, char** argv)
{
	fprintf(stderr, "WITH GUI\n");

	QApplication app(argc, argv);

    VesselGraph *newVG = new VesselGraph();


    // SET UP SIMULATION THREAD
    SimulationThread *simThread = new SimulationThread();
	setlocale(LC_ALL,"en_US");
    simThread->argc = argc;
    simThread->argv = argv;
    simThread->vg = newVG;
    simThread->start();

    // SET UP GUI
    Window *gui = new Window();
    gui->vg = newVG;
    gui->show();

    fprintf(stderr, "START APP\n");
    return app.exec();

    //perfusion(argc, argv);
}

#else


#include "Perfusion.hpp"
int main(int argc, char** argv)
{
	fprintf(stderr, "NO GUI\n");
    perfusion(argc, argv,0);
}

#endif
