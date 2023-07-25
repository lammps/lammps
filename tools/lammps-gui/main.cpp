#include "lammpsgui.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    const char *infile = nullptr;
    if (argc > 1) infile = argv[1];

    LammpsGui w(nullptr, infile);
    w.show();
    return a.exec();
}

// Local Variables:
// c-basic-offset: 4
// End:
