#include "lammpsgui.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    LammpsGui w;
    w.show();
    return a.exec();
}
