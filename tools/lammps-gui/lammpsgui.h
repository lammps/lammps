#ifndef LAMMPSGUI_H
#define LAMMPSGUI_H

#include <QMainWindow>
#include <QFont>
#include <QString>

QT_BEGIN_NAMESPACE
namespace Ui { class LammpsGui; }
QT_END_NAMESPACE

class LammpsGui : public QMainWindow
{
    Q_OBJECT

public:
    LammpsGui(QWidget *parent = nullptr);
    ~LammpsGui();

private slots:
    void new_document();
    void open();
    void save();
    void save_as();
    void quit();
    void copy();
    void cut();
    void paste();
    void undo();
    void redo();
    void clear();
    void run_buffer();
    void run_line();
    void about();

private:
    Ui::LammpsGui *ui;
    QString current_file;
    QFont text_font;
    int current_line;
    void *lammps_handle;
};
#endif // LAMMPSGUI_H
