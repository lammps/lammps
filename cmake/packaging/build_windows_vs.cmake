# CMake script to be run post installation to build zipped package

# clean up old zipfile and deployment tree
file(REMOVE LAMMPS_GUI-Win10-amd64.zip)
file(REMOVE_RECURSE LAMMPS_GUI)
file(RENAME ${INSTNAME} LAMMPS_GUI)

# move all executables and dlls to main folder and delete bin folder
file(GLOB BINFILES LIST_DIRECTORIES FALSE LAMMPS_GUI/bin/*.exe LAMMPS_GUI/bin/*.dll)
foreach(bin ${BINFILES})
  get_filename_component(exe ${bin} NAME)
  file(RENAME ${bin} LAMMPS_GUI/${exe})
endforeach()
file(REMOVE_RECURSE LAMMPS_GUI/bin)

# create qt.conf so Qt will find its plugins
file(WRITE LAMMPS_GUI/qt.conf "[Paths]\r\nPlugins = qt5plugins\r\n")

# initialize environment and then run windeployqt to populate folder with missing dependencies and Qt plugins
file(WRITE qtdeploy.bat "@ECHO OFF\r\nset VSCMD_DEBUG=0\r\nCALL ${VC_INIT} x64\r\nset PATH=${QT5_BIN_DIR};%PATH%\r\nwindeployqt --plugindir LAMMPS_GUI/qt5plugins --release  LAMMPS_GUI/lammps-gui.exe --no-quick-import --no-webkit2 --no-translations --no-system-d3d-compiler --no-angle --no-opengl-sw\r\n")
execute_process(COMMAND cmd.exe /c qtdeploy.bat COMMAND_ECHO STDERR)
file(REMOVE qtdeploy.bat)

# download and uncompress static FFMpeg and gzip binaries
file(DOWNLOAD "https://download.lammps.org/thirdparty/ffmpeg-gzip.zip" ffmpeg-gzip.zip)
file(WRITE unpackzip.ps1 "Expand-Archive -Path ffmpeg-gzip.zip -DestinationPath LAMMPS_GUI")
execute_process(COMMAND powershell -ExecutionPolicy Bypass -File unpackzip.ps1)
file(REMOVE unpackzip.ps1)
file(REMOVE ffmpeg-gzip.zip)

# create zip archive
file(WRITE makearchive.ps1 "Compress-Archive -Path LAMMPS_GUI -CompressionLevel Optimal -DestinationPath LAMMPS_GUI-Win10-amd64.zip")
execute_process(COMMAND powershell -ExecutionPolicy Bypass -File makearchive.ps1)
file(REMOVE makearchive.ps1)
file(REMOVE_RECURSE LAMMPS_GUI)
