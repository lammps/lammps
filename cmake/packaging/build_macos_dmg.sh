#!/bin/bash

APP_NAME=lammps-gui
VERSION="$1"
LAMMPS_GUI_APP="$2"
BUILD_DIR="${PWD}"
PACKAGING_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE_DIR="${BUILD_DIR}/dmg-staging"
DMG_FILE="LAMMPS-macOS-multiarch-GUI-${VERSION}.dmg"

# install/upgrade dmgbuild helper script
python3 -m pip install --upgrade --user pip
python3 -m pip install --upgrade --user dmgbuild

if ! python3 -c 'import dmgbuild' > /dev/null 2>&1
then
    echo "ERROR: dmgbuild is required. Install with: python3 -m pip install --user dmgbuild"
    exit 1
fi

rm -rv ${APP_NAME}.app
mv -v ${LAMMPS_GUI_APP} .

echo "Delete old files, if they exist"
rm -f ${APP_NAME}.dmg ${APP_NAME}-rw.dmg LAMMPS-macOS-multiarch-GUI-*.dmg
rm -rf "${STAGE_DIR}"

echo "Force ad hoc signing of dynamic LAMMPS library and LAMMPS-GUI"
codesign --force -s - ${BUILD_DIR}/liblammps.0.dylib
codesign --force -s - ${BUILD_DIR}/lammps-gui.app/Contents/MacOS/lammps-gui
codesign --force -s - ${BUILD_DIR}/lammps-gui.app/Contents/Frameworks/liblammps.0.dylib

echo "Bundle Qt frameworks and plugins with macdeployqt"
macdeployqt ${APP_NAME}.app

echo "Stage a copy of the app bundle plus README and background image"
mkdir -p "${STAGE_DIR}"
ditto ${APP_NAME}.app "${STAGE_DIR}/LAMMPS-GUI.app"
pushd "${STAGE_DIR}"
mv LAMMPS-GUI.app/Contents/Resources/README.txt .
mv LAMMPS-GUI.app/Contents/Resources/LAMMPS_DMG_Background.png background.png
cd LAMMPS-GUI.app/Contents

echo "Update rpath for LAMMPS to link to the bundled liblammps.0.dylib copy"
install_name_tool -delete_rpath ${BUILD_DIR} bin/lmp
install_name_tool -add_rpath '@executable_path/../Frameworks' bin/lmp

echo "Codesign bundled plugins"
codesign --force -s - PlugIns/*/*.dylib

echo "Codesign bundled frameworks"
codesign --force -s - Frameworks/Qt*.framework/Versions/A/Qt*

echo "Codesign bundled executables"
for s in bin/*
do \
    test "$s" = "bin/ffmpeg" && continue
    test "$s" = "bin/lammps-gui" && continue
    test -f $s && codesign --force -s - $s
done
codesign --force -s - MacOS/lammps-gui

echo "Codesign status for LAMMPS-GUI and bundle"
codesign -v --verbose=4 MacOS/lammps-gui

echo "Attach icons to LAMMPS console and GUI executables and lib"
echo "read 'icns' (-16455) \"Resources/lammps.icns\";" > icon.rsrc
Rez -a icon.rsrc -o bin/lmp
SetFile -a C bin/lmp
echo "read 'icns' (-16455) \"Resources/lammps-gui.icns\";" > icon.rsrc
Rez -a icon.rsrc -o MacOS/lammps-gui
SetFile -a C MacOS/lammps-gui
if [ -f Frameworks/liblammps.0.dylib ]; then
    Rez -a icon.rsrc -o Frameworks/liblammps.0.dylib
    SetFile -a C Frameworks/liblammps.0.dylib
fi
rm icon.rsrc
popd

echo "Create compressed disk image using dmgbuild"
python3 -m dmgbuild -s "${PACKAGING_DIR}/dmg_settings.py" \
    -D app="${STAGE_DIR}/LAMMPS-GUI.app" \
    -D readme="${STAGE_DIR}/README.txt" \
    -D background="${STAGE_DIR}/background.png" \
    "${APP_NAME}" "${DMG_FILE}"

echo "Attach icon to .dmg file"
echo "read 'icns' (-16455) \"${APP_NAME}.app/Contents/Resources/lammps.icns\";" > icon.rsrc
Rez -a icon.rsrc -o ${DMG_FILE}
SetFile -a C ${DMG_FILE}
rm icon.rsrc

echo "Delete staging directory"
rm -rf "${STAGE_DIR}"

exit 0
