# dmgbuild settings for the LAMMPS-GUI drag-and-drop installer disk image.
#
# Used by build_macos_dmg.sh. Replaces the former Finder AppleScript, so the
# image can be built without scripting the Finder (e.g. over ssh or in CI).
# dmgbuild writes the window layout directly into the image's .DS_Store.
#
# Paths are passed from the build script with -D:
#   app=...         staged LAMMPS-GUI.app bundle
#   readme=...      README.txt to place next to the app
#   background=...  background image (1024x768 px at 96 dpi = 768x576 pt)

import os.path

for _key in ("app", "readme", "background"):
    if _key not in defines:
        raise ValueError(f"dmg_settings.py: missing -D {_key}=<path>")

_app = defines["app"]
_readme = defines["readme"]

# compressed read-only image, same as the former 'hdiutil convert -format UDZO'
format = "UDZO"
filesystem = "HFS+"

files = [_app, _readme]
symlinks = {"Applications": "/Applications"}
background = defines["background"]

# Same window geometry as the former AppleScript 'bounds {100, 40, 868, 640}',
# with sidebar, toolbar and status bar hidden.
window_rect = ((100, 100), (768, 600))
default_view = "icon-view"
show_status_bar = False
show_tab_view = False
show_toolbar = False
show_pathbar = False
show_sidebar = False
sidebar_width = 0

# Icon view: free placement, 64pt icons, Finder's default 12pt labels
arrange_by = None
icon_size = 64
text_size = 12
show_icon_preview = False

# Icon centers, same as the former AppleScript 'set position of item ...'
icon_locations = {
    os.path.basename(_app.rstrip("/")): (190, 216),
    "Applications": (576, 216),
    os.path.basename(_readme): (190, 400),
}
