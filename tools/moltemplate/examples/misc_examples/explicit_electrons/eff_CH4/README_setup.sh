cd moltemplate_files/

  moltemplate.sh -atomstyle electron system_ch4.lt

  # This will generate 3 files:
  #  "system_ch4.data", "system_ch4.in.init", "system_ch4.in.settings"

  mv -f system_ch4.data system_ch4.in.init system_ch4.in.settings ../

  # optional: delete temporary files

  rm -rf output_ttree

cd ../
