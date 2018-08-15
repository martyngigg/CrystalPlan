from PyInstaller.utils.hooks import copy_metadata, collect_submodules, collect_data_files

datas = copy_metadata('mayavi')
datas += collect_data_files('mayavi')
hiddenimports = collect_submodules('mayavi')
