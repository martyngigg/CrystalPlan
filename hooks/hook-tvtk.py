from PyInstaller.utils.hooks import copy_metadata, collect_submodules, collect_data_files

# datas = copy_metadata('tvtk')
datas = collect_data_files('tvtk')
hiddenimports = collect_submodules('tvtk')
