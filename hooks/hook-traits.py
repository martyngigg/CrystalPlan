from PyInstaller.utils.hooks import copy_metadata, collect_submodules, collect_data_files

datas = copy_metadata('traits')
datas += collect_data_files('traits', include_py_files=True)
hiddenimports = collect_submodules('traits')
