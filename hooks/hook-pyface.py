from PyInstaller.utils.hooks import copy_metadata, collect_submodules, collect_data_files

datas = copy_metadata('pyface')
datas += collect_data_files('pyface', include_py_files=True)
hiddenimports = collect_submodules('pyface')
