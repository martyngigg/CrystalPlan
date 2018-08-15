from PyInstaller.utils.hooks import copy_metadata, collect_submodules, collect_data_files

datas = copy_metadata('traitsui')
datas += collect_data_files('traitsui', include_py_files=True)
hiddenimports = collect_submodules('traitsui')
hiddenimports += [
    'traitsui.wx.range_editor',
    'traitsui.wx.null_editor',
]
