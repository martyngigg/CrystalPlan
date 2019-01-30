# -*- mode: python -*-

block_cipher = None


a = Analysis(['crystalplan.py'],
             pathex=['C:\\Users\mantid\git\CrystalPlan'],
             binaries=[],
             datas=[("instruments", "instruments"), ("gui\icons", "gui\icons")],
             hiddenimports=['_sysconfigdata', 'scipy._lib.messagestream', 'scipy.linalg.cython_blas', 'scipy.linalg.cython_lapack', 'scipy.linalg._fblas', 'tvtk'],
             hookspath=['hooks'],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='crystalplan',
          debug=True,
          strip=False,
          upx=True,
          console=True,
          icon='gui\icons\CrystalPlan_icon.ico'
         )

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               console=True,
               debug=True,
               name='crystalplan')