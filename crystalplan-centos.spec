# -*- mode: python -*-

block_cipher = None


a = Analysis(['crystalplan.py'],
             pathex=['/CrystalPlan'],
             binaries=[("/usr/local/lib/libiomp5.so", ".")],
             datas=[("instruments", "instruments"), ("gui/icons", "gui/icons")],
             hiddenimports=['scipy._lib.messagestream', 'pubsub.core.publisherbase', 'tvtk'],
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
          debug=False,
          strip=False,
          upx=True,
          console=False,
          icon='gui/icons/CrystalPlan_icon.ico'
         )

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='crystalplan')
