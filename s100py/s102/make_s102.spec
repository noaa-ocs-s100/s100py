# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['utils.py'],
             pathex=['C:\\PydroXL_19_Dev\\envs\\nomkl2\\Lib\\site-packages\\s100py\\s102'],
             binaries=[],
             datas=[('C:\\PydroXL_19_Dev\\envs\\nomkl2\\Library\\share\\proj', 'Library\\share\\proj')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='make_s102',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
