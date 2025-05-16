# -*- mode: python ; coding: utf-8 -*-
import os
from pathlib import Path
from PyInstaller.utils.hooks import collect_submodules
import plotly

# Get the absolute path to the src/faster directory
faster_dir = Path('src/faster')

# Collect all config files
config_files = [
    ('src/faster/config/marker_info.json', 'config'),
    ('src/faster/config/variant_catalog.thermofisher_24markers.json', 'config')
]

# Collect all binary files
bin_files = [
    ('src/faster/bin/ExpansionHunter', 'bin')
]

# Collect all data files
datas = config_files + bin_files

plotly_package_data = (str(Path(plotly.__file__).parent / 'package_data'), 'plotly/package_data')
datas.append(plotly_package_data)

hidden_plotly_graph_objs = collect_submodules('plotly.graph_objs')
hidden_plotly_validators = collect_submodules('plotly.validators')

a = Analysis(
    ['src/faster/gui.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[
        'pandas',
        'plotly',
        'tkinter',
        'numpy',
        'json',
        'logging',
        'threading',
        'pathlib',
    ] + hidden_plotly_graph_objs + hidden_plotly_validators,
    hookspath=['.'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='FASTER-GUI',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='src/faster/assets/icon.ico' if os.path.exists('src/faster/assets/icon.ico') else None,
)
