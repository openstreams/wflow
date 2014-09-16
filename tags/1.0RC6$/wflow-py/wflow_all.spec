# -*- mode: python -*-
adapt = Analysis(['wflow\\wflow_adapt.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
			 
wave = Analysis(['wflow\\wflow_wave.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
sbm = Analysis(['wflow\\wflow_sbm.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
hbv = Analysis(['wflow\\wflow_hbv.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
gr4 = Analysis(['wflow\\wflow_gr4.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
floodmap = Analysis(['wflow\\wflow_floodmap.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
plottss = Analysis(['wflow\\plottss.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

prep1 = Analysis(['Scripts\\wflow_prepare_step1.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
prep2 = Analysis(['Scripts\\wflow_prepare_step2.py'],
             pathex=['d:\\repos\\Hydrology\\OpenStreams\\src\\wflow-py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)			 

			 
MERGE ((plottss,'plottss','plottss'),(adapt,'wflow_adapt','wflow_adapt'),(wave,'wflow_wave','wflow_wave'),(sbm,'wflow_sbm','wflow_sbm'),(hbv,'wflow_hbv','wflow_hbv'),(gr4,'wflow_gr4','wflow_gr4'),(floodmap,'wflow_floodmap','wflow_floodmap'),(prep1,'wflow_prepare_step1','wflow_prepare_step1'),(prep2,'wflow_prepare_step2','wflow_prepare_step2'))			 


adapt_pyz = PYZ(adapt.pure)
wave_pyz = PYZ(wave.pure)
sbm_pyz = PYZ(sbm.pure)
hbv_pyz = PYZ(hbv.pure)
gr4_pyz = PYZ(gr4.pure)
floodmap_pyz = PYZ(floodmap.pure)
plottss_pyz = PYZ(plottss.pure)
prep1_pyz = PYZ(prep1.pure)
prep2_pyz = PYZ(prep2.pure)


adapt_exe = EXE(adapt_pyz,  adapt.scripts, exclude_binaries=True, name='wflow_adapt.exe',debug=False,strip=None, upx=True, console=True )
adapt_coll = COLLECT(adapt_exe, adapt.binaries,  adapt.zipfiles,    adapt.datas,   strip=None,  upx=True, name='wflow_adapt')

wave_exe = EXE(wave_pyz,  wave.scripts, exclude_binaries=True, name='wflow_wave.exe',debug=False,strip=None, upx=True, console=True )
wave_coll = COLLECT(wave_exe, wave.binaries,  wave.zipfiles,    wave.datas,   strip=None,  upx=True, name='wflow_wave')

sbm_exe = EXE(sbm_pyz,  sbm.scripts, exclude_binaries=True, name='wflow_sbm.exe',debug=False,strip=None, upx=True, console=True )
sbm_coll = COLLECT(sbm_exe, sbm.binaries,  sbm.zipfiles,    sbm.datas,   strip=None,  upx=True, name='wflow_sbm')

hbv_exe = EXE(hbv_pyz,  hbv.scripts, exclude_binaries=True, name='wflow_hbv.exe',debug=False,strip=None, upx=True, console=True )
hbv_coll = COLLECT(hbv_exe, hbv.binaries,  hbv.zipfiles,    hbv.datas,   strip=None,  upx=True, name='wflow_hbv')

gr4_exe = EXE(gr4_pyz,  gr4.scripts, exclude_binaries=True, name='wflow_gr4.exe',debug=False,strip=None, upx=True, console=True )
gr4_coll = COLLECT(gr4_exe, gr4.binaries,  gr4.zipfiles,    a=gr4.datas,   strip=None,  upx=True, name='wflow_gr4')

floodmap_exe = EXE(floodmap_pyz,  floodmap.scripts, exclude_binaries=True, name='wflow_floodmap.exe',debug=False,strip=None, upx=True, console=True )
floodmap_coll = COLLECT(floodmap_exe, floodmap.binaries,  floodmap.zipfiles,    floodmap.datas,   strip=None,  upx=True, name='wflow_floodmap')

plottss_exe = EXE(plottss_pyz,  plottss.scripts, exclude_binaries=True, name='plottss.exe',debug=False,strip=None, upx=True, console=True )
plottss_coll = COLLECT(plottss_exe, plottss.binaries,  plottss.zipfiles,    plottss.datas,   strip=None,  upx=True, name='plottss')

prep1_exe = EXE(prep1_pyz,  prep1.scripts, exclude_binaries=True, name='wflow_prepare_step1.exe',debug=False,strip=None, upx=True, console=True )
prep1_coll = COLLECT(prep1_exe, prep1.binaries,  prep1.zipfiles,    prep1.datas,   strip=None,  upx=True, name='prep1')

prep2_exe = EXE(prep2_pyz,  prep2.scripts, exclude_binaries=True, name='wflow_prepare_step2.exe',debug=False,strip=None, upx=True, console=True )
prep2_coll = COLLECT(prep2_exe, prep2.binaries,  prep2.zipfiles,    prep2.datas,   strip=None,  upx=True, name='prep2')

