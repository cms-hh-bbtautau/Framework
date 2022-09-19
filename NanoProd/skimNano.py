import os
if os.path.exists(os.path.join(os.path.dirname(__file__), "BaselineSelection.py")):
  import BaselineSelection as Baseline
else:
  import Common.BaselineSelection as Baseline

def _def_met_cuts(met_thr, met_collections):
  cut = ' || '.join([ f'{v}_pt > {met_thr}' for v in met_collections ])
  return f"( {cut} )"

met_cuts = _def_met_cuts(80, ["MET", "DeepMETResolutionTune", "DeepMETResponseTune", "PuppiMET" ])

def skim_B0(df):
  Baseline.Initialize()
  df = Baseline.ApplyRecoBaseline0(df)
  return df

def skim_B0B1(df):
  Baseline.Initialize()
  df = Baseline.ApplyRecoBaseline0(df)
  df = Baseline.ApplyRecoBaseline1(df)
  return df

def skim_failed_B0B1(df):
  Baseline.Initialize()
  df, b0_filter = Baseline.ApplyRecoBaseline0(df, apply_filter=False)
  df, b1_filter = Baseline.ApplyRecoBaseline1(df, apply_filter=False)
  df = df.Filter(f'!(({b0_filter}) && ({b1_filter}))')
  return df

def skim_B0B1_MET(df):
  Baseline.Initialize()
  df, b0_filter = Baseline.ApplyRecoBaseline0(df, apply_filter=False)
  df, b1_filter = Baseline.ApplyRecoBaseline1(df, apply_filter=False)
  df = df.Filter(f'(({b0_filter}) && ({b1_filter})) || {met_cuts}')
  return df
def skim_failed_B0B1_MET(df):
  Baseline.Initialize()
  df, b0_filter = Baseline.ApplyRecoBaseline0(df, apply_filter=False)
  df, b1_filter = Baseline.ApplyRecoBaseline1(df, apply_filter=False)
  df = df.Filter(f'!( (({b0_filter}) && ({b1_filter})) || {met_cuts} )')
  return df
