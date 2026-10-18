#!/usr/bin/env python3
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import zipfile
import requests
from shapely.geometry import box

ROOT = Path(os.environ.get('CISC_ROOT', '/proj/sebastianlab/users/hquintal/cisc'))
DATA_OUT = Path(os.environ.get('CISC_DATA_OUT', str(ROOT / 'data/out')))
CLUSTER_BASE = Path(os.environ.get('CISC_CLUSTER_DIR', str(DATA_OUT / '06_clusters')))
WWA_DIR = Path(
    os.environ.get(
        "CISC_WWA_DIR",
        str(ROOT / "figures/fig_02_wwa"),
    )
)
ZONE_SHP = Path(os.environ.get('CISC_ZONE_SHP', str(ROOT / 'data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp')))
COUNTY_SHP_TEXT = os.environ.get('CISC_COUNTY_SHP', '').strip()
COUNTY_SHP = Path(COUNTY_SHP_TEXT) if COUNTY_SHP_TEXT else None
PERF_BASE = Path(os.environ.get('CISC_PERF_DIR', str(DATA_OUT / '13_performance')))
CONTROL = PERF_BASE / '_control'
LOOKUP_XY = PERF_BASE / 'lookup_xy_to_county'
LOOKUP_ZONE = PERF_BASE / 'lookup_zone_to_county'
OBS_DIR = PERF_BASE / 'wwa_county_days_definitive'
MODELED_DIR = PERF_BASE / 'modeled_county_daily'
MET_DIR = PERF_BASE / 'metrics_county'
PART_DIR = MET_DIR / 'parts'
QA_DIR = PERF_BASE / 'qa_wwa_eval'
TMP_DIR = PERF_BASE / '_tmp'

REGISTRY = CONTROL / 'dataset_registry.csv'
LOOKUP_MANIFEST = CONTROL / 'lookup_manifest.csv'
MODELED_MANIFEST = CONTROL / 'modeled_manifest.csv'
EVAL_MANIFEST = CONTROL / 'eval_manifest.csv'
COUNTY_CACHE = CONTROL / 'counties_2022_cb.gpkg'
ZONE_COUNTY = LOOKUP_ZONE / 'lookup_zone_to_county_bbox.parquet'
OBS_LOOKUP = OBS_DIR / 'wwa_county_day_lookup.parquet'

COUNTY_ZIP_CACHE = (
    CONTROL
    / "cb_2022_us_county_500k.zip"
)

COUNTY_URL = os.environ.get(
    "CISC_COUNTY_URL",
    (
        "https://www2.census.gov/geo/tiger/"
        "GENZ2022/shp/cb_2022_us_county_500k.zip"
    ),
)

MODE = os.environ.get('CISC_WWA_EVAL_MODE', 'setup').strip().lower()
WORKER_INDEX = int(os.environ.get('SLURM_ARRAY_TASK_ID', os.environ.get('CISC_WWA_WORKER_INDEX', '0')))
WORKER_COUNT = int(os.environ.get('CISC_WWA_WORKER_COUNT', '4'))
OVERWRITE_LOOKUPS = os.environ.get('CISC_OVERWRITE_LOOKUPS', 'false').lower() in {'1','true','yes'}
OVERWRITE_OBSERVED = os.environ.get('CISC_OVERWRITE_OBSERVED', 'false').lower() in {'1','true','yes'}
OVERWRITE_MODELED = os.environ.get('CISC_OVERWRITE_MODELED', 'false').lower() in {'1','true','yes'}
OVERWRITE_METRICS = os.environ.get('CISC_OVERWRITE_METRICS', 'true').lower() in {'1','true','yes'}
EVAL_INDIVIDUAL = os.environ.get('CISC_EVALUATE_INDIVIDUAL_WWA_TAGS', 'false').lower() in {'1','true','yes'}
YEAR_START = int(os.environ.get('CISC_EVAL_YEAR_START', '2002'))
YEAR_END = int(os.environ.get('CISC_EVAL_YEAR_END', '2023'))
START = pd.Timestamp(f'{YEAR_START}-01-01')
END = pd.Timestamp(f'{YEAR_END}-12-31')
BBOX = tuple(float(v) for v in os.environ.get('CISC_EVAL_BBOX', '-95,22,-75,40').split(','))
SCENARIOS = tuple(
    value.strip()
    for value in os.environ.get(
        "CISC_EVAL_SCENARIOS",
        "filtered_rp5,filtered_rp10,filtered_rp25,filtered_rp50",
    ).split(",")
    if value.strip()
)

EXPECTED_SCENARIOS = {
    "filtered_rp5",
    "filtered_rp10",
    "filtered_rp25",
    "filtered_rp50",
}

if set(SCENARIOS) != EXPECTED_SCENARIOS:
    raise ValueError(
        "CISC_EVAL_SCENARIOS must contain exactly "
        "filtered_rp5, filtered_rp10, "
        "filtered_rp25, and filtered_rp50"
    )
R_LIBS_USER = os.environ.get('R_LIBS_USER', '/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5')
RSCRIPT = os.environ.get('CISC_RSCRIPT', 'Rscript')

WWA_GROUPS = {
    'all_heat': ('EH_W',),
    'all_precip': ('SV_W','FF_W','HI_W','HU_W','TI_W','TR_W'),
}
HEAT_PHEN = {'EH'}
# PRECIP_PHEN = {'SV','FF','HI','HU','TI','TR'}
# VALID_KEYS = {
#     'eh_a','eh_w','ht_a','ht_w','bz_a','bz_w','cf_a','cf_w','df_a','df_w',
#     'fa_a','fa_w','ff_a','ff_w','fl_a','fl_w','hi_a','hi_w','hs_a','hs_w',
#     'hu_a','hu_w','sr_a','sr_w','sv_a','sv_w','tr_a','tr_w','ty_a','ty_w',
#     'ws_a','ws_w',
# }
PRECIP_PHEN = {'HI','HU','TI','TR'}
VALID_KEYS = {
    'eh_w','hi_w','hu_w','ti_w','tr_w',
}

FIXED_WINDOWS = (
    ('20020101_20231231', pd.Timestamp('2002-01-01'), pd.Timestamp('2023-12-31')),
    ('20140101_20231231', pd.Timestamp('2014-01-01'), pd.Timestamp('2023-12-31')),
    ('20190101_20231231', pd.Timestamp('2019-01-01'), pd.Timestamp('2023-12-31')),
)


def log(msg: str) -> None:
    print(f'[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}', flush=True)

def mkdirs() -> None:
    for p in (CONTROL, LOOKUP_XY, LOOKUP_ZONE, OBS_DIR, MODELED_DIR, MET_DIR, PART_DIR, QA_DIR, TMP_DIR):
        p.mkdir(parents=True, exist_ok=True)


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.csv')
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_parquet(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f'{path.stem}.', suffix='.tmp.parquet', dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        df.to_parquet(tmp, index=False, engine='pyarrow', compression='zstd')
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def geoid5(values) -> pd.Series:
    s = pd.Series(values, dtype='string').str.extract(r'(\d+)', expand=False)
    return s.where(s.str.len().le(5)).str.zfill(5)


def var_group(tag: str) -> str:
    if tag.lower().startswith('heatindex'):
        return 'heatindex'
    if tag.lower().startswith('precipitation'):
        return 'precipitation'
    raise ValueError(f'Cannot infer variable from dataset tag: {tag}')


def scenario_rp(name: str) -> float:
    m = re.fullmatch(r'filtered_rp(\d+(?:p\d+)?)', name)
    if not m:
        raise ValueError(f'Unsupported scenario: {name}')
    return float(m.group(1).replace('p','.'))


def key_var_group(key: str) -> str | None:
    k = key.lower()
    if k == 'all_heat': return 'heatindex'
    if k == 'all_precip': return 'precipitation'
    m = re.fullmatch(r'([a-z]{2})_[a-z]', k)
    if not m: return None
    p = m.group(1).upper()
    if p in HEAT_PHEN: return 'heatindex'
    if p in PRECIP_PHEN: return 'precipitation'
    return None


def assigned(df: pd.DataFrame) -> pd.DataFrame:
    if WORKER_COUNT < 1 or not 0 <= WORKER_INDEX < WORKER_COUNT:
        raise ValueError(f'Invalid worker {WORKER_INDEX}/{WORKER_COUNT}')
    return df.iloc[WORKER_INDEX::WORKER_COUNT].copy()


def parse_wwa_tag(path: Path) -> str | None:
    m = re.match(r'^wwa_freq_([A-Z]{2})_([A-Z])_.+\.qs$', path.name)
    return f'{m.group(1)}_{m.group(2)}'.lower() if m else None


def qs_to_parquet(qs_path: Path, parquet_path: Path) -> Path:
    if parquet_path.exists() and not OVERWRITE_OBSERVED:
        return parquet_path
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    code = ('suppressPackageStartupMessages({library(qs);library(data.table);library(arrow)});'
            'x<-as.data.table(qread(commandArgs(TRUE)[1]));'
            'write_parquet(x,commandArgs(TRUE)[2],compression="zstd")')
    env = os.environ.copy(); env['R_LIBS_USER'] = R_LIBS_USER
    cp = subprocess.run([RSCRIPT,'-e',code,str(qs_path),str(parquet_path)], env=env, capture_output=True, text=True)
    if cp.returncode:
        raise RuntimeError(f'Failed converting {qs_path}\n{cp.stdout}\n{cp.stderr}')
    return parquet_path


def discover_registry() -> pd.DataFrame:
    rows = []
    for scenario in SCENARIOS:
        sdir = CLUSTER_BASE / scenario
        if not sdir.exists():
            log(f'Missing scenario; skipping: {sdir}')
            continue
        for ddir in sorted(p for p in sdir.iterdir() if p.is_dir()):
            files = sorted(ddir.glob('cluster_rp*_*.parquet'))
            if files:
                rows.append(dict(scenario=scenario, scenario_rp=scenario_rp(scenario), dataset_tag=ddir.name,
                                 var_group=var_group(ddir.name), dataset_dir=str(ddir), n_cluster_files=len(files)))
    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError(f'No cluster Parquet files found under {CLUSTER_BASE}')
    out = (
        out.sort_values(
            [
                "scenario_rp",
                "var_group",
                "dataset_tag",
            ]
        )
        .reset_index(drop=True)
    )

    counts = (
        out.groupby(
            [
                "scenario",
                "var_group",
            ]
        )["dataset_tag"]
        .nunique()
    )

    invalid = counts[counts.ne(5)]

    if not invalid.empty:
        raise RuntimeError(
            "Expected five resolution/resampling datasets "
            "for each variable and return period; found:\n"
            + invalid.to_string()
        )

    return out


def download_county_shapefile() -> Path:
    if COUNTY_SHP is not None:
        if not COUNTY_SHP.exists():
            raise FileNotFoundError(
                f"Missing county shapefile: {COUNTY_SHP}"
            )

        return COUNTY_SHP

    if COUNTY_CACHE.exists():
        return COUNTY_CACHE

    COUNTY_ZIP_CACHE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    log(
        "County cache missing; downloading Census "
        f"county boundaries from {COUNTY_URL}"
    )

    response = requests.get(
        COUNTY_URL,
        timeout=300,
    )

    response.raise_for_status()

    temporary_zip = COUNTY_ZIP_CACHE.with_suffix(
        ".tmp.zip"
    )

    temporary_zip.write_bytes(
        response.content
    )

    os.replace(
        temporary_zip,
        COUNTY_ZIP_CACHE,
    )

    extract_dir = (
        CONTROL
        / "cb_2022_us_county_500k"
    )

    extract_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    with zipfile.ZipFile(
        COUNTY_ZIP_CACHE,
        "r",
    ) as archive:
        archive.extractall(
            extract_dir
        )

    shapefiles = sorted(
        extract_dir.glob(
            "*.shp"
        )
    )

    if not shapefiles:
        raise RuntimeError(
            "County ZIP contained no shapefile"
        )

    source_shp = shapefiles[0]

    counties = gpd.read_file(
        source_shp
    )

    counties.to_file(
        COUNTY_CACHE,
        driver="GPKG",
    )

    log(
        f"Wrote county cache: {COUNTY_CACHE}"
    )

    return COUNTY_CACHE


def prepare_counties() -> gpd.GeoDataFrame:
    county_path = download_county_shapefile()

    counties = gpd.read_file(
        county_path
    )

    required = {
        "GEOID",
        "STATEFP",
    }

    missing = required.difference(
        counties.columns
    )

    if missing:
        raise ValueError(
            "County data missing columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    counties = counties[
        ~counties[
            "STATEFP"
        ].astype(str).isin(
            {
                "02",
                "15",
                "60",
                "66",
                "69",
                "72",
                "78",
            }
        )
    ].copy()

    counties["GEOID"] = geoid5(
        counties["GEOID"]
    )

    counties = counties[
        counties["GEOID"].notna()
    ].copy()

    counties = counties.to_crs(
        4326
    )

    counties = counties[
        counties.geometry.intersects(
            box(*BBOX)
        )
    ].copy()

    counties = counties[
        [
            "GEOID",
            "geometry",
        ]
    ].copy()

    counties.geometry = (
        counties.geometry.make_valid()
    )

    return counties


def prepare_zones() -> gpd.GeoDataFrame:
    if not ZONE_SHP.exists(): raise FileNotFoundError(ZONE_SHP)
    g = gpd.read_file(ZONE_SHP)
    if 'STATE_ZONE' not in g.columns: raise ValueError('ZONE_SHP missing STATE_ZONE')
    g['STATE_ZONE'] = g['STATE_ZONE'].astype(str).str.strip().str.upper()
    g = g[g['STATE_ZONE'].ne('')].to_crs(4326)
    g.geometry = g.geometry.make_valid()
    return g[g.geometry.intersects(box(*BBOX))][['STATE_ZONE','geometry']].copy()


def build_zone_county(counties, zones) -> pd.DataFrame:
    if ZONE_COUNTY.exists() and not OVERWRITE_LOOKUPS:
        return pd.read_parquet(ZONE_COUNTY)
    j = gpd.sjoin(zones.to_crs(5070), counties.to_crs(5070), how='inner', predicate='intersects')
    out = j[['STATE_ZONE','GEOID']].drop_duplicates().sort_values(['STATE_ZONE','GEOID']).reset_index(drop=True)
    atomic_parquet(out, ZONE_COUNTY); atomic_csv(out, ZONE_COUNTY.with_suffix('.csv'))
    return out


def build_observed(z2c: pd.DataFrame) -> pd.DataFrame:
    files = []
    for p in sorted(WWA_DIR.glob('wwa_freq_*.qs')):
        tag = parse_wwa_tag(p)
        if tag: files.append((tag,p))
    if not files: raise RuntimeError(f'No WWA .qs files in {WWA_DIR}')
    tag_paths = {}
    for tag, qs_path in files:
        out_path = OBS_DIR / f'wwa_county_days_{tag}.parquet'; tag_paths[tag] = out_path
        if out_path.exists() and not OVERWRITE_OBSERVED: continue
        cache = TMP_DIR / 'wwa_qs_cache' / f'{qs_path.stem}.parquet'
        qs_to_parquet(qs_path, cache)
        x = pd.read_parquet(cache)
        if not {'STATE_ZONE','date'} <= set(x.columns):
            atomic_parquet(pd.DataFrame(columns=['event_tag','GEOID','date','observed']), out_path); continue
        x = x[['STATE_ZONE','date']].copy()
        x['STATE_ZONE'] = x['STATE_ZONE'].astype(str).str.strip().str.upper()
        x['date'] = pd.to_datetime(x['date'], errors='coerce').dt.normalize()
        x = x[x['date'].between(START,END)].merge(z2c,on='STATE_ZONE',how='inner')
        x = x[['GEOID','date']].drop_duplicates().assign(event_tag=tag, observed=np.int8(1))
        x = x[['event_tag','GEOID','date','observed']].sort_values(['date','GEOID'])
        atomic_parquet(x,out_path); atomic_csv(x,out_path.with_suffix('.csv'))
    parts = []
    if EVAL_INDIVIDUAL:
        for tag,p in tag_paths.items():
            if tag in VALID_KEYS and p.exists():
                x = pd.read_parquet(p)
                if not x.empty:
                    parts.append(x[['GEOID','date']].drop_duplicates().assign(key=tag,event_tag=tag,event_group=pd.NA,observed=np.int8(1)))
    for group,tags in WWA_GROUPS.items():
        sub = []
        for tag in tags:
            p = tag_paths.get(tag.lower())
            if p and p.exists():
                x = pd.read_parquet(p)
                if not x.empty: sub.append(x[['GEOID','date']])
        if sub:
            g = pd.concat(sub,ignore_index=True).drop_duplicates().assign(key=group,event_tag=pd.NA,event_group=group,observed=np.int8(1))
        else:
            g = pd.DataFrame(columns=['GEOID','date','key','event_tag','event_group','observed'])
        gp = OBS_DIR / f'wwa_county_days__{group}.parquet'
        atomic_parquet(g,gp); atomic_csv(g,gp.with_suffix('.csv'))
        if not g.empty: parts.append(g)
    lookup = pd.concat(parts,ignore_index=True).drop_duplicates() if parts else pd.DataFrame(columns=['key','event_tag','event_group','GEOID','date','observed'])
    lookup = lookup.sort_values(['key','GEOID','date']) if not lookup.empty else lookup
    atomic_parquet(lookup,OBS_LOOKUP); atomic_csv(lookup,OBS_LOOKUP.with_suffix('.csv'))
    return lookup


def write_manifests(reg: pd.DataFrame, obs: pd.DataFrame) -> None:
    atomic_csv(reg,REGISTRY)
    lm = reg[['dataset_tag','var_group']].drop_duplicates().sort_values(['var_group','dataset_tag']).reset_index(drop=True)
    lm.insert(0,'task_index',range(len(lm))); atomic_csv(lm,LOOKUP_MANIFEST)
    mm = reg.copy(); mm.insert(0,'task_index',range(len(mm)))
    mm['lookup_path'] = mm['dataset_tag'].map(lambda t: str(LOOKUP_XY/f'lookup_xy_to_county_{t}.parquet'))
    mm['output_path'] = mm.apply(lambda r: str(MODELED_DIR/r['scenario']/f"modeled_county_daily_{r['dataset_tag']}.parquet"),axis=1)
    atomic_csv(mm,MODELED_MANIFEST)
    keys = sorted(obs['key'].dropna().astype(str).unique())
    if not EVAL_INDIVIDUAL: keys = [k for k in keys if k in {'all_heat','all_precip'}]
    rows=[]
    for key in keys:
        vg=key_var_group(key)
        if not vg: continue
        for r in mm[mm['var_group'].eq(vg)].itertuples(index=False):
            rows.append(dict(event=key,var_group=vg,scenario=r.scenario,scenario_rp=r.scenario_rp,dataset_tag=r.dataset_tag,
                             modeled_path=r.output_path,part_path=str(PART_DIR/f'part_{vg}_{key}__{r.scenario}__{r.dataset_tag}.parquet')))
    em=pd.DataFrame(rows); em.insert(0,'task_index',range(len(em))); atomic_csv(em,EVAL_MANIFEST)
    (QA_DIR/'preflight_summary.txt').write_text(f'Registry rows: {len(reg)}\nLookup tasks: {len(lm)}\nModeled tasks: {len(mm)}\nEvaluation tasks: {len(em)}\nKeys: {keys}\n')
    log(f'lookup={len(lm)} modeled={len(mm)} eval={len(em)}')


def build_xy_lookup(tag: str) -> Path:
    out = LOOKUP_XY / f'lookup_xy_to_county_{tag}.parquet'
    if out.exists() and not OVERWRITE_LOOKUPS:
        return out

    src = CLUSTER_BASE / "clustered_rp5" / tag
    parts = []

    for p in sorted(src.glob("cluster_rp5_*.parquet")):
        m = re.search(r'(\d{4})$', p.stem)
        if not m or not YEAR_START <= int(m.group(1)) <= YEAR_END:
            continue

        x = pd.read_parquet(
            p,
            columns=[
                "x",
                "y",
                "cluster",
            ],
        )

        x["cluster"] = pd.to_numeric(
            x["cluster"],
            errors="coerce",
        )

        x = x[
            x["cluster"].gt(0)
        ][
            [
                "x",
                "y",
            ]
        ].copy()

        x = (
            x.apply(
                pd.to_numeric,
                errors="coerce",
            )
            .dropna()
            .drop_duplicates()
        )

        parts.append(x)

    if not parts:
        raise RuntimeError(
            f"No positive clustered RP5 coordinates found for {tag}"
        )

    xy=pd.concat(parts,ignore_index=True).drop_duplicates().sort_values(['y','x']).reset_index(drop=True)
    xy['x_key']=xy['x'].round(8); xy['y_key']=xy['y'].round(8)
    pts=gpd.GeoDataFrame(xy,geometry=gpd.points_from_xy(xy.x,xy.y),crs=4326)
    counties=prepare_counties()
    j=gpd.sjoin(pts,counties,how='left',predicate='within')
    miss=j['GEOID'].isna()
    if miss.any():
        b=gpd.sjoin(pts.loc[miss],counties,how='left',predicate='intersects').dropna(subset=['GEOID']).reset_index().drop_duplicates('index')
        j.loc[b['index'],'GEOID']=b['GEOID'].to_numpy()
    o=j[['x','y','x_key','y_key','GEOID']].dropna(subset=['GEOID']).drop_duplicates(['x_key','y_key']).sort_values(['y','x']).reset_index(drop=True)
    o['GEOID']=geoid5(o['GEOID']); atomic_parquet(o,out); atomic_csv(o,out.with_suffix('.csv')); return out


def run_lookup() -> None:
    m=pd.read_csv(LOOKUP_MANIFEST); jobs=assigned(m); log(f'lookup worker {WORKER_INDEX}: {len(jobs)}/{len(m)}')
    for i,r in enumerate(jobs.itertuples(index=False),1):
        log(f'lookup {i}/{len(jobs)} {r.dataset_tag}'); build_xy_lookup(r.dataset_tag)


def build_modeled(row: pd.Series) -> Path:
    out=Path(row['output_path'])
    if out.exists() and not OVERWRITE_MODELED: return out
    lk=pd.read_parquet(row['lookup_path'])[['x_key','y_key','GEOID']].drop_duplicates()
    parts=[]
    for p in sorted(Path(row['dataset_dir']).glob('cluster_rp*_*.parquet')):
        m=re.search(r'(\d{4})$',p.stem)
        if not m or not YEAR_START <= int(m.group(1)) <= YEAR_END: continue
        x=pd.read_parquet(p,columns=['x','y','time_utc','cluster','cluster_id'])
        x['cluster']=pd.to_numeric(x['cluster'],errors='coerce'); x=x[x['cluster'].gt(0)&x['cluster_id'].notna()].copy()
        if x.empty: continue
        x['date']=pd.to_datetime(x['time_utc'],errors='coerce',utc=True).dt.tz_convert(None).dt.normalize(); x=x[x['date'].between(START,END)]
        x['x_key']=pd.to_numeric(x['x'],errors='coerce').round(8); x['y_key']=pd.to_numeric(x['y'],errors='coerce').round(8)
        x=x.merge(lk,on=['x_key','y_key'],how='inner')
        if not x.empty: parts.append(x[['GEOID','date']].drop_duplicates())
    o=pd.concat(parts,ignore_index=True).drop_duplicates().assign(modeled=np.int8(1)).sort_values(['date','GEOID']) if parts else pd.DataFrame(columns=['GEOID','date','modeled'])
    atomic_parquet(o,out); return out


def run_modeled() -> None:
    m=pd.read_csv(MODELED_MANIFEST); jobs=assigned(m); log(f'modeled worker {WORKER_INDEX}: {len(jobs)}/{len(m)}')
    for i,(_,r) in enumerate(jobs.iterrows(),1):
        log(f"modeled {i}/{len(jobs)} {r['scenario']} {r['dataset_tag']}"); build_modeled(r)


def window_table() -> pd.DataFrame:
    rows=[]
    for y in range(YEAR_START,YEAR_END+1):
        s=pd.Timestamp(f'{y}-01-01'); e=pd.Timestamp(f'{y}-12-31'); rows.append(dict(window_tag=str(y),window_type='year',year=y,start=s,end=e,total_days=(e-s).days+1))
    for tag,s,e in FIXED_WINDOWS: rows.append(dict(window_tag=tag,window_type='fixed',year=pd.NA,start=s,end=e,total_days=(e-s).days+1))
    return pd.DataFrame(rows)


def attach_windows(x: pd.DataFrame) -> pd.DataFrame:
    if x.empty: return pd.DataFrame(columns=['GEOID','date','window_tag'])
    parts=[]; a=x.copy(); a['window_tag']=a['date'].dt.year.astype(str); parts.append(a[['GEOID','date','window_tag']])
    for tag,s,e in FIXED_WINDOWS:
        z=x[x['date'].between(s,e)].copy(); z['window_tag']=tag; parts.append(z[['GEOID','date','window_tag']])
    return pd.concat(parts,ignore_index=True)


def evaluate(row: pd.Series) -> Path:
    out=Path(row['part_path'])
    if out.exists() and not OVERWRITE_METRICS: return out
    obs=pd.read_parquet(OBS_LOOKUP); obs=obs[obs['key'].eq(row['event'])][['GEOID','date']].drop_duplicates()
    mp=Path(row['modeled_path']); mod=pd.read_parquet(mp)[['GEOID','date']].drop_duplicates() if mp.exists() else pd.DataFrame(columns=['GEOID','date'])
    for x in (obs,mod): x['GEOID']=geoid5(x['GEOID']); x['date']=pd.to_datetime(x['date'],errors='coerce').dt.normalize()
    obs=obs.dropna().drop_duplicates(); mod=mod.dropna().drop_duplicates(); counties=sorted(prepare_counties()['GEOID'].dropna().unique()); wins=window_table()
    base=pd.MultiIndex.from_product([counties,wins['window_tag']],names=['GEOID','window_tag']).to_frame(index=False)
    ow=attach_windows(obs); mw=attach_windows(mod)
    on=ow.groupby(['GEOID','window_tag']).size().rename('obs_n').reset_index(); mn=mw.groupby(['GEOID','window_tag']).size().rename('mod_n').reset_index()
    inter=obs.merge(mod,on=['GEOID','date'],how='inner').drop_duplicates() if not obs.empty and not mod.empty else pd.DataFrame(columns=['GEOID','date'])
    tw=attach_windows(inter); tn=tw.groupby(['GEOID','window_tag']).size().rename('tp_n').reset_index()
    c=base.merge(on,how='left').merge(mn,how='left').merge(tn,how='left').fillna({'obs_n':0,'mod_n':0,'tp_n':0}).merge(wins,on='window_tag',how='left')
    for q in ('obs_n','mod_n','tp_n'): c[q]=c[q].astype(np.int64)
    TP=c['tp_n']; FP=(c['mod_n']-TP).clip(lower=0); FN=(c['obs_n']-TP).clip(lower=0); TN=(c['total_days']-TP-FP-FN).clip(lower=0); N=TP+FP+FN+TN
    div=lambda a,b: a/b.replace(0,np.nan)
    precision=div(TP,TP+FP); recall=div(TP,TP+FN); specificity=div(TN,TN+FP); accuracy=div(TP+TN,N); fpr=div(FP,FP+TN); fnr=div(FN,FN+TP)
    f1=div(2*precision*recall,precision+recall); f2=div(5*precision*recall,4*precision+recall); auc=(recall+specificity)/2
    met=pd.DataFrame(dict(TP=TP,FP=FP,FN=FN,TN=TN,N=N,precision=precision,recall=recall,specificity=specificity,accuracy=accuracy,f1=f1,f2=f2,fpr=fpr,fnr=fnr,auc=auc))
    o=pd.concat([c[['GEOID','window_tag','window_type','year','start','end']].reset_index(drop=True),met],axis=1)
    for name,val in reversed([('scenario',row['scenario']),('scenario_rp',row['scenario_rp']),('dataset_tag',row['dataset_tag']),('var_group',row['var_group']),('event',row['event'])]): o.insert(0,name,val)
    atomic_parquet(o,out); return out


def run_eval() -> None:
    m=pd.read_csv(EVAL_MANIFEST); jobs=assigned(m); log(f'eval worker {WORKER_INDEX}: {len(jobs)}/{len(m)}')
    for i,(_,r) in enumerate(jobs.iterrows(),1):
        log(f"eval {i}/{len(jobs)} {r['event']} {r['scenario']} {r['dataset_tag']}"); evaluate(r)


def summarize() -> None:
    m=pd.read_csv(EVAL_MANIFEST); paths=[Path(p) for p in m['part_path']]; missing=[p for p in paths if not p.exists()]
    if missing: raise RuntimeError(f'Missing {len(missing)} metric parts; first={missing[:10]}')
    x=pd.concat([pd.read_parquet(p) for p in paths],ignore_index=True).sort_values(['var_group','event','scenario_rp','dataset_tag','window_tag','GEOID'])
    for (vg,event),g in x.groupby(['var_group','event'],dropna=False):
        p=MET_DIR/f'metrics_{vg}_{event}_multiwindow.parquet'; atomic_parquet(g,p); atomic_csv(g,p.with_suffix('.csv'))
    s=x.groupby(['var_group','event','scenario','scenario_rp','dataset_tag','window_tag'],dropna=False).agg(n=('GEOID','size'),mean_f1=('f1','mean'),mean_auc=('auc','mean'),mean_precision=('precision','mean'),mean_recall=('recall','mean')).reset_index()
    atomic_parquet(s,MET_DIR/'metrics_summary_all.parquet'); atomic_csv(s,MET_DIR/'metrics_summary_all.csv'); log('summary written')


def setup() -> None:
    reg=discover_registry(); counties=prepare_counties(); zones=prepare_zones(); z2c=build_zone_county(counties,zones); obs=build_observed(z2c); write_manifests(reg,obs)


def main() -> None:
    mkdirs(); log(f'MODE={MODE} ROOT={ROOT} SCENARIOS={SCENARIOS} YEARS={YEAR_START}-{YEAR_END}')
    if MODE=='setup': setup()
    elif MODE=='lookup_worker': run_lookup()
    elif MODE=='modeled_worker': run_modeled()
    elif MODE=='eval_worker': run_eval()
    elif MODE=='summarize': summarize()
    else: raise ValueError('Mode must be setup, lookup_worker, modeled_worker, eval_worker, or summarize')
    log('Done.')

if __name__ == '__main__':
    main()
