"""Reduce the v2 per-cell bleaching output (all reef locations) to per-country event counts.

For each threshold, country (cells.csv `admin`), year, scenario and sample:
  factual = number of cells with pred_factual > thr
  cf      = number of those same cells that also have pred_cf > thr   (conditional counting,
            identical to get_batch_counts_corrected_v2 in 06_Figure3.R / 10_TablesSupp7-10.R)
Counts are additive over cells, so any country grouping (Global, AOSIS, USA + territories ...)
is a sum over countries and can be done downstream in R.
Three cell sets: "reefs" = all modelled reef cells, "sites" = only cells holding survey sites (cells.csv in_sites),
"sitesweighted" = survey-site cells counted once per survey site (bleaching_survey_data/attribution_sites.csv);
sites in one cell share DHW and province model, so a cell's prediction is each of its sites' prediction.

usage: reef_counts_v2.py <seed> [<seed> ...]     -> writes OUT_DIR/counts_s<seed>.npz
       reef_counts_v2.py --combine               -> writes OUT_DIR/reef_counts_v2.nc
"""
import os, sys, glob
import numpy as np
import pandas as pd
import xarray as xr

IN_DIR = "/esarchive/Coracle/bleaching_attribution_v2/reef_locations"
CELLS = "/esarchive/Coracle/climate_data/reef_locations/cells.csv"
SITES = "/esarchive/Coracle/bleaching_survey_data/attribution_sites.csv"
OUT_DIR = os.environ.get("OUT_DIR", "/gpfs/scratch/bsc32/bsc963510/coracle_fig3/v3_sets")
CELLSETS = ["reefs", "sites", "sitesweighted"]
THRESHOLDS = [10, 30]
N_CHUNKS = 16

cells = pd.read_csv(CELLS)
cells = cells[cells.chunk >= 0]
countries = np.array(sorted(cells.admin.unique()))
admin_of = cells.set_index("coord_id").admin
in_sites_of = cells.set_index("coord_id").in_sites.astype(bool)
# number of survey sites per cell (sites in unmodelled cells, e.g. Tropical Eastern Pacific, drop out)
n_sites_of = pd.read_csv(SITES).groupby("coord_id").size().reindex(cells.coord_id, fill_value=0)


def onehot(cell_ids):
    idx = np.searchsorted(countries, admin_of.loc[cell_ids].values)
    m = np.zeros((len(countries), len(cell_ids)), dtype=np.float32)
    m[idx, np.arange(len(cell_ids))] = 1
    sites = in_sites_of.loc[cell_ids].values.astype(np.float32)
    w = n_sites_of.loc[cell_ids].values.astype(np.float32)
    return np.stack([m, m * sites[None, :], m * w[None, :]])   # cellset, country, cell


def do_seed(seed):
    out = os.path.join(OUT_DIR, f"counts_s{seed}.npz")
    if os.path.exists(out):
        return
    fact = cf = None
    n_cells = 0
    for c in range(N_CHUNKS):
        with xr.open_dataset(os.path.join(IN_DIR, f"bleach_s{seed}_c{c}.nc")) as ds:
            ds = ds.transpose("cell", "year", "scenario", "sample")
            pf = ds.pred_factual.values            # cell, year, sample
            pcf = ds.pred_cf.values                # cell, year, scenario, sample
            cell_ids = ds.cell.values
            years, scen, samp = ds.year.values, ds.scenario.values.astype(str), ds["sample"].values
        oh = onehot(cell_ids)
        n_cells += len(cell_ids)
        if fact is None:
            ny, ns, nsm = len(years), len(scen), len(samp)
            fact = np.zeros((len(THRESHOLDS), len(CELLSETS), len(countries), ny, nsm), np.int32)
            cf = np.zeros((len(THRESHOLDS), len(CELLSETS), len(countries), ny, ns, nsm), np.int32)
        for t, thr in enumerate(THRESHOLDS):
            F = pf > thr                            # NaN -> False (= na.rm in R)
            C = F[:, :, None, :] & (pcf > thr)
            fact[t] += np.rint(oh @ F.reshape(len(cell_ids), -1).astype(np.float32)).astype(np.int32).reshape(len(CELLSETS), -1, ny, nsm)
            cf[t] += np.rint(oh @ C.reshape(len(cell_ids), -1).astype(np.float32)).astype(np.int32).reshape(len(CELLSETS), -1, ny, ns, nsm)
        print(f"seed {seed} chunk {c} done", flush=True)
    np.savez_compressed(out + ".tmp.npz", fact=fact, cf=cf, years=years, scenario=scen,
                        sample=samp + 10 * seed, n_cells=n_cells)
    os.replace(out + ".tmp.npz", out)


def combine():
    files = sorted(glob.glob(os.path.join(OUT_DIR, "counts_s*.npz")),
                   key=lambda f: int(f.split("_s")[-1].split(".")[0]))
    assert len(files) == 84, f"only {len(files)} seed files"
    parts = [np.load(f) for f in files]
    fact = np.concatenate([p["fact"] for p in parts], axis=-1)
    cf = np.concatenate([p["cf"] for p in parts], axis=-1)
    sample = np.concatenate([p["sample"] for p in parts])
    p0 = parts[0]
    n_cells_country = np.stack([cells.admin.value_counts().reindex(countries, fill_value=0).values,
                                cells[cells.in_sites.astype(bool)].admin.value_counts().reindex(countries, fill_value=0).values,
                                n_sites_of.groupby(cells.set_index("coord_id").admin).sum().reindex(countries, fill_value=0).values])
    ds = xr.Dataset(
        {"factual_count": (("threshold", "cellset", "country", "year", "sample"), fact),
         "cf_count": (("threshold", "cellset", "country", "year", "scenario", "sample"), cf),
         "n_cells": (("cellset", "country"), n_cells_country)},
        coords={"threshold": THRESHOLDS, "cellset": CELLSETS, "country": countries, "year": p0["years"],
                "scenario": p0["scenario"], "sample": sample},
        attrs={"description": "v2_mixed bleaching event counts per country (cells.csv admin) for all reef cells (reefs) "
                              "survey-site cells only (sites) and survey-site cells weighted by number of sites (sitesweighted; "
                              "n_cells = number of sites); "
                              "cf_count = factual events that also exceed the threshold in the counterfactual",
               "source": IN_DIR})
    enc = {v: {"zlib": True, "complevel": 4} for v in ["factual_count", "cf_count"]}
    ds.to_netcdf(os.path.join(OUT_DIR, "reef_counts_v2.nc"), encoding=enc)
    print(ds)


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    if sys.argv[1] == "--combine":
        combine()
    else:
        for s in sys.argv[1:]:
            do_seed(int(s))
