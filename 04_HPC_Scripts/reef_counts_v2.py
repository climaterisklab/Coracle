"""Reduce the v2 per-cell bleaching output (all reef locations) to per-country event counts.

For each threshold, cell set, country (cells.csv `admin`), year, scenario and sample:
  factual = sum of cell weights over cells with pred_factual > thr
  cf      = sum of cell weights over those same cells that also have pred_cf > thr   (conditional counting,
            identical to get_batch_counts_corrected_v2 in 06_Figure3.R / 10_TablesSupp7-10.R)
Weighted counts are additive over cells, so any country grouping (Global, AOSIS, USA + territories ...)
is a sum over countries and can be done downstream in R.
Cell sets (cells: all modelled reef cells, or only those holding survey sites = cells.csv in_sites):
  "reefs_area" / "sites_area": weight cos(latitude) = proportional to the 5 km cell's area, i.e. coral reef area
                               (agreed with Chris and Max for the figures and tables; same as Max's
                               summary/attribution_{reefs,sites}.csv)
  "reefs" / "sites":           weight 1 (plain cell counts; Max's attribution_{reefs,sites}_cells.csv)

usage: reef_counts_v2.py <seed> [<seed> ...]     -> writes OUT_DIR/counts_s<seed>.npz
       reef_counts_v2.py --combine               -> writes OUT_DIR/reef_counts_v2.nc
"""
import os, sys, glob
import numpy as np
import pandas as pd
import xarray as xr

IN_DIR = "/esarchive/Coracle/bleaching_attribution_v2/reef_locations"
CELLS = "/esarchive/Coracle/climate_data/reef_locations/cells.csv"
OUT_DIR = os.environ.get("OUT_DIR", "/gpfs/scratch/bsc32/bsc963510/coracle_fig3/v4_area")
CELLSETS = ["reefs_area", "sites_area", "reefs", "sites"]
THRESHOLDS = [10, 30]
N_CHUNKS = 16

cells = pd.read_csv(CELLS)
cells = cells[cells.chunk >= 0]
countries = np.array(sorted(cells.admin.unique()))
admin_of = cells.set_index("coord_id").admin
in_sites = cells.in_sites.astype(bool).values
coslat = np.cos(np.deg2rad(cells.lat.values))
# per-cell weight of each cell set (rows follow CELLSETS)
weights = pd.DataFrame({"reefs_area": coslat, "sites_area": coslat * in_sites,
                        "reefs": 1.0, "sites": in_sites * 1.0}, index=cells.coord_id)[CELLSETS]


def onehot(cell_ids):
    idx = np.searchsorted(countries, admin_of.loc[cell_ids].values)
    m = np.zeros((len(countries), len(cell_ids)))
    m[idx, np.arange(len(cell_ids))] = 1
    w = weights.loc[cell_ids].values.T                       # cellset, cell
    return (m[None, :, :] * w[:, None, :]).reshape(-1, len(cell_ids))   # (cellset x country), cell


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
            fact = np.zeros((len(THRESHOLDS), len(CELLSETS), len(countries), ny, nsm))
            cf = np.zeros((len(THRESHOLDS), len(CELLSETS), len(countries), ny, ns, nsm))
        for t, thr in enumerate(THRESHOLDS):
            F = pf > thr                            # NaN -> False (= na.rm in R)
            C = F[:, :, None, :] & (pcf > thr)
            fact[t] += (oh @ F.reshape(len(cell_ids), -1)).reshape(len(CELLSETS), -1, ny, nsm)
            cf[t] += (oh @ C.reshape(len(cell_ids), -1)).reshape(len(CELLSETS), -1, ny, ns, nsm)
        print(f"seed {seed} chunk {c} done", flush=True)
    np.savez_compressed(out + ".tmp.npz", fact=fact, cf=cf, years=years, scenario=scen,
                        sample=samp + 10 * seed, n_cells=n_cells)
    os.replace(out + ".tmp.npz", out)


def combine():
    files = sorted(glob.glob(os.path.join(OUT_DIR, "counts_s*.npz")),
                   key=lambda f: int(f.split("_s")[-1].split(".")[0]))
    assert len(files) == 84, f"only {len(files)} seed files"
    # fill one preallocated float32 array seed by seed (concatenating float64 copies needs > 64 GB)
    p0 = np.load(files[0])
    nsm = p0["fact"].shape[-1]
    fact = np.empty(p0["fact"].shape[:-1] + (nsm * len(files),), np.float32)
    cf = np.empty(p0["cf"].shape[:-1] + (nsm * len(files),), np.float32)
    sample = np.empty(nsm * len(files), np.int64)
    for i, f in enumerate(files):
        with np.load(f) as p:
            sl = slice(i * nsm, (i + 1) * nsm)
            fact[..., sl] = p["fact"]
            cf[..., sl] = p["cf"]
            sample[sl] = p["sample"]
    # N in the tables = number of grid cells; area = sum of the weights (cos-latitude cell equivalents)
    n_cells_country = np.stack([(weights[cs] > 0).groupby(cells.admin.values).sum().reindex(countries, fill_value=0).values
                                for cs in CELLSETS])
    area_country = np.stack([weights[cs].groupby(cells.admin.values).sum().reindex(countries, fill_value=0).values
                             for cs in CELLSETS])
    ds = xr.Dataset(
        {"factual_count": (("threshold", "cellset", "country", "year", "sample"), fact),
         "cf_count": (("threshold", "cellset", "country", "year", "scenario", "sample"), cf),
         "n_cells": (("cellset", "country"), n_cells_country),
         "area": (("cellset", "country"), area_country)},
        coords={"threshold": THRESHOLDS, "cellset": CELLSETS, "country": countries, "year": p0["years"],
                "scenario": p0["scenario"], "sample": sample},
        attrs={"description": "v2_mixed bleaching events per country (cells.csv admin), summed cell weights: "
                              "*_area = cos(latitude) (coral reef area), plain = 1 per cell; reefs = all modelled reef cells, "
                              "sites = cells holding survey sites. cf_count = factual events that also exceed the "
                              "threshold in the counterfactual. n_cells = number of grid cells, area = sum of weights",
               "source": IN_DIR})
    enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in ["factual_count", "cf_count"]}
    ds.to_netcdf(os.path.join(OUT_DIR, "reef_counts_v2.nc"), encoding=enc)
    print(ds)


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    if sys.argv[1] == "--combine":
        combine()
    else:
        for s in sys.argv[1:]:
            do_seed(int(s))
