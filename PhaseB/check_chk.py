"""Check that the last two AMReX checkpoints of each run are complete (Header, Level_*/ *_H, listed FAB files non-empty).
  python check_chk.py RUN [RUN ...]   -> OK/BAD per chk; move BAD ones aside before restarting (e.g. after a Lustre outage)."""
# A checkpoint is complete if Header exists and every Level_*/ has a Cell_H whose listed data files all exist & are non-empty.
for run in sys.argv[1:]:
    chks = sorted(glob.glob(os.path.join(run, "chk?????")))
    if not chks: print(f"{run}: no chk"); continue
    res = []
    for c in chks[-2:]:
        ok = os.path.exists(os.path.join(c, "Header"))
        levs = sorted(glob.glob(os.path.join(c, "Level_*")))
        ok = ok and len(levs) > 0
        for L in levs:
            for h in glob.glob(os.path.join(L, "*_H")):
                txt = open(h).read().split()
                files = {t.split()[0] for t in open(h).read().splitlines() if t.startswith("FabOnDisk:")}
                for line in open(h).read().splitlines():
                    if line.startswith("FabOnDisk:"):
                        fn = os.path.join(L, line.split()[1])
                        if not os.path.exists(fn) or os.path.getsize(fn) == 0: ok = False
            if not glob.glob(os.path.join(L, "*_H")): ok = False
        res.append(f"{os.path.basename(c)}:{'OK' if ok else 'BAD'}")
    print(f"{os.path.basename(run):48s} {' '.join(res)}")
