# One-control-donor 5hmC validation

bcontrol1 is held out; the other three donors provide the ONT reference.
The paper FANS benchmark catalog is the sole modifiedC backbone. Both sorted
FANS populations provide reference labels for precision and recall.

Read [the blog update](BLOG_UPDATE.md) and [implementation details](IMPLEMENTATION.md).
The genome-wide fixed-control identity audit and regional historical-code
compatibility checks passed. This is a single-donor proof; other folds,
bulk benchmarking, and single-channel ablations remain deferred.

Reproduce the three bar-plot figures from the included aggregates:

```bash
python -m pip install -r requirements.txt
python scripts/plot_benchmarks.py
```

PNG and PDF outputs are written to `results/`. Precision and recall are
separate; all controls, including no-evidence and ambiguous reads, remain
in the recall denominator. Error bars are not shown because this comparison
contains one biological donor.

The full blog can be refreshed from an authorized local result tree using
`python scripts/render_update.py /path/to/validation-results`. It expects the
private fixed-control and genome-control summaries; individual reads and
server paths are not exported.
