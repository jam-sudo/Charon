"""Layer 3 FIH dose benchmark — DEFERRED.

This benchmark script is intentionally empty. FIH dose validation requires
a dataset of compounds with known published MRSD (Maximum Recommended
Starting Dose) from IND/NDA filings — the ``tier2_drugs_at_fda`` dataset.

That dataset has not been built. The directories
``validation/data/tier2_drugs_at_fda/`` and ``validation/data/tier3_literature/``
exist as scaffolds but are currently empty.

**Why not use the Obach panel?** The Obach 1999 panel provides IV PK
parameters (CL, Vss, t1/2) for Layer 2 validation. These old drugs predate
modern FIH dose-finding methodology (ICH M3, FDA 2005 guidance). Their
clinical starting doses are not documented in the MRSD sense.

**When will this be implemented?** After the tier2_drugs_at_fda dataset is
curated from publicly available FDA NDA reviews (Drugs@FDA database) or
published Phase I dose-escalation papers.

See also: README.md Known Limitations, Validation Status.
"""
