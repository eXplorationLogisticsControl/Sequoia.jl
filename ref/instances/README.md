# Mixed-integer SDP test instances

The usual MISDP benchmark is the set maintained with [SCIP-SDP](https://www.opt.tu-darmstadt.de/scipsdp/) at TU Darmstadt. It comes from Matter and Pfetsch, [Presolving for mixed-integer semidefinite optimization](https://www2.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving.html) (2021), and is a subset of the 194-instance set in Gally, Pfetsch, and Ulbrich, [A framework for solving mixed-integer semidefinite programs](https://optimization-online.org/wp-content/uploads/2016/04/5394.pdf) (2018).

Files are gzip-compressed sparse SDPA (`.dat-s`). Integer variable indices are listed after an `*INTEGERS*` marker, which a plain SDPA reader treats as a comment. SCIP-SDP also reads Conic Benchmark Format.

| Set | Instances | Archive | Contents |
| --- | ---: | ---: | --- |
| [Truss topology](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/TrussTopology.tgz) | 38 | 196 KB | Discrete bar sizes for a truss that stays stable under an ellipsoidal load set. One SDP block of order about 10–44, plus a diagonal linear block. Binary counts run from 27 to a few hundred. Stored in `TrussTopology/`. |
| [Minimum k-partition](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/minKPartition.tgz) | 32 | 588 KB | Eisenblätter's formulation. The matrix entries are the binary variables, so a 34-node VLSI instance already has 561 of them. The data are very sparse. |
| [RIP](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/RIP.tgz) | 46 | 1.5 MB | Restricted isometry constants, also sparse PCA. About 30–70 binary support variables and one continuous PSD matrix. |
| [Cardinality least squares](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/cardLeastSquares.tgz) | 43 | 639 MB | Sparse regression, including colon-cancer instances with 100 binaries and one dense SDP block of order 63. |
| [Random](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/randomMISDP.tgz) | 26 | 275 MB | Kobayashi–Takano style. SDP order and the binary count are each in {60, 90, 120}. |

The five archives together are [MISDP-Presolving-Instances.tgz](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-presolving/MISDP-Presolving-Instances.tgz) (916 MB), 185 instances in total.

Hojny and Pfetsch added two stable-set families on the [symmetry page](https://www2.mathematik.tu-darmstadt.de/~pfetsch/MISDP-symmetries.html): [Color02 graphs](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-symmetry/StableSetColor02.tgz) (6 MB) and [Flower Snarks](https://wwwopt.mathematik.tu-darmstadt.de/~pfetsch/MISDP-symmetry/StableSetFlowerSnarks.tgz) (2 MB).
