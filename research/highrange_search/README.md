# High-dimensional search diagnostics

These numerical angular and radial search tools generated candidates used during the 45–47D investigation. They are research programs, not proof verifiers. Install requirements.txt; the three scripts expose --help. The angular LP is dimension-specific; radial programs start from the supplied frozen 48D radial/angular inputs. Generated coefficient files, sampled sign tests and optimizer status never establish a density bound.

The retained numerical outputs are historical checkpoints and may fail global signs. They are not substitutes for the final rational certificates in ../../verifications/dim45, dim46 and dim47. In particular, the successful 45D witness required a new rational angular polynomial and additional exact-root-guided separation rounds, as explained in its proof. This package preserves the generic candidate-generation code and selected checkpoints, not a claim to reproduce the entire adaptive discovery chain from an empty initialization. The full exact verification of final witnesses is self-contained in the registered proof packages.

Run searches in a writable copy: they create output JSON and optional numerical cache files in this directory. No cache files are proof inputs.
