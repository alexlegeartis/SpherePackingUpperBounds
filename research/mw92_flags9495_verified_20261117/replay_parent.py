from pathlib import Path
from fractions import Fraction as F
import json,numpy as np
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
D=Path(__file__).resolve().parent;P=D.parent/'mw92_sections_20261115'
src=(P/'verify.py').read_text();exec(src[src.index('def parent():'):src.index('def rogers_comparison_intervals():')])
G,result=parent();assert G.det()==47**14
(D/'PARENT_RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
