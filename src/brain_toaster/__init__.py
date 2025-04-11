"""Python port of g_ci_spm. It calculates the estimation of effect size g and its confidence interval from SPM t maps."""

from .utils import prepare_contrast, compute_effect_size, load_design_matrix  
from .parallel import compute_ci_for_all
from .core import es_ci_spm_parallel

__version__ = "0.1.0"

__all__ = ['es_ci_spm_parallel']