try:
    import veros  
except ImportError:
    raise RuntimeError(
        'Plugin needs Veros to be installed (try `pip install veros`)'
    )

__version__ ="0.1.0" 

from veros_global_realistic.real import real
from veros_global_realistic.set_inits import set_inits
from veros_global_realistic.variables import VARIABLES
from veros_global_realistic.settings import SETTINGS


__VEROS_INTERFACE__ = dict(
    name = 'veros_global_realistic',
    setup_entrypoint = set_inits,
    run_entrypoint = real,
    settings = SETTINGS,
    variables = VARIABLES,
)
